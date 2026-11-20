! ═══════════════════════════════════════════════════════════════════════
!  tycho_fit — Sun–Earth–Moon orbit determination from Tycho Brahe's
!  raw naked-eye readings of 1586 and 1595 (Hven / Uraniborg).
!
!  Input (written by export_fit_data.py in the tycho repository):
!    <obs file>    one line per reading: moment (TT), the rotation ICRS ->
!                  true equator of date, local apparent sidereal time, the
!                  kind of reading (declination, altitude, meridian
!                  altitude, hour angle, azimuth, distance to a star or to
!                  the Sun, RA difference, diameter), the limb, the value,
!                  and JPL DE441's prediction of it (for validation only).
!    <state file>  JPL DE441 barycentric Sun/Earth/Moon states at the
!                  epoch (used only as the "truth" to compare with, and to
!                  build the synthetic data).
!
!  Model.  14 orbital parameters at the epoch (1590-07-01 TDB):
!    Earth (heliocentric, ICRS): e, i, Omega, omega, M0, n, a
!    Moon  (geocentric,   ICRS): e, i, Omega, omega, M0, n, a
!  with GM_sun + GM_earth = n_E^2 a_E^3 and GM_earth + GM_moon =
!  n_M^2 a_M^3 (moon/earth mass ratio fixed at 0.0123000371).  The three
!  bodies are integrated as an N-body problem (Yoshida 4th order, 1 h
!  steps), so the Moon's solar perturbations (evection, variation, annual
!  equation, and the parallactic inequality, which depends on a_M/a_E)
!  come from the dynamics.  The Earth–Sun distance a_E is therefore
!  measurable in principle through the parallactic inequality (amplitude
!  about 2') and the Sun's diurnal parallax (9").
!  Three nuisance parameters absorb instrument zero points: an altitude
!  offset (quadrants) and a declination offset (armillae) per year.
!
!  Each reading is predicted topocentrically from Uraniborg: annual
!  aberration (first order), rotation to the true equator of date,
!  diurnal parallax, limb offsets, refraction (Saemundsson, 10 C,
!  1010 mbar — the same model as compare_jpl.py).
!
!  Fits: batch least squares (Levenberg–Marquardt, robust per-group
!  weights) and the iterated extended Kalman filter (sequential
!  measurement updates; with constant parameters, F = I and Q = 0, it is
!  the same estimator, which the program checks).  A Monte-Carlo run on
!  synthetic data (JPL + noise with Tycho's scatter, same moments and
!  kinds of readings) shows how well each parameter can be determined.
!
!  Mode "blind" (tycho_obs_blind.dat) uses nothing from JPL: every reading,
!  times from the star clock or solved with the model's own Sun (local
!  apparent time, each daytime session's clock) and Moon (transits), signs
!  and limbs chosen by the model, no planets, the Earth-Moon barycentre's
!  elements, and orbits built up from ellipses fitted to the readings
!  (see blind_fit below).  JPL's state is read only for the comparison.
!  Optionally the Moon's radius x(18), the strength of refraction x(19)
!  and Uraniborg's latitude x(20) are fitted too, per year if asked
!  (x(21:24) for 1595), and the final fit can be restarted from an earlier
!  solution on a subset of the readings (arguments 9, 11, 12).
!  With "tycho" in argument 11 (readings from export_tycho_only.py) the
!  model uses nothing but Tycho's values: his frame and stars, no
!  aberration or J2, a massless Moon, a spherical Earth, and Cassini's
!  homogeneous atmosphere for refraction; "huber", "tukey" or "noclip"
!  choose the weighting of the final fit.
!
!  Build and run:
!    gfortran -O3 -ffree-line-length-none -J build/tycho -o build/tycho/tycho_fit \
!        src/constants_mod.f90 app/tycho_fit.f90
!    build/tycho/tycho_fit tycho_obs_blind.dat jpl_state.dat blind
! ═══════════════════════════════════════════════════════════════════════
module tycho_fit_mod
  use constants_mod
  implicit none

  ! NMON 30-day windows from JD 2300000.5 (1585-1595): a declination zero
  ! (x(96:220)), an altitude zero (x(221:345)), a zero of the Moon's
  ! distances (x(346:470)) and an excess of the distances read to the far
  ! limb (x(471:595)) of each window
  integer, parameter :: NMON = 125
  ! reference stars of the distances: an offset (east, north; arcmin) of
  ! each from the catalogue, x(596:...), "rstar" (random effects)
  integer, parameter :: NSTAR = 80, IST0 = 95 + 4*NMON
  integer, parameter :: NP_ORB = 14, NP_NUIS = 81 + 4*NMON + 2*NSTAR, NP = NP_ORB + NP_NUIS
  integer, parameter :: MAXOBS = 20000
  real(dp), parameter :: LAT0_DEG = 55.0_dp + 54.0_dp/60 + 27.0_dp/3600
  ! latitude used by the predictions: LAT0_DEG + x(20) (arcmin); refraction
  ! scaled by x(19); the Moon's radius MOON_RADIUS*x(18).  Set from x at the
  ! start of predict_all.  x(18:20) are fitted only where param_on allows.
  real(dp) :: lat_deg = LAT0_DEG, ref_scale = 1.0_dp
  ! x(26:29): 1594's altitude offset, Moon radius, refraction, latitude
  ! (peryear); x(30): 1594's declination offset
  ! x(43:87): one block of 9 for each further year (1587, 1588, 1590, 1591,
  ! 1592): altitude offset, declination offset, Moon radius, refraction,
  ! latitude (peryear), armillary tilts (armil), quadrant tilt and azimuth
  ! zero (ewinst); see ip_year
  integer :: ii_
  logical :: param_on(NP) = [(.true., ii_ = 1, 17), (.false., ii_ = 18, 29), .true., (.false., ii_ = 31, NP)]
  ! x(37:39): tilt of the altitude quadrants, an altitude error growing with
  ! sin(azimuth) (east minus west), arcmin; x(40:42): zero of the azimuths,
  ! arcmin; each for 1586, 1594, 1595 ("ewinst")
  logical :: use_ew = .false.
  ! altitude of the body read (the lower of it and the Sun for readings
  ! against the Sun), deg, under the last predicted orbits
  real(dp) :: alt_read(MAXOBS) = 90.0_dp
  ! hour angle of the body read (deg, west +), under the last predicted orbits
  real(dp) :: ha_read(MAXOBS) = 0.0_dp
  ! x(31:36): tilt of the armillary sphere's polar axis (arcmin) along the
  ! meridian and east-west, for 1586, 1594, 1595 ("armil"); applied to every
  ! reading on the armillae (declinations, hour angles, RA differences and
  ! the clock stars)
  logical :: use_armil = .false.
  ! "phase": a limb reading of the Moon is the extreme point of its LIT part
  ! in the direction read (a horn of the crescent, the edge of the
  ! terminator of a gibbous Moon), from the model's own Sun
  logical :: use_phase = .false.
  ! "irr" (x(88), arcmin): irradiation, the Moon's apparent semidiameter
  ! grows by x(88) times its lit fraction (1 - cos elongation)/2;
  ! "kew" (x(89)): the east and west limbs (and distances' near/far lit
  ! limbs) with their own radius scale instead of k_moon
  logical :: use_irr = .false., use_kew = .false.
  ! "dfar" (x(90)), "dnear" (x(91)), arcmin: offsets of the Moon-star and
  ! Moon-Sun distances read to the Moon's limb on the far side from the
  ! other body (across the disk) and on the near side
  logical :: use_dfar = .false., use_dnear = .false.
  ! "limbew" (x(92:93)) / "limbo" (x(92:95)), arcmin: each named limb of the
  ! Moon (W, E, U, L) read outward by its own amount, in every kind of
  ! reading (the bright limb read differently on the two sides)
  logical :: use_limbo = .false.
  ! "decns": the upper/lower limb of a declination reading is the north/south one
  logical :: dec_ns = .false.
  ! "mdec" / "malt": the zero of the declinations / of the altitudes in each
  ! 30-day window (Sun and Moon together; the year's own offset is then off)
  logical :: use_mdec = .false., use_malt = .false., use_mdist = .false., use_mfar = .false.
  logical :: auto_sel = .false.
  ! mode "alien": the procedure of data/fit/ALIEN_PROCEDURE.md, its choices
  ! fixed here (arguments 7 and 9-12 are not read)
  logical :: alien = .false.
  character(len=*), parameter :: ALIEN_MODEL = 'tycho+aberr,moonr,refr,lat,phase,nosunpar,tgroup,rdec,ralt,rdist,rfar'
  character(len=*), parameter :: ALIEN_SUBSET = 'moonstar'
  real(dp), parameter :: ALIEN_GAP = 5.0_dp
  ! random effects ("rdec", "ralt", "rdist", "rfar"): the campaigns' values
  ! of a correction scatter about its mean (the year's offset; for the far
  ! limb one excess, x(90)) with a spread tau (arcmin) found from the data;
  ! kinds 1 declination, 2 altitude, 3 distance, 4 far limb
  logical :: reff(5) = .false.
  real(dp) :: tau_re(5) = 1.0_dp
  ! derivatives of a reading with respect to its reference star's east and
  ! north offsets (arcmin per arcmin), set in predict_pt
  real(dp) :: star_dp(2, MAXOBS) = 0.0_dp
  logical :: use_rstar = .false.
  ! reading being predicted, and whether a distance was read to the Moon's
  ! far limb (set in predict_pt)
  integer :: cur_i = 0
  ! length of the windows, days (environment TYCHO_WIN, default 30)
  real(dp) :: win_days = 30.0_dp
  ! campaigns (environment TYCHO_GAP, days): a new window wherever the Moon
  ! readings pause for more than camp_gap days; 0 = fixed windows
  real(dp) :: camp_gap = 0.0_dp
  logical :: farlimb(MAXOBS) = .false.
  ! "nosunpar": the Sun seen from the Earth's centre (its 9" diurnal
  ! parallax, which also carries the Earth-Sun distance, left out), so
  ! that only the Moon's motion measures the distance
  logical :: no_sun_par = .false.
  ! "tycho" model (argument 11): nothing but Tycho's values and physics --
  ! no aberration, no J2, a massless Moon, a spherical Earth (distances in
  ! Earth radii; km only as the unit 6378.137), Cassini's homogeneous
  ! atmosphere for refraction (refractive excess 2.9e-4*x(19), height
  ! 1.3e-3*x(25) Earth radii, both fitted), the Sun's radius 16' at the
  ! mean distance a_E (only two Sun readings are of a limb)
  logical :: tycho_model = .false., use_aberr = .true.
  ! diagnostics: put one modern piece back into the tycho model
  logical :: use_saem = .false., use_ell = .false.
  ! DIAGNOSTIC 'sunfix': the Sun's apparent radius 16' at the mean distance (not 696000 km)
  logical :: sun_sd_fixed = .false.
  real(dp) :: moon_ratio = 0.0123000371_dp, h_scale = 1.0_dp
  real(dp), parameter :: SUN_SD_RAD = 16.0_dp/60.0_dp*3.14159265358979324_dp/180.0_dp
  ! weighting in the final fit: 0 reject beyond clip*sigma (default),
  ! 1 Huber (k = 1.5), 2 Tukey biweight (c = 4.685), 3 none (plain least squares)
  integer :: robust_kind = 0
  ! peryear: the 1595 readings have their own altitude offset, Moon radius,
  ! refraction and latitude, x(21:24) (otherwise x(15), x(18:20) for both)
  logical :: peryear = .false.
  ! elongation of the Moon from the Sun (deg, 0-360, eastward, in right
  ! ascension) at each reading, under the last predicted orbits
  real(dp) :: elong(MAXOBS) = 0.0_dp
  real(dp), parameter :: ELEV_KM = 0.030_dp
  real(dp), parameter :: TEMP_C = 10.0_dp, PRES_MBAR = 1010.0_dp
  real(dp), parameter :: MOON_EARTH_RATIO = 0.0123000371_dp
  real(dp), parameter :: SUN_RADIUS = 696000.0_dp, MOON_RADIUS = 1737.4_dp
  real(dp), parameter :: STEP_S = 3600.0_dp
  ! bodies: 1 Sun, 2 Earth, 3 Moon, 4 Mercury, 5 Venus, 6 Mars, 7 Jupiter,
  ! 8 Saturn (planet barycentres; perturbers only, started from DE441)
  integer, parameter :: NBODY = 8
  real(dp), parameter :: GM_PLANETS(5) = [22031.868551_dp, 324858.592_dp, 42828.375816_dp, &
                                          126712764.1_dp, 37940584.8418_dp]
  real(dp) :: sun_state(6), planet_state(6,5)
  ! the planets move about the Sun with DE441's GM_sun, whatever the fit
  ! says about the Earth's orbit (they are known perturbers, not fitted)
  real(dp), parameter :: GM_SUN_PLANETS = 132712440041.279419_dp
  logical :: use_planets = .true.
  ! Earth's oblateness acting on the Moon (and back on the Earth)
  real(dp), parameter :: J2_EARTH = 1.08263d-3, RE_KM = 6378.137_dp
  real(dp) :: earth_pole(3) = [0.0_dp, 0.0_dp, 1.0_dp]
  logical :: use_j2 = .true.

  character(len=6), parameter :: PNAME0(95) = [character(len=6) :: &
      'e_E', 'i_E', 'Om_E', 'w_E', 'M0_E', 'n_E', 'a_E', &
      'e_M', 'i_M', 'Om_M', 'w_M', 'M0_M', 'n_M', 'a_M', &
      'b_alt', 'b_d86', 'b_d95', 'k_moon', 'k_refr', 'd_lat', 'b_al95', 'k_mo95', 'k_re95', 'd_la95', 'h_atm', 'b_al94', 'k_mo94', 'k_re94', 'd_la94', 'b_d94', &
      'ar_m86', 'ar_e86', 'ar_m94', 'ar_e94', 'ar_m95', 'ar_e95', &
      'qt86', 'qt94', 'qt95', 'az86', 'az94', 'az95', &
      'b_al87', 'b_d87', 'k_mo87', 'k_re87', 'd_la87', 'ar_m87', 'ar_e87', 'qt87', 'az87', &
      'b_al88', 'b_d88', 'k_mo88', 'k_re88', 'd_la88', 'ar_m88', 'ar_e88', 'qt88', 'az88', &
      'b_al90', 'b_d90', 'k_mo90', 'k_re90', 'd_la90', 'ar_m90', 'ar_e90', 'qt90', 'az90', &
      'b_al91', 'b_d91', 'k_mo91', 'k_re91', 'd_la91', 'ar_m91', 'ar_e91', 'qt91', 'az91', &
      'b_al92', 'b_d92', 'k_mo92', 'k_re92', 'd_la92', 'ar_m92', 'ar_e92', 'qt92', 'az92', &
      'irr', 'k_ew', 'd_far', 'd_near', 'o_W', 'o_E', 'o_U', 'o_L']
  character(len=6) :: PNAME(NP)

  type :: reading
    character(len=12) :: id
    integer :: year, body          ! body 1 = Sun, 2 = Moon
    character(len=6) :: qty
    character(len=1) :: point, refpoint
    character(len=4) :: refkind
    real(dp) :: jd, last, M(3,3), sra, sdec, value, jpl, group_err
    integer :: group
    ! blind files (tycho_obs_blind.dat): how the moment is found
    !   S true moment (star clock) | A local apparent time lat_h, turned
    !   into a moment with the model's own Sun | M meridian transit near jd,
    !   found with the model's own Sun or Moon
    character(len=1) :: tkind = 'S', sgn = '+'
    real(dp) :: lat_h = 0.0_dp
    ! observing session of an A reading (one clock); the clock correction
    ! (minutes, + = true time later than written) found with the model's
    ! own Sun; clk = 'x' if the logbook says not to use it for the clock
    integer :: sess = 0
    character(len=1) :: clk = '-'
    real(dp) :: dclock = 0.0_dp
    ! window (observing campaign) of the reading, for the window zeros
    integer :: win = 1
    ! reference star (1..nstar_used), 0 if none
    integer :: star = 0
  end type

  type(reading) :: obs(MAXOBS), obs_all(MAXOBS)
  ! integrated states at each reading for the orbits of c_key (predict_all)
  integer, parameter :: KEYN = 24
  real(dp), allocatable :: c_pos(:,:,:), c_vel(:,:,:)
  real(dp) :: c_key(KEYN) = -1.0d300
  integer :: nobs = 0, nall = 0
  ! limb chosen for a reading whose limb is not named (point n or v); kept
  ! fixed while the Jacobian is formed
  character(len=1) :: lchoice(MAXOBS) = 'C'
  logical :: fix_choice = .false.
  real(dp) :: t_used(MAXOBS) = 0.0_dp ! moment each reading was predicted at (JD TT)
  integer :: nb_used = 8              ! bodies integrated (3 = Sun, Earth, Moon only)
  ! Elements 1-7 describe the Earth-Moon barycentre about the Sun (instead
  ! of the Earth); its orbit is a clean ellipse whatever the Earth-Sun
  ! distance, while the Earth's own osculating elements carry its monthly
  ! wobble about the barycentre
  logical :: emb_param = .false.
  integer :: order(MAXOBS)          ! obs sorted by time
  real(dp) :: epoch_jd
  real(dp) :: obs_xyz(3)            ! observer, geocentric (km), in the frame of date at LAST = 0
  real(dp) :: coslat_rho, sinlat_rho

  ! Groups (year x body x quantity) for robust weighting
  ! groups: year x body x quantity (NGB), times three timing classes of the
  ! Moon's readings when tgroup is set (star clock, local time, transit)
  integer, parameter :: NQ = 8, NY = 8, NGB = NY*2*NQ, NG = 3*NGB
  integer, parameter :: YEARS(NY) = [1586, 1594, 1595, 1587, 1588, 1590, 1591, 1592]
  character(len=6), parameter :: QNAMES(NQ) = [character(len=6) :: &
      'DEC', 'ALT', 'ALTM', 'HA', 'AZ', 'DIST', 'RADIFF', 'DIAM']
  real(dp) :: sigma_g(NG)
  logical :: freeze_sigma = .false.   ! keep sigma_g fixed during a fit (profiles)
  ! DIAGNOSTIC: group sigmas read from the file in TYCHO_SIGMAS and kept for
  ! the whole run (the same weights for simulated and real readings)
  logical :: sigma_fixed = .false.
  logical :: freeze_use = .false.     ! keep the set of used readings fixed (profiles)
  logical :: kepler_moon = .false.    ! Moon on an unperturbed Kepler ellipse about the Earth
  real(dp) :: last_chi2 = 0.0_dp

contains

  ! ── small helpers ─────────────────────────────────────────────────
  pure function wrap180(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = modulo(x + 180.0_dp, 360.0_dp) - 180.0_dp
  end function

  ! index into YEARS (1586, 1594, 1595, then the further years); anything
  ! else counts as 1595
  integer function yidx(year)
    integer, intent(in) :: year
    integer :: k
    yidx = 3
    do k = 1, NY
      if (YEARS(k) == year) yidx = k
    end do
  end function

  ! index of the further year's parameter j (0 altitude offset, 1
  ! declination offset, 2 Moon radius, 3 refraction, 4 latitude, 5-6
  ! armillary tilts, 7 quadrant tilt, 8 azimuth zero); 0 for 1586/94/95
  integer function ip_year(year, j)
    integer, intent(in) :: year, j
    ip_year = 0
    if (yidx(year) >= 4) ip_year = 43 + 9*(yidx(year) - 4) + j
  end function

  ! the reference stars: the same star has the same place in the frame of
  ! the readings file (the date's place rotated back with the reading's
  ! matrix M)
  subroutine set_stars()
    real(dp) :: v(3), u(3), list(3, NSTAR)
    integer :: i, k, ns
    ns = 0
    do i = 1, nobs
      obs(i)%star = 0
      if (trim(obs(i)%refkind) /= 'STAR') cycle
      u = [cos(obs(i)%sdec*DEG2RAD)*cos(obs(i)%sra*DEG2RAD), cos(obs(i)%sdec*DEG2RAD)*sin(obs(i)%sra*DEG2RAD), &
           sin(obs(i)%sdec*DEG2RAD)]
      v = matmul(transpose(obs(i)%M), u)
      do k = 1, ns
        if (dot_product(v, list(:, k)) > cos(5.0_dp/60.0_dp*DEG2RAD)) exit
      end do
      if (k > ns) then
        if (ns == NSTAR) stop 'set_stars: more than NSTAR reference stars'
        ns = ns + 1; list(:, ns) = v; k = ns
      end if
      obs(i)%star = k
    end do
    print '(A,I4,A)', '   reference stars of the distances:', ns, ' (an offset each with rstar)'
  end subroutine

  ! window of every reading: fixed win_days windows, or campaigns of Moon
  ! readings separated by pauses longer than camp_gap days (the Sun's
  ! readings join the campaign nearest in time)
  subroutine set_windows()
    real(dp) :: t(MAXOBS), edge(NMON + 1), tk
    integer :: i, n, ne, k
    if (camp_gap <= 0.0_dp) then
      do i = 1, nobs
        obs(i)%win = mon_of(obs(i)%jd)
      end do
      return
    end if
    n = 0
    do i = 1, nobs
      if (obs(i)%body == 2) then
        n = n + 1; t(n) = obs(i)%jd
      end if
    end do
    do i = 2, n                                ! insertion sort
      tk = t(i); k = i - 1
      do while (k >= 1)
        if (t(k) <= tk) exit
        t(k + 1) = t(k); k = k - 1
      end do
      t(k + 1) = tk
    end do
    ne = 0
    do i = 2, n
      if (t(i) - t(i - 1) > camp_gap) then
        ne = ne + 1
        if (ne > NMON - 1) stop 'set_windows: too many campaigns; raise NMON or TYCHO_GAP'
        edge(ne) = 0.5_dp*(t(i) + t(i - 1))
      end if
    end do
    do i = 1, nobs
      k = 1
      do while (k <= ne)
        if (obs(i)%jd < edge(k)) exit
        k = k + 1
      end do
      obs(i)%win = k
    end do
    print '(A,F6.1,A,I4,A)', '   campaigns: pauses of more than ', camp_gap, ' days give ', ne + 1, ' windows'
  end subroutine

  ! Forward selection of the per-campaign corrections by the Bayesian
  ! information criterion: a candidate with k parameters (campaigns that
  ! have such readings) is adopted when it lowers chi2 by more than
  ! k ln N, N the readings used; the best candidate first, until none is
  ! worth its parameters.  Weights and the set of readings are frozen
  ! during the tests, so the chi2 are comparable.
  subroutine select_corrections(xb, v, act, cv, u, rms)
    real(dp), intent(inout) :: xb(NP), cv(NP,NP), rms
    real(dp), intent(in) :: v(nobs)
    logical, intent(in) :: act(NP)
    logical, intent(inout) :: u(nobs)
    character(len=6), parameter :: CAND(4) = [character(len=6) :: 'mdec', 'mfar', 'malt', 'mdist']
    logical :: taken(4)
    real(dp) :: chi0, chi1, best, xt(NP), xbest(NP), sgf(NG), crit
    integer :: c, k, kbest, nused, i, kpar
    logical :: haswin(NMON)
    taken = .false.
    sgf = sigma_g
    freeze_sigma = .true.; freeze_use = .true.
    call fit_lm(xb, v, act, 4.0_dp, cv, u, rms, .false.)
    chi0 = last_chi2; nused = count(u)
    print '(/A,F12.2,A,I6,A)', '   choosing the corrections per campaign (chi2 ', chi0, ', ', nused, ' readings)'
    do
      best = 0.0_dp; kbest = 0
      do c = 1, 4
        if (taken(c)) cycle
        call switch_corr(CAND(c), .true.)
        ! parameters: campaigns with readings of that kind
        haswin = .false.
        do i = 1, nobs
          if (.not. u(i)) cycle
          select case (c)
          case (1); if (trim(obs(i)%qty) == 'DEC') haswin(obs(i)%win) = .true.
          case (2); if (trim(obs(i)%qty) == 'DIST' .and. obs(i)%body == 2 .and. farlimb(i)) haswin(obs(i)%win) = .true.
          case (3); if (obs(i)%qty(1:3) == 'ALT') haswin(obs(i)%win) = .true.
          case (4); if (trim(obs(i)%qty) == 'DIST' .and. obs(i)%body == 2) haswin(obs(i)%win) = .true.
          end select
        end do
        kpar = count(haswin)
        xt = xb
        call fit_lm(xt, v, act, 4.0_dp, cv, u, rms, .false.)
        chi1 = last_chi2
        crit = (chi0 - chi1) - kpar*log(real(nused, dp))
        print '(A,A6,A,I4,A,F12.2,A,F10.2,A,F14.1)', '      + ', CAND(c), ': ', kpar, ' parameters, chi2 ', chi1, &
            '   gain - penalty ', crit, '   a_E ', xt(7)
        if (crit > best) then
          best = crit; kbest = c; xbest = xt
        end if
        call switch_corr(CAND(c), .false.)
      end do
      if (kbest == 0) exit
      taken(kbest) = .true.
      call switch_corr(CAND(kbest), .true.)
      xb = xbest
      call fit_lm(xb, v, act, 4.0_dp, cv, u, rms, .false.)
      chi0 = last_chi2
      print '(A,A6,A,F12.2)', '      adopted ', CAND(kbest), ';  chi2 ', chi0
    end do
    freeze_sigma = .false.; freeze_use = .false.
    sigma_g = sgf
    call fit_lm(xb, v, act, 4.0_dp, cv, u, rms, .false.)
    print '(A,4(1X,A6,L2))', '   corrections adopted:', (CAND(c), taken(c), c = 1, 4)
  end subroutine

  ! command argument n, or the procedure's own choice in mode alien
  subroutine get_arg(n, arg)
    integer, intent(in) :: n
    character(len=*), intent(out) :: arg
    if (.not. alien) then
      call get_command_argument(n, arg); return
    end if
    select case (n)
    case (7); arg = 'blind'
    case (9); arg = 'scratch'
    case (10); arg = 'clock'
    case (11)
      arg = ALIEN_MODEL
      ! (a variant of the model for comparison, environment TYCHO_ALIEN_MODEL)
      block
        integer :: st
        character(len=256) :: mv
        call get_environment_variable('TYCHO_ALIEN_MODEL', mv, status=st)
        if (st == 0 .and. len_trim(mv) > 0) arg = mv
      end block
    case (12); arg = ALIEN_SUBSET
    case default; call get_command_argument(n, arg)
    end select
  end subroutine

  ! spread of the prior of parameter k (0: no prior)
  real(dp) function prior_tau(k)
    integer, intent(in) :: k
    integer :: kind
    prior_tau = 0.0_dp
    if (k <= 95) return
    kind = (k - 96)/NMON + 1
    if (k > IST0) kind = 5
    if (kind >= 1 .and. kind <= 5) then
      if (reff(kind)) prior_tau = tau_re(kind)
    end if
  end function

  ! the spreads of the random effects from the campaigns' values and their
  ! uncertainties (one EM step): tau^2 = mean(x^2 + var) over campaigns
  ! that have readings of the kind
  subroutine update_tau(x, cv)
    real(dp), intent(in) :: x(NP), cv(NP,NP)
    integer :: kind, w, k, n
    real(dp) :: s
    do kind = 1, 5
      if (.not. reff(kind)) cycle
      s = 0.0_dp; n = 0
      do w = 1, merge(2*NSTAR, NMON, kind == 5)
        k = 95 + (kind - 1)*NMON + w
        if (kind == 5) k = IST0 + w
        if (cv(k, k) <= 0.0_dp) cycle
        s = s + x(k)**2 + cv(k, k); n = n + 1
      end do
      if (n > 0) tau_re(kind) = max(sqrt(s/n), 0.05_dp)
    end do
  end subroutine

  ! switch one per-campaign correction on or off (and with it the year's
  ! own offset of the same kind, the other way)
  subroutine switch_corr(name, on)
    character(len=*), intent(in) :: name
    logical, intent(in) :: on
    integer :: k1
    select case (trim(name))
    case ('mdec')
      use_mdec = on; param_on(96:95 + NMON) = on
      param_on([16, 17, 30]) = .not. on
      do k1 = 4, NY
        param_on(ip_year(YEARS(k1), 1)) = .not. on
      end do
    case ('malt')
      use_malt = on; param_on(96 + NMON:95 + 2*NMON) = on
      param_on(15) = .not. on
      do k1 = 4, NY
        if (peryear) param_on(ip_year(YEARS(k1), 0)) = .not. on
      end do
    case ('mdist')
      use_mdist = on; param_on(96 + 2*NMON:95 + 3*NMON) = on
    case ('mfar')
      use_mfar = on; param_on(96 + 3*NMON:IST0) = on
    end select
  end subroutine

  ! "tgroup": the Moon's readings weighted in three timing classes, each
  ! with its own scatter: timed by stars (S), by local time and the Sun (A),
  ! meridian transits (M)
  subroutine set_time_groups()
    integer :: i, c
    do i = 1, nall
      if (obs_all(i)%body /= 2 .or. obs_all(i)%group > NGB) cycle
      c = index('SAM', obs_all(i)%tkind) - 1
      obs_all(i)%group = obs_all(i)%group + NGB*max(c, 0)
    end do
    do i = 1, nobs
      if (obs(i)%body /= 2 .or. obs(i)%group > NGB) cycle
      c = index('SAM', obs(i)%tkind) - 1
      obs(i)%group = obs(i)%group + NGB*max(c, 0)
    end do
  end subroutine

  ! 30-day window of a moment (1..NMON)
  integer function mon_of(jd)
    real(dp), intent(in) :: jd
    mon_of = max(1, min(NMON, int((jd - 2300000.5_dp)/win_days) + 1))
  end function

  subroutine init_names()
    integer :: k
    PNAME(1:95) = PNAME0
    do k = 1, NMON
      write(PNAME(95 + k), '(A,I3.3)') 'dd', k
      write(PNAME(95 + NMON + k), '(A,I3.3)') 'da', k
      write(PNAME(95 + 2*NMON + k), '(A,I3.3)') 'di', k
      write(PNAME(95 + 3*NMON + k), '(A,I3.3)') 'df', k
    end do
    do k = 1, NSTAR
      write(PNAME(IST0 + 2*k - 1), '(A,I2.2)') 'stE', k
      write(PNAME(IST0 + 2*k), '(A,I2.2)') 'stN', k
    end do
  end subroutine

  ! declination offset of a year
  integer function ip_dec(year)
    integer, intent(in) :: year
    ip_dec = 17
    if (year == 1586) ip_dec = 16
    if (year == 1594) ip_dec = 30
    if (ip_year(year, 1) > 0) ip_dec = ip_year(year, 1)
  end function

  ! quadrant tilt (j = 0) and azimuth zero (j = 1) of a year
  integer function ip_ew(year, j)
    integer, intent(in) :: year, j
    if (yidx(year) >= 4) then
      ip_ew = ip_year(year, 7 + j)
    else
      ip_ew = 37 + 3*j + yidx(year) - 1
    end if
  end function

  integer function qindex(q)
    character(len=*), intent(in) :: q
    integer :: k
    qindex = 0
    do k = 1, NQ
      if (trim(q) == trim(QNAMES(k))) qindex = k
    end do
  end function

  subroutine eq_to_hor(H, d, a, A_az)
    real(dp), intent(in) :: H, d
    real(dp), intent(out) :: a, A_az
    real(dp) :: hr, dr, p, sa
    hr = H*DEG2RAD; dr = d*DEG2RAD; p = lat_deg*DEG2RAD
    sa = sin(p)*sin(dr) + cos(p)*cos(dr)*cos(hr)
    a = asin(max(-1.0_dp, min(1.0_dp, sa)))*RAD2DEG
    A_az = modulo(atan2(-cos(dr)*sin(hr), sin(dr)*cos(p) - cos(dr)*cos(hr)*sin(p))*RAD2DEG, 360.0_dp)
  end subroutine

  subroutine hor_to_eq(a, A_az, H, d)
    real(dp), intent(in) :: a, A_az
    real(dp), intent(out) :: H, d
    real(dp) :: alr, azr, p, sd
    alr = a*DEG2RAD; azr = A_az*DEG2RAD; p = lat_deg*DEG2RAD
    sd = sin(p)*sin(alr) + cos(p)*cos(alr)*cos(azr)
    d = asin(max(-1.0_dp, min(1.0_dp, sd)))*RAD2DEG
    H = atan2(-sin(azr)*cos(alr), cos(p)*sin(alr) - sin(p)*cos(alr)*cos(azr))*RAD2DEG
  end subroutine

  function refraction_deg(alt) result(r)
    real(dp), intent(in) :: alt
    real(dp) :: r
    if (alt < -1.5_dp) then
      r = 0.0_dp; return
    end if
    if (tycho_model .and. .not. use_saem) then
      r = cassini_deg(alt); return
    end if
    r = 1.02_dp / tan((alt + 10.3_dp/(alt + 5.11_dp))*DEG2RAD)
    r = r * (PRES_MBAR/1010.0_dp) * (283.0_dp/(273.0_dp + TEMP_C))
    r = ref_scale*max(r, 0.0_dp) / 60.0_dp
  end function

  ! Cassini: a homogeneous shell of height h (Earth radii), index 1 + n1; the
  ! ray at apparent zenith distance z crosses it straight and bends once at
  ! the top.  Input the true altitude; the apparent one is found by iteration.
  function cassini_deg(alt) result(r)
    real(dp), intent(in) :: alt
    real(dp) :: r, n1, h, zt, za, si, so
    integer :: it
    n1 = 2.9d-4*ref_scale; h = 1.3d-3*h_scale
    zt = (90.0_dp - alt)*DEG2RAD; za = zt
    do it = 1, 8
      si = sin(za)/(1.0_dp + h); so = min(1.0_dp, (1.0_dp + n1)*si)
      za = zt - (asin(so) - asin(si))
    end do
    r = (zt - za)*RAD2DEG
  end function

  subroutine refract(H, d, H2, d2)
    real(dp), intent(in) :: H, d
    real(dp), intent(out) :: H2, d2
    real(dp) :: a, A_az
    call eq_to_hor(H, d, a, A_az)
    call hor_to_eq(a + refraction_deg(a), A_az, H2, d2)
  end subroutine

  subroutine limb_point(H, d, sd, pt, H2, d2)
    real(dp), intent(in) :: H, d, sd
    character(len=1), intent(in) :: pt
    real(dp), intent(out) :: H2, d2
    real(dp) :: a, A_az
    H2 = H; d2 = d
    select case (pt)
    case ('U', 'L')
      call eq_to_hor(H, d, a, A_az)
      if (pt == 'U') then
        call hor_to_eq(a + sd, A_az, H2, d2)
      else
        call hor_to_eq(a - sd, A_az, H2, d2)
      end if
    case ('E')
      H2 = H - sd/max(cos(d*DEG2RAD), 1.0d-6)
    case ('W')
      H2 = H + sd/max(cos(d*DEG2RAD), 1.0d-6)
    case ('A')                                ! the limb toward the north pole
      d2 = d + sd
    case ('B')                                ! toward the south pole
      d2 = d - sd
    end select
  end subroutine

  ! Hour angle and declination as read on an armillary sphere whose polar
  ! axis is tilted by am (along the meridian) and ae (east-west), arcmin
  subroutine armil(H, d, x, year, H2, d2)
    real(dp), intent(in) :: H, d, x(NP)
    integer, intent(in) :: year
    real(dp), intent(out) :: H2, d2
    real(dp) :: u(3), v(3), am, ae
    integer :: k
    if (.not. use_armil) then
      H2 = H; d2 = d; return
    end if
    if (yidx(year) >= 4) then
      k = ip_year(year, 5) - 31
    else
      k = 2*(yidx(year) - 1)
    end if
    am = x(31 + k)/60.0_dp*DEG2RAD; ae = x(32 + k)/60.0_dp*DEG2RAD
    u = [cos(d*DEG2RAD)*cos(H*DEG2RAD), cos(d*DEG2RAD)*sin(H*DEG2RAD), sin(d*DEG2RAD)]
    v = [u(1)*cos(am) - u(3)*sin(am), u(2), u(1)*sin(am) + u(3)*cos(am)]
    u = [v(1), v(2)*cos(ae) + v(3)*sin(ae), -v(2)*sin(ae) + v(3)*cos(ae)]
    H2 = atan2(u(2), u(1))*RAD2DEG
    d2 = asin(max(-1.0_dp, min(1.0_dp, u(3))))*RAD2DEG
  end subroutine

  ! Unit direction (east, north) on the sky from (H1, d1) toward (H2, d2)
  subroutine sky_dir(H1, d1, H2, d2, ue, un)
    real(dp), intent(in) :: H1, d1, H2, d2
    real(dp), intent(out) :: ue, un
    real(dp) :: da, n
    da = -(H2 - H1)*DEG2RAD                ! right ascension grows eastward, hour angle westward
    ue = sin(da)*cos(d2*DEG2RAD)
    un = cos(d1*DEG2RAD)*sin(d2*DEG2RAD) - sin(d1*DEG2RAD)*cos(d2*DEG2RAD)*cos(da)
    n = max(sqrt(ue*ue + un*un), 1.0d-12); ue = ue/n; un = un/n
  end subroutine

  ! Apparent semidiameter of the lit Moon in the direction of limb point pt
  ! (U, L in altitude; E, W; N toward, F away from (Hx, dx)); sd for the full
  ! disk.  Sun at (Hs, ds).  On the dark side the lit figure reaches
  ! sd*sqrt(cos^2 phi + max(cos i, 0)^2 sin^2 phi), phi from the line of the
  ! horns, i the phase angle (about 180 deg - elongation).
  function lit_sd(pt, H, d, sd, Hs, ds, Hx, dx) result(sdl)
    character(len=1), intent(in) :: pt
    real(dp), intent(in) :: H, d, sd, Hs, ds, Hx, dx
    real(dp) :: sdl, se, sn, ue, un, a, A_az, H2, d2, ct, cp, ci, n
    sdl = sd
    if (.not. use_phase) return
    select case (pt)
    case ('U', 'L')
      call eq_to_hor(H, d, a, A_az)
      call hor_to_eq(a + 0.05_dp, A_az, H2, d2)
      ue = -wrap180(H2 - H)*cos(d*DEG2RAD); un = d2 - d
      n = max(sqrt(ue*ue + un*un), 1.0d-12); ue = ue/n; un = un/n
      if (pt == 'L') then
        ue = -ue; un = -un
      end if
    case ('E'); ue = 1.0_dp; un = 0.0_dp
    case ('W'); ue = -1.0_dp; un = 0.0_dp
    case ('A'); ue = 0.0_dp; un = 1.0_dp
    case ('B'); ue = 0.0_dp; un = -1.0_dp
    case ('N', 'F')
      call sky_dir(H, d, Hx, dx, ue, un)
      if (pt == 'F') then
        ue = -ue; un = -un
      end if
    case default
      return
    end select
    call sky_dir(H, d, Hs, ds, se, sn)
    ct = ue*se + un*sn
    if (ct >= 0.0_dp) return                 ! the lit side: the limb itself
    cp = abs(-ue*sn + un*se)                 ! cos of the angle to the line of the horns
    ci = max(-cos(sep(H, d, Hs, ds)*DEG2RAD), 0.0_dp)   ! cos(phase angle) on a gibbous Moon
    sdl = sd*sqrt(cp*cp + ci*ci*ct*ct)
  end function

  function sep(H1, d1, H2, d2) result(s)
    real(dp), intent(in) :: H1, d1, H2, d2
    real(dp) :: s, c
    c = sin(d1*DEG2RAD)*sin(d2*DEG2RAD) + cos(d1*DEG2RAD)*cos(d2*DEG2RAD)*cos((H1-H2)*DEG2RAD)
    s = acos(max(-1.0_dp, min(1.0_dp, c)))*RAD2DEG
  end function

  function side(H1, d1, H2, d2, pt) result(s)
    ! positive when (H2, d2) lies on the named side of (H1, d1)
    real(dp), intent(in) :: H1, d1, H2, d2
    character(len=1), intent(in) :: pt
    real(dp) :: s, a1, a2, A_az, dx
    s = 0.0_dp
    select case (pt)
    case ('E', 'W')
      dx = wrap180(H2 - H1)*cos(d1*DEG2RAD)
      s = dx
      if (pt == 'E') s = -dx
    case ('U', 'L')
      call eq_to_hor(H1, d1, a1, A_az)
      call eq_to_hor(H2, d2, a2, A_az)
      s = a2 - a1
      if (pt == 'L') s = a1 - a2
    end select
  end function

  ! ── Kepler ──────────────────────────────────────────────────────────
  function solve_kepler(M_in, ecc) result(E)
    real(dp), intent(in) :: M_in, ecc
    real(dp) :: E, M, dE
    integer :: it
    M = modulo(M_in, TAU); E = M
    if (ecc > 0.8_dp) E = PI
    do it = 1, 50
      dE = (E - ecc*sin(E) - M) / (1.0_dp - ecc*cos(E))
      E = E - dE
      if (abs(dE) < 1.0d-15) exit
    end do
  end function

  subroutine kepler_to_cart(a, ecc, inc, raan, argp, M, mu, r, v)
    real(dp), intent(in) :: a, ecc, inc, raan, argp, M, mu
    real(dp), intent(out) :: r(3), v(3)
    real(dp) :: E, nu, rm, hm, rp(2), vp(2), cO, sO, ci, si, cw, sw
    E = solve_kepler(M, ecc)
    nu = atan2(sqrt(1.0_dp - ecc**2)*sin(E), cos(E) - ecc)
    rm = a*(1.0_dp - ecc*cos(E))
    hm = sqrt(mu*a*(1.0_dp - ecc**2))
    rp = [rm*cos(nu), rm*sin(nu)]
    vp = [-(mu/hm)*sin(nu), (mu/hm)*(ecc + cos(nu))]
    cO = cos(raan); sO = sin(raan); ci = cos(inc); si = sin(inc); cw = cos(argp); sw = sin(argp)
    r(1) = (cO*cw - sO*sw*ci)*rp(1) + (-cO*sw - sO*cw*ci)*rp(2)
    r(2) = (sO*cw + cO*sw*ci)*rp(1) + (-sO*sw + cO*cw*ci)*rp(2)
    r(3) = (sw*si)*rp(1) + (cw*si)*rp(2)
    v(1) = (cO*cw - sO*sw*ci)*vp(1) + (-cO*sw - sO*cw*ci)*vp(2)
    v(2) = (sO*cw + cO*sw*ci)*vp(1) + (-sO*sw + cO*cw*ci)*vp(2)
    v(3) = (sw*si)*vp(1) + (cw*si)*vp(2)
  end subroutine

  subroutine cart_to_kepler(r, v, mu, a, ecc, inc, raan, argp, M0)
    real(dp), intent(in) :: r(3), v(3), mu
    real(dp), intent(out) :: a, ecc, inc, raan, argp, M0
    real(dp) :: h(3), n(3), ev(3), rm, hm, nm, nu, E, cnu
    rm = norm2(r)
    h = [r(2)*v(3) - r(3)*v(2), r(3)*v(1) - r(1)*v(3), r(1)*v(2) - r(2)*v(1)]
    hm = norm2(h)
    n = [-h(2), h(1), 0.0_dp]; nm = norm2(n)
    ev(1) = (v(2)*h(3) - v(3)*h(2))/mu - r(1)/rm
    ev(2) = (v(3)*h(1) - v(1)*h(3))/mu - r(2)/rm
    ev(3) = (v(1)*h(2) - v(2)*h(1))/mu - r(3)/rm
    ecc = norm2(ev)
    a = -mu/(2.0_dp*(0.5_dp*dot_product(v, v) - mu/rm))
    inc = acos(max(-1.0_dp, min(1.0_dp, h(3)/hm)))
    raan = modulo(atan2(n(2), n(1)), TAU)
    cnu = max(-1.0_dp, min(1.0_dp, dot_product(n, ev)/(nm*ecc)))
    argp = acos(cnu); if (ev(3) < 0.0_dp) argp = TAU - argp
    cnu = max(-1.0_dp, min(1.0_dp, dot_product(ev, r)/(ecc*rm)))
    nu = acos(cnu); if (dot_product(r, v) < 0.0_dp) nu = TAU - nu
    E = atan2(sqrt(1.0_dp - ecc**2)*sin(nu), ecc + cos(nu))
    M0 = modulo(E - ecc*sin(E), TAU)
  end subroutine

  ! ── N-body (Sun, Earth, Moon) ───────────────────────────────────────
  subroutine accel(pos, gm, acc)
    real(dp), intent(in) :: pos(3,NBODY), gm(NBODY)
    real(dp), intent(out) :: acc(3,NBODY)
    real(dp) :: rij(3), r2, f
    integer :: i, j
    acc = 0.0_dp
    do i = 1, nb_used
      do j = i+1, nb_used
        rij = pos(:,j) - pos(:,i)
        r2 = dot_product(rij, rij)
        f = 1.0_dp/(r2*sqrt(r2))
        acc(:,i) = acc(:,i) + gm(j)*rij*f
        if (i == 1 .and. j > 3) then
          acc(:,j) = acc(:,j) - GM_SUN_PLANETS*rij*f
        else
          acc(:,j) = acc(:,j) - gm(i)*rij*f
        end if
      end do
    end do
    if (use_j2) then
      ! J2 of the Earth on the Moon; the reaction (mass ratio) on the Earth
      rij = pos(:,3) - pos(:,2); r2 = dot_product(rij, rij)
      block
        real(dp) :: z, r5, a(3), rr
        rr = sqrt(r2); z = dot_product(rij, earth_pole)
        r5 = r2*r2*rr
        a = -1.5_dp*J2_EARTH*gm(2)*RE_KM**2/r5*((1.0_dp - 5.0_dp*z*z/r2)*rij + 2.0_dp*z*earth_pole)
        acc(:,3) = acc(:,3) + a
        acc(:,2) = acc(:,2) - a*gm(3)/gm(2)
      end block
    end if
    if (kepler_moon) then
      ! the Sun pulls the Moon exactly as it pulls the Earth: no solar
      ! perturbation of the geocentric orbit (and none of the Moon on the Sun)
      rij = pos(:,1) - pos(:,2); r2 = dot_product(rij, rij); f = 1.0_dp/(r2*sqrt(r2))
      acc(:,3) = acc(:,3) - gm(1)*(pos(:,1) - pos(:,3))/norm2(pos(:,1) - pos(:,3))**3 + gm(1)*rij*f
    end if
  end subroutine

  subroutine step4(pos, vel, gm, dt)
    real(dp), intent(inout) :: pos(3,NBODY), vel(3,NBODY)
    real(dp), intent(in) :: gm(NBODY), dt
    real(dp), parameter :: CBRT2 = 2.0_dp**(1.0_dp/3.0_dp)
    real(dp), parameter :: W1 = 1.0_dp/(2.0_dp - CBRT2), W0 = -CBRT2/(2.0_dp - CBRT2)
    real(dp), parameter :: C1 = 0.5_dp*W1, C2 = 0.5_dp*(W0 + W1)
    real(dp) :: acc(3,NBODY)
    pos = pos + C1*dt*vel; call accel(pos, gm, acc); vel = vel + W1*dt*acc
    pos = pos + C2*dt*vel; call accel(pos, gm, acc); vel = vel + W0*dt*acc
    pos = pos + C2*dt*vel; call accel(pos, gm, acc); vel = vel + W1*dt*acc
    pos = pos + C1*dt*vel
  end subroutine

  subroutine propagate(pos, vel, gm, dt_total)
    real(dp), intent(inout) :: pos(3,NBODY), vel(3,NBODY)
    real(dp), intent(in) :: gm(NBODY), dt_total
    integer :: n, k
    real(dp) :: dt
    if (dt_total == 0.0_dp) return
    n = max(1, ceiling(abs(dt_total)/STEP_S))
    dt = dt_total/n
    do k = 1, n
      call step4(pos, vel, gm, dt)
    end do
  end subroutine

  subroutine initial_state(x, pos, vel, gm)
    ! columns: 1 Sun, 2 Earth, 3 Moon, 4.. planets (barycentric, km, km/s).
    ! The Sun and the planets start from DE441 (the Sun's barycentric motion
    ! is due to the planets); the Earth and the Moon from the fitted elements.
    real(dp), intent(in) :: x(NP)
    real(dp), intent(out) :: pos(3,NBODY), vel(3,NBODY), gm(NBODY)
    real(dp) :: mu_se, mu_em, gm_e, gm_m, gm_s, rE(3), vE(3), rM(3), vM(3)
    integer :: k
    mu_se = (x(6)/DAY_S)**2 * x(7)**3
    mu_em = (x(13)/DAY_S)**2 * x(14)**3
    gm_m = mu_em*moon_ratio/(1.0_dp + moon_ratio)
    gm_e = mu_em - gm_m
    gm_s = mu_se - gm_e
    if (emb_param) gm_s = mu_se - mu_em        ! mu_se = GM_sun + GM_earth + GM_moon
    gm(1:3) = [gm_s, gm_e, gm_m]
    call kepler_to_cart(x(7), x(1), x(2), x(3), x(4), x(5), mu_se, rE, vE)
    call kepler_to_cart(x(14), x(8), x(9), x(10), x(11), x(12), mu_em, rM, vM)
    pos(:,1) = sun_state(1:3);    vel(:,1) = sun_state(4:6)
    if (emb_param) then
      ! rE, vE are the Earth-Moon barycentre's
      rE = rE - (gm_m/mu_em)*rM; vE = vE - (gm_m/mu_em)*vM
    end if
    pos(:,2) = pos(:,1) + rE;     vel(:,2) = vel(:,1) + vE
    pos(:,3) = pos(:,2) + rM;     vel(:,3) = vel(:,2) + vM
    do k = 1, 5
      pos(:,3+k) = planet_state(1:3,k); vel(:,3+k) = planet_state(4:6,k)
      gm(3+k) = merge(GM_PLANETS(k), 0.0_dp, use_planets)
    end do
  end subroutine

  ! ── Topocentric apparent (H, dec, semidiameter) of a body ───────────
  subroutine topo(rb_icrs, v_earth, o, H, d, sd, radius, geo)
    real(dp), intent(in) :: rb_icrs(3), v_earth(3), radius
    logical, intent(in), optional :: geo   ! no diurnal parallax
    type(reading), intent(in) :: o
    real(dp), intent(out) :: H, d, sd
    real(dp) :: r(3), u(3), dn, rd(3), ob(3), lst, ra
    r = rb_icrs
    dn = norm2(r); u = r/dn
    if (use_aberr) r = r + dn*(v_earth/C_LIGHT_KMS - u*dot_product(u, v_earth/C_LIGHT_KMS))   ! annual aberration
    rd = matmul(o%M, r)
    lst = o%last*DEG2RAD
    ob = [coslat_rho*cos(lst), coslat_rho*sin(lst), sinlat_rho]
    if (present(geo)) then
      if (geo) ob = 0.0_dp
    end if
    rd = rd - ob
    dn = norm2(rd)
    ra = atan2(rd(2), rd(1))*RAD2DEG
    d = asin(rd(3)/dn)*RAD2DEG
    H = wrap180(o%last - ra)
    sd = asin(radius/dn)*RAD2DEG
  end subroutine

  ! ── Prediction of one reading (deg) ─────────────────────────────────
  ! A limb not named in the logbook (point n: nearer or farther, v: upper
  ! or lower) is the one closer to the reading under this model; the
  ! choice is recorded in lchoice(i) and kept while fix_choice is set.
  function predict_one(o, x, rsun, rmoon, vearth, i) result(p)
    type(reading), intent(in) :: o
    real(dp), intent(in) :: x(NP), rsun(3), rmoon(3), vearth(3)
    integer, intent(in) :: i
    real(dp) :: p, p1, p2
    type(reading) :: o2
    character(len=2) :: cand
    if (o%point /= 'n' .and. o%point /= 'v') then
      p = predict_pt(o, x, rsun, rmoon, vearth); return
    end if
    cand = merge('NF', 'UL', o%point == 'n')
    o2 = o
    if (fix_choice) then
      o2%point = lchoice(i)
      p = predict_pt(o2, x, rsun, rmoon, vearth); return
    end if
    o2%point = cand(1:1); p1 = predict_pt(o2, x, rsun, rmoon, vearth)
    o2%point = cand(2:2); p2 = predict_pt(o2, x, rsun, rmoon, vearth)
    if (abs(dresid(o, o%value, p1)) <= abs(dresid(o, o%value, p2))) then
      p = p1; lchoice(i) = cand(1:1)
    else
      p = p2; lchoice(i) = cand(2:2)
    end if
  end function

  ! Residual of a reading (arcmin): a sign not printed (sgn ?) is the one
  ! that fits the prediction
  function dresid(o, value, p) result(r)
    type(reading), intent(in) :: o
    real(dp), intent(in) :: value, p
    real(dp) :: r, r2
    r = resid_arcmin(o, value, p)
    if (o%sgn == '?') then
      r2 = resid_arcmin(o, -value, p)
      if (abs(r2) < abs(r)) r = r2
    end if
  end function

  function predict_pt(o, x, rsun, rmoon, vearth) result(p)
    type(reading), intent(in) :: o
    real(dp), intent(in) :: x(NP), rsun(3), rmoon(3), vearth(3)
    real(dp) :: p, H, d, sd, Hl, dl, Hr, dr, a, A_az, Hs, ds, sds, Hc, dc, Ho, dov, Hsun, dsun, sdx, h1, h2
    character(len=1) :: ptd
    if (o%body == 3) then
      ! a clock star's hour angle on the armillae
      call refract(wrap180(o%last - o%sra), o%sdec, Hr, dr)
      call armil(Hr, dr, x, o%year, p, dc)
      return
    end if
    if (o%body == 1) then
      call topo(rsun, vearth, o, H, d, sd, merge(x(7)*sin(SUN_SD_RAD), SUN_RADIUS, tycho_model .or. sun_sd_fixed), no_sun_par)
    else
      call topo(rmoon, vearth, o, H, d, sd, MOON_RADIUS*x(iy(18, o%year)))
    end if
    Hsun = 0.0_dp; dsun = 0.0_dp
    if ((use_phase .or. use_irr) .and. o%body == 2) call topo(rsun, vearth, o, Hsun, dsun, sdx, SUN_RADIUS, no_sun_par)
    p = 0.0_dp
    select case (trim(o%qty))
    case ('ALT', 'ALTM')
      if (o%point == 'U' .or. o%point == 'L') then
        call limb_point(H, d, msd(o%point), o%point, Hl, dl)
      else
        Hl = H; dl = d
      end if
      call eq_to_hor(Hl, dl, a, A_az)
      p = a + refraction_deg(a) + x(iy(15, o%year))/60.0_dp
      if (use_malt) p = p + x(95 + NMON + o%win)/60.0_dp
      if (use_ew) p = p + x(ip_ew(o%year, 0))/60.0_dp*sin(A_az*DEG2RAD)
    case ('DEC')
      ! on the armillae the upper and lower limbs of a declination are the
      ! limbs toward the north and the south pole (dec_ns), not in altitude
      ptd = o%point
      if (dec_ns .and. ptd == 'U') ptd = 'A'
      if (dec_ns .and. ptd == 'L') ptd = 'B'
      call limb_point(H, d, msd(ptd), ptd, Hl, dl)
      call refract(Hl, dl, Hr, dr)
      call armil(Hr, dr, x, o%year, Hc, dr)
      p = dr + x(ip_dec(o%year))/60.0_dp
      if (use_mdec) p = p + x(95 + o%win)/60.0_dp
    case ('HA')
      call limb_point(H, d, msd(o%point), o%point, Hl, dl)
      call refract(Hl, dl, Hr, dr)
      call armil(Hr, dr, x, o%year, p, dc)
    case ('AZ')
      call limb_point(H, d, msd(o%point), o%point, Hl, dl)
      call refract(Hl, dl, Hr, dr)
      call eq_to_hor(Hr, dr, a, A_az)
      p = wrap180(A_az - 180.0_dp)
      if (use_ew) p = p + x(ip_ew(o%year, 1))/60.0_dp
    case ('DIAM')
      p = msd('U') + msd('L')
    case ('DIST', 'RADIFF')
      if (trim(o%refkind) == 'STAR') then
        Hs = wrap180(o%last - o%sra); ds = o%sdec; sds = 0.0_dp
        if (use_rstar .and. o%star > 0) then
          Hs = Hs - x(IST0 + 2*o%star - 1)/60.0_dp/max(cos(ds*DEG2RAD), 1.0d-6)
          ds = ds + x(IST0 + 2*o%star)/60.0_dp
        end if
      else
        call topo(rsun, vearth, o, Hs, ds, sds, SUN_RADIUS, no_sun_par)
      end if
      if (trim(o%qty) == 'RADIFF') then
        call limb_point(H, d, msd(o%point), o%point, Hl, dl)
        call refract(Hl, dl, Hc, dc)
        call armil(Hc, dc, x, o%year, Hr, dr)
        call limb_point(Hs, ds, sds, o%refpoint, Hl, dl)
        call refract(Hl, dl, Hc, dc)
        call armil(Hc, dc, x, o%year, Ho, dov)
        p = wrap180(Ho - Hr)
        if (cur_i > 0 .and. o%star > 0) star_dp(:, cur_i) = [-1.0_dp/max(cos(ds*DEG2RAD), 1.0d-6), 0.0_dp]
      else
        call refract(H, d, Hc, dc)
        call refract(Hs, ds, Ho, dov)
        p = sep(Hc, dc, Ho, dov)
        if (cur_i > 0 .and. o%star > 0) then
          ! the separation's change for a star moved 1' east / north
          star_dp(1, cur_i) = (sep(Hc, dc, Ho - 1.0_dp/60.0_dp/max(cos(dov*DEG2RAD), 1.0d-6), dov) - p)*60.0_dp
          star_dp(2, cur_i) = (sep(Hc, dc, Ho, dov + 1.0_dp/60.0_dp) - p)*60.0_dp
        end if
        if (use_mdist .and. o%body == 2) p = p + x(95 + 2*NMON + o%win)/60.0_dp
        select case (o%point)
        case ('N'); p = p - lit_sd('N', H, d, sd, Hsun, dsun, Hs, ds)
        case ('F'); p = p + lit_sd('F', H, d, sd, Hsun, dsun, Hs, ds)
        case ('U', 'L', 'E', 'W')
          if (side(Hc, dc, Ho, dov, o%point) > 0.0_dp) then
            p = p - msd(o%point)
            if (use_dnear) p = p + x(91)/60.0_dp
            if (cur_i > 0) farlimb(cur_i) = .false.
          else
            p = p + msd(o%point)
            if (use_dfar) p = p + x(90)/60.0_dp
            if (cur_i > 0) farlimb(cur_i) = .true.
            if (use_mfar .and. o%body == 2) p = p + x(95 + 3*NMON + o%win)/60.0_dp
          end if
        end select
        if (o%refpoint /= 'C' .and. sds > 0.0_dp) then
          if (side(Ho, dov, Hc, dc, o%refpoint) > 0.0_dp) then
            p = p - sds
          else
            p = p + sds
          end if
        end if
      end if
    end select
  contains
    ! the Moon's (lit) semidiameter toward limb point pt; the full sd for the Sun
    function msd(pt) result(sdl)
      character(len=1), intent(in) :: pt
      real(dp) :: sdl
      real(dp) :: sd2
      sdl = sd
      if (o%body /= 2) return
      sd2 = sd
      if (use_kew .and. index('EWNF', pt) > 0) sd2 = sd*x(89)/x(iy(18, o%year))
      if (use_irr) sd2 = sd2 + x(88)/60.0_dp*0.5_dp*(1.0_dp - cos(sep(H, d, Hsun, dsun)*DEG2RAD))
      if (use_limbo) then
        select case (pt)
        case ('W'); sd2 = sd2 + x(92)/60.0_dp
        case ('E'); sd2 = sd2 + x(93)/60.0_dp
        case ('U', 'A'); sd2 = sd2 + x(94)/60.0_dp
        case ('L', 'B'); sd2 = sd2 + x(95)/60.0_dp
        end select
      end if
      sdl = lit_sd(pt, H, d, sd2, Hsun, dsun, 0.0_dp, 0.0_dp)
    end function
  end function

  ! ── Predictions for all readings for a parameter vector ─────────────
  subroutine predict_all(x, pred)
    real(dp), intent(in) :: x(NP)
    real(dp), intent(out) :: pred(nobs)
    real(dp) :: pos0(3,NBODY), vel0(3,NBODY), gm(NBODY), pos(3,NBODY), vel(3,NBODY), t
    real(dp) :: key(KEYN)
    integer :: k, i
    call initial_state(x, pos0, vel0, gm)
    ! the same orbits for the same readings (only instrument parameters
    ! changed, as in most columns of the Jacobian): reuse the integration
    key(1:14) = x(1:14); key(15) = epoch_jd; key(16) = nobs; key(17) = sum(obs(1:nobs)%jd)
    key(18) = moon_ratio; key(19) = sum(planet_state); key(20) = sum(sun_state); key(21) = nb_used
    key(22) = merge(1, 0, use_planets) + 2*merge(1, 0, use_j2) + 4*merge(1, 0, kepler_moon) + &
              8*merge(1, 0, emb_param) + 16*merge(1, 0, tycho_model)
    key(23) = sum(earth_pole*[1.0_dp, 2.0_dp, 3.0_dp]); key(24) = sum(obs(1:nobs)%jd*order(1:nobs))
    if (allocated(c_pos)) then
      if (size(c_pos, 3) >= nobs .and. all(key == c_key)) then
        do i = 1, nobs
          pred(i) = predict_at(i, x, c_pos(:,:,i), c_vel(:,:,i), gm)
        end do
        return
      end if
    end if
    if (.not. allocated(c_pos)) allocate(c_pos(3, NBODY, MAXOBS), c_vel(3, NBODY, MAXOBS))
    ! forward from the epoch
    pos = pos0; vel = vel0; t = epoch_jd
    do k = 1, nobs
      i = order(k)
      if (obs(i)%jd < epoch_jd) cycle
      call propagate(pos, vel, gm, (obs(i)%jd - t)*DAY_S); t = obs(i)%jd
      c_pos(:,:,i) = pos; c_vel(:,:,i) = vel
      pred(i) = predict_at(i, x, pos, vel, gm)
    end do
    ! backward from the epoch
    pos = pos0; vel = vel0; t = epoch_jd
    do k = nobs, 1, -1
      i = order(k)
      if (obs(i)%jd >= epoch_jd) cycle
      call propagate(pos, vel, gm, (obs(i)%jd - t)*DAY_S); t = obs(i)%jd
      c_pos(:,:,i) = pos; c_vel(:,:,i) = vel
      pred(i) = predict_at(i, x, pos, vel, gm)
    end do
    c_key = key
  end subroutine

  ! Prediction of reading i from the state at its (guessed) moment.  For
  ! tkind A the moment is where the model's own Sun has the hour angle of
  ! the written local apparent time; for tkind M where the body crosses
  ! the meridian.  The Sun and Moon are moved to that moment by a
  ! second-order Taylor step (error < 1 km over an hour).
  function predict_at(i, x, pos, vel, gm) result(p)
    integer, intent(in) :: i
    real(dp), intent(in) :: x(NP), pos(3,NBODY), vel(3,NBODY), gm(NBODY)
    real(dp) :: p, acc(3,NBODY), rs0(3), vs0(3), as0(3), rm0(3), vm0(3), am0(3), rs(3), rm(3)
    real(dp) :: dt, H, d, sd, target, rate
    type(reading) :: o2
    integer :: it
    real(dp), parameter :: SID_DEG_PER_S = 360.98564736629_dp/86400.0_dp
    o2 = obs(i)
    cur_i = i
    call set_site(x, o2%year)
    rs0 = pos(:,1) - pos(:,2); rm0 = pos(:,3) - pos(:,2)
    elong(i) = modulo(atan2(rm0(2), rm0(1)) - atan2(rs0(2), rs0(1)), TAU)*RAD2DEG
    if (o2%tkind == 'S') then
      t_used(i) = o2%jd
      call body_alt(o2, rs0, rm0, vel(:,2), i)
      p = predict_one(o2, x, rs0, rm0, vel(:,2), i); return
    end if
    vs0 = vel(:,1) - vel(:,2); vm0 = vel(:,3) - vel(:,2)
    call accel(pos, gm, acc)
    as0 = acc(:,1) - acc(:,2); am0 = acc(:,3) - acc(:,2)
    dt = 0.0_dp
    do it = 1, 4
      rs = rs0 + vs0*dt + 0.5_dp*as0*dt**2
      rm = rm0 + vm0*dt + 0.5_dp*am0*dt**2
      o2%last = obs(i)%last + SID_DEG_PER_S*dt
      if (o2%tkind == 'A') then
        call topo(rs, vel(:,2), o2, H, d, sd, SUN_RADIUS, no_sun_par)
        target = (o2%lat_h + o2%dclock/60.0_dp - 12.0_dp)*15.0_dp; rate = 360.0_dp
      else if (o2%body == 1) then
        call topo(rs, vel(:,2), o2, H, d, sd, SUN_RADIUS, no_sun_par)
        target = 0.0_dp; rate = 360.0_dp
      else
        call topo(rm, vel(:,2), o2, H, d, sd, MOON_RADIUS)
        target = 0.0_dp; rate = 347.8_dp
      end if
      dt = dt + wrap180(target - H)/rate*86400.0_dp
    end do
    rs = rs0 + vs0*dt + 0.5_dp*as0*dt**2
    rm = rm0 + vm0*dt + 0.5_dp*am0*dt**2
    o2%last = obs(i)%last + SID_DEG_PER_S*dt
    o2%jd = obs(i)%jd + dt/86400.0_dp
    t_used(i) = o2%jd
    call body_alt(o2, rs, rm, vel(:,2), i)
    p = predict_one(o2, x, rs, rm, vel(:,2), i)
  end function

  subroutine body_alt(o, rs, rm, ve, i)
    type(reading), intent(in) :: o
    real(dp), intent(in) :: rs(3), rm(3), ve(3)
    integer, intent(in) :: i
    real(dp) :: H, d, sd, a, A_az
    if (o%body == 3) return
    if (o%body == 1) then
      call topo(rs, ve, o, H, d, sd, SUN_RADIUS)
    else
      call topo(rm, ve, o, H, d, sd, MOON_RADIUS)
    end if
    call eq_to_hor(H, d, alt_read(i), A_az)
    ha_read(i) = H
    if (trim(o%refkind) == 'SUN') then
      call topo(rs, ve, o, H, d, sd, SUN_RADIUS)
      call eq_to_hor(H, d, a, A_az)
      alt_read(i) = min(alt_read(i), a)
    end if
  end subroutine

  ! Moon radius and refraction scales of the further years start at 1
  subroutine one_scales(x)
    real(dp), intent(inout) :: x(NP)
    integer :: k
    do k = 4, NY
      x(ip_year(YEARS(k), 2):ip_year(YEARS(k), 3)) = 1.0_dp
    end do
    x(89) = 1.0_dp
  end subroutine

  ! index of a per-year parameter for a reading of the given year
  integer function iy(k, year)
    integer, intent(in) :: k, year
    iy = k
    if (peryear .and. year == 1595) iy = merge(21, k + 4, k == 15)
    if (peryear .and. year == 1594) iy = merge(26, k + 9, k == 15)
    if (peryear .and. ip_year(year, 0) > 0) iy = ip_year(year, merge(0, k - 16, k == 15))
  end function

  subroutine set_site(x, year)
    real(dp), intent(in) :: x(NP)
    integer, intent(in) :: year
    if (lat_deg /= LAT0_DEG + x(iy(20, year))/60.0_dp) then
      lat_deg = LAT0_DEG + x(iy(20, year))/60.0_dp
      call setup_observer()
    end if
    ref_scale = x(iy(19, year))
    h_scale = x(25)
  end subroutine

  function resid_arcmin(o, value, p) result(r)
    type(reading), intent(in) :: o
    real(dp), intent(in) :: value, p
    real(dp) :: r
    r = (value - p)*60.0_dp
    if (trim(o%qty) == 'HA' .or. trim(o%qty) == 'AZ' .or. trim(o%qty) == 'RADIFF') &
      r = wrap180(value - p)*60.0_dp
  end function

  ! ── Input ───────────────────────────────────────────────────────────
  subroutine read_obs(fname)
    character(len=*), intent(in) :: fname
    integer :: u, ios, k
    character(len=1024) :: line
    character(len=1) :: b
    type(reading) :: o
    real(dp) :: m9(9), err
    open(newunit=u, file=fname, status='old', action='read')
    nobs = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
      read(line, *, iostat=ios) o%id, o%year, b, o%qty, o%point, o%refkind, o%refpoint, &
          o%jd, o%last, m9, o%sra, o%sdec, o%value, o%jpl, err
      if (ios /= 0) then
        print *, 'cannot parse: ', trim(line(1:80)); cycle
      end if
      o%body = merge(1, 2, b == 'S')
      if (b == 'X') o%body = 3          ! a clock star (hour angle), used only for the clocks
      o%M = transpose(reshape(m9, [3, 3]))
      o%group_err = err
      o%group = (yidx(o%year) - 1)*2*NQ + (min(o%body, 2) - 1)*NQ + qindex(o%qty)
      ! blind files carry three more columns: tkind lat_h sgn
      o%tkind = 'S'; o%lat_h = 0.0_dp; o%sgn = '+'; o%sess = 0; o%clk = '-'; o%dclock = 0.0_dp
      read(line, *, iostat=ios) o%id, o%year, b, o%qty, o%point, o%refkind, o%refpoint, &
          o%jd, o%last, m9, o%sra, o%sdec, o%value, o%jpl, err, o%tkind, o%lat_h, o%sgn, o%sess, o%clk
      if (ios /= 0) then
        o%sess = 0; o%clk = '-'
        read(line, *, iostat=ios) o%id, o%year, b, o%qty, o%point, o%refkind, o%refpoint, &
            o%jd, o%last, m9, o%sra, o%sdec, o%value, o%jpl, err, o%tkind, o%lat_h, o%sgn
        if (ios /= 0) then
          o%tkind = 'S'; o%lat_h = 0.0_dp; o%sgn = '+'
        end if
      end if
      nobs = nobs + 1
      obs(nobs) = o
    end do
    close(u)
    nall = nobs
    call set_windows()
    call set_stars()
    obs_all(1:nall) = obs(1:nobs)
    call sort_order()
    k = 0
  end subroutine

  subroutine sort_order()
    ! sort by time (insertion sort on an index)
    integer :: i, j, tmp
    do i = 1, nobs
      order(i) = i
    end do
    do i = 2, nobs
      tmp = order(i); j = i - 1
      do while (j >= 1)
        if (obs(order(j))%jd <= obs(tmp)%jd) exit
        order(j+1) = order(j); j = j - 1
      end do
      order(j+1) = tmp
    end do
  end subroutine

  ! Work on the readings of a time window only (body 0 = both)
  subroutine set_window(t1, t2, body)
    real(dp), intent(in) :: t1, t2
    integer, intent(in) :: body
    integer :: i
    nobs = 0
    do i = 1, nall
      if (obs_all(i)%jd < t1 .or. obs_all(i)%jd > t2) cycle
      if (body /= 0 .and. obs_all(i)%body /= body) cycle
      nobs = nobs + 1
      obs(nobs) = obs_all(i)
    end do
    call sort_order()
  end subroutine

  ! Move the parameter vector to another epoch (osculating elements there)
  subroutine transfer_epoch(x, tnew)
    real(dp), intent(inout) :: x(NP)
    real(dp), intent(in) :: tnew
    real(dp) :: pos(3,NBODY), vel(3,NBODY), gm(NBODY), mu_se, mu_em, a, e, inc, O, w, M
    call initial_state(x, pos, vel, gm)
    call propagate(pos, vel, gm, (tnew - epoch_jd)*DAY_S)
    call state_to_elements(pos, vel, gm, x)
    epoch_jd = tnew
  end subroutine

  ! Osculating elements (in the current parametrization) of a state
  subroutine state_to_elements(pos, vel, gm, x)
    real(dp), intent(in) :: pos(3,NBODY), vel(3,NBODY), gm(NBODY)
    real(dp), intent(inout) :: x(NP)
    real(dp) :: mu_se, mu_em, a, e, inc, O, w, M, rb(3), vb(3)
    mu_se = gm(1) + gm(2); mu_em = gm(2) + gm(3)
    rb = pos(:,2); vb = vel(:,2)
    if (emb_param) then
      mu_se = gm(1) + gm(2) + gm(3)
      rb = (gm(2)*pos(:,2) + gm(3)*pos(:,3))/mu_em
      vb = (gm(2)*vel(:,2) + gm(3)*vel(:,3))/mu_em
    end if
    call cart_to_kepler(rb - pos(:,1), vb - vel(:,1), mu_se, a, e, inc, O, w, M)
    x(1:7) = [e, inc, O, w, M, sqrt(mu_se/a**3)*DAY_S, a]
    call cart_to_kepler(pos(:,3) - pos(:,2), vel(:,3) - vel(:,2), mu_em, a, e, inc, O, w, M)
    x(8:14) = [e, inc, O, w, M, sqrt(mu_em/a**3)*DAY_S, a]
  end subroutine

  subroutine read_state(fname, x_true)
    character(len=*), intent(in) :: fname
    real(dp), intent(out) :: x_true(NP)
    integer :: u, ios, k
    character(len=512) :: line
    character(len=8) :: nm
    real(dp) :: s(6), sun(6), earth(6), moon(6), mu_se, mu_em, gm_s, gm_e, gm_m
    real(dp) :: a, e, i, O, w, M
    gm_s = 132712440041.279419_dp; gm_e = 398600.435507_dp; gm_m = 4902.800118_dp
    open(newunit=u, file=fname, status='old', action='read')
    k = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle
      if (k == 0) then
        read(line, *) epoch_jd; k = 1; cycle
      end if
      read(line, *) nm, s
      select case (trim(nm))
      case ('sun'); sun = s
      case ('earth'); earth = s
      case ('moon'); moon = s
      case ('mercury'); planet_state(:,1) = s
      case ('venus'); planet_state(:,2) = s
      case ('mars'); planet_state(:,3) = s
      case ('jupiter'); planet_state(:,4) = s
      case ('saturn'); planet_state(:,5) = s
      end select
    end do
    close(u)
    sun_state = sun
    mu_se = gm_s + gm_e; mu_em = gm_e + gm_m
    call cart_to_kepler(earth(1:3) - sun(1:3), earth(4:6) - sun(4:6), mu_se, a, e, i, O, w, M)
    x_true(1:7) = [e, i, O, w, M, sqrt(mu_se/a**3)*DAY_S, a]
    call cart_to_kepler(moon(1:3) - earth(1:3), moon(4:6) - earth(4:6), mu_em, a, e, i, O, w, M)
    x_true(8:14) = [e, i, O, w, M, sqrt(mu_em/a**3)*DAY_S, a]
    x_true(15:17) = 0.0_dp
    x_true(18:19) = 1.0_dp; x_true(20) = 0.0_dp
    x_true(21) = 0.0_dp; x_true(22:23) = 1.0_dp; x_true(24) = 0.0_dp; x_true(25) = 1.0_dp
    x_true(26) = 0.0_dp; x_true(27:28) = 1.0_dp; x_true(29:30) = 0.0_dp; x_true(31:NP) = 0.0_dp
    call one_scales(x_true)
  end subroutine

  subroutine setup_observer()
    ! geocentric position of Uraniborg (WGS84), km
    real(dp) :: f, e2, phi, N
    f = 1.0_dp/WGS84_INVF; e2 = f*(2.0_dp - f)
    phi = lat_deg*DEG2RAD
    if (tycho_model .and. .not. use_ell) then
      coslat_rho = RE_KM*cos(phi); sinlat_rho = RE_KM*sin(phi)
      return
    end if
    N = (WGS84_RADIUS/1000.0_dp)/sqrt(1.0_dp - e2*sin(phi)**2)
    coslat_rho = (N + ELEV_KM)*cos(phi)
    sinlat_rho = (N*(1.0_dp - e2) + ELEV_KM)*sin(phi)
  end subroutine

  ! ── Robust per-group sigma (arcmin) ─────────────────────────────────
  subroutine update_sigmas(r, use)
    real(dp), intent(in) :: r(nobs)
    logical, intent(in) :: use(nobs)
    real(dp) :: v(MAXOBS)
    integer :: g, i, n
    if (freeze_sigma .or. sigma_fixed) return
    do g = 1, NG
      n = 0
      do i = 1, nobs
        if (obs(i)%group == g .and. use(i) .and. obs(i)%body /= 3) then
          n = n + 1; v(n) = abs(r(i))
        end if
      end do
      if (n >= 5) sigma_g(g) = max(0.2_dp, 1.4826_dp*median(v(1:n)))
    end do
  end subroutine

  function median(v) result(m)
    real(dp), intent(in) :: v(:)
    real(dp) :: m, w(size(v)), t
    integer :: i, j, n
    n = size(v); w = v
    do i = 2, n
      t = w(i); j = i - 1
      do while (j >= 1)
        if (w(j) <= t) exit
        w(j+1) = w(j); j = j - 1
      end do
      w(j+1) = t
    end do
    if (mod(n, 2) == 1) then
      m = w((n+1)/2)
    else
      m = 0.5_dp*(w(n/2) + w(n/2+1))
    end if
  end function

  ! ── Finite-difference steps ─────────────────────────────────────────
  ! Natural scale of each parameter: the Jacobian is formed per unit of
  ! this scale, so the normal equations are well conditioned.
  function pscale(x, k) result(s)
    real(dp), intent(in) :: x(NP)
    integer, intent(in) :: k
    real(dp) :: s
    select case (k)
    case (1, 8);           s = 1.0d-2
    case (2:5, 9:12);      s = 1.0d-2
    case (6, 13);          s = abs(x(k))*1.0d-5
    case (7, 14);          s = abs(x(k))*1.0d-2
    case (18, 19, 22, 23, 25, 27, 28, 45, 46, 54, 55, 63, 64, 72, 73, 81, 82, 89); s = 1.0d-2
    case default;          s = 1.0_dp
    end select
  end function

  function fd_step(x, k) result(h)
    real(dp), intent(in) :: x(NP)
    integer, intent(in) :: k
    real(dp) :: h
    select case (k)
    case (1, 8);           h = 1.0d-6
    case (2:5, 9:12);      h = 1.0d-6
    case (6, 13);          h = abs(x(k))*1.0d-7
    case (7, 14);          h = abs(x(k))*1.0d-5
    case (18, 19, 22, 23, 25, 27, 28, 45, 46, 54, 55, 63, 64, 72, 73, 81, 82, 89); h = 1.0d-4
    case default;          h = 0.05_dp
    end select
  end function

  subroutine jacobian(x, values, J, r0, active)
    real(dp), intent(in) :: x(NP), values(nobs)
    real(dp), intent(out) :: J(nobs, NP), r0(nobs)
    logical, intent(in) :: active(NP)
    real(dp) :: xp(NP), xm(NP), pp(nobs), pm(nobs), p0(nobs), h
    integer :: k, i
    call predict_all(x, p0)
    do i = 1, nobs
      r0(i) = dresid(obs(i), values(i), p0(i))
    end do
    J = 0.0_dp
    fix_choice = .true.          ! the same limbs for the perturbed orbits
    do k = 1, NP
      if (.not. active(k) .or. .not. param_on(k)) cycle
      if (k > IST0) then
        ! a reference star's offset
        do i = 1, nobs
          if (obs(i)%star == (k - IST0 + 1)/2) J(i, k) = star_dp(2 - mod(k - IST0, 2), i)
        end do
        cycle
      end if
      if (k > 95) then
        ! a window's zero: the prediction moves one-for-one (arcmin)
        do i = 1, nobs
          if (obs(i)%body == 3) cycle
          if (k <= 95 + NMON) then
            if (trim(obs(i)%qty) == 'DEC' .and. obs(i)%win == k - 95) J(i, k) = 1.0_dp
          else if (k <= 95 + 2*NMON) then
            if (obs(i)%qty(1:3) == 'ALT' .and. obs(i)%win == k - 95 - NMON) J(i, k) = 1.0_dp
          else if (k <= 95 + 3*NMON) then
            if (trim(obs(i)%qty) == 'DIST' .and. obs(i)%body == 2 .and. &
                obs(i)%win == k - 95 - 2*NMON) J(i, k) = 1.0_dp
          else
            if (trim(obs(i)%qty) == 'DIST' .and. obs(i)%body == 2 .and. farlimb(i) .and. &
                obs(i)%win == k - 95 - 3*NMON) J(i, k) = 1.0_dp
          end if
        end do
        cycle
      end if
      h = fd_step(x, k)
      xp = x; xp(k) = x(k) + h
      xm = x; xm(k) = x(k) - h
      call predict_all(xp, pp)
      call predict_all(xm, pm)
      do i = 1, nobs
        ! derivative of the prediction, in arcmin per unit parameter
        J(i, k) = resid_arcmin(obs(i), pp(i), pm(i))/(2.0_dp*h)*pscale(x, k)   ! d(prediction)/d(scaled param), arcmin
      end do
    end do
    fix_choice = .false.
  end subroutine

  ! ── Symmetric positive-definite solve (Cholesky) ────────────────────
  subroutine chol_solve(A, b, x, ok)
    real(dp), intent(in) :: A(:,:), b(:)
    real(dp), intent(out) :: x(size(b))
    logical, intent(out) :: ok
    real(dp) :: L(size(b), size(b)), s, y(size(b))
    integer :: n, i, j, k
    n = size(b); L = 0.0_dp; ok = .true.
    do j = 1, n
      s = A(j,j) - sum(L(j,1:j-1)**2)
      if (s <= 0.0_dp) then
        ok = .false.; x = 0.0_dp; return
      end if
      L(j,j) = sqrt(s)
      do i = j+1, n
        L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1)))/L(j,j)
      end do
    end do
    do i = 1, n
      y(i) = (b(i) - sum(L(i,1:i-1)*y(1:i-1)))/L(i,i)
    end do
    do i = n, 1, -1
      x(i) = (y(i) - sum(L(i+1:n,i)*x(i+1:n)))/L(i,i)
    end do
    k = 0
  end subroutine

  ! ── Batch least squares (Levenberg–Marquardt, iteratively reweighted,
  !    with optional rejection beyond clip*sigma) ─────────────────────
  subroutine fit_lm(x, values, active, clip, cov, use, rms, verbose)
    real(dp), intent(inout) :: x(NP)
    real(dp), intent(in) :: values(nobs), clip
    logical, intent(in) :: active(NP), verbose
    real(dp), intent(out) :: cov(NP, NP), rms
    logical, intent(inout) :: use(nobs)
    real(dp), allocatable :: J(:,:), r(:)
    real(dp) :: N(NP,NP), g(NP), dx(NP), lam, chi2, chi2n, xn(NP), w, pn(nobs), rn(nobs), Ninv(NP,NP), e(NP)
    real(dp) :: wr(nobs), uu
    integer :: it, k, i, na, idx(NP), m, nstall
    logical :: ok
    allocate(J(nobs, NP), r(nobs))
    nstall = 0
    if (.not. freeze_use) use = .true.
    if (.not. (freeze_sigma .or. sigma_fixed)) sigma_g = 1.0_dp
    wr = 1.0_dp
    lam = 1.0d-3
    na = 0
    do k = 1, NP
      if (active(k) .and. param_on(k)) then
        na = na + 1; idx(na) = k
      end if
    end do
    do it = 1, 80
      call jacobian(x, values, J, r, active)
      call update_sigmas(r, use)
      if (robust_kind == 0 .and. clip > 0.0_dp .and. .not. freeze_use) then
        do i = 1, nobs
          use(i) = abs(r(i)) <= clip*sigma_g(obs(i)%group)
        end do
        call update_sigmas(r, use)
      end if
      ! iteratively reweighted: the weights of this iteration from its residuals
      do i = 1, nobs
        uu = abs(r(i))/sigma_g(obs(i)%group)
        select case (robust_kind)
        case (1); wr(i) = merge(1.0_dp, 1.5_dp/max(uu, 1.0d-12), uu <= 1.5_dp)
        case (2); wr(i) = merge((1.0_dp - (uu/4.685_dp)**2)**2, 0.0_dp, uu < 4.685_dp)
        case default; wr(i) = 1.0_dp
        end select
        if (robust_kind == 2 .and. .not. freeze_use) use(i) = wr(i) > 0.0_dp
        if (robust_kind == 2 .and. freeze_use .and. use(i)) wr(i) = max(wr(i), 1.0d-3)
        if (obs(i)%body == 3) use(i) = .false.     ! clock stars are not fitted
      end do
      chi2 = 0.0_dp; N = 0.0_dp; g = 0.0_dp
      do i = 1, nobs
        if (.not. use(i)) cycle
        w = wr(i)/sigma_g(obs(i)%group)**2
        chi2 = chi2 + w*r(i)**2
        do k = 1, na
          if (J(i, idx(k)) == 0.0_dp) cycle
          g(idx(k)) = g(idx(k)) + w*J(i, idx(k))*r(i)
          N(idx(k), idx(1:na)) = N(idx(k), idx(1:na)) + w*J(i, idx(k))*J(i, idx(1:na))
        end do
      end do
      ! random effects: a campaign's value is drawn about zero with spread tau
      do k = 1, na
        uu = prior_tau(idx(k))
        if (uu > 0.0_dp .and. N(idx(k), idx(k)) > 0.0_dp) then
          N(idx(k), idx(k)) = N(idx(k), idx(k)) + 1.0_dp/uu**2
          g(idx(k)) = g(idx(k)) - x(idx(k))/uu**2
          chi2 = chi2 + (x(idx(k))/uu)**2
        end if
      end do
      ! a parameter no used reading depends on stays where it is
      do k = 1, na
        if (N(idx(k), idx(k)) <= 0.0_dp) then
          N(idx(k), idx(k)) = 1.0_dp; g(idx(k)) = 0.0_dp
        end if
      end do
      ! Marquardt loop
      do m = 1, 12
        block
          real(dp) :: A(na,na), b(na), s(na)
          A = N(idx(1:na), idx(1:na)); b = g(idx(1:na))
          do k = 1, na
            A(k,k) = A(k,k)*(1.0_dp + lam)
          end do
          call chol_solve(A, b, s, ok)
          if (.not. ok .and. verbose) print '(A,17ES10.2)', '   singular; diag(N) =', (N(idx(k),idx(k)), k = 1, na)
          dx = 0.0_dp
          ! at most 20 natural units per step (0.2 in e, 11 deg, 20 % in a)
          if (ok .and. maxval(abs(s)) > 20.0_dp) s = s*20.0_dp/maxval(abs(s))
          if (ok) dx(idx(1:na)) = s
          do k = 1, NP
            dx(k) = dx(k)*pscale(x, k)
          end do
        end block
        xn = x + dx
        call normalize_elements(xn)
        call predict_all(xn, pn)
        chi2n = 0.0_dp
        do i = 1, nobs
          if (.not. use(i)) cycle
          rn(i) = dresid(obs(i), values(i), pn(i))
          chi2n = chi2n + wr(i)*rn(i)**2/sigma_g(obs(i)%group)**2
        end do
        do k = 1, na
          uu = prior_tau(idx(k))
          if (uu > 0.0_dp .and. xn(idx(k)) /= 0.0_dp) chi2n = chi2n + (xn(idx(k))/uu)**2
        end do
        if (ok .and. chi2n < chi2) then
          lam = max(lam/5.0_dp, 1.0d-9); exit
        end if
        lam = lam*10.0_dp
      end do
      if (verbose) print '(A,I3,A,ES12.5,A,ES12.5,A,ES9.2,A,I6)', '   iter', it, '  chi2', chi2, ' -> ', chi2n, &
          '  lambda', lam, '  used', count(use)
      if (chi2n < chi2) then
        x = xn
        nstall = 0
      else
        nstall = nstall + 1
      end if
      if (abs(chi2 - chi2n) < 1.0d-7*chi2) nstall = nstall + 1
      if (nstall >= 3 .and. it > 5) exit
    end do
    ! covariance of the active parameters
    cov = 0.0_dp
    do k = 1, na
      e = 0.0_dp; e(k) = 1.0_dp
      block
        real(dp) :: A(na,na), s(na)
        A = N(idx(1:na), idx(1:na))
        call chol_solve(A, e(1:na), s, ok)
        if (ok) cov(idx(1:na), idx(k)) = s
      end block
    end do
    do k = 1, NP
      do i = 1, NP
        cov(i, k) = cov(i, k)*pscale(x, i)*pscale(x, k)
      end do
    end do
    call predict_all(x, pn)
    rms = 0.0_dp; k = 0
    do i = 1, nobs
      if (.not. use(i)) cycle
      rms = rms + dresid(obs(i), values(i), pn(i))**2; k = k + 1
    end do
    rms = sqrt(rms/max(k, 1))
    last_chi2 = 0.0_dp
    do i = 1, nobs
      if (use(i)) last_chi2 = last_chi2 + wr(i)*(dresid(obs(i), values(i), pn(i))/sigma_g(obs(i)%group))**2
    end do
    Ninv = 0.0_dp
    deallocate(J, r)
  end subroutine

  ! ── Clock of each daytime session from the model's own Sun ──────────
  ! Tycho's local apparent time comes from his clock.  Where no star fixes
  ! the true time, the Sun's hour angles, altitudes and azimuths read in
  ! the same session do (as Tycho himself corrected his clocks): the clock
  ! error that cancels each one is residual / (rate of change with time);
  ! a reading takes the median of those within 45 minutes (at least
  ! three), else of its session's within 1.5 hours.  The correction is
  ! added to what is already applied, so calling this repeatedly between
  ! fits converges.
  subroutine update_clocks(x, nsess, medcorr)
    real(dp), intent(in) :: x(NP)
    integer, intent(out) :: nsess
    real(dp), intent(out) :: medcorr
    real(dp) :: p0(nobs), p1(nobs), rate, off(nobs), corr(nobs), tmp(nobs), cand(nobs), cj(nobs)
    logical :: isclk(nobs)
    integer :: i, j, nc, ncand, ns
    call predict_all(x, p0)
    do i = 1, nobs
      if (obs(i)%tkind == 'A') obs(i)%dclock = obs(i)%dclock + 1.0_dp
    end do
    call predict_all(x, p1)
    do i = 1, nobs
      if (obs(i)%tkind == 'A') obs(i)%dclock = obs(i)%dclock - 1.0_dp
    end do
    isclk = .false.
    do i = 1, nobs
      if (obs(i)%tkind /= 'A' .or. obs(i)%sess == 0 .or. obs(i)%body == 2 .or. obs(i)%clk == 'x') cycle
      if (trim(obs(i)%qty) /= 'HA' .and. trim(obs(i)%qty) /= 'ALT' .and. trim(obs(i)%qty) /= 'AZ') cycle
      if (trim(obs(i)%qty) == 'ALT' .and. obs(i)%value < 5.0_dp) cycle
      rate = resid_arcmin(obs(i), p1(i), p0(i))          ! arcmin per minute of time
      if (abs(rate) < 3.0_dp) cycle
      off(i) = dresid(obs(i), obs(i)%value, p0(i))/rate
      if (abs(off(i)) < 90.0_dp) isclk(i) = .true.
    end do
    ! as in compare_jpl.py: in a session with hour angles (read on the
    ! armillae, a direct measure of time) only those set the clock;
    ! altitudes and azimuths depend on the quadrant's zero point and on
    ! refraction
    do i = 1, nobs
      if (.not. isclk(i) .or. trim(obs(i)%qty) == 'HA') cycle
      do j = 1, nobs
        if (isclk(j) .and. obs(j)%sess == obs(i)%sess .and. trim(obs(j)%qty) == 'HA') then
          isclk(i) = .false.; exit
        end if
      end do
    end do
    corr = 0.0_dp; nc = 0
    do j = 1, nobs
      if (obs(j)%tkind /= 'A' .or. obs(j)%sess == 0) cycle
      ncand = 0
      do i = 1, nobs
        if (isclk(i) .and. obs(i)%sess == obs(j)%sess .and. abs(obs(i)%jd - obs(j)%jd)*24.0_dp <= 0.75_dp) then
          ncand = ncand + 1; cand(ncand) = off(i)
        end if
      end do
      if (ncand < 3) then
        ncand = 0
        do i = 1, nobs
          if (isclk(i) .and. obs(i)%sess == obs(j)%sess .and. abs(obs(i)%jd - obs(j)%jd)*24.0_dp <= 1.5_dp) then
            ncand = ncand + 1; cand(ncand) = off(i)
          end if
        end do
      end if
      if (ncand == 0) cycle
      corr(j) = median(cand(1:ncand))
      nc = nc + 1; cj(nc) = abs(corr(j))
    end do
    do j = 1, nobs
      obs(j)%dclock = obs(j)%dclock + corr(j)
    end do
    ns = 0
    do i = 1, nobs
      if (obs(i)%sess > ns) ns = obs(i)%sess
    end do
    nsess = nc
    medcorr = 0.0_dp
    if (nc > 0) medcorr = median(cj(1:nc))
    tmp = 0.0_dp
  end subroutine

  ! A step through e = 0 is the same orbit with the perihelion on the other
  ! side: (-e, w, M) = (e, w + pi, M + pi).  A negative inclination likewise
  ! is (|i|, Omega + pi, w + pi).
  subroutine normalize_elements(x)
    real(dp), intent(inout) :: x(NP)
    integer :: b
    do b = 0, 7, 7
      if (x(1+b) < 0.0_dp) then
        x(1+b) = -x(1+b); x(4+b) = x(4+b) + PI; x(5+b) = x(5+b) + PI
      end if
      x(1+b) = max(1.0d-7, min(0.5_dp, x(1+b)))
      if (x(2+b) < 0.0_dp) then
        x(2+b) = -x(2+b); x(3+b) = x(3+b) + PI; x(4+b) = x(4+b) + PI
      end if
      x(3+b) = modulo(x(3+b), TAU); x(4+b) = modulo(x(4+b), TAU); x(5+b) = modulo(x(5+b), TAU)
    end do
  end subroutine

  ! ── Iterated extended Kalman filter: sequential measurement updates
  !    of the constant parameter vector (F = I, Q = 0) ────────────────
  subroutine fit_iekf(x, values, active, use, sigma_prior, n_outer, P)
    real(dp), intent(inout) :: x(NP)
    real(dp), intent(in) :: values(nobs), sigma_prior(NP)
    logical, intent(in) :: active(NP), use(nobs)
    integer, intent(in) :: n_outer
    real(dp), intent(out) :: P(NP,NP)
    real(dp), allocatable :: J(:,:), r(:)
    real(dp) :: xr(NP), xs(NP), Hrow(NP), PH(NP), S, gain(NP), y
    integer :: it, i, k, kk
    allocate(J(nobs, NP), r(nobs))
    do it = 1, n_outer
      xr = x
      call jacobian(xr, values, J, r, active)
      P = 0.0_dp
      do k = 1, NP
        if (active(k) .and. param_on(k)) P(k,k) = sigma_prior(k)**2
      end do
      do k = 1, NP
        if (active(k) .and. param_on(k)) P(k,k) = P(k,k)/pscale(xr, k)**2
      end do
      ! xs, P in scaled units relative to xr
      xs = 0.0_dp
      do kk = 1, nobs
        i = order(kk)
        if (.not. use(i)) cycle
        Hrow = J(i, :)
        ! innovation of the linearised measurement about the reference xr
        y = r(i) - dot_product(Hrow, xs)
        PH = matmul(P, Hrow)
        S = dot_product(Hrow, PH) + sigma_g(obs(i)%group)**2
        gain = PH/S
        xs = xs + gain*y
        P = P - spread(gain, 2, NP)*spread(PH, 1, NP)
      end do
      do k = 1, NP
        x(k) = xr(k) + xs(k)*pscale(xr, k)
      end do
      call normalize_elements(x)
      print '(A,I3,A,ES10.3,A,F14.1)', '   EKF pass', it, '   largest scaled update', maxval(abs(xs)), &
          '   a_E', x(7)
      do k = 1, NP
        do kk = 1, NP
          P(k, kk) = P(k, kk)*pscale(xr, k)*pscale(xr, kk)
        end do
      end do
      if (maxval(abs(xs)) < 1.0d-4) exit
    end do
    deallocate(J, r)
  end subroutine

end module tycho_fit_mod

! ═══════════════════════════════════════════════════════════════════════
program tycho_fit
  use constants_mod
  use tycho_fit_mod
  implicit none
  character(len=512) :: obs_file, state_file, mode, arg
  real(dp) :: x_true(NP), x(NP), x0(NP), cov(NP,NP), P(NP,NP), rms, clip
  real(dp), allocatable :: pred(:), values(:), synth(:)
  real(dp) :: sig(NP), r
  logical :: active(NP)
  logical, allocatable :: use(:)
  integer :: i, k, g, n, n_mc, mc, seed0
  real(dp) :: sigma_prior(NP), gn(2), best_sigma(NG), chi_best
  real(dp), parameter :: FACTORS(15) = [0.3_dp, 0.5_dp, 0.7_dp, 0.8_dp, 0.9_dp, 0.95_dp, 1.0_dp, 1.05_dp, &
                                        1.1_dp, 1.2_dp, 1.4_dp, 2.0_dp, 3.0_dp, 5.0_dp, 20.0_dp]

  ! usage: tycho_fit <obs> <state> <mode> [clip] [n] [seed0]
  !   mode = check | fit | mc | profile | kepler
  call get_command_argument(1, obs_file)
  call get_command_argument(2, state_file)
  call get_command_argument(3, mode)
  if (trim(mode) == 'alien') then
    alien = .true.; mode = 'blind'; camp_gap = ALIEN_GAP
    print '(A)', ' THE PROCEDURE (data/fit/ALIEN_PROCEDURE.md): model ' // ALIEN_MODEL // ', subset ' // ALIEN_SUBSET
  end if
  if (len_trim(mode) == 0) mode = 'fit'
  clip = 0.0_dp; n_mc = 10; seed0 = 1
  call get_command_argument(4, arg); if (len_trim(arg) > 0) read(arg, *) clip
  call get_command_argument(5, arg); if (len_trim(arg) > 0) read(arg, *) n_mc
  call get_command_argument(6, arg); if (len_trim(arg) > 0) read(arg, *) seed0
  call get_arg(7, arg); if (trim(arg) == 'noplanets') use_planets = .false.

  call get_arg(11, arg)
  if (index(arg, 'nosunpar') > 0) no_sun_par = .true.
  if (index(arg, 'sunfix') > 0) sun_sd_fixed = .true.
  if (index(arg, 'tycho') > 0) then
    tycho_model = .true.; use_aberr = .false.; use_j2 = .false.; moon_ratio = 0.0_dp
    ! DIAGNOSTIC switches: modern pieces back, one at a time
    if (index(arg, '+saem') > 0) use_saem = .true.          ! Saemundsson refraction (scaled)
    if (index(arg, '+aberr') > 0) use_aberr = .true.        ! annual aberration
    if (index(arg, '+j2') > 0) use_j2 = .true.              ! Earth's J2
    if (index(arg, '+mass') > 0) moon_ratio = 0.0123000371_dp   ! Moon/Earth mass ratio
    if (index(arg, '+ell') > 0) use_ell = .true.            ! WGS84 ellipsoid for the observer
  end if
  call init_names()
  block
    character(len=32) :: wv
    integer :: st
    call get_environment_variable('TYCHO_WIN', wv, status=st)
    if (st == 0 .and. len_trim(wv) > 0) read(wv, *) win_days
    if (win_days /= 30.0_dp) print '(A,F6.1,A)', '   windows of ', win_days, ' days'
    call get_environment_variable('TYCHO_GAP', wv, status=st)
    if (st == 0 .and. len_trim(wv) > 0) read(wv, *) camp_gap
  end block
  call setup_observer()
  call read_obs(trim(obs_file))
  if (alien) then
    ! the procedure reads nothing from the state file but the epoch (a
    ! date chosen in the middle of the readings); JPL's orbits are not
    ! loaded at all (compare_orbits.py compares afterwards)
    block
      integer :: un, ios
      character(len=512) :: ln
      open(newunit=un, file=trim(state_file), status='old', action='read')
      do
        read(un, '(A)', iostat=ios) ln
        if (ios /= 0) stop 'no epoch in the state file'
        if (ln(1:1) == '#' .or. len_trim(ln) == 0) cycle
        read(ln, *) epoch_jd; exit
      end do
      close(un)
    end block
    x_true = 0.0_dp; sun_state = 0.0_dp; planet_state = 0.0_dp
  else
    call read_state(trim(state_file), x_true)
  end if
  ! Earth's pole (ICRS) at the epoch: third row of the rotation to the
  ! equator of date of the reading nearest to the epoch
  k = minloc([(abs(obs(i)%jd - epoch_jd), i = 1, nobs)], 1)
  earth_pole = obs(k)%M(3,:)
  call get_command_argument(8, arg); if (trim(arg) == 'noj2') use_j2 = .false.
  call get_arg(7, arg)
  if (trim(mode) == 'blind' .or. trim(arg) == 'blind') then
    ! nothing from JPL in the model: Sun, Earth and Moon only, the Sun
    ! starting at rest at the origin
    nb_used = 3; use_planets = .false.; sun_state = 0.0_dp
    ! DIAGNOSTIC (argument 11 contains 'planets'): the planets as perturbers,
    ! started from DE441
    call get_arg(11, arg)
    if (index(arg, 'planets') > 0 .and. .not. alien) then
      nb_used = 8; use_planets = .true.
      call read_state(trim(state_file), x_true)
      print '(A)', '   DIAGNOSTIC: planets from DE441 as perturbers'
    end if
  end if
  allocate(pred(nobs), values(nobs), synth(nobs), use(nobs))
  block
    character(len=512) :: sf
    integer :: st, un, g
    call get_environment_variable('TYCHO_SIGMAS', sf, status=st)
    if (st == 0 .and. len_trim(sf) > 0) then
      open(newunit=un, file=trim(sf), status='old', action='read')
      do g = 1, NG
        read(un, *) sigma_g(g)
      end do
      close(un)
      sigma_fixed = .true.
      print '(A,A)', '   DIAGNOSTIC: group sigmas fixed from ', trim(sf)
    end if
  end block
  print '(A,I6,A,A,A,A)', ' Readings: ', nobs, '  from ', trim(obs_file), '   mode ', trim(mode)
  print '(A,F12.3,A,F5.1,A,L2)', ' Epoch (JD TDB): ', epoch_jd, '   clip (sigma, 0 = none): ', clip, &
      '   planets:', use_planets
  values = [(obs(i)%value, i = 1, nobs)]
  active = .true.
  if (.not. alien) call predict_all(x_true, pred)

  select case (trim(mode))

  case ('check')
    call get_arg(10, arg)
    if (trim(arg) == 'clock') then
      do k = 1, 3
        call update_clocks(x_true, n, r)
      end do
      call predict_all(x_true, pred)
      print '(A,I6,A,F8.3,A)', ' session clocks from the (JPL) Sun:', n, ' readings, last median change', r, ' min'
    end if
    block
      integer :: un
      open(newunit=un, file='check_times.dat', status='replace')
      write(un, '(A)') '# id tkind jd_used dclock_min model_deg jpl_deg value_deg'
      do i = 1, nobs
        write(un, '(A12,1X,A1,F18.8,F10.3,3F14.7)') obs(i)%id, obs(i)%tkind, t_used(i), obs(i)%dclock, pred(i), &
            obs(i)%jpl, obs(i)%value
      end do
      close(un)
    end block
    print '(/A)', ' Model check: this program with the JPL orbits vs. compare_jpl.py (JPL DE441)'
    print '(A)', '    group               n    rms(model-JPL)   rms(Tycho-JPL)'
    do g = 1, NG
      n = 0; gn = 0.0_dp
      do i = 1, nobs
        if (obs(i)%group /= g) cycle
        n = n + 1
        gn(1) = gn(1) + resid_arcmin(obs(i), pred(i), obs(i)%jpl)**2
        gn(2) = gn(2) + dresid(obs(i), obs(i)%value, obs(i)%jpl)**2
      end do
      if (n > 0) print '(4X,A,I6,F12.2,A,F12.2,A)', group_name(g), n, sqrt(gn(1)/n)*60, '"', sqrt(gn(2)/n), "'"
    end do

  case ('mc')
    ! synthetic data: the model with JPL's orbits + Gaussian noise with
    ! the scatter of Tycho's readings about JPL in each group
    call tycho_sigmas()
    do mc = seed0, seed0 + n_mc - 1
      call add_noise(pred, synth, mc)
      x = start_guess(x_true)
      call fit_lm(x, synth, active, 0.0_dp, cov, use, rms, .false.)
      write(*, '(A,I5,42ES20.11)') ' MC ', mc, x - x_true
    end do

  case ('fit')
    print '(/A)', " Tycho's readings: batch least squares (Levenberg-Marquardt)"
    x = start_guess(x_true); x0 = x
    call fit_lm(x, values, active, clip, cov, use, rms, .true.)
    do k = 1, NP
      sig(k) = sqrt(max(cov(k,k), 0.0_dp))
    end do
    print '(A,F8.3,A,I6,A,I6,A,F12.2)', '   rms residual ', rms, "'   readings used ", count(use), ' of ', nobs, &
        '   chi2 ', last_chi2
    call report(" Fitted orbits from Tycho's readings (formal 1-sigma)", x, sig, x_true)
    call print_sigmas()
    call write_residuals(x, use)
    write(*, '(A,755ES20.11)') ' FIT ', x
    write(*, '(A,755ES20.11)') ' SIG ', sig
    print '(/A)', ' Iterated extended Kalman filter (sequential updates, F = I, Q = 0), same readings and weights'
    sigma_prior = [0.05_dp, 0.05_dp, 0.05_dp, 0.5_dp, 0.5_dp, 1.0d-4, 1.0d7, &
                   0.05_dp, 0.05_dp, 0.05_dp, 0.5_dp, 0.5_dp, 1.0d-2, 4.0d4, 5.0_dp, 5.0_dp, 5.0_dp, &
                   0.1_dp, 0.3_dp, 5.0_dp, 5.0_dp, 0.1_dp, 0.3_dp, 5.0_dp, 0.3_dp, &
                   5.0_dp, 0.1_dp, 0.3_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, &
                   5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, &
                   (5.0_dp, 5.0_dp, 0.1_dp, 0.3_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, ii_ = 1, 5), &
                   1.0_dp, 0.1_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, &
                   (5.0_dp, ii_ = 1, 4*NMON), (1.0_dp, ii_ = 1, 2*NSTAR)]
    x = x0
    call fit_iekf(x, values, active, use, sigma_prior, 30, P)
    do k = 1, NP
      sig(k) = sqrt(max(P(k,k), 0.0_dp))
    end do
    call predict_all(x, pred)
    rms = 0.0_dp; n = 0; r = 0.0_dp
    do i = 1, nobs
      if (.not. use(i)) cycle
      gn(1) = dresid(obs(i), values(i), pred(i))
      rms = rms + gn(1)**2; r = r + (gn(1)/sigma_g(obs(i)%group))**2; n = n + 1
    end do
    print '(A,F8.3,A,F12.2)', '   rms residual ', sqrt(rms/n), "'   chi2 ", r
    call report(' Iterated EKF result (1-sigma from the filter covariance)', x, sig, x_true)
    write(*, '(A,755ES20.11)') ' EKF ', x

  case ('profile')
    ! chi^2 as a function of the Earth-Sun distance a_E, all other
    ! parameters refitted, with the weights of the free fit
    x = start_guess(x_true)
    call fit_lm(x, values, active, clip, cov, use, rms, .false.)
    best_sigma = sigma_g; chi_best = last_chi2
    freeze_sigma = .true.; freeze_use = .true.
    print '(A,F14.1,A,F12.3)', ' free fit: a_E = ', x(7), ' km   chi2 = ', last_chi2
    print '(A)', '  a_E/JPL      a_E (km)        chi2     dchi2    rms(arcmin)'
    active(7) = .false.
    do k = 1, size(FACTORS)
      if (.not. sigma_fixed) sigma_g = best_sigma
      x0 = x; x0(7) = x_true(7)*FACTORS(k)       ! warm start from the free fit
      call fit_lm(x0, values, active, clip, cov, use, rms, .false.)
      print '(F9.3,F16.1,F12.3,F10.3,F10.3,I7)', FACTORS(k), x0(7), last_chi2, last_chi2 - chi_best, rms, count(use)
    end do

  case ('jplfit', 'jplfit0')
    ! noise-free data: JPL DE441's own predictions of the same readings.
    ! The fitted orbits then show the bias of this model itself.
    values = [(obs(i)%jpl, i = 1, nobs)]
    x = start_guess(x_true)
    if (trim(mode) == 'jplfit0') x = x_true
    call fit_lm(x, values, active, 0.0_dp, cov, use, rms, .false.)
    do k = 1, NP
      sig(k) = 0.0_dp
    end do
    print '(A,F8.3,A)', '   rms residual (model vs JPL DE441) ', rms*60, '"'
    call report(' Fit to JPL DE441''s own predictions (bias of the 3-body model)', x, sig, x_true)

  case ('blind')
    call blind_fit()

  case ('kepler')
    ! the Moon on a fixed Kepler ellipse about the Earth (no Sun)
    kepler_moon = .true.
    x = start_guess(x_true)
    call fit_lm(x, values, active, clip, cov, use, rms, .true.)
    print '(A,F8.3,A)', '   rms residual with a fixed Kepler ellipse for the Moon: ', rms, "'"
    call print_sigmas()

  end select

contains

  ! ═══ Blind fit: Tycho's readings only, nothing from JPL ═══════════════
  ! The orbits are built up from the readings themselves:
  !   1. the Sun: a Kepler ellipse fitted to all Sun readings, from a
  !      circle (eccentricity 0.02, 8 starting longitudes); the year of the
  !      Julian calendar and the obliquity from the extreme declinations
  !      only start it;
  !   2. the Moon: a Kepler ellipse fitted to the ten days of 1586 with the
  !      most Moon readings, from 72 starts (node and longitude every 60 and
  !      30 deg, 5 deg to the fitted ecliptic, the sidereal month and the
  !      distance of 60 Earth radii known since antiquity);
  !   3. the Moon with the Sun's pull (N-body), the window widened step by
  !      step to all of 1586;
  !   4. the Earth-Sun distance: scanned from Tycho's 1150 Earth radii up,
  !      everything else refitted at each step;
  !   5. all readings of 1586 and 1595 together, every parameter free,
  !      at the epoch 1590-07-01; outliers rejected from the fit's own
  !      residuals (4 sigma, robust scatter per kind of reading).
  ! JPL's orbits (x_true) are read only to print the comparison at the end.
  subroutine blind_fit()
    real(dp) :: xb(NP), xt(NP), xbest(NP), xk(NP), score, best, tc, epoch0, half, inc, Om, ob
    real(dp) :: nE(3), Xv(3), Yv(3), nm(3), sc_scan(40), cv(NP,NP), rr, sgv(NP)
    real(dp), allocatable :: dec(:)
    logical :: act(NP)
    integer :: k1, k2, nd, nw, nbest, kk
    real(dp), parameter :: RE = 6378.137_dp
    real(dp), parameter :: HALVES(7) = [5.0_dp, 10.0_dp, 20.0_dp, 40.0_dp, 80.0_dp, 160.0_dp, 400.0_dp]
    real(dp), parameter :: AE_RE(14) = [1150.0_dp, 2000.0_dp, 4000.0_dp, 7000.0_dp, 10000.0_dp, 13000.0_dp, &
        16000.0_dp, 19000.0_dp, 22000.0_dp, 25000.0_dp, 30000.0_dp, 40000.0_dp, 60000.0_dp, 100000.0_dp]
    real(dp), parameter :: PROFW(11) = [0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp, 0.9_dp, 1.0_dp, 1.1_dp, 1.25_dp, &
                                        1.5_dp, 1.75_dp, 2.0_dp]
    real(dp), parameter :: PROF(11) = [0.80_dp, 0.85_dp, 0.90_dp, 0.95_dp, 0.98_dp, 1.0_dp, 1.02_dp, 1.05_dp, &
                                       1.10_dp, 1.15_dp, 1.20_dp]

    epoch0 = epoch_jd
    print '(/A)', ' BLIND FIT: Tycho''s readings only.  No JPL orbit, planet, time, sign, limb or'
    print '(A)',  ' selection; JPL is printed only for comparison at the end.'
    print '(A)',  ' Elements 1-7 are those of the Earth-Moon barycentre about the Sun.'
    ! JPL's orbits in the same parametrization, for the comparison only
    ! (not in the procedure, which has no JPL orbits)
    xt = 0.0_dp
    if (.not. alien) then
    block
      real(dp) :: pos(3,NBODY), vel(3,NBODY), gm(NBODY)
      real(dp), allocatable :: p1(:), p2(:)
      emb_param = .false.
      call initial_state(x_true, pos, vel, gm)
      emb_param = .true.
      xt = x_true
      call state_to_elements(pos, vel, gm, xt)
      allocate(p1(nobs), p2(nobs))
      emb_param = .false.; call predict_all(x_true, p1)
      emb_param = .true.;  call predict_all(xt, p2)
      print '(A,ES10.2,A)', ' (check of the barycentre elements: largest change of a prediction ', &
          maxval(abs(p1 - p2))*3600.0_dp, '")'
    end block
    end if

    ! Diagnostic only (argument 9 = fromjpl): skip the build-up and start the
    ! final fit from JPL's orbits, to see whether the blind fit found the
    ! lowest chi2 of these readings
    call get_arg(9, arg)
    if (trim(arg) == 'fromjpl') then
      print '(A)', ' DIAGNOSTIC: final fit started from JPL''s orbits (optimizer check, not a result)'
      xb = xt; xk = xt
      go to 500
    end if
    ! argument 9 = fromfile: start the final fit from the vector in
    ! blind_start.dat (a BLIND line of an earlier blind fit, all of Tycho's
    ! own making), to refit subsets or with more parameters free
    if (trim(arg) == 'fromfile') then
      block
        integer :: un, ios
        character(len=4096) :: ln
        character(len=16) :: tag
        open(newunit=un, file='blind_start.dat', status='old', action='read')
        read(un, '(A)') ln
        close(un)
        xb = 0.0_dp; xb(18:19) = 1.0_dp; xb(22:23) = 1.0_dp; xb(25) = 1.0_dp; xb(27:28) = 1.0_dp; call one_scales(xb)
        read(ln, *, iostat=ios) tag, xb
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:595)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:345)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:95)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:91)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:89)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:87)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:42)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:36)
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:30)
        if (ios /= 0) then
          read(ln, *, iostat=ios) tag, xb(1:25)
          ! a file without 1594: start 1594 from 1595's constants
          xb(26:29) = [xb(21), xb(22), xb(23), xb(24)]; xb(30) = xb(17)
        end if
        if (ios /= 0) read(ln, *, iostat=ios) tag, xb(1:24)
        if (ios /= 0) then
          xb(21:24) = [xb(15), xb(18), xb(19), xb(20)]
          read(ln, *, iostat=ios) tag, xb(1:20)
          if (ios /= 0) then
            xb(18:20) = [1.0_dp, 1.0_dp, 0.0_dp]
            read(ln, *) tag, xb(1:17)
          end if
          xb(21:24) = [xb(15), xb(18), xb(19), xb(20)]
        end if
      end block
      print '(A)', ' Final fit started from blind_start.dat (an earlier blind solution)'
      xk = xb
      go to 500
    end if

    ! ── 1. the Sun ──
    allocate(dec(nall)); nd = 0
    do k1 = 1, nall
      ! the Sun's meridian altitudes: half their yearly range is the obliquity
      if (obs_all(k1)%body == 1 .and. trim(obs_all(k1)%qty) == 'ALTM') then
        nd = nd + 1; dec(nd) = obs_all(k1)%value
      end if
    end do
    call sort_vals(dec(1:nd))
    ob = 0.5_dp*(dec(max(1, nint(0.98_dp*nd))) - dec(max(1, nint(0.02_dp*nd))))
    print '(A,I5,A,F8.3,A)', '   Sun meridian altitudes:', nd, ';  half their range (2-98 %):', ob, ' deg'
    xb = 0.0_dp; xb(18:19) = 1.0_dp; xb(22:23) = 1.0_dp; xb(25) = 1.0_dp; xb(27:28) = 1.0_dp; call one_scales(xb)
    xb(1:7) = [0.02_dp, ob*DEG2RAD, 0.0_dp, 0.0_dp, 0.0_dp, TAU/365.25_dp, 1150.0_dp*RE]
    xb(8:14) = [0.05_dp, ob*DEG2RAD, 0.0_dp, 0.0_dp, 0.0_dp, TAU/27.3217_dp, 60.0_dp*RE]
    call set_window(-1.0d30, 1.0d30, 1)
    print '(A,I6,A)', '   1. Sun:', nobs, ' Sun readings.  (a) a circle from 12 starting longitudes'
    act = .false.; act(2:3) = .true.; act(5:6) = .true.; act(15:17) = .true.
    best = huge(1.0_dp)
    do kk = 0, 11
      xk = xb; xk(1) = 1.0d-7; xk(5) = kk*30.0_dp*DEG2RAD
      call run_fit(xk, act, 5.0_dp, score, nw)
      print '(A,I4,A,F10.3,A,I6,A,F9.4,A,F9.3)', '      start M0 =', kk*30, ' deg:  median |resid| ', score, &
          "'  used", nw, '   i', xk(2)*RAD2DEG, '   longitude', modulo(xk(3) + xk(4) + xk(5), TAU)*RAD2DEG
      if (score < best) then
        best = score; xbest = xk
      end if
    end do
    xb = xbest
    print '(A)', '      (b) the ellipse: eccentricity released from 4 perihelion directions'
    act(1) = .true.; act(4) = .true.
    best = huge(1.0_dp)
    do kk = 0, 3
      xk = xb; xk(1) = 0.02_dp; xk(4) = kk*90.0_dp*DEG2RAD; xk(5) = xb(4) + xb(5) - xk(4)
      call run_fit(xk, act, 5.0_dp, score, nw)
      print '(A,I4,A,F10.3,A,I6,A,F8.5,A,F9.4)', '      start w =', kk*90, ' deg:  median |resid| ', score, &
          "'  used", nw, '   e', xk(1), '   i', xk(2)*RAD2DEG
      if (score < best) then
        best = score; xbest = xk
      end if
    end do
    xb = xbest
    print '(A,F8.3,A,F10.5,A,F9.5,A,F9.4,A)', '   Sun fitted: median |resid| ', best, "'   e ", xb(1), &
        '   i ', xb(2)*RAD2DEG, ' deg   year ', TAU/xb(6), ' d'

    ! ── 2. the Moon: Kepler ellipse to the densest ten days of 1586 ──
    nbest = 0; tc = 0.0_dp
    do k1 = 1, nall
      if (obs_all(k1)%body /= 2 .or. obs_all(k1)%year /= 1586) cycle
      nw = count([(obs_all(k2)%body == 2 .and. obs_all(k2)%jd >= obs_all(k1)%jd .and. &
                   obs_all(k2)%jd <= obs_all(k1)%jd + 10.0_dp, k2 = 1, nall)])
      if (nw > nbest) then
        nbest = nw; tc = obs_all(k1)%jd + 5.0_dp
      end if
    end do
    tc = anint(tc)
    call transfer_epoch(xb, tc)
    inc = xb(2); Om = xb(3)
    nE = [sin(inc)*sin(Om), -sin(inc)*cos(Om), cos(inc)]
    Xv = [cos(Om), sin(Om), 0.0_dp]
    Yv = [nE(2)*Xv(3) - nE(3)*Xv(2), nE(3)*Xv(1) - nE(1)*Xv(3), nE(1)*Xv(2) - nE(2)*Xv(1)]
    kepler_moon = .true.
    call set_window(tc - 5.0_dp, tc + 5.0_dp, 2)
    print '(/A,F12.2,A,I5,A)', '   2. Moon: Kepler ellipse, 10 days around JD', tc, ' (', nobs, &
        ' Moon readings): (a) a circle from 72 starts'
    act = .false.; act(9:10) = .true.; act(12) = .true.
    best = huge(1.0_dp)
    do k1 = 0, 5
      do k2 = 0, 11
        nm = cos(5.0_dp*DEG2RAD)*nE + sin(5.0_dp*DEG2RAD)*(sin(k1*60.0_dp*DEG2RAD)*Xv - cos(k1*60.0_dp*DEG2RAD)*Yv)
        xk = xb
        xk(8:14) = [1.0d-7, acos(nm(3)), modulo(atan2(nm(1), -nm(2)), TAU), 0.0_dp, k2*30.0_dp*DEG2RAD, &
                    TAU/27.3217_dp, 60.0_dp*RE]
        call run_fit(xk, act, 5.0_dp, score, nw)
        if (score < best) then
          best = score; xbest = xk
        end if
      end do
    end do
    xb = xbest
    print '(A,F9.3,A)', '      best circle: median |resid| ', best, "'"
    print '(A)', '      (b) the ellipse from 4 perigee directions'
    act(8) = .true.; act(11) = .true.
    best = huge(1.0_dp)
    do kk = 0, 3
      xk = xb; xk(8) = 0.05_dp; xk(11) = kk*90.0_dp*DEG2RAD; xk(12) = xb(11) + xb(12) - xk(11)
      call run_fit(xk, act, 5.0_dp, score, nw)
      if (score < best) then
        best = score; xbest = xk
      end if
    end do
    xb = xbest
    print '(A,F9.3,A,F8.4,A,F8.3,A)', '      best start: median |resid| ', best, "'   e ", xb(8), &
        '   i(equator) ', xb(9)*RAD2DEG, ' deg'

    ! ── 3. the Moon with the Sun's pull, widening the window over 1586 ──
    kepler_moon = .false.
    print '(/A)', '   3. Moon, N-body, window widened over 1586'
    do kk = 1, size(HALVES)
      half = HALVES(kk)
      act = .false.; act(8:12) = .true.
      if (half >= 20.0_dp) act(13) = .true.
      if (half >= 40.0_dp) act(14) = .true.
      call set_window(tc - half, tc + half, 2)
      call run_fit(xb, act, 5.0_dp, score, nw)
      print '(A,F6.0,A,I5,A,F9.3,A,I5,A,F9.5,A,F10.0,A)', '      +-', half, ' d:', nobs, &
          ' readings  median |resid| ', score, "'  used", nw, '  month ', TAU/xb(13), ' d  a_M ', xb(14), ' km'
    end do

    ! ── 4. the Earth-Sun distance, from all of 1586 ──
    call set_window(-1.0d30, epoch0, 0)
    act = .true.; act(7) = .false.
    print '(/A,I5,A)', '   4. Earth-Sun distance scanned on all of 1586 (', nobs, &
        ' readings; everything else refitted)'
    best = huge(1.0_dp)
    do kk = 1, size(AE_RE)
      xb(7) = AE_RE(kk)*RE
      call run_fit(xb, act, 5.0_dp, score, nw)
      sc_scan(kk) = score
      print '(A,F9.0,A,F7.2,A,F9.4,A,I6)', '      a_E =', AE_RE(kk), ' Earth radii (', xb(7)/1.0d6, &
          ' M km):  median |resid| ', score, "'  used", nw
      if (score < best) then
        best = score; xbest = xb
      end if
    end do
    xb = xbest
    act = .true.
    call run_fit(xb, act, 5.0_dp, score, nw)
    print '(A,F8.2,A,F9.4,A)', '      free: a_E = ', xb(7)/1.0d6, ' M km   median |resid| ', score, "'"
    xk = xb                                   ! the 1586 solution (start of the Kalman filter below)

    ! ── 5. 1586 and 1595 together, epoch 1590-07-01 ──
500 continue
    call transfer_epoch(xb, epoch0)
    call transfer_epoch(xk, epoch0)
    call set_window(-1.0d30, 1.0d30, 0)
    ! argument 11: what else is fitted here (comma list):
    !   moonr  the Moon's radius (otherwise 1737.4 km)
    !   refr   the strength of refraction (scale of the Saemundsson formula)
    !   lat    Uraniborg's latitude (otherwise 55 54' 27")
    call get_arg(11, arg)
    param_on(18) = index(arg, 'moonr') > 0
    param_on(19) = index(arg, 'refr') > 0
    param_on(20) = index(arg, 'lat') > 0
    param_on(25) = param_on(19) .and. tycho_model   ! Cassini's height with the refraction
    !   armil    the tilts of the armillary sphere's axis, each year
    !   ewinst   quadrant tilt (east-west) and azimuth zero, each year
    if (index(arg, 'ewinst') > 0) then
      use_ew = .true.; param_on(37:42) = .true.
      do k1 = 4, NY
        param_on(ip_year(YEARS(k1), 7):ip_year(YEARS(k1), 8)) = .true.
      end do
    end if
    if (index(arg, 'phase') > 0) use_phase = .true.
    if (index(arg, 'irr') > 0) then
      use_irr = .true.; param_on(88) = .true.
    end if
    if (index(arg, 'kew') > 0) then
      use_kew = .true.; param_on(89) = .true.
    end if
    if (index(arg, 'dfar') > 0) then
      use_dfar = .true.; param_on(90) = .true.
    end if
    if (index(arg, 'decns') > 0) dec_ns = .true.
    if (index(arg, 'tgroup') > 0) call set_time_groups()
    auto_sel = index(arg, 'auto') > 0
    if (index(arg, 'rdec') > 0) then
      reff(1) = .true.; use_mdec = .true.; param_on(96:95 + NMON) = .true.
    end if
    if (index(arg, 'ralt') > 0) then
      reff(2) = .true.; use_malt = .true.; param_on(96 + NMON:95 + 2*NMON) = .true.
    end if
    if (index(arg, 'rdist') > 0) then
      reff(3) = .true.; use_mdist = .true.; param_on(96 + 2*NMON:95 + 3*NMON) = .true.
    end if
    if (index(arg, 'rfar') > 0) then
      reff(4) = .true.; use_mfar = .true.; param_on(96 + 3*NMON:IST0) = .true.
      use_dfar = .true.; param_on(90) = .true.
    end if
    if (index(arg, 'rstar') > 0) then
      reff(5) = .true.; use_rstar = .true.; param_on(IST0 + 1:NP) = .true.
    end if
    if (index(arg, 'mdec') > 0) then
      use_mdec = .true.; param_on(96:95 + NMON) = .true.
      param_on([16, 17, 30]) = .false.
      do k1 = 4, NY
        param_on(ip_year(YEARS(k1), 1)) = .false.
      end do
    end if
    if (index(arg, 'mdist') > 0) then
      use_mdist = .true.; param_on(96 + 2*NMON:95 + 3*NMON) = .true.
    end if
    if (index(arg, 'mfar') > 0) then
      use_mfar = .true.; param_on(96 + 3*NMON:IST0) = .true.
    end if
    if (index(arg, 'malt') > 0) then
      use_malt = .true.; param_on(96 + NMON:95 + 2*NMON) = .true.
      param_on([15, 21, 26]) = .false.
      do k1 = 4, NY
        param_on(ip_year(YEARS(k1), 0)) = .false.
      end do
    end if
    if (index(arg, 'limbew') > 0) then
      use_limbo = .true.; param_on(92:93) = .true.
    end if
    if (index(arg, 'limbo') > 0) then
      use_limbo = .true.; param_on(92:95) = .true.
    end if
    if (index(arg, 'dnear') > 0) then
      use_dnear = .true.; param_on(91) = .true.
    end if
    if (index(arg, 'armil') > 0) then
      use_armil = .true.; param_on(31:36) = .true.
      do k1 = 4, NY
        param_on(ip_year(YEARS(k1), 5):ip_year(YEARS(k1), 6)) = .true.
      end do
    end if
    !   huber / tukey / noclip  weighting of the final fit (default: reject
    !            beyond 4 sigma)
    if (index(arg, 'huber') > 0) robust_kind = 1
    if (index(arg, 'tukey') > 0) robust_kind = 2
    if (index(arg, 'noclip') > 0) robust_kind = 3
    if (tycho_model) print '(A)', '   TYCHO MODEL: Tycho''s frame and stars, no aberration/J2, massless Moon, sphere, Cassini'
    if (robust_kind > 0) print '(A,I2)', '   weighting of the final fit (1 Huber, 2 Tukey, 3 none):', robust_kind
    !   peryear  the 1595 readings get their own altitude offset and their
    !            own values of whichever of the three above are fitted
    peryear = index(arg, 'peryear') > 0
    if (peryear) then
      param_on(21) = .true.; param_on(22:24) = param_on(18:20)
      param_on(26) = .true.; param_on(27:29) = param_on(18:20)
      do k1 = 4, NY
        param_on(ip_year(YEARS(k1), 0)) = .true.
        param_on(ip_year(YEARS(k1), 2):ip_year(YEARS(k1), 4)) = param_on(18:20)
      end do
      print '(A)', '   per-year instrument parameters for every year'
    end if
    if (any(param_on(18:20))) print '(A,3L2)', '   also fitted (Moon radius, refraction, latitude):', param_on(18:20)
    ! argument 12: subset of the readings for the final fit
    !   y1586 / y1595  one year;  moonstar  the Moon only where the stars
    !   give the time (all Sun readings kept)
    call get_arg(12, arg)
    call quarter_table(xb)
    if (len_trim(arg) > 0 .and. trim(arg) /= 'all') call apply_subset(trim(arg))
    print '(/A,I6,A)', '   5. All readings (', nobs, '), all years together'
    block
      real(dp), allocatable :: v(:), pr(:)
      logical, allocatable :: u(:)
      real(dp) :: r1586, r1595
      allocate(v(nobs), pr(nobs), u(nobs))
      v = obs(1:nobs)%value
      call predict_all(xb, pr)
      call med_resid_year(v, pr, r1586, r1595)
      print '(A,F9.3,A,F9.3,A)', '      1586 solution predicts 1586 to ', r1586, "' and 1595 to ", r1595, &
          "' (median |resid|)"
      act = .true.; act(7) = .false.
      call run_fit(xb, act, 5.0_dp, score, nw)
      print '(A,F9.4,A,I6)', '      a_E fixed:  median |resid| ', score, "'  used", nw
      act = .true.
      call fit_lm(xb, v, act, 4.0_dp, cv, u, rms, .true.)
      print '(A,F14.1,A,F8.3,A,F12.2)', '      without clock correction of the daytime sessions: a_E ', xb(7), &
          ' km   rms ', rms, "'   chi2 ", last_chi2
      ! the clocks of the daytime sessions from the model's own Sun
      ! (argument 10 = noclock leaves Tycho's clock times as written)
      call get_arg(10, arg)
      if (trim(arg) /= 'noclock') then
        do kk = 1, 4
          call update_clocks(xb, nw, score)
          call fit_lm(xb, v, act, 4.0_dp, cv, u, rms, .false.)
          print '(A,I2,A,I5,A,F8.2,A,F14.1,A,F8.3,A,F12.2)', '      session clocks, pass', kk, ':', nw, &
              ' readings corrected, median |correction| ', score, ' min;  a_E ', xb(7), ' km  rms ', rms, &
              "'  chi2 ", last_chi2
        end do
      end if
      ! "auto": the corrections of the instruments per campaign, each taken
      ! only if the readings ask for it (Bayesian information criterion,
      ! with the weights and readings of the fit so far held fixed)
      if (auto_sel) call select_corrections(xb, v, act, cv, u, rms)
      if (any(reff)) then
        do kk = 1, 6
          call update_tau(xb, cv)
          call fit_lm(xb, v, act, 4.0_dp, cv, u, rms, .false.)
          print '(A,I2,A,5F7.3,A,F14.1)', '      random effects, pass', kk, ':  tau (dec, alt, dist, far, star) ', tau_re, &
              "'   a_E ", xb(7)
        end do
      end if
      do kk = 1, NP
        sgv(kk) = sqrt(max(cv(kk,kk), 0.0_dp))
      end do
      print '(A,F8.3,A,I6,A,I6,A,F12.2)', '   rms residual ', rms, "'   readings used ", count(u), ' of ', nobs, &
          '   chi2 ', last_chi2
      call report(" BLIND FIT (nothing from JPL; JPL column for comparison only)", xb, sgv, xt)
      call print_sigmas()
      call write_blind_residuals(xb, u, v)
      write(*, '(A,755ES20.11)') ' BLIND ', xb
      write(*, '(A,755ES20.11)') ' BLINDSIG ', sgv
      block
        integer :: un, g
        open(newunit=un, file='sigmas_out.dat', status='replace', action='write')
        do g = 1, NG
          write(un, '(ES16.8,2X,A)') sigma_g(g), group_name(g)
        end do
        close(un)
      end block

      ! the procedure's uncertainty: leave out one year at a time and refit
      ! (jackknife over the years); then stop
      if (alien) then
        block
          type(reading), allocatable :: keep_obs(:)
          real(dp), allocatable :: vj(:)
          logical, allocatable :: uj(:)
          real(dp) :: aj(NY), xj(NP), rj, mean, var
          real(dp), allocatable :: cvj(:,:)
          integer :: nkeep, y, nj, j2, m2
          allocate(cvj(NP, NP))
          nkeep = nobs; allocate(keep_obs(nkeep)); keep_obs = obs(1:nkeep)
          nj = 0
          print '(/A)', '   jackknife: each year left out in turn'
          do y = 1, NY
            if (count(keep_obs%year == YEARS(y)) == 0) cycle
            m2 = 0
            do j2 = 1, nkeep
              if (keep_obs(j2)%year /= YEARS(y)) then
                m2 = m2 + 1; obs(m2) = keep_obs(j2)
              end if
            end do
            nobs = m2
            call sort_order()
            allocate(vj(nobs), uj(nobs)); vj = obs(1:nobs)%value
            xj = xb
            call fit_lm(xj, vj, act, 4.0_dp, cvj, uj, rj, .false.)
            if (any(reff)) then
              do kk = 1, 3
                call update_tau(xj, cvj)
                call fit_lm(xj, vj, act, 4.0_dp, cvj, uj, rj, .false.)
              end do
            end if
            nj = nj + 1; aj(nj) = xj(7)
            print '(A,I5,A,F14.1,A)', '      without ', YEARS(y), ':  a_E ', xj(7), ' km'
            deallocate(vj, uj)
          end do
          obs(1:nkeep) = keep_obs; nobs = nkeep
          call sort_order()
          mean = sum(aj(1:nj))/nj
          var = (nj - 1.0_dp)/nj*sum((aj(1:nj) - mean)**2)
          print '(/A,F14.1,A,F12.1,A,F12.1,A)', ' RESULT  Earth-Sun distance ', xb(7), ' km  +- ', sqrt(var), &
              ' (jackknife over years; formal ', sgv(7), ')'
          print '(A,F14.4,A,F10.4,A)', '         in Earth radii ', xb(7)/RE_KM, '   Earth-Moon ', xb(14)/RE_KM, ' Earth radii'
        end block
        return
      end if

      ! DIAGNOSTIC (environment TYCHO_PULL): each reading's share of the
      ! pull on a_E -- w J r / sum(w J^2), with J the derivative in a_E
      ! after the other fitted parameters are projected out -- written to
      ! pull.dat.  TYCHO_PULL=<km> first moves a_E there and refits the
      ! rest (weights and readings of the blind fit kept); then the shares
      ! add up to the step from there toward the blind solution.
      block
        character(len=32) :: pv
        integer :: st, i2, k3, na2, idx2(NP), un2
        real(dp), allocatable :: Jm(:,:), rr0(:), jp(:)
        real(dp), allocatable :: Noo(:,:), bo(:), so(:)
        real(dp) :: den, w2, aeval, sgf2(NG)
        logical :: ok2
        call get_environment_variable('TYCHO_PULL', pv, status=st)
        if (st == 0) then
          allocate(Jm(nobs, NP), rr0(nobs), jp(nobs))
          aeval = 0.0_dp
          if (len_trim(pv) > 0) read(pv, *) aeval
          if (aeval > 0.0_dp) then
            sgf2 = sigma_g
            freeze_sigma = .true.; freeze_use = .true.
            xb(7) = aeval; act(7) = .false.
            call fit_lm(xb, v, act, 4.0_dp, cv, u, rr, .false.)
            act(7) = .true.; sigma_g = sgf2
            freeze_sigma = .false.; freeze_use = .false.
            call write_blind_residuals(xb, u, v)     ! the residuals at the given a_E
          end if
          call jacobian(xb, v, Jm, rr0, act)
          na2 = 0
          do k3 = 1, NP
            if (k3 /= 7 .and. act(k3) .and. param_on(k3)) then
              na2 = na2 + 1; idx2(na2) = k3
            end if
          end do
          allocate(Noo(na2, na2), bo(na2), so(na2))
          Noo = 0.0_dp; bo = 0.0_dp
          do i2 = 1, nobs
            if (.not. u(i2)) cycle
            w2 = 1.0_dp/sigma_g(obs(i2)%group)**2
            do k3 = 1, na2
              Noo(k3, :) = Noo(k3, :) + w2*Jm(i2, idx2(k3))*Jm(i2, idx2(1:na2))
            end do
            bo = bo + w2*Jm(i2, idx2(1:na2))*Jm(i2, 7)
          end do
          do k3 = 1, na2
            if (Noo(k3, k3) <= 0.0_dp) then
              Noo(k3, k3) = 1.0_dp; bo(k3) = 0.0_dp
            end if
          end do
          call chol_solve(Noo, bo, so, ok2)
          if (.not. ok2) print '(A)', '   PULL: projection singular'
          jp = Jm(:, 7) - matmul(Jm(:, idx2(1:na2)), so)
          den = 0.0_dp
          do i2 = 1, nobs
            if (u(i2)) den = den + jp(i2)**2/sigma_g(obs(i2)%group)**2
          end do
          open(newunit=un2, file='pull.dat', status='replace', action='write')
          write(un2, '(A)') '# id year body qty point refkind tkind used elong alt ha resid_arcmin sigma jperp share_km'
          do i2 = 1, nobs
            w2 = merge(1.0_dp/sigma_g(obs(i2)%group)**2, 0.0_dp, u(i2))
            write(un2, '(A12,I5,I2,1X,A6,1X,A1,1X,A4,1X,A1,L2,3F9.3,F10.3,F9.3,ES12.4,ES13.5)') obs(i2)%id, &
                obs(i2)%year, obs(i2)%body, obs(i2)%qty, obs(i2)%point, obs(i2)%refkind, obs(i2)%tkind, u(i2), &
                elong(i2), alt_read(i2), ha_read(i2), rr0(i2), sigma_g(obs(i2)%group), jp(i2), &
                w2*jp(i2)*rr0(i2)/den*pscale(xb, 7)
          end do
          close(un2)
          print '(A,F14.1,A,F14.1,A,ES10.3)', '   PULL at a_E ', xb(7), ' km: sum of shares ', &
              sum(merge(jp*rr0/sigma_g(obs(1:nobs)%group)**2, 0.0_dp, u))/den*pscale(xb, 7), &
              ' km;  1-sigma ', pscale(xb, 7)/sqrt(den)
          return
        end if
      end block

      ! χ² profile in a_E around the blind solution (weights and readings fixed)
      block
        real(dp) :: sgf(NG), chi0, xp(NP)
        logical :: actp(NP)
        sgf = sigma_g; chi0 = last_chi2
        freeze_sigma = .true.; freeze_use = .true.
        actp = .true.; actp(7) = .false.
        print '(/A)', '   chi2 against a_E (all else refitted; weights of the blind fit)'
        print '(A)', '     a_E/fit      a_E (km)        chi2     dchi2    rms(arcmin)'
        do kk = 1, 11
          sigma_g = sgf
          xp = xb; xp(7) = xb(7)*PROF(kk)
          call get_arg(12, arg)
          if (len_trim(arg) > 0 .and. trim(arg) /= 'all') xp(7) = xb(7)*PROFW(kk)
          call fit_lm(xp, v, actp, 4.0_dp, cv, u, rr, .false.)
          print '(F9.3,F16.1,F12.3,F10.3,F10.3,I7)', xp(7)/xb(7), xp(7), last_chi2, last_chi2 - chi0, rr, count(u)
        end do
        sigma_g = sgf
        freeze_sigma = .false.; freeze_use = .false.

        ! Kalman filter started from the 1586-only solution, same weights
        call get_arg(9, arg)
        if (trim(arg) == 'fromfile') return
        print '(/A)', ' Iterated EKF (F = I, Q = 0) on all readings, started from the 1586 solution'
        call fit_lm(xb, v, act, 4.0_dp, cv, u, rr, .false.)   ! the readings and weights to use
        xp = xk
        sigma_prior = [0.05_dp, 0.05_dp, 0.05_dp, 0.5_dp, 0.5_dp, 1.0d-4, 1.0d7, &
                       0.05_dp, 0.05_dp, 0.05_dp, 0.5_dp, 0.5_dp, 1.0d-2, 4.0d4, 5.0_dp, 5.0_dp, 5.0_dp, &
                   0.1_dp, 0.3_dp, 5.0_dp, 5.0_dp, 0.1_dp, 0.3_dp, 5.0_dp, 0.3_dp, &
                   5.0_dp, 0.1_dp, 0.3_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, &
                   5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, &
                   (5.0_dp, 5.0_dp, 0.1_dp, 0.3_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, ii_ = 1, 5), &
                   1.0_dp, 0.1_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, &
                   (5.0_dp, ii_ = 1, 4*NMON), (1.0_dp, ii_ = 1, 2*NSTAR)]
        call fit_iekf(xp, v, act, u, sigma_prior, 30, P)
        call predict_all(xp, pr)
        rr = 0.0_dp
        do kk = 1, nobs
          if (u(kk)) rr = rr + (dresid(obs(kk), v(kk), pr(kk))/sigma_g(obs(kk)%group))**2
        end do
        do kk = 1, NP
          sgv(kk) = sqrt(max(P(kk,kk), 0.0_dp))
        end do
        print '(A,F12.2,A,F12.2)', '   chi2 ', rr, '   (batch least squares ', last_chi2, ')'
        call report(' Iterated EKF from the 1586 solution (filter 1-sigma)', xp, sgv, xt)
        write(*, '(A,755ES20.11)') ' BLINDEKF ', xp
      end block
    end block
  end subroutine

  ! How many Moon readings fall near the quarters, and how each is timed
  subroutine quarter_table(xr)
    real(dp), intent(in) :: xr(NP)
    real(dp) :: pr(nobs)
    integer :: cnt(NY, 3, 3), y, ph, tk, j
    call predict_all(xr, pr)
    cnt = 0
    do j = 1, nobs
      if (obs(j)%body /= 2) cycle
      y = yidx(obs(j)%year)
      if (elong(j) >= 45.0_dp .and. elong(j) <= 135.0_dp) then
        ph = 1
      else if (elong(j) >= 225.0_dp .and. elong(j) <= 315.0_dp) then
        ph = 2
      else
        ph = 3
      end if
      tk = index('SAM', obs(j)%tkind)
      if (tk == 0) tk = 1
      cnt(y, ph, tk) = cnt(y, ph, tk) + 1
    end do
    print '(/A)', '   Moon readings by phase (elongation 45-135 deg: first quarter, 225-315: last) and timing'
    print '(A)',  '     year  phase           star clock  local time  transit'
    do y = 1, NY
      if (sum(cnt(y, :, :)) == 0) cycle
      print '(5X,I4,2X,A,3I12)', YEARS(y), 'first quarter', cnt(y, 1, :)
      print '(5X,I4,2X,A,3I12)', YEARS(y), 'last quarter ', cnt(y, 2, :)
      print '(5X,I4,2X,A,3I12)', YEARS(y), 'other        ', cnt(y, 3, :)
    end do
  end subroutine

  subroutine apply_subset(name)
    ! name: tokens joined by '+': y1586 / y1594 / y1595 (keep that year),
    ! no1594, moonstar (the Moon only where the stars give the time),
    ! old1586 (1586 without the readings of printed pp. 30-47)
    character(len=*), intent(in) :: name
    logical :: keep(nobs)
    integer :: j, m
    character(len=6) :: ytok
    do j = 1, nobs
      keep(j) = .true.
      ! yYYYY (one or more: keep those years), noYYYY (drop that year)
      if (index(name, 'y15') > 0) then
        write(ytok, '(A,I4)') 'y', obs(j)%year
        keep(j) = keep(j) .and. index(name, trim(ytok)) > 0
      end if
      write(ytok, '(A,I4)') 'no', obs(j)%year
      if (index(name, trim(ytok)) > 0) keep(j) = .false.
      if (index(name, 'moonstar') > 0) keep(j) = keep(j) .and. (obs(j)%body == 1 .or. obs(j)%tkind == 'S')
      ! moonsm: the Moon timed by stars or read at the meridian (no clock needed)
      if (index(name, 'moonsm') > 0) keep(j) = keep(j) .and. (obs(j)%body == 1 .or. obs(j)%tkind == 'S' .or. &
          obs(j)%tkind == 'M')
      if (index(name, 'starsess') > 0 .and. obs(j)%body == 2) keep(j) = keep(j) .and. obs(j)%tkind == 'A' .and. &
          obs(j)%sess > 0 .and. any(obs(1:nobs)%body == 3 .and. obs(1:nobs)%sess == obs(j)%sess)
      ! notime: the Moon only by readings that hardly depend on the clock
      ! (meridian altitudes, declinations, diameters; +dist: and distances)
      if (index(name, 'notime') > 0 .and. obs(j)%body == 2) keep(j) = keep(j) .and. &
          (trim(obs(j)%qty) == 'ALTM' .or. trim(obs(j)%qty) == 'DEC' .or. trim(obs(j)%qty) == 'DIAM' .or. &
           (index(name, 'dist') > 0 .and. trim(obs(j)%qty) == 'DIST'))
      ! alt15: the body read (and the Sun, for readings against it) at least
      ! 15 deg high; nosunref: no Moon readings against the Sun
      if (index(name, 'alt15') > 0) keep(j) = keep(j) .and. (alt_read(j) >= 15.0_dp .or. obs(j)%body == 3)
      if (index(name, 'nosunref') > 0) keep(j) = keep(j) .and. trim(obs(j)%refkind) /= 'SUN'
      ! xDEC, xDIST, xRADIFF, xHA, xALT (with ALTM), xAZ: leave that kind of
      ! Moon reading out
      if (obs(j)%body == 2) then
        if (index(name, 'xDEC') > 0 .and. trim(obs(j)%qty) == 'DEC') keep(j) = .false.
        if (index(name, 'xDIST') > 0 .and. trim(obs(j)%qty) == 'DIST') keep(j) = .false.
        if (index(name, 'xRADIFF') > 0 .and. trim(obs(j)%qty) == 'RADIFF') keep(j) = .false.
        if (index(name, 'xHA') > 0 .and. trim(obs(j)%qty) == 'HA') keep(j) = .false.
        if (index(name, 'xALT') > 0 .and. obs(j)%qty(1:3) == 'ALT') keep(j) = .false.
        if (index(name, 'xAZ') > 0 .and. trim(obs(j)%qty) == 'AZ') keep(j) = .false.
      end if
      ! noupdec: no declinations of the Moon's upper limb; nowdist: no
      ! distances from its west limb
      if (obs(j)%body == 2 .and. index(name, 'noupdec') > 0 .and. trim(obs(j)%qty) == 'DEC' .and. &
          obs(j)%point == 'U') keep(j) = .false.
      if (obs(j)%body == 2 .and. index(name, 'nowdist') > 0 .and. trim(obs(j)%qty) == 'DIST' .and. &
          obs(j)%point == 'W') keep(j) = .false.
      ! nocres: no Moon readings at a crescent (within 70 deg of the Sun)
      if (index(name, 'nocres') > 0 .and. obs(j)%body == 2) keep(j) = keep(j) .and. &
          elong(j) > 70.0_dp .and. elong(j) < 290.0_dp
      if (index(name, 'old1586') > 0) keep(j) = keep(j) .and. &
          .not. (obs(j)%id(1:5) == 'm86c:' .or. obs(j)%id(1:5) == 'm86d:')
    end do
    m = 0
    do j = 1, nobs
      if (keep(j)) then
        m = m + 1; obs(m) = obs(j)
      end if
    end do
    nobs = m
    call sort_order()
    print '(A,A,A,I6,A)', '   subset ', name, ': ', nobs, ' readings'
  end subroutine

  subroutine run_fit(xr, act, clp, score, nused)
    ! fit on the current window; score = median |residual| over all its readings
    real(dp), intent(inout) :: xr(NP)
    logical, intent(in) :: act(NP)
    real(dp), intent(in) :: clp
    real(dp), intent(out) :: score
    integer, intent(out) :: nused
    real(dp), allocatable :: v(:), pr(:), ar(:)
    logical, allocatable :: u(:)
    real(dp) :: cv(NP,NP), rr
    integer :: j
    allocate(v(nobs), pr(nobs), u(nobs), ar(nobs))
    v = obs(1:nobs)%value
    call fit_lm(xr, v, act, clp, cv, u, rr, .false.)
    call predict_all(xr, pr)
    nused = 0
    do j = 1, nobs
      if (obs(j)%body == 3) cycle
      nused = nused + 1
      ar(nused) = abs(dresid(obs(j), v(j), pr(j)))
    end do
    score = median(ar(1:nused))
    nused = count(u)
  end subroutine

  subroutine med_resid_year(v, pr, r86, r95)
    real(dp), intent(in) :: v(nobs), pr(nobs)
    real(dp), intent(out) :: r86, r95
    real(dp) :: a(nobs), b(nobs)
    integer :: j, na, nb
    na = 0; nb = 0
    do j = 1, nobs
      if (obs(j)%body == 3) cycle
      if (obs(j)%year == 1586) then
        na = na + 1; a(na) = abs(dresid(obs(j), v(j), pr(j)))
      else
        nb = nb + 1; b(nb) = abs(dresid(obs(j), v(j), pr(j)))
      end if
    end do
    r86 = median(a(1:na)); r95 = median(b(1:nb))
  end subroutine

  subroutine sort_vals(a)
    real(dp), intent(inout) :: a(:)
    real(dp) :: t
    integer :: i2, j2
    do i2 = 2, size(a)
      t = a(i2); j2 = i2 - 1
      do while (j2 >= 1)
        if (a(j2) <= t) exit
        a(j2+1) = a(j2); j2 = j2 - 1
      end do
      a(j2+1) = t
    end do
  end subroutine

  subroutine write_blind_residuals(xf, u, v)
    real(dp), intent(in) :: xf(NP), v(nobs)
    logical, intent(in) :: u(nobs)
    real(dp) :: pr(nobs)
    integer :: un, j
    call predict_all(xf, pr)
    open(newunit=un, file='tycho_blind_residuals.dat', status='replace')
    write(un, '(A)') '# id year body qty tkind jd_tt used value_deg fit_deg jpl_deg resid_fit_arcmin ' // &
        'resid_jpl_arcmin jd_used dclock_min sess'
    do j = 1, nobs
      write(un, '(A12,I5,1X,A1,1X,A6,1X,A1,F18.8,L3,3F14.7,2F10.3,F18.8,F10.3,I6)') obs(j)%id, obs(j)%year, &
          merge('S', 'M', obs(j)%body == 1), obs(j)%qty, obs(j)%tkind, obs(j)%jd, u(j), v(j), pr(j), obs(j)%jpl, &
          dresid(obs(j), v(j), pr(j)), dresid(obs(j), v(j), obs(j)%jpl), t_used(j), obs(j)%dclock, obs(j)%sess
    end do
    close(un)
    print '(A)', '   residuals written to tycho_blind_residuals.dat'
  end subroutine

  function group_name(g) result(s)
    integer, intent(in) :: g
    character(len=18) :: s
    integer :: y, b, q, c
    c = (g - 1)/NGB
    y = mod(g - 1, NGB)/(2*NQ); b = mod((g - 1)/NQ, 2); q = mod(g - 1, NQ) + 1
    write(s, '(I4,1X,A4,1X,A6,1X,A1)') YEARS(y + 1), merge('Sun ', 'Moon', b == 0), QNAMES(q), 'SAM'(c+1:c+1)
  end function

  subroutine tycho_sigmas()
    ! scatter of Tycho's readings about JPL, per group, for the synthetic noise
    real(dp) :: rr(nobs)
    logical :: u(nobs)
    do i = 1, nobs
      rr(i) = dresid(obs(i), obs(i)%value, obs(i)%jpl)
    end do
    u = .true.
    sigma_g = 1.0_dp
    call update_sigmas(rr, u)
  end subroutine

  subroutine add_noise(p, s, seed)
    real(dp), intent(in) :: p(nobs)
    real(dp), intent(out) :: s(nobs)
    integer, intent(in) :: seed
    integer :: nseed
    integer, allocatable :: sd(:)
    real(dp) :: u1, u2
    call random_seed(size=nseed)
    allocate(sd(nseed)); sd = 12345 + 7919*seed
    call random_seed(put=sd)
    do i = 1, nobs
      call random_number(u1); call random_number(u2)
      s(i) = p(i) + sigma_g(obs(i)%group)/60.0_dp*sqrt(-2.0_dp*log(max(u1, 1.0d-12)))*cos(TAU*u2)
    end do
  end subroutine

  function start_guess(xt) result(xs)
    ! a deliberately wrong start: angles off by 0.3-1 deg, eccentricities by
    ! 30 %, periods by 1e-4, distances by 10-20 %
    real(dp), intent(in) :: xt(NP)
    real(dp) :: xs(NP)
    xs = xt
    xs(1) = xt(1)*1.3_dp;  xs(8) = xt(8)*0.7_dp
    xs(2) = xt(2) + 0.3_dp*DEG2RAD; xs(9) = xt(9) - 0.5_dp*DEG2RAD
    xs(3) = xt(3) - 0.5_dp*DEG2RAD; xs(10) = xt(10) + 1.0_dp*DEG2RAD
    xs(4) = xt(4) + 1.0_dp*DEG2RAD; xs(11) = xt(11) - 1.0_dp*DEG2RAD
    xs(5) = xt(5) - 0.8_dp*DEG2RAD; xs(12) = xt(12) + 0.8_dp*DEG2RAD
    xs(6) = xt(6)*(1.0_dp + 1.0d-4); xs(13) = xt(13)*(1.0_dp - 1.0d-4)
    xs(7) = xt(7)*1.2_dp; xs(14) = xt(14)*0.9_dp
    xs(15:17) = 0.0_dp
  end function

  subroutine report(title, xf, sg, xt)
    character(len=*), intent(in) :: title
    real(dp), intent(in) :: xf(NP), sg(NP), xt(NP)
    real(dp) :: f
    integer :: k2
    character(len=12) :: unit
    print '(/A)', trim(title)
    print '(A)', '   param        fitted            1-sigma            JPL             fit-JPL   (fit-JPL)/sigma'
    do k2 = 1, NP
      select case (k2)
      case (2:5, 9:12); f = RAD2DEG; unit = 'deg'
      case (6, 13);     f = RAD2DEG; unit = 'deg/day'
      case (7, 14);     f = 1.0_dp;  unit = 'km'
      case (15:17, 20, 21, 24, 26, 29:44, 47:53, 56:62, 65:71, 74:80, 83:88, 90:NP); f = 1.0_dp;  unit = 'arcmin'
      case default;     f = 1.0_dp;  unit = ''
      end select
      print '(3X,A6,ES17.9,ES14.4,ES17.9,ES12.3,F9.2,2X,A)', PNAME(k2), xf(k2)*f, sg(k2)*f, xt(k2)*f, &
          (xf(k2) - xt(k2))*f, (xf(k2) - xt(k2))/max(sg(k2), 1.0d-30), trim(unit)
    end do
    print '(A,F14.1,A,F10.1,A,F14.1,A,F8.4)', '   a_E  = ', xf(7), ' km  +-', sg(7), '   JPL ', xt(7), &
        ' km    ratio ', xf(7)/xt(7)
    print '(A,F14.1,A,F10.1,A)', '        = ', xf(7)/6378.137_dp, ' Earth radii  +-', sg(7)/6378.137_dp, &
        '   (Tycho himself assumed about 1150)'
    print '(A,F14.3,A,F10.3,A,F14.3,A)', '   a_M  = ', xf(14), ' km  +-', sg(14), '   JPL ', xt(14), ' km'
    print '(A,F12.4,A,F8.4,A)', '   a_M  = ', xf(14)/6378.137_dp, ' Earth radii  +-', sg(14)/6378.137_dp, &
        '   (JPL 60.6347)'
    print '(A,F10.2,A,F8.2,A,F7.4,A)', '   Moon radius = ', 1737.4_dp*xf(18), ' km  +-', 1737.4_dp*sg(18), &
        '  (', 1737.4_dp*xf(18)/6378.137_dp, ' Earth radii; true 1737.4 km, 0.2724)'
    print '(A,F8.4,A,F7.4,A,F8.3,A,F7.3,A)', '   refraction scale = ', xf(19), ' +-', sg(19), &
        '    latitude offset = ', xf(20), ' +-', sg(20), "'"
    if (tycho_model) print '(A,ES11.4,A,ES9.2,A,ES11.4,A,ES9.2,A,F7.2,A)', '   Cassini: refractive excess ', &
        2.9d-4*xf(19), ' +-', 2.9d-4*sg(19), '   height ', 1.3d-3*xf(25), ' +-', 1.3d-3*sg(25), &
        ' Earth radii (', 1.3d-3*xf(25)*6378.137_dp, ' km)'
    if (peryear) then
      print '(A,F10.2,A,F8.2,A,F8.4,A,F7.4,A,F8.3,A,F7.3,A)', '   1595: Moon radius = ', 1737.4_dp*xf(22), &
          ' km  +-', 1737.4_dp*sg(22), '   refraction ', xf(23), ' +-', sg(23), '   latitude ', xf(24), ' +-', sg(24), "'"
      print '(A,F10.2,A,F8.2,A,F8.4,A,F7.4,A,F8.3,A,F7.3,A)', '   1594: Moon radius = ', 1737.4_dp*xf(27), &
          ' km  +-', 1737.4_dp*sg(27), '   refraction ', xf(28), ' +-', sg(28), '   latitude ', xf(29), ' +-', sg(29), "'"
    end if
    print '(A,F8.3,A,F8.3)', '   a_E/a_M = ', xf(7)/xf(14), ' +-', xf(7)/xf(14)*sg(7)/xf(7)
    print '(A,F14.4,A,F14.4,A)', '   sidereal year  = ', TAU/xf(6), ' d   (JPL osculating ', TAU/xt(6), ' d)'
    print '(A,F14.5,A,F14.5,A)', '   sidereal month = ', TAU/xf(13), ' d   (JPL osculating ', TAU/xt(13), ' d)'
    print '(A,ES14.6,A,ES14.6)', '   GM_sun+GM_earth  = ', (xf(6)/DAY_S)**2*xf(7)**3, '   JPL ', (xt(6)/DAY_S)**2*xt(7)**3
    print '(A,ES14.6,A,ES14.6)', '   GM_earth+GM_moon = ', (xf(13)/DAY_S)**2*xf(14)**3, '   JPL ', (xt(13)/DAY_S)**2*xt(14)**3
  end subroutine

  subroutine print_sigmas()
    print '(/A)', '   Robust scatter of the residuals per group (arcmin):'
    do g = 1, NG
      n = count([(obs(i)%group == g, i = 1, nobs)])
      if (n > 0) print '(5X,A,I6,F8.2)', group_name(g), n, sigma_g(g)
    end do
  end subroutine

  subroutine write_residuals(xf, u)
    real(dp), intent(in) :: xf(NP)
    logical, intent(in) :: u(nobs)
    integer :: un
    call predict_all(xf, pred)
    open(newunit=un, file='tycho_fit_residuals.dat', status='replace')
    write(un, '(A)') '# id year body qty jd_tt used value_deg fit_deg jpl_deg resid_fit_arcmin resid_jpl_arcmin'
    do i = 1, nobs
      write(un, '(A12,I5,1X,A1,1X,A6,F18.8,L3,3F14.7,2F10.3)') obs(i)%id, obs(i)%year, &
          merge('S', 'M', obs(i)%body == 1), obs(i)%qty, obs(i)%jd, u(i), obs(i)%value, pred(i), obs(i)%jpl, &
          dresid(obs(i), obs(i)%value, pred(i)), dresid(obs(i), obs(i)%value, obs(i)%jpl)
    end do
    close(un)
    print '(A)', '   residuals written to tycho_fit_residuals.dat'
  end subroutine

end program tycho_fit
