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
!  Build and run:
!    gfortran -O2 -o build/tycho_fit src/constants_mod.f90 app/tycho_fit.f90
!    build/tycho_fit tycho_obs_2arcmin.dat jpl_state.dat
! ═══════════════════════════════════════════════════════════════════════
module tycho_fit_mod
  use constants_mod
  implicit none

  integer, parameter :: NP_ORB = 14, NP_NUIS = 3, NP = NP_ORB + NP_NUIS
  integer, parameter :: MAXOBS = 20000
  real(dp), parameter :: LAT_DEG = 55.0_dp + 54.0_dp/60 + 27.0_dp/3600
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

  character(len=6), parameter :: PNAME(NP) = [character(len=6) :: &
      'e_E', 'i_E', 'Om_E', 'w_E', 'M0_E', 'n_E', 'a_E', &
      'e_M', 'i_M', 'Om_M', 'w_M', 'M0_M', 'n_M', 'a_M', &
      'b_alt', 'b_d86', 'b_d95']

  type :: reading
    character(len=12) :: id
    integer :: year, body          ! body 1 = Sun, 2 = Moon
    character(len=6) :: qty
    character(len=1) :: point, refpoint
    character(len=4) :: refkind
    real(dp) :: jd, last, M(3,3), sra, sdec, value, jpl, group_err
    integer :: group
  end type

  type(reading) :: obs(MAXOBS)
  integer :: nobs = 0
  integer :: order(MAXOBS)          ! obs sorted by time
  real(dp) :: epoch_jd
  real(dp) :: obs_xyz(3)            ! observer, geocentric (km), in the frame of date at LAST = 0
  real(dp) :: coslat_rho, sinlat_rho

  ! Groups (year x body x quantity) for robust weighting
  integer, parameter :: NQ = 8, NG = 2*2*NQ
  character(len=6), parameter :: QNAMES(NQ) = [character(len=6) :: &
      'DEC', 'ALT', 'ALTM', 'HA', 'AZ', 'DIST', 'RADIFF', 'DIAM']
  real(dp) :: sigma_g(NG)
  logical :: freeze_sigma = .false.   ! keep sigma_g fixed during a fit (profiles)
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
    hr = H*DEG2RAD; dr = d*DEG2RAD; p = LAT_DEG*DEG2RAD
    sa = sin(p)*sin(dr) + cos(p)*cos(dr)*cos(hr)
    a = asin(max(-1.0_dp, min(1.0_dp, sa)))*RAD2DEG
    A_az = modulo(atan2(-cos(dr)*sin(hr), sin(dr)*cos(p) - cos(dr)*cos(hr)*sin(p))*RAD2DEG, 360.0_dp)
  end subroutine

  subroutine hor_to_eq(a, A_az, H, d)
    real(dp), intent(in) :: a, A_az
    real(dp), intent(out) :: H, d
    real(dp) :: alr, azr, p, sd
    alr = a*DEG2RAD; azr = A_az*DEG2RAD; p = LAT_DEG*DEG2RAD
    sd = sin(p)*sin(alr) + cos(p)*cos(alr)*cos(azr)
    d = asin(max(-1.0_dp, min(1.0_dp, sd)))*RAD2DEG
    H = atan2(-sin(azr)*cos(alr), cos(p)*sin(alr) - sin(p)*cos(alr)*cos(azr))*RAD2DEG
  end subroutine

  pure function refraction_deg(alt) result(r)
    real(dp), intent(in) :: alt
    real(dp) :: r
    if (alt < -1.5_dp) then
      r = 0.0_dp; return
    end if
    r = 1.02_dp / tan((alt + 10.3_dp/(alt + 5.11_dp))*DEG2RAD)
    r = r * (PRES_MBAR/1010.0_dp) * (283.0_dp/(273.0_dp + TEMP_C))
    r = max(r, 0.0_dp) / 60.0_dp
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
    end select
  end subroutine

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
    do i = 1, NBODY
      do j = i+1, NBODY
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
    gm_m = mu_em*MOON_EARTH_RATIO/(1.0_dp + MOON_EARTH_RATIO)
    gm_e = mu_em - gm_m
    gm_s = mu_se - gm_e
    gm(1:3) = [gm_s, gm_e, gm_m]
    call kepler_to_cart(x(7), x(1), x(2), x(3), x(4), x(5), mu_se, rE, vE)
    call kepler_to_cart(x(14), x(8), x(9), x(10), x(11), x(12), mu_em, rM, vM)
    pos(:,1) = sun_state(1:3);    vel(:,1) = sun_state(4:6)
    pos(:,2) = pos(:,1) + rE;     vel(:,2) = vel(:,1) + vE
    pos(:,3) = pos(:,2) + rM;     vel(:,3) = vel(:,2) + vM
    do k = 1, 5
      pos(:,3+k) = planet_state(1:3,k); vel(:,3+k) = planet_state(4:6,k)
      gm(3+k) = merge(GM_PLANETS(k), 0.0_dp, use_planets)
    end do
  end subroutine

  ! ── Topocentric apparent (H, dec, semidiameter) of a body ───────────
  subroutine topo(rb_icrs, v_earth, o, H, d, sd, radius)
    real(dp), intent(in) :: rb_icrs(3), v_earth(3), radius
    type(reading), intent(in) :: o
    real(dp), intent(out) :: H, d, sd
    real(dp) :: r(3), u(3), dn, rd(3), ob(3), lst, ra
    r = rb_icrs
    dn = norm2(r); u = r/dn
    r = r + dn*(v_earth/C_LIGHT_KMS - u*dot_product(u, v_earth/C_LIGHT_KMS))   ! annual aberration
    rd = matmul(o%M, r)
    lst = o%last*DEG2RAD
    ob = [coslat_rho*cos(lst), coslat_rho*sin(lst), sinlat_rho]
    rd = rd - ob
    dn = norm2(rd)
    ra = atan2(rd(2), rd(1))*RAD2DEG
    d = asin(rd(3)/dn)*RAD2DEG
    H = wrap180(o%last - ra)
    sd = asin(radius/dn)*RAD2DEG
  end subroutine

  ! ── Prediction of one reading (deg) ─────────────────────────────────
  function predict_one(o, x, rsun, rmoon, vearth) result(p)
    type(reading), intent(in) :: o
    real(dp), intent(in) :: x(NP), rsun(3), rmoon(3), vearth(3)
    real(dp) :: p, H, d, sd, Hl, dl, Hr, dr, a, A_az, Hs, ds, sds, Hc, dc, Ho, dov
    if (o%body == 1) then
      call topo(rsun, vearth, o, H, d, sd, SUN_RADIUS)
    else
      call topo(rmoon, vearth, o, H, d, sd, MOON_RADIUS)
    end if
    p = 0.0_dp
    select case (trim(o%qty))
    case ('ALT', 'ALTM')
      if (o%point == 'U' .or. o%point == 'L') then
        call limb_point(H, d, sd, o%point, Hl, dl)
      else
        Hl = H; dl = d
      end if
      call eq_to_hor(Hl, dl, a, A_az)
      p = a + refraction_deg(a) + x(15)/60.0_dp
    case ('DEC')
      call limb_point(H, d, sd, o%point, Hl, dl)
      call refract(Hl, dl, Hr, dr)
      p = dr + merge(x(16), x(17), o%year == 1586)/60.0_dp
    case ('HA')
      call limb_point(H, d, sd, o%point, Hl, dl)
      call refract(Hl, dl, Hr, dr)
      p = Hr
    case ('AZ')
      call limb_point(H, d, sd, o%point, Hl, dl)
      call refract(Hl, dl, Hr, dr)
      call eq_to_hor(Hr, dr, a, A_az)
      p = wrap180(A_az - 180.0_dp)
    case ('DIAM')
      p = 2.0_dp*sd
    case ('DIST', 'RADIFF')
      if (trim(o%refkind) == 'STAR') then
        Hs = wrap180(o%last - o%sra); ds = o%sdec; sds = 0.0_dp
      else
        call topo(rsun, vearth, o, Hs, ds, sds, SUN_RADIUS)
      end if
      if (trim(o%qty) == 'RADIFF') then
        call limb_point(H, d, sd, o%point, Hl, dl)
        call refract(Hl, dl, Hr, dr)
        call limb_point(Hs, ds, sds, o%refpoint, Hl, dl)
        call refract(Hl, dl, Ho, dov)
        p = wrap180(Ho - Hr)
      else
        call refract(H, d, Hc, dc)
        call refract(Hs, ds, Ho, dov)
        p = sep(Hc, dc, Ho, dov)
        select case (o%point)
        case ('N'); p = p - sd
        case ('F'); p = p + sd
        case ('U', 'L', 'E', 'W')
          if (side(Hc, dc, Ho, dov, o%point) > 0.0_dp) then
            p = p - sd
          else
            p = p + sd
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
  end function

  ! ── Predictions for all readings for a parameter vector ─────────────
  subroutine predict_all(x, pred)
    real(dp), intent(in) :: x(NP)
    real(dp), intent(out) :: pred(nobs)
    real(dp) :: pos0(3,NBODY), vel0(3,NBODY), gm(NBODY), pos(3,NBODY), vel(3,NBODY), t
    integer :: k, i
    call initial_state(x, pos0, vel0, gm)
    ! forward from the epoch
    pos = pos0; vel = vel0; t = epoch_jd
    do k = 1, nobs
      i = order(k)
      if (obs(i)%jd < epoch_jd) cycle
      call propagate(pos, vel, gm, (obs(i)%jd - t)*DAY_S); t = obs(i)%jd
      pred(i) = predict_one(obs(i), x, pos(:,1) - pos(:,2), pos(:,3) - pos(:,2), vel(:,2))
    end do
    ! backward from the epoch
    pos = pos0; vel = vel0; t = epoch_jd
    do k = nobs, 1, -1
      i = order(k)
      if (obs(i)%jd >= epoch_jd) cycle
      call propagate(pos, vel, gm, (obs(i)%jd - t)*DAY_S); t = obs(i)%jd
      pred(i) = predict_one(obs(i), x, pos(:,1) - pos(:,2), pos(:,3) - pos(:,2), vel(:,2))
    end do
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
    integer :: u, ios, i, j, k, tmp
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
      o%M = transpose(reshape(m9, [3, 3]))
      o%group_err = err
      o%group = (merge(0, 1, o%year == 1586))*2*NQ + (o%body - 1)*NQ + qindex(o%qty)
      nobs = nobs + 1
      obs(nobs) = o
    end do
    close(u)
    ! sort by time (insertion sort on an index)
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
    k = 0
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
  end subroutine

  subroutine setup_observer()
    ! geocentric position of Uraniborg (WGS84), km
    real(dp) :: f, e2, phi, N
    f = 1.0_dp/WGS84_INVF; e2 = f*(2.0_dp - f)
    phi = LAT_DEG*DEG2RAD
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
    if (freeze_sigma) return
    do g = 1, NG
      n = 0
      do i = 1, nobs
        if (obs(i)%group == g .and. use(i)) then
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
      r0(i) = resid_arcmin(obs(i), values(i), p0(i))
    end do
    J = 0.0_dp
    do k = 1, NP
      if (.not. active(k)) cycle
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
    integer :: it, k, i, na, idx(NP), m, nstall
    logical :: ok
    allocate(J(nobs, NP), r(nobs))
    nstall = 0
    if (.not. freeze_use) use = .true.
    if (.not. freeze_sigma) sigma_g = 1.0_dp
    lam = 1.0d-3
    na = 0
    do k = 1, NP
      if (active(k)) then
        na = na + 1; idx(na) = k
      end if
    end do
    do it = 1, 80
      call jacobian(x, values, J, r, active)
      call update_sigmas(r, use)
      if (clip > 0.0_dp .and. .not. freeze_use) then
        do i = 1, nobs
          use(i) = abs(r(i)) <= clip*sigma_g(obs(i)%group)
        end do
        call update_sigmas(r, use)
      end if
      chi2 = 0.0_dp; N = 0.0_dp; g = 0.0_dp
      do i = 1, nobs
        if (.not. use(i)) cycle
        w = 1.0_dp/sigma_g(obs(i)%group)**2
        chi2 = chi2 + w*r(i)**2
        do k = 1, na
          g(idx(k)) = g(idx(k)) + w*J(i, idx(k))*r(i)
          N(idx(k), idx(1:na)) = N(idx(k), idx(1:na)) + w*J(i, idx(k))*J(i, idx(1:na))
        end do
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
          if (ok) dx(idx(1:na)) = s
          do k = 1, NP
            dx(k) = dx(k)*pscale(x, k)
          end do
        end block
        xn = x + dx
        xn(1) = max(1.0d-6, min(0.5_dp, xn(1))); xn(8) = max(1.0d-6, min(0.5_dp, xn(8)))
        call predict_all(xn, pn)
        chi2n = 0.0_dp
        do i = 1, nobs
          if (.not. use(i)) cycle
          rn(i) = resid_arcmin(obs(i), values(i), pn(i))
          chi2n = chi2n + rn(i)**2/sigma_g(obs(i)%group)**2
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
      rms = rms + resid_arcmin(obs(i), values(i), pn(i))**2; k = k + 1
    end do
    rms = sqrt(rms/max(k, 1))
    last_chi2 = 0.0_dp
    do i = 1, nobs
      if (use(i)) last_chi2 = last_chi2 + (resid_arcmin(obs(i), values(i), pn(i))/sigma_g(obs(i)%group))**2
    end do
    Ninv = 0.0_dp
    deallocate(J, r)
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
        if (active(k)) P(k,k) = sigma_prior(k)**2
      end do
      do k = 1, NP
        if (active(k)) P(k,k) = P(k,k)/pscale(xr, k)**2
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
      x(1) = max(1.0d-6, min(0.5_dp, x(1))); x(8) = max(1.0d-6, min(0.5_dp, x(8)))
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
  if (len_trim(mode) == 0) mode = 'fit'
  clip = 0.0_dp; n_mc = 10; seed0 = 1
  call get_command_argument(4, arg); if (len_trim(arg) > 0) read(arg, *) clip
  call get_command_argument(5, arg); if (len_trim(arg) > 0) read(arg, *) n_mc
  call get_command_argument(6, arg); if (len_trim(arg) > 0) read(arg, *) seed0
  call get_command_argument(7, arg); if (trim(arg) == 'noplanets') use_planets = .false.

  call setup_observer()
  call read_obs(trim(obs_file))
  call read_state(trim(state_file), x_true)
  ! Earth's pole (ICRS) at the epoch: third row of the rotation to the
  ! equator of date of the reading nearest to the epoch
  k = minloc([(abs(obs(i)%jd - epoch_jd), i = 1, nobs)], 1)
  earth_pole = obs(k)%M(3,:)
  call get_command_argument(8, arg); if (trim(arg) == 'noj2') use_j2 = .false.
  allocate(pred(nobs), values(nobs), synth(nobs), use(nobs))
  print '(A,I6,A,A,A,A)', ' Readings: ', nobs, '  from ', trim(obs_file), '   mode ', trim(mode)
  print '(A,F12.3,A,F5.1,A,L2)', ' Epoch (JD TDB): ', epoch_jd, '   clip (sigma, 0 = none): ', clip, &
      '   planets:', use_planets
  values = [(obs(i)%value, i = 1, nobs)]
  active = .true.
  call predict_all(x_true, pred)

  select case (trim(mode))

  case ('check')
    print '(/A)', ' Model check: this program with the JPL orbits vs. compare_jpl.py (JPL DE441)'
    print '(A)', '    group               n    rms(model-JPL)   rms(Tycho-JPL)'
    do g = 1, NG
      n = 0; gn = 0.0_dp
      do i = 1, nobs
        if (obs(i)%group /= g) cycle
        n = n + 1
        gn(1) = gn(1) + resid_arcmin(obs(i), pred(i), obs(i)%jpl)**2
        gn(2) = gn(2) + resid_arcmin(obs(i), obs(i)%value, obs(i)%jpl)**2
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
      write(*, '(A,I5,17ES20.11)') ' MC ', mc, x - x_true
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
    write(*, '(A,17ES20.11)') ' FIT ', x
    write(*, '(A,17ES20.11)') ' SIG ', sig
    print '(/A)', ' Iterated extended Kalman filter (sequential updates, F = I, Q = 0), same readings and weights'
    sigma_prior = [0.05_dp, 0.05_dp, 0.05_dp, 0.5_dp, 0.5_dp, 1.0d-4, 1.0d7, &
                   0.05_dp, 0.05_dp, 0.05_dp, 0.5_dp, 0.5_dp, 1.0d-2, 4.0d4, 5.0_dp, 5.0_dp, 5.0_dp]
    x = x0
    call fit_iekf(x, values, active, use, sigma_prior, 30, P)
    do k = 1, NP
      sig(k) = sqrt(max(P(k,k), 0.0_dp))
    end do
    call predict_all(x, pred)
    rms = 0.0_dp; n = 0; r = 0.0_dp
    do i = 1, nobs
      if (.not. use(i)) cycle
      gn(1) = resid_arcmin(obs(i), values(i), pred(i))
      rms = rms + gn(1)**2; r = r + (gn(1)/sigma_g(obs(i)%group))**2; n = n + 1
    end do
    print '(A,F8.3,A,F12.2)', '   rms residual ', sqrt(rms/n), "'   chi2 ", r
    call report(' Iterated EKF result (1-sigma from the filter covariance)', x, sig, x_true)
    write(*, '(A,17ES20.11)') ' EKF ', x

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
      sigma_g = best_sigma
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

  case ('kepler')
    ! the Moon on a fixed Kepler ellipse about the Earth (no Sun)
    kepler_moon = .true.
    x = start_guess(x_true)
    call fit_lm(x, values, active, clip, cov, use, rms, .true.)
    print '(A,F8.3,A)', '   rms residual with a fixed Kepler ellipse for the Moon: ', rms, "'"
    call print_sigmas()

  end select

contains

  function group_name(g) result(s)
    integer, intent(in) :: g
    character(len=18) :: s
    integer :: y, b, q
    y = (g - 1)/(2*NQ); b = mod((g - 1)/NQ, 2); q = mod(g - 1, NQ) + 1
    write(s, '(I4,1X,A4,1X,A6)') merge(1586, 1595, y == 0), merge('Sun ', 'Moon', b == 0), QNAMES(q)
  end function

  subroutine tycho_sigmas()
    ! scatter of Tycho's readings about JPL, per group, for the synthetic noise
    real(dp) :: rr(nobs)
    logical :: u(nobs)
    do i = 1, nobs
      rr(i) = resid_arcmin(obs(i), obs(i)%value, obs(i)%jpl)
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
      case (15:17);     f = 1.0_dp;  unit = 'arcmin'
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
          resid_arcmin(obs(i), obs(i)%value, pred(i)), resid_arcmin(obs(i), obs(i)%value, obs(i)%jpl)
    end do
    close(un)
    print '(A)', '   residuals written to tycho_fit_residuals.dat'
  end subroutine

end program tycho_fit
