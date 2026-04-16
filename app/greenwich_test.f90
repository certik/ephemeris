program greenwich_test
  ! Compare topocentric vs geocentric Moon Dec for a Greenwich
  ! transit circle observation: 1948 Jan 18, ~17h UT.
  ! Observed: RA = 1h 13m 27.56s, Dec = +4d 52' 52.65"
  !
  ! A transit circle measures the moment an object crosses the
  ! meridian, so the RA IS the precise time measurement.
  ! We find the UT when computed Moon RA = observed RA (transit),
  ! then compare the computed Dec with the observed Dec.
  use constants_mod
  use linalg_mod
  use spk_reader_mod
  use nutation_mod
  use astro_mod
  implicit none

  type(spk_kernel) :: kernel
  real(dp) :: lat_deg, lon_deg, elev_m
  real(dp) :: delta_t
  integer  :: leap_sec
  real(dp) :: obs_ra_deg, obs_dec_deg
  real(dp) :: topo_ut, topo_dec, geo_ut, geo_dec
  character(len=256) :: bsp_file

  call get_command_argument(1, bsp_file)
  if (len_trim(bsp_file) == 0) bsp_file = 'de440s.bsp'
  call load_nutation('nutation.dat')
  call spk_open(trim(bsp_file), kernel)

  ! ── Greenwich Royal Observatory ──
  lat_deg = 51.4772_dp
  lon_deg = -0.0015_dp
  elev_m  = 30.0_dp

  ! delta_T(1948) ~ 29.15 s; no leap seconds before 1972
  delta_t = 29.15_dp
  leap_sec = 0

  ! Observed values from Greenwich 1948 PDF, page A12
  obs_ra_deg  = (1.0_dp + 13.0_dp/60.0_dp + 27.56_dp/3600.0_dp) * 15.0_dp
  obs_dec_deg = 4.0_dp + 52.0_dp/60.0_dp + 52.65_dp/3600.0_dp

  print '(A)', '================================================================='
  print '(A)', ' Greenwich Transit Circle - Moon - 1948 Jan 18'
  print '(A)', ' Observed: RA = 1h 13m 27.56s, Dec = +4d 52'' 52.65"'
  print '(A)', '================================================================='
  print '(A)', ''
  print '(A)', ' Finding UT when computed Moon RA = observed RA ...'

  ! Find transit time & Dec for TOPOCENTRIC
  call find_transit_dec(kernel, 1948, 1, 18, lat_deg, lon_deg, elev_m, &
       delta_t, leap_sec, obs_ra_deg, .true., topo_ut, topo_dec)

  ! Find transit time & Dec for GEOCENTRIC
  call find_transit_dec(kernel, 1948, 1, 18, lat_deg, lon_deg, elev_m, &
       delta_t, leap_sec, obs_ra_deg, .false., geo_ut, geo_dec)

  print '(A)', ''
  print '(A)',        '  OBSERVED (from PDF page A12):'
  call print_ra_dec('    ', obs_ra_deg, obs_dec_deg)
  print '(A)', ''

  print '(A)',        '  TOPOCENTRIC (from Greenwich surface):'
  print '(A,F14.8,A)', '    Transit UT = ', topo_ut, ' h'
  print '(A)',         '    RA  = (matched to observed by construction)'
  call print_dec('    Dec = ', topo_dec)
  print '(A,F10.2,A)', '    Dec residual (obs - comp): ', &
       (obs_dec_deg - topo_dec) * 3600.0_dp, '"'
  print '(A)', ''

  print '(A)',        '  GEOCENTRIC (from Earth center):'
  print '(A,F14.8,A)', '    Transit UT = ', geo_ut, ' h'
  print '(A)',         '    RA  = (matched to observed by construction)'
  call print_dec('    Dec = ', geo_dec)
  print '(A,F10.2,A)', '    Dec residual (obs - comp): ', &
       (obs_dec_deg - geo_dec) * 3600.0_dp, '"'
  print '(A)', ''

  print '(A)', ' ========================================='
  print '(A,F10.2,A)', '  Topocentric Dec residual: ', &
       (obs_dec_deg - topo_dec) * 3600.0_dp, '"'
  print '(A,F10.2,A)', '  Geocentric  Dec residual: ', &
       (obs_dec_deg - geo_dec) * 3600.0_dp, '"'
  print '(A,F10.2,A)', '  Topo - Geo (parallax):    ', &
       (topo_dec - geo_dec) * 3600.0_dp, '"'
  print '(A)', ' ========================================='

  call spk_close(kernel)

contains

  subroutine compute_moon_radec(kernel, jd_whole, utc_frac, &
       lat_deg, lon_deg, elev_m, delta_t, leap_sec, topocentric, &
       ra_deg, dec_deg)
    type(spk_kernel), intent(inout) :: kernel
    real(dp), intent(in) :: jd_whole, utc_frac
    real(dp), intent(in) :: lat_deg, lon_deg, elev_m, delta_t
    integer, intent(in)  :: leap_sec
    logical, intent(in)  :: topocentric
    real(dp), intent(out) :: ra_deg, dec_deg

    real(dp) :: tt_frac, tdb_frac, ut1_frac
    real(dp) :: M(3,3), d_psi, d_eps, mean_ob
    real(dp) :: gmst_h, gast_h
    real(dp) :: R_itrs(3,3), RT(3,3)
    real(dp) :: itrs_pos(3), itrs_vel(3)
    real(dp) :: obs_gcrs(3), obs_vel_gcrs(3)
    real(dp) :: earth_pos(3), earth_vel(3)
    real(dp) :: bcrs_pos(3), bcrs_vel(3)
    real(dp) :: moon_dir(3), moon_vel(3), moon_lt
    real(dp) :: app(3), zero3(3)
    real(dp) :: jd_tt, jd_tdb

    call utc_to_tt(jd_whole, utc_frac, leap_sec, tt_frac)
    call tt_to_tdb(jd_whole, tt_frac, tdb_frac)
    call tt_to_ut1(jd_whole, tt_frac, delta_t, ut1_frac)

    jd_tt  = jd_whole + tt_frac
    jd_tdb = jd_whole + tdb_frac

    call compute_M(jd_tt, jd_tdb, M, d_psi, d_eps, mean_ob)
    gmst_h = greenwich_mean_sidereal_time(jd_whole, ut1_frac, jd_tdb)
    gast_h = greenwich_apparent_sidereal_time(gmst_h, d_psi, mean_ob, jd_tt)
    R_itrs = itrs_rotation(gast_h, M)

    call earth_position_au(kernel, jd_whole, tdb_frac, earth_pos, earth_vel)

    if (topocentric) then
      call wgs84_to_itrs_au(lat_deg, lon_deg, elev_m, itrs_pos)
      call itrs_velocity_au_per_day(itrs_pos, itrs_vel)
      RT = mat33_T(R_itrs)
      obs_gcrs = mat33_vec(RT, itrs_pos)
      obs_vel_gcrs = mat33_vec(RT, itrs_vel)
      bcrs_pos = earth_pos + obs_gcrs
      bcrs_vel = earth_vel + obs_vel_gcrs
    else
      obs_gcrs = 0.0_dp
      obs_vel_gcrs = 0.0_dp
      bcrs_pos = earth_pos
      bcrs_vel = earth_vel
    end if

    call correct_light_travel_time(bcrs_pos, bcrs_vel, kernel, &
         jd_whole, tdb_frac, 2, moon_dir, moon_vel, moon_lt)
    zero3 = 0.0_dp
    call add_deflection(moon_dir, bcrs_pos, obs_gcrs, kernel, &
         jd_whole, tdb_frac)
    call add_aberration(moon_dir, bcrs_vel, moon_lt)

    ! GCRS -> true equator & equinox of date
    app = mat33_vec(M, moon_dir)
    ra_deg = atan2(app(2), app(1)) * 180.0_dp / PI
    if (ra_deg < 0.0_dp) ra_deg = ra_deg + 360.0_dp
    dec_deg = atan2(app(3), sqrt(app(1)**2 + app(2)**2)) * 180.0_dp / PI
  end subroutine

  subroutine find_transit_dec(kernel, year, month, day, &
       lat_deg, lon_deg, elev_m, delta_t, leap_sec, &
       target_ra, topocentric, transit_ut_h, transit_dec)
    ! Bisect to find UT when Moon RA = target_ra
    type(spk_kernel), intent(inout) :: kernel
    integer, intent(in)  :: year, month, day, leap_sec
    real(dp), intent(in) :: lat_deg, lon_deg, elev_m, delta_t
    real(dp), intent(in) :: target_ra
    logical, intent(in)  :: topocentric
    real(dp), intent(out) :: transit_ut_h, transit_dec

    real(dp) :: jd_whole, ut_lo, ut_hi, ut_mid
    real(dp) :: ra_lo, dec_lo, ra_hi, dec_hi, ra_mid, dec_mid
    real(dp) :: dra_lo, dra_mid
    integer :: jd_int, iter

    jd_int = julian_day(year, month, day)
    jd_whole = real(jd_int, dp)

    ! Search around 17h UT (+/- 4h)
    ut_lo = 13.0_dp / 24.0_dp - 0.5_dp
    ut_hi = 21.0_dp / 24.0_dp - 0.5_dp

    call compute_moon_radec(kernel, jd_whole, ut_lo, lat_deg, lon_deg, &
         elev_m, delta_t, leap_sec, topocentric, ra_lo, dec_lo)
    call compute_moon_radec(kernel, jd_whole, ut_hi, lat_deg, lon_deg, &
         elev_m, delta_t, leap_sec, topocentric, ra_hi, dec_hi)

    ! Bisection on RA difference
    do iter = 1, 60
      ut_mid = 0.5_dp * (ut_lo + ut_hi)
      call compute_moon_radec(kernel, jd_whole, ut_mid, lat_deg, lon_deg, &
           elev_m, delta_t, leap_sec, topocentric, ra_mid, dec_mid)

      dra_lo  = ra_lo  - target_ra
      dra_mid = ra_mid - target_ra
      ! Handle wrap-around
      if (dra_lo > 180.0_dp) dra_lo = dra_lo - 360.0_dp
      if (dra_lo < -180.0_dp) dra_lo = dra_lo + 360.0_dp
      if (dra_mid > 180.0_dp) dra_mid = dra_mid - 360.0_dp
      if (dra_mid < -180.0_dp) dra_mid = dra_mid + 360.0_dp

      if (dra_lo * dra_mid <= 0.0_dp) then
        ut_hi = ut_mid
        ra_hi = ra_mid
        dec_hi = dec_mid
      else
        ut_lo = ut_mid
        ra_lo = ra_mid
        dec_lo = dec_mid
      end if

      if (abs(ut_hi - ut_lo) * 86400.0_dp < 1.0d-6) exit
    end do

    transit_ut_h = (ut_mid + 0.5_dp) * 24.0_dp
    transit_dec  = dec_mid
  end subroutine

  subroutine print_ra_dec(prefix, ra_deg, dec_deg)
    character(len=*), intent(in) :: prefix
    real(dp), intent(in) :: ra_deg, dec_deg
    integer :: ra_h, ra_m, dec_d, dec_m
    real(dp) :: ra_s, dec_s, ra_rem, dec_rem
    character(len=1) :: dec_sign

    ra_rem = ra_deg / 15.0_dp
    ra_h = int(ra_rem)
    ra_rem = (ra_rem - ra_h) * 60.0_dp
    ra_m = int(ra_rem)
    ra_s = (ra_rem - ra_m) * 60.0_dp

    if (dec_deg < 0.0_dp) then
      dec_sign = '-'; dec_rem = -dec_deg
    else
      dec_sign = '+'; dec_rem = dec_deg
    end if
    dec_d = int(dec_rem)
    dec_rem = (dec_rem - dec_d) * 60.0_dp
    dec_m = int(dec_rem)
    dec_s = (dec_rem - dec_m) * 60.0_dp

    print '(A,A,I2.2,A,I2.2,A,F6.2,A)', prefix, 'RA  = ', ra_h, 'h ', &
         ra_m, 'm ', ra_s, 's'
    call print_dec(prefix // 'Dec = ', dec_deg)
  end subroutine

  subroutine print_dec(prefix, dec_deg)
    character(len=*), intent(in) :: prefix
    real(dp), intent(in) :: dec_deg
    integer :: dec_d, dec_m
    real(dp) :: dec_s, dec_rem
    character(len=1) :: dec_sign

    if (dec_deg < 0.0_dp) then
      dec_sign = '-'; dec_rem = -dec_deg
    else
      dec_sign = '+'; dec_rem = dec_deg
    end if
    dec_d = int(dec_rem)
    dec_rem = (dec_rem - dec_d) * 60.0_dp
    dec_m = int(dec_rem)
    dec_s = (dec_rem - dec_m) * 60.0_dp

    print '(A,A1,I2.2,A,I2.2,A,F5.2,A)', prefix, dec_sign, &
         dec_d, 'd ', dec_m, "' ", dec_s, '"'
  end subroutine

end program greenwich_test
