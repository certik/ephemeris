! ─────────────────────────────────────────────────────────────────────
!  obs_reader_mod — reads transit circle RA/Dec observations
! ─────────────────────────────────────────────────────────────────────
module obs_reader_mod
  use constants_mod
  implicit none

  integer, parameter :: MAX_OBS = 10000

  type :: observation
    real(dp) :: jd
    integer  :: body       ! 1=Sun, 2=Moon
    real(dp) :: ra_deg     ! right ascension (degrees, apparent of date)
    real(dp) :: dec_deg    ! declination (degrees, apparent of date)
    real(dp) :: dist       ! AU (0 if unavailable)
  end type

contains

  subroutine read_observations(filename, obs, n_obs, jd_max)
    character(len=*), intent(in)  :: filename
    type(observation), intent(out) :: obs(MAX_OBS)
    integer, intent(out) :: n_obs
    real(dp), intent(in) :: jd_max

    integer :: u, ios
    character(len=512) :: line
    character(1) :: body_c
    real(dp) :: jd, ra, dec, dist
    integer :: bi

    open(newunit=u, file=filename, status='old', action='read')
    n_obs = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle
      read(line, *, iostat=ios) jd, body_c, ra, dec, dist
      if (ios /= 0) cycle
      if (jd > jd_max) cycle
      if (body_c == 'S') then; bi = 1; else; bi = 2; end if
      n_obs = n_obs + 1
      if (n_obs > MAX_OBS) then
        n_obs = MAX_OBS; exit
      end if
      obs(n_obs)%jd      = jd
      obs(n_obs)%body    = bi
      obs(n_obs)%ra_deg  = ra
      obs(n_obs)%dec_deg = dec
      obs(n_obs)%dist    = dist
    end do
    close(u)
  end subroutine

end module obs_reader_mod
