program exp_ph_misfit
  use exponentiated_phase_misfit
  use sacio

  implicit none

  character(len=MAX_STRING_LEN) :: fsyn, fobs, fadj, fsac, sshort_p, slong_p, &
                                   stime_min, stime_max
  real(kind=dp) :: short_p, long_p, time_min, time_max, dt
  real(kind=dp), dimension(:), allocatable :: syn, dat
  real(kind=dp), dimension(:,:), allocatable :: windows
  integer :: ier
  type(sachead) :: header
  type(ExponentiatedPhaseMisfit) :: epm

  if (command_argument_count() /= 7) then
    write(*,*) 'Usage: exp_ph_misfit fobs fsyn short_p long_p time_min time_max fadj'
    write(*,*) '  fobs: observed seismogram in SAC format'
    write(*,*) '  fsyn: synthetic seismogram in SAC format'
    write(*,*) '  short_p: short period of bandpass filter (s)'
    write(*,*) '  long_p: long period of bandpass filter (s)'
    write(*,*) '  time_min: start time of time window (s)'
    write(*,*) '  time_max: end time of time window (s)'
    write(*,*) '  fadj: output path to adjoint source'
    stop
  end if

  ! get command line arguments
  call get_command_argument(1, fobs)
  call get_command_argument(2, fsyn)
  call get_command_argument(3, sshort_p)
  call get_command_argument(4, slong_p)
  call get_command_argument(5, stime_min)
  call get_command_argument(6, stime_max)
  call get_command_argument(7, fadj)
  read(sshort_p, *) short_p
  read(slong_p, *) long_p
  read(stime_min, *) time_min
  read(stime_max, *) time_max

  ! read observed and synthetic seismograms
  call sacio_readsac(fobs, header, dat, ier)
  call sacio_readsac(fsyn, header, syn, ier)
  dt = dble(header%delta)

  ! filter seismograms
  cfg%min_period = short_p
  cfg%max_period = long_p
  call bandpass_dp(syn, dt, real(1/long_p), real(1/short_p), IORD)
  call bandpass_dp(dat, dt, real(1/long_p), real(1/short_p), IORD)

  ! set up time windows
  allocate(windows(1, 2))
  windows(1, 1) = time_min - header%b
  windows(1, 2) = time_max - header%b

  ! calculate exponentiated phase misfit and adjoint source
  call epm%calc_adjoint_source(dat, syn, dt, windows)

  ! display results
  write(*,'(a,F8.5)') 'Real part residual: ', epm%residuals_real(1)
  write(*,'(a,F8.5)') 'Imaginary part residual: ', epm%residuals_imag(1)
  write(*,'(a,F8.5)') 'Real part misfit: ', epm%misfits_real(1)
  write(*,'(a,F8.5)') 'Imaginary part misfit: ', epm%misfits_imag(1)
  write(*,'(a,F8.5)') 'Total misfit: ', epm%total_misfit

  ! write adjoint source to SAC file
  fsac = trim(fadj)//'/'//trim(header%knetwk)//'.'//trim(header%kstnm) &
         //'.BXZ.expphase.sac'
  call sacio_writesac(fsac, header, epm%adj_src, ier)

  write(*,'(a)') 'Adjoint source written to: '//trim(fsac)

end program exp_ph_misfit