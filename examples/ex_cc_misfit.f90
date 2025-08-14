program cc_misfit
  use cross_correlate
  use sacio

  implicit none

  real(kind=dp), parameter :: dt_sigma_min = 1.0_dp ! minimum error in time shift (s)
  real(kind=dp), parameter :: dlna_sigma_min = 0.5_dp ! minimum error in ln(amplitude)
  character(len=MAX_STRING_LEN) :: fsyn, fobs, fadj, fsac, sshort_p, slong_p, &
                                   stime_min, stime_max
  real(kind=dp) :: short_p, long_p, time_min, time_max, tshift, dlna, sigma_dt, sigma_dlna, &
                   misfit_p, misfit_q, dt
  real(kind=dp), dimension(:), allocatable :: syn, dat, syn_tw, dat_tw, adj_tw_p, adj_tw_q, &
                                              adj_p, adj_q
  integer :: ier, nb, ne
  type(sachead) :: header

  if (command_argument_count() /= 7) then
    write(*,*) 'Usage: cc_misfit fobs fsyn short_p long_p time_min time_max fadj'
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
  call bandpass_dp(syn, dt, real(1/long_p), real(1/short_p), IORD)
  call bandpass_dp(dat, dt, real(1/long_p), real(1/short_p), IORD)

  ! cut time window
  nb = int((time_min - header%b) / dt) + 1
  ne = int((time_max - header%b) / dt) + 1
  syn_tw = syn(nb:ne)
  dat_tw = dat(nb:ne)

  ! calculate cross-correlation shift
  call calc_cc_shift(dat_tw, syn_tw, dt, dt_sigma_min, dlna_sigma_min, &
                tshift, dlna, sigma_dt, sigma_dlna)

  write(*,'(a,F7.5,a,F7.5)') 'Time shift (s): ', tshift, '+/-', dt_sigma_min
  write(*,'(a,F7.5,a,F7.5)') 'Amplitude anomaly (ln): ', dlna, '+/-', dlna_sigma_min

  call calc_cc_adjsrc(syn_tw, tshift, dlna, dt, sigma_dt, sigma_dlna, &
                      misfit_p, misfit_q, adj_tw_p, adj_tw_q)
  write(*,'(a,F7.5)') 'Time shift misfit: ', misfit_p
  write(*,'(a,F7.5)') 'Amplitude anomaly misfit: ', misfit_q

  allocate(adj_p(header%npts), adj_q(header%npts))
  adj_p = 0.0_dp
  adj_q = 0.0_dp
  adj_p(nb:ne) = adj_tw_p
  adj_q(nb:ne) = adj_tw_q

  ! write adjoint source to SAC file
  fsac = trim(fadj)//'/'//trim(header%knetwk)//'.'//trim(header%kstnm) &
         //'BXZ.ccdt.sac'
  call sacio_writesac(fsac, header, adj_p, ier)

  fsac = trim(fadj)//'/'//trim(header%knetwk)//'.'//trim(header%kstnm) &
         //'BXZ.ccdlna.sac'
  call sacio_writesac(fsac, header, adj_q, ier)

  
end program cc_misfit
