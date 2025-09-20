program ex_rf_misfit
  use rf_misfit
  use sacio

  implicit none

  character(len=MAX_STRING_LEN) :: fsyn, fobs, fadj, fsac, sshort_p, slong_p,&
                                   stime_min, stime_max, fsynr, fsynz
  real(kind=dp) :: time_min, time_max, dt, tp, short_p, long_p
  real(kind=dp), dimension(:), allocatable :: syn, dat, synr, synz
  real(kind=dp), dimension(2) :: window
  integer :: ier
  real(kind=dp) :: f0 = 2.5
  type(sachead) :: header
  type(RFMisfit) :: rfm

  if (command_argument_count() /= 9) then
    write(*,*) 'Usage: rf_misfit rfobs rfsyn synr synz time_min time_max fadj'
    write(*,*) '  fobs: observed RF in SAC format'
    write(*,*) '  fsyn: synthetic RF in SAC format'
    write(*,*) '  synr: Synthetic seismogram radial component in SAC format'
    write(*,*) '  synz: Synthetic seismogram vertical component in SAC format'
    write(*,*) '  short_p: Short period for filtering (s)'
    write(*,*) '  long_p: Long period for filtering (s)'
    write(*,*) '  time_min: start time of time window (s)'
    write(*,*) '  time_max: end time of time window (s)'
    write(*,*) '  fadj: output path to adjoint source'
    stop
  end if

  ! get command line arguments
  call get_command_argument(1, fobs)
  call get_command_argument(2, fsyn)
  call get_command_argument(3, fsynr)
  call get_command_argument(4, fsynz)
  call get_command_argument(5, sshort_p)
  call get_command_argument(6, slong_p)
  call get_command_argument(7, stime_min)
  call get_command_argument(8, stime_max)
  call get_command_argument(9, fadj)
  read(sshort_p, *) short_p
  read(slong_p, *) long_p
  read(stime_min, *) time_min
  read(stime_max, *) time_max

  ! read observed and synthetic seismograms
  call sacio_readsac(fobs, header, dat, ier)
  if (ier /= 0) then
    print *, 'Error reading observed SAC file:', trim(fobs)
    stop
  end if
  call sacio_readsac(fsyn, header, syn, ier)
  if (ier /= 0) then
    print *, 'Error reading synthetic SAC file:', trim(fsyn)
    stop
  end if
  call sacio_readsac(fsynr, header, synr, ier)
  if (ier /= 0) then
    print *, 'Error reading synthetic radial SAC file:', trim(fsynr)
    stop
  end if
  call sacio_readsac(fsynz, header, synz, ier)
  if (ier /= 0) then
    print *, 'Error reading synthetic vertical SAC file:', trim(fsynz)
    stop
  end if
  dt = dble(header%delta)

  ! filter seismograms
  cfg%min_period = short_p
  cfg%max_period = long_p
  call bandpass_dp(synr, dt, real(1/long_p), real(1/short_p), IORD)
  call bandpass_dp(synz, dt, real(1/long_p), real(1/short_p), IORD)

  tp = maxloc(synz, dim=1) * dt
  window = [time_min, time_max]
  call rfm%calc_adjoint_source(dat, syn, synr, synz, dt, tp, window, f0)

  print *, 'RF misfit:', rfm%total_misfit
  ! write adjoint source to SAC file
  fsac = trim(fadj)//trim(header%knetwk)//'.'//trim(header%kstnm) &
         //'.BXR.rf.adj.sac'
  call sacio_writesac(fsac, header, rfm%adj_src_r, ier)
  fsac = trim(fadj)//trim(header%knetwk)//'.'//trim(header%kstnm) &
         //'.BXZ.rf.adj.sac'
  call sacio_writesac(fsac, header, rfm%adj_src_z, ier)

end program