program ex_wcc_misfit
  use waveform_conv_misfit
  use sacio
  use signal
  use config, only: MAX_STRING_LEN, IORD, dp
  use adj_config, only: cfg => adj_config_global

  implicit none

  character(len=MAX_STRING_LEN) :: fobs_r, fobs_z, fsyn_r, fsyn_z, fadj, &
                                   stime_min, stime_max
  real(kind=dp) :: time_min, time_max, dt, tp
  real(kind=dp), dimension(:), allocatable :: syn_r, syn_z, dat_r, dat_z
  real(kind=dp), dimension(2) :: window
  integer :: ier
  type(sachead) :: header_r, header_z
  type(WaveformConvMisfit) :: wccm

  if (command_argument_count() /= 7) then
    write(*,*) 'Usage: ex_wcc_misfit fobs_r fobs_z fsyn_r fsyn_z short_p long_p time_min time_max fadj'
    write(*,*) '  fobs_r: observed radial seismogram in SAC format'
    write(*,*) '  fobs_z: observed vertical seismogram in SAC format' 
    write(*,*) '  fsyn_r: synthetic radial seismogram in SAC format'
    write(*,*) '  fsyn_z: synthetic vertical seismogram in SAC format'
    write(*,*) '  time_min: start time of time window (s)'
    write(*,*) '  time_max: end time of time window (s)'
    write(*,*) '  fadj: output directory for adjoint sources'
    stop
  end if

  ! get command line arguments
  call get_command_argument(1, fobs_r)
  call get_command_argument(2, fobs_z)
  call get_command_argument(3, fsyn_r)
  call get_command_argument(4, fsyn_z)
  call get_command_argument(5, stime_min)
  call get_command_argument(6, stime_max)
  call get_command_argument(7, fadj)
  read(stime_min, *) time_min
  read(stime_max, *) time_max

  ! read observed and synthetic seismograms
  call sacio_readsac(fobs_r, header_r, dat_r, ier)
  if (ier /= 0) then
    write(*,*) 'Error reading observed radial seismogram'
    stop
  end if
  call sacio_readsac(fobs_z, header_z, dat_z, ier)
  if (ier /= 0) then
    write(*,*) 'Error reading observed vertical seismogram'
    stop
  end if
  call sacio_readsac(fsyn_r, header_r, syn_r, ier)
  if (ier /= 0) then
    write(*,*) 'Error reading synthetic radial seismogram'
    stop
  end if
  call sacio_readsac(fsyn_z, header_z, syn_z, ier)
  if (ier /= 0) then
    write(*,*) 'Error reading synthetic vertical seismogram'
    stop
  end if
  dt = dble(header_r%delta)

  ! filter seismograms
  ! cfg%min_period = short_p
  ! cfg%max_period = long_p
  ! call bandpass_dp(syn_r, dt, real(1/long_p), real(1/short_p), IORD)
  ! call bandpass_dp(dat_r, dt, real(1/long_p), real(1/short_p), IORD)
  ! call bandpass_dp(syn_z, dt, real(1/long_p), real(1/short_p), IORD)
  ! call bandpass_dp(dat_z, dt, real(1/long_p), real(1/short_p), IORD)

  ! set measurement window
  tp = maxloc(syn_z, dim=1) * dt
  window(1) = time_min + tp
  window(2) = time_max + tp

  ! calculate adjoint source using waveform convolution misfit
  call wccm%calc_adjoint_source(dat_r, dat_z, syn_r, syn_z, dt, window)

  ! display results
  write(*,'(a,F12.6)') 'Waveform convolution misfit: ', wccm%misfits(1)
  write(*,'(a,F12.6)') 'Waveform convolution residual: ', wccm%residuals(1)

  ! write adjoint sources to SAC files
  ! radial component adjoint source
  call sacio_writesac(trim(fadj)//'/'//trim(header_r%knetwk)//'.'//trim(header_r%kstnm) &
                      //'.BXR.wcc.adj.sac', header_r, wccm%adj_src_r, ier)
  if (ier /= 0) then
    write(*,*) 'Error writing radial adjoint source'
  else
    write(*,*) 'Radial adjoint source written to: ', &
               trim(fadj)//'/'//trim(header_r%knetwk)//'.'//trim(header_r%kstnm)//'.BXR.wcc.adj.sac'
  end if

  ! vertical component adjoint source  
  call sacio_writesac(trim(fadj)//'/'//trim(header_z%knetwk)//'.'//trim(header_z%kstnm) &
                      //'.BXZ.wcc.adj.sac', header_z, wccm%adj_src_z, ier)
  if (ier /= 0) then
    write(*,*) 'Error writing vertical adjoint source'
  else
    write(*,*) 'Vertical adjoint source written to: ', &
               trim(fadj)//'/'//trim(header_z%knetwk)//'.'//trim(header_z%kstnm)//'.BXZ.wcc.adj.sac'
  end if

end program ex_wcc_misfit