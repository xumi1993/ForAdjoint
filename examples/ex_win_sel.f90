program ex_win_sel
  use config
  use sacio
  use signal
  use win_sel
  use travel_times_mod
  implicit none

  type(sachead) :: obs_head, syn_head
  type(win_sel_type) :: win
  character(len=256) :: obs_file, syn_file, output_file
  character(len=256) :: smin_period, smax_period, sis_split_phases
  real(kind=dp), allocatable :: obs_data(:), syn_data(:)
  real(kind=dp) :: dt, t0, tp, dis, evdp, gcarc
  real(kind=dp) :: min_period, max_period
  integer :: i, npts, ier, funit
  logical :: is_split_phases

  ! Check command line arguments
  if (command_argument_count() /= 5) then
    print *, 'Usage: xex_win_sel obs_file syn_file min_period max_period is_split_phases'
    print *, '  obs_file: Observed seismogram in SAC format'
    print *, '  syn_file: Synthetic seismogram in SAC format'
    print *, '  min_period: Minimum period for bandpass filter (s)'
    print *, '  max_period: Maximum period for bandpass filter (s)'
    print *, '  is_split_phases: Split windows by phases (true/false)'
    print *, ''
    print *, 'Example:'
    print *, '  ./bin/xex_win_sel example_data/obs.II.ABKT.LHZ.sac \'
    print *, '                    example_data/syn.II.ABKT.LHZ.sac 50 150'
    stop
  end if

  ! Get command line arguments
  call get_command_argument(1, obs_file)
  call get_command_argument(2, syn_file)
  call get_command_argument(3, smin_period)
  call get_command_argument(4, smax_period)
  read(smin_period, *) min_period
  read(smax_period, *) max_period
  call get_command_argument(5, sis_split_phases)
  read(sis_split_phases, *) is_split_phases

  ! Configure window selection parameters
  win_config_global%threshold_corr = 0.8_dp        ! cc_threshold
  win_config_global%threshold_shift_fac = 0.3_dp   ! time_shift_threshold = 15s = 0.3*min_period
  win_config_global%jump_fac = 0.1_dp              ! jump_buffer = 0.1 * min_period
  win_config_global%min_velocity = 2.4_dp          ! min_velocity (km/s)
  win_config_global%min_win_len_fac = 1.5_dp       ! min_length_period = 1.5 * min_period
  win_config_global%min_peaks_troughs = 3          ! min_peaks_troughs
  win_config_global%max_noise_window = 5.0_dp      ! max_noise_window
  win_config_global%max_energy_ratio = 10.0_dp     ! max_energy_ratio
  win_config_global%is_split_phases = is_split_phases ! is_split_phases

  print *, '=================================='
  print *, 'Window Selection Test Program'
  print *, '=================================='
  print *, ''
  print *, 'Parameters:'
  print *, '  Min period: ', min_period, ' s'
  print *, '  Max period: ', max_period, ' s'
  print *, ''
  print *, 'Reading SAC files...'
  print *, '  Observed: ', trim(obs_file)
  print *, '  Synthetic: ', trim(syn_file)

  call sacio_readsac(obs_file, obs_head, obs_data, ier)
  if (ier /= 0) then
    print *, 'Error reading observed SAC file:', trim(obs_file)
    stop
  end if

  call sacio_readsac(syn_file, syn_head, syn_data, ier)
  if (ier /= 0) then
    print *, 'Error reading synthetic SAC file:', trim(syn_file)
    stop
  end if

  ! Get basic parameters
  npts = obs_head%npts
  dt = dble(obs_head%delta)
  t0 = -dble(obs_head%b)
  evdp = dble(obs_head%evdp)
  gcarc = dble(obs_head%gcarc)
  dis = dble(obs_head%dist)

  print *, ''
  print *, 'SAC file information:'
  print *, '  npts  = ', npts
  print *, '  dt    = ', dt, ' s'
  print *, '  t0 (b)= ', t0, ' s'
  print *, '  evdp  = ', evdp, ' km'
  print *, '  gcarc = ', gcarc, ' deg'
  print *, '  dist  = ', dis, ' km'

  ! Allocate data arrays (already allocated by sacio_readsac, just assign)
  ! obs_data and syn_data are already allocated by sacio_readsac


  ! Bandpass filter
  print *, ''
  print *, 'Applying bandpass filter...'
  print *, '  Period range: ', min_period, ' - ', max_period, ' s'
  print *, '  Frequency range: ', 1.0_dp/max_period, ' - ', 1.0_dp/min_period, ' Hz'
  
  call bandpass_dp(obs_data, dt, real(1.0_dp/max_period, kind=cr), &
                   real(1.0_dp/min_period, kind=cr), 2)
  call bandpass_dp(syn_data, dt, real(1.0_dp/max_period, kind=cr), &
                   real(1.0_dp/min_period, kind=cr), 2)

  ! Calculate first arrival time using velocity model
  print *, ''
  print *, 'Calculating first arrival time using AK135 model...'
  
  ! Initialize velocity model (AK135)
!   call vmodel%init('AK135')
  
  ! Calculate travel time from source to receiver
  ! get_travel_time(dsrc, drec, gcarc_rec) returns travel time array
  print*, '  Source depth (evdp): ', evdp, ' km'
  print*, '  Great circle distance (gcarc): ', gcarc, ' deg'
!   ttime = vmodel%get_travel_time(evdp, [0.0_dp], [gcarc])
  block 
    integer :: nphases
    character(len=8), dimension(:), allocatable :: names
    real(kind=dp), dimension(:), allocatable :: times
    call ttimes(gcarc,evdp,nphases,names,times)
    tp = times(1)
  end block
  
  print *, '  Source depth: ', evdp, ' km'
  print *, '  Great circle distance: ', gcarc, ' deg'
  print *, '  Calculated first arrival time (tp): ', tp, ' s'

  ! Initialize window selector
  print *, ''
  print *, 'Initializing window selector...'
  print *, '  Window length: ', 2.0_dp * min_period, ' s'
  
  call win%init(obs_data, syn_data, dt, t0, tp, dis, min_period)

  print *, '  tstart = ', win%tstart, ' s'
  print *, '  tend   = ', win%tend, ' s'
  print *, '  Noise level = ', win%noise_level

  ! Generate good windows
  print *, ''
  print *, 'Generating good windows...'
  print *, 'Window selection criteria:'
  print *, '  CC threshold        : ', win_config_global%threshold_corr
  print *, '  Time shift threshold: ', win_config_global%threshold_shift_fac * min_period, ' s'
  print *, '  Jump buffer         : ', win%jump_buffer, ' s'
  print *, '  Min window length   : ', win_config_global%min_win_len_fac * min_period, ' s'
  print *, '  Min peaks/troughs   : ', win_config_global%min_peaks_troughs
  print *, '  Max noise ratio     : ', win_config_global%max_noise_window
  print *, '  Max energy ratio    : ', win_config_global%max_energy_ratio
  print *, ''
  
  call win%gen_good_windows()

  ! Print results
  print *, '=================================='
  print *, 'Window Selection Results'
  print *, '=================================='
  print *, ''
  print *, 'Number of good windows found: ', win%n_win
  print *, ''

  if (win%n_win > 0) then
    print *, 'Window details:'
    print *, '  Win#    Start(s)      End(s)   Duration(s)   Avg CC   Avg Shift(s)'
    print *, '  ----------------------------------------------------------------'
    
    do i = 1, win%n_win
      print '(I5, 2F12.2, F12.2, F10.3, F12.3)', i, &
            win%twin(i, 1), win%twin(i, 2), &
            win%twin(i, 2) - win%twin(i, 1), &
            sum(win%cc_coe(win%win_samp(i,1):win%win_samp(i,2))) / &
              real(win%win_samp(i,2) - win%win_samp(i,1) + 1, kind=dp), &
            sum(win%time_shift(win%win_samp(i,1):win%win_samp(i,2))) / &
              real(win%win_samp(i,2) - win%win_samp(i,1) + 1, kind=dp)
    end do
    
    print *, ''
    print *, 'Summary statistics:'
    print *, '  Total window duration: ', sum(win%twin(:, 2) - win%twin(:, 1)), ' s'
    print *, '  Coverage: ', sum(win%twin(:, 2) - win%twin(:, 1)) / (win%tend - win%tstart) * 100.0_dp, ' %'
  else
    print *, 'No good windows found!'
    print *, 'Possible reasons:'
    print *, '  - CC threshold too high'
    print *, '  - Time shift threshold too low'
    print *, '  - Signal-to-noise ratio too low'
    print *, '  - Energy ratio outside acceptable range'
  end if

  print *, ''
  print *, '=================================='
  print *, 'Test completed successfully!'
  print *, '=================================='

  ! Write window information to text files for visualization
  output_file = 'window_results.txt'
  print *, ''
  print *, 'Writing window results to: ', trim(output_file)
  
  open(newunit=funit, file=output_file, status='replace', action='write')
  write(funit, '(A)') '# Window Selection Results'
  write(funit, '(A)') '# Format: window_number  start_time(s)  end_time(s)  duration(s)  avg_cc  avg_shift(s)'
  write(funit, '(A,I0)') '# Number of windows: ', win%n_win
  write(funit, '(A,F10.2)') '# tstart: ', win%tstart
  write(funit, '(A,F10.2)') '# tend: ', win%tend
  write(funit, '(A,F10.4)') '# noise_level: ', win%noise_level
  write(funit, '(A)')
  
  if (win%n_win > 0) then
    do i = 1, win%n_win
      write(funit, '(I5, 2F12.2, F12.2, F10.3, F12.3)') i, &
            win%twin(i, 1), win%twin(i, 2), &
            win%twin(i, 2) - win%twin(i, 1), &
            sum(win%cc_coe(win%win_samp(i,1):win%win_samp(i,2))) / &
              real(win%win_samp(i,2) - win%win_samp(i,1) + 1, kind=dp), &
            sum(win%time_shift(win%win_samp(i,1):win%win_samp(i,2))) / &
              real(win%win_samp(i,2) - win%win_samp(i,1) + 1, kind=dp)
    end do
  end if
  close(funit)

  ! Write detailed CC and time shift data for visualization
  output_file = 'window_cc_shift.txt'
  print *, 'Writing CC and shift data to: ', trim(output_file)
  
  open(newunit=funit, file=output_file, status='replace', action='write')
  write(funit, '(A)') '# Sliding window cross-correlation and time shift results'
  write(funit, '(A)') '# Format: time(s)  cc_coefficient  time_shift(s)'
  do i = 1, size(win%times_cc)
    write(funit, '(F12.4, F10.4, F12.4)') win%times_cc(i), win%cc_coe(i), win%time_shift(i)
  end do
  close(funit)
  
  print *, ''
  print *, 'Output files created successfully!'
  print *, 'You can now visualize the results using the Jupyter notebook.'

  ! Clean up
  deallocate(obs_data, syn_data)

end program ex_win_sel
