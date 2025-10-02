module win_sel
  use config
  use signal
  use cross_correlate
  implicit none

  integer, private, parameter :: max_good_win = 1000

  type :: win_config
    real(kind=dp) :: min_cc = 0.4, max_noise = 0.2, jump_fac = 0.1, min_velocity = 2.4, &
                     threshold_shift_fac = 0.3, threshold_corr = 0.7, min_win_len_fac = 1.5, &
                     max_noise_window = 5.0, max_energy_ratio = 5.0
    integer :: min_peaks_troughs = 3
  end type win_config

  type :: win_sel_type
    integer :: n_win = 0
    real(kind=dp), allocatable :: dat(:)  ! data array for windowing
    real(kind=dp), allocatable :: syn(:)  ! synthetic array for windowing
    real(kind=dp), allocatable :: twin(:,:)   ! window start time (s)
    real(kind=dp), allocatable :: wt(:)   ! window weight
    real(kind=dp), allocatable :: cc_coe(:)   ! cross-correlation coefficient
    real(kind=dp), allocatable :: time_shift(:)   ! time shift
    real(kind=dp), allocatable :: times_cc(:) ! logical array for good windows
    real(kind=dp) :: noise_level, tstart, tend, min_period, max_period, jump_buffer, dt, t0
    integer :: nstart, nend, npts
    integer, allocatable :: win_samp(:,:)
    ! real(kind=dp), allocatable :: snr(:)  ! signal to noise ratio
  contains
    procedure :: init => initialize, gen_good_windows
    procedure, private :: sliding_cc
  end type win_sel_type

  type(win_config) :: win_config_global

contains
  subroutine initialize(this, dat, syn, dt, t0, tp, dis, min_period, max_period)
    class(win_sel_type), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt, t0, tp, dis, min_period, max_period

    this%dat = dat / maxval(abs(syn))
    this%syn = syn / maxval(abs(syn))
    this%dt = dt
    this%t0 = t0 ! t0 is the same as specfem t0
    this%npts = size(dat)
    this%min_period = min_period
    this%max_period = max_period
    this%tstart = tp - min_period / 2.0_dp + t0
    this%tend = dis / win_config_global%min_velocity + min_period / 2.0_dp + t0
    this%jump_buffer = win_config_global%jump_fac * this%min_period
    this%nstart = int(this%tstart / dt) + 1
    this%nend = min(size(dat), int(this%tend / dt) + 1)
    this%noise_level = sqrt(sum(dat(1:this%nstart)**2) / real(this%nstart, kind=dp))

  end subroutine initialize

  subroutine sliding_cc(this)
    class(win_sel_type), intent(inout) :: this
    integer :: nlen, i, nb, ne
    real(kind=dp), allocatable :: dat_win(:), syn_win(:)

    nlen = int(2.0_dp * this%min_period / this%dt) + 1
    allocate(dat_win(nlen), syn_win(nlen))
    allocate(this%time_shift(this%npts - nlen))
    allocate(this%cc_coe(this%npts - nlen))
    allocate(this%times_cc(this%npts - nlen))
    this%time_shift = 0.0_dp
    this%cc_coe = 0.0_dp
    this%times_cc = 0.0_dp
    do i = 1, this%npts - nlen
      nb = i
      ne = i + nlen - 1
      dat_win(:) = this%dat(nb:ne)
      syn_win(:) = this%syn(nb:ne)
      call window_taper(dat_win, 1.0_dp, 1)
      call window_taper(syn_win, 1.0_dp, 1)
      call xcorr_shift_coe(dat_win, syn_win, this%dt, this%time_shift(i), this%cc_coe(i))
      this%times_cc(i) = (dble(i + nlen)/2.0_dp - 1.0_dp) * this%dt
    enddo
    
  end subroutine sliding_cc

  subroutine gen_good_windows(this)
    class(win_sel_type), intent(inout) :: this
    integer :: i, n_jump, n_groups, group_start, group_end, ib, ie, iwin, n_peaks, &
               groups(max_good_win, 2)
    integer, allocatable, dimension(:) :: max_peaks, min_peaks, all_peaks
    logical, allocatable, dimension(:) :: good_windows, good_cc, good_shift, good_time
    real(kind=dp), allocatable :: time_shift_grad(:)
    logical :: in_group, good_groups(max_good_win)
    real(kind=dp) :: max_amp, eng_dat, eng_syn
    

    good_cc = (this%cc_coe > win_config_global%threshold_corr)
    good_shift = (abs(this%time_shift) < win_config_global%threshold_shift_fac * this%min_period)
    good_time = (this%times_cc >= this%tstart) .and. (this%times_cc <= this%tend)
    good_windows = good_cc .and. good_shift .and. good_time

    allocate(time_shift_grad(size(this%time_shift)))
    time_shift_grad = 0.0_dp
    do concurrent (i = 2:size(this%time_shift))
      time_shift_grad(i) = this%time_shift(i) - this%time_shift(i-1)
    end do

    ! Mark windows around large jumps as bad
    do i = 1, size(this%time_shift)
      if (abs(time_shift_grad(i)) > this%jump_buffer) then
        n_jump = nint(this%jump_buffer / this%dt)
        ib = max(1, i - n_jump)
        ie = min(size(this%time_shift), i + n_jump)
        good_windows(ib:ie) = .false.
      end if
    end do

    ! Identify continuous good window segments
    n_groups = 0
    in_group = .false.
    good_groups = .false.
    group_start = 0
    do i = 1, size(this%time_shift)
      if (good_windows(i) .and. .not. in_group) then
        in_group = .true.
        group_start = i
      else if (.not. good_windows(i) .and. in_group) then
        in_group = .false.
        group_end = i - 1
        n_groups = n_groups + 1
        good_groups(n_groups) = .true.
        groups(n_groups, 1) = group_start
        groups(n_groups, 2) = group_end
      end if
    end do

    if (in_group) then
      group_end = size(good_windows)
      n_groups = n_groups + 1
      good_groups(n_groups) = .true.
      groups(n_groups, 1) = group_start
      groups(n_groups, 2) = group_end
    end if

    ! Reject windows
    do i = 1, n_groups
      ! time length too short
      if ((groups(i, 2) - groups(i, 1)) * this%dt < &
           win_config_global%min_win_len_fac * this%min_period) then
        good_groups(i) = .false.
        cycle
      end if

      ! number of peaks too few
      max_peaks = find_maxima(this%syn(groups(i, 1):groups(i, 2)))
      min_peaks = find_maxima(-this%syn(groups(i, 1):groups(i, 2)))
      n_peaks = size(max_peaks) + size(min_peaks)
      max_peaks = find_maxima(this%dat(groups(i, 1):groups(i, 2)))
      min_peaks = find_maxima(-this%dat(groups(i, 1):groups(i, 2)))
      n_peaks = min(n_peaks, size(max_peaks) + size(min_peaks))
      if (n_peaks < win_config_global%min_peaks_troughs) then
        good_groups(i) = .false.
        cycle
      end if

      ! peaks too close
      all_peaks = merge_sorted_arrays(min_peaks, max_peaks)
      do ib = 2, size(all_peaks)
        if ((all_peaks(ib) - all_peaks(ib-1)) * this%dt < 0.5_dp * this%min_period) then
          good_groups(i) = .false.
          exit
        end if
      end do
      if (.not. good_groups(i)) cycle

      ! signal to noise ratio too low
      max_amp = maxval(abs(this%dat(groups(i, 1):groups(i, 2))))
      if (max_amp / this%noise_level < win_config_global%max_noise_window) then
        good_groups(i) = .false.
        cycle
      end if

      ! Energy ratio too low
      eng_dat = sum(this%dat(groups(i, 1):groups(i, 2))**2) / real(groups(i, 2) - groups(i, 1) + 1, kind=dp)
      eng_syn = sum(this%syn(groups(i, 1):groups(i, 2))**2) / real(groups(i, 2) - groups(i, 1) + 1, kind=dp)
      if (eng_dat / eng_syn > win_config_global%max_energy_ratio .or. &
          eng_dat / eng_syn < 1/win_config_global%max_energy_ratio) then
        good_groups(i) = .false.
        cycle
      end if

    end do

    this%n_win = count(good_groups(1:n_groups))
    if (this%n_win == 0) then
      print *, "Warning: No good windows found."
      return
    end if
    allocate(this%twin(this%n_win, 2))
    allocate(this%win_samp(this%n_win, 2))
    iwin = 0
    do i = 1, n_groups
      if (good_groups(i)) then
        iwin = iwin + 1
        this%twin(iwin, :) = groups(i, :) * this%dt
        this%win_samp(iwin, :) = groups(i, :)
      end if
    end do

  end subroutine gen_good_windows

  subroutine xcorr_shift_coe(d, s, dt, tshift, cc_max_coef)
    real(kind=dp), dimension(:), intent(in) :: s, d
    real(kind=dp), intent(in) :: dt
    real(kind=dp), intent(out) :: cc_max_coef, tshift
    real(kind=dp), dimension(:), allocatable :: s_shift, d_shift
    integer :: nlen, ishift, index, index_shift

    nlen = size(s)
    allocate(s_shift(nlen))
    s_shift = 0.0_dp
    d_shift = d(:)

    ishift = xcorr_shift(d, s)
    do index = 1, nlen
      index_shift = index - ishift
      if (index_shift >= 1 .and. index_shift <= nlen) then
        s_shift(index) = s(index_shift)
      end if
    end do
    tshift = real(ishift, kind=dp) * dt

    ! taper the 0 edges on d
    if (ishift > 0) then
      d_shift(1:ishift) = 0.0_dp
    else if (ishift < 0) then
      d_shift(nlen+ishift+1:nlen) = 0.0_dp
    end if

    cc_max_coef = sum(d_shift * s_shift) / sqrt(sum(d_shift**2) * sum(s_shift**2))

  end subroutine xcorr_shift_coe

  function merge_sorted_arrays(arr1, arr2) result(merged)
    integer, intent(in) :: arr1(:), arr2(:)
    integer, allocatable :: merged(:)
    integer :: n1, n2, n_total
    integer :: i1, i2, idx
    
    n1 = size(arr1)
    n2 = size(arr2)
    n_total = n1 + n2
    
    allocate(merged(n_total))
    
    i1 = 1
    i2 = 1
    idx = 1
    
    ! merge the two arrays
    do while (i1 <= n1 .and. i2 <= n2)
      if (arr1(i1) <= arr2(i2)) then
        merged(idx) = arr1(i1)
        i1 = i1 + 1
      else
        merged(idx) = arr2(i2)
        i2 = i2 + 1
      end if
      idx = idx + 1
    end do
    
    ! append remaining elements
    do while (i1 <= n1)
      merged(idx) = arr1(i1)
      i1 = i1 + 1
      idx = idx + 1
    end do
    
    do while (i2 <= n2)
      merged(idx) = arr2(i2)
      i2 = i2 + 1
      idx = idx + 1
    end do
    
  end function merge_sorted_arrays

end module win_sel