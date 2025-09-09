module mt_tt_misfit
  use config
  use signal
  use adj_config
  use fftpack

  implicit none

  type, extends(AdjointMeasurement) :: MTTTMisfit
    real(kind=dp), private, dimension(:), allocatable :: tshift, dlna, sigma_dt, sigma_dlna, &
                                            misfit_p, misfit_q, cc_max
    real(kind=dp) :: dt, tlen
    integer :: nlen_f, nlen
    logical :: is_mtm
  contains
    procedure :: calc_adjoint_source
    procedure, private :: initialize, check_time_series_acceptability, &
                          prepare_data_for_mtm, calculate_freq_limits, &
                          calculate_multitaper
  end type MTTTMisfit

  type(fft_cls), private :: fft_obj
contains

  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:,:), allocatable :: tapers
    real(kind=dp), dimension(:), allocatable :: s, d, adj_tw_p, adj_tw_q, freq, wvec, &
                                                ey1, ey2, phi_mtm, abs_mtm, dtau_mtm, dlna_mtm
    real(kind=dp) :: df
    integer :: iwin, nb, ne, nlen_win, nfreq_min, nfreq_max
    logical :: is_mtm
    type(fft_cls) :: fftins


    ! num of measurements
    if (size(windows, 2) /= 2) then
      write(*,*) 'Error: windows must have two columns (start and end times)'
      error stop
    end if
    ! allocate measurement arrays
    call this%initialize()

    this%nwin = size(windows, 1)

    this%nlen_f = 2 ** lnpt
    this%dt = dt
    this%nlen = size(dat)
    this%tlen = this%nlen * dt
    if (allocated(this%adj_src)) deallocate(this%adj_src)
    allocate(this%adj_src(this%nlen))
    this%adj_src = 0.0_dp

    !loop over windows
    do iwin = 1, this%nwin
      is_mtm = .true.

      ! trim the windows to the length of the data
      call get_window_info(windows(iwin,:), dt, nb, ne, nlen_win)
      s = syn(nb:ne)
      d = dat(nb:ne)

      ! taper the windows
      call window_taper(s, taper_percentage, itaper_type)
      call window_taper(d, taper_percentage, itaper_type)

      ! calculate cross-correlation shift
      call calc_cc_shift(d, s, dt, dt_sigma_min, dlna_sigma_min, &
                         this%tshift(iwin), this%dlna(iwin), &
                         this%sigma_dt(iwin), this%sigma_dlna(iwin))

      do while (is_mtm)
        ! prepare data for MTM
        is_mtm = this%check_time_series_acceptability(iwin, nlen_win)
        if (.not. is_mtm) exit

        call this%prepare_data_for_mtm(d, iwin, windows(iwin,:), is_mtm)
        if (.not. is_mtm) exit

        freq = fftins%fftfreq(this%nlen_f, dt)
        df = freq(2) - freq(1)
        wvec = 2.0_dp * PI * freq

        ! calculate frequency limits
        call this%calculate_freq_limits(df, nfreq_min, nfreq_max, is_mtm)
        if (.not. is_mtm) exit

        ! calculate multitaper transfer function and measurements
        call staper(nlen_win, num_taper, tapers, this%nlen_f, ey1, ey2)
        tapers = tapers * dsqrt(dble(nlen_win))
        call this%calculate_multitaper(d, s, tapers, wvec, nfreq_min, nfreq_max, &
                                      this%tshift(iwin), this%dlna(iwin), &
                                      phi_mtm, abs_mtm, dtau_mtm, dlna_mtm)

      end do
    end do
  end subroutine calc_adjoint_source

  subroutine calculate_multitaper(this, d, s, tapers, wvec, nfreq_min, nfreq_max, cc_tshift, cc_dlna,&
                                  phi_w, abs_w, dtau_w, dlna_w)
    ! Calculate the multitaper transfer function and measurements
    ! Arguments:
    ! d, s: time series of data and synthetics
    ! tapers: array of tapers (nlen_t, ntaper)
    ! wvec: angular frequency vector (nlen_f)
    ! nfreq_min, nfreq_max: frequency indices for analysis
    ! cc_tshift, cc_dlna: cross-correlation time shift and log amplitude ratio
    ! Outputs:
    ! phi_w: phase of transfer function
    ! abs_w: amplitude of transfer function
    ! dtau_w: derivative of phase with respect to frequency
    ! dlna_w: derivative of log amplitude with respect to frequency
    
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: d, s, wvec
    real(kind=dp), dimension(:,:), intent(in) :: tapers
    integer, intent(in) :: nfreq_min, nfreq_max
    real(kind=dp), intent(in) :: cc_tshift, cc_dlna
    real(kind=dp), dimension(:), allocatable, intent(out) :: phi_w, abs_w, dtau_w, dlna_w
    complex(kind=dp), dimension(:), allocatable :: top_tf, bot_tf, d_tw, s_tw, trans_func
    real(kind=dp), dimension(:), allocatable :: taper, d_t, s_t
    real(kind=dp) :: wlevel
    integer :: nlen_t, ntaper, fnum, itaper, i
    type(fft_cls) :: fftins
    
    nlen_t = size(d)
    ntaper = size(tapers, 2)
    fnum = int(this%nlen_f / 2) + 1

    allocate(top_tf(this%nlen_f), bot_tf(this%nlen_f))
    top_tf = (0.0_dp, 0.0_dp)
    bot_tf = (0.0_dp, 0.0_dp)

    do itaper = 1, ntaper
      taper = tapers(1:nlen_t, itaper)

      ! Apply taper to data and synthetics
      d_t = d * taper
      s_t = s * taper

      d_tw = fftins%fft(d_t, this%nlen_f) * this%dt
      s_tw = fftins%fft(s_t, this%nlen_f) * this%dt

      ! Accumulate the top and bottom terms for multitaper
      top_tf = top_tf + d_tw * conjg(s_tw)
      bot_tf = bot_tf + s_tw * conjg(s_tw)
    enddo

    ! Calculate the transfer function with water level stabilization
    wlevel = transfunc_waterlevel**2 * maxval(abs(bot_tf(1:fnum)))

    ! Calculate the transfer function
    do i = nfreq_min, nfreq_max
      if (abs(bot_tf(i)) >= wlevel) then
        trans_func(i) = top_tf(i) / bot_tf(i)
      else
        trans_func(i) = top_tf(i) / (bot_tf(i) + wlevel)
      end if
    enddo

    allocate(phi_w(this%nlen_f), abs_w(this%nlen_f), dtau_w(this%nlen_f), dlna_w(this%nlen_f))
    phi_w = 0.0_dp
    abs_w = 0.0_dp
    dtau_w = 0.0_dp
    dlna_w = 0.0_dp

    ! Calculate phase, amplitude, and their derivatives
    phi_w(nfreq_min:nfreq_max) = atan2(aimag(trans_func(nfreq_min:nfreq_max)), real(trans_func(nfreq_min:nfreq_max)))
    call process_cycle_skipping(phi_w, nfreq_min, nfreq_max, wvec(1:this%nlen_f), this%nlen_f)

    ! Calculate amplitude
    abs_w(nfreq_min:nfreq_max) = abs(trans_func(nfreq_min:nfreq_max))

    ! Add CC measurement to transfer function
    dtau_w(1) = cc_tshift
    dtau_w(max(2, nfreq_min):nfreq_max) = &
      - 1.0_dp / wvec(max(2, nfreq_min):nfreq_max) * &
      phi_w(max(2, nfreq_min):nfreq_max) + cc_tshift

    dlna_w(nfreq_min:nfreq_max) = log(abs_w(nfreq_min:nfreq_max)) + cc_dlna

  end subroutine calculate_multitaper
    
  function check_time_series_acceptability(this, iwin, nlen_w) result(is_acceptable)
    class(MTTTMisfit), intent(in) :: this
    integer, intent(in) :: nlen_w, iwin
    logical :: is_acceptable

    ! Check if the time shift is within acceptable limits
    if (this%tshift(iwin) <= this%dt) then
      is_acceptable = .false.
    elseif (min_cycle_in_window * min_period > nlen_w) then
      is_acceptable = .false.
    else
      is_acceptable = .true.
    end if
  end function check_time_series_acceptability

  subroutine prepare_data_for_mtm(this, d, iwin, window, is_acceptable)
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(inout) :: d
    real(kind=dp), dimension(:), intent(in) :: window
    integer, intent(in) :: iwin
    integer :: nb, ne, ishift, nb_d, ne_d, nlen_d, nlen_w
    logical, intent(out) :: is_acceptable

    ! Prepare data for MTM analysis
    call get_window_info(window, this%dt, nb, ne, nlen_w)
    ishift = int(this%tshift(iwin) / this%dt)

    nb_d = max(1, nb + ishift)
    ne_d = min(this%nlen, ne + ishift)
    nlen_d = ne_d - nb_d + 1

    if (nlen_d == nlen_w) then
      ! If the data length matches the window length, use it directly
      d(1:nlen_w) = d(nb_d:ne_d)
      d = d * exp(-this%dlna(iwin))
      call window_taper(d, taper_percentage, itaper_type)
      is_acceptable = .true.
    else
      is_acceptable = .false.
    end if
  end subroutine prepare_data_for_mtm

  subroutine calculate_freq_limits(this, df, nfreq_min, nfreq_max, is_acceptable)
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), intent(in) :: df
    integer, intent(out) :: nfreq_min, nfreq_max
    logical, intent(out) :: is_acceptable
    complex(kind=dp), dimension(:), allocatable :: s_spec
    real(kind=dp) :: ampmax, scaled_wl, half_taper_bandwidth, &
                     chosen_bandwidth
    integer :: fnum, i_ampmax, ifreq_min, ifreq_max, nlen_win, iw
    logical :: is_search

    ! calculate frequency limits for MTM analysis
    fnum = int(this%nlen_f / 2) + 1
    s_spec = fft_obj%fft(this%adj_src, this%nlen_f) * this%dt

    ! Calculate the frequency limits based on the sampling rate and window length
    ampmax = maxval(abs(s_spec))
    i_ampmax = maxloc(abs(s_spec), dim=1)

    ! Scale the maximum amplitude to the water threshold
    scaled_wl = water_threshold * ampmax

    ! Find the frequency index corresponding to the maximum amplitude
    ifreq_min = int(1.0_dp / (max_period * df))
    ifreq_max = int(1.0_dp / (min_period * df))

    ! get the frequency limits
    nfreq_max = fnum
    is_search = .true.
    do iw = 1, fnum
      if (iw <= i_ampmax) cycle
      call search_frequency_limit(is_search, iw, nfreq_max, s_spec, scaled_wl, 10)
    end do
    ! Make sure `nfreq_max` does not go beyond the Nyquist frequency
    nfreq_max = min(nfreq_max, ifreq_max, int(1.0_dp / (2.0_dp * this%dt) / df))

    nfreq_min = 1
    is_search = .true.
    do iw = fnum, 1, -1
      if (iw >= i_ampmax) cycle
      call search_frequency_limit(is_search, iw, nfreq_min, s_spec, scaled_wl, 10)
    end do
    ! Make sure `nfreq_min` does not go below the minimum frequency
    nfreq_min = max(nfreq_min, ifreq_min, int(min_cycle_in_window / this%tlen / df))

    half_taper_bandwidth = mt_nw / (4.0_dp * this%tlen)
    chosen_bandwidth = (nfreq_max - nfreq_min) * df
    if (chosen_bandwidth < half_taper_bandwidth) then
      is_acceptable = .false.
    else
      is_acceptable = .true.
    end if

  end subroutine calculate_freq_limits

  subroutine search_frequency_limit(is_search, index, nfreq_limit, spectra, water_threshold, c)
    integer, intent(in) :: index          ! index of the frequency to check
    integer, intent(in) :: c              ! scaling factor
    integer, intent(inout) :: nfreq_limit ! frequency limit index
    real(kind=dp), intent(in) :: water_threshold ! water threshold value
    complex(kind=dp), intent(in) :: spectra(:)     ! spectrum data
    logical, intent(inout) :: is_search   ! search flag

    if (abs(spectra(index)) < water_threshold .and. is_search) then
      is_search   = .false.
      nfreq_limit = index
    end if

    if (abs(spectra(index)) > c * water_threshold .and. .not. is_search) then
      is_search   = .true.
      nfreq_limit = index
    end if

  end subroutine search_frequency_limit

  subroutine initialize(this)
    class(MTTTMisfit), intent(inout) :: this

    ! deallocate arrays
    if (allocated(this%adj_src)) deallocate(this%adj_src)
    if (allocated(this%misfits)) deallocate(this%misfits)
    if (allocated(this%tshift)) deallocate(this%tshift)
    if (allocated(this%dlna)) deallocate(this%dlna)
    if (allocated(this%sigma_dt)) deallocate(this%sigma_dt)
    if (allocated(this%sigma_dlna)) deallocate(this%sigma_dlna)
    if (allocated(this%misfit_p)) deallocate(this%misfit_p)
    if (allocated(this%misfit_q)) deallocate(this%misfit_q)
    if (allocated(this%cc_max)) deallocate(this%cc_max)
    if (allocated(this%select_meas)) deallocate(this%select_meas)
    if (allocated(this%imeas)) deallocate(this%imeas)
    if (allocated(this%residuals)) deallocate(this%residuals)
    if (allocated(this%errors)) deallocate(this%errors)

    allocate(this%misfits(this%nwin))
    allocate(this%tshift(this%nwin))
    allocate(this%dlna(this%nwin))
    allocate(this%sigma_dt(this%nwin))
    allocate(this%sigma_dlna(this%nwin))
    allocate(this%misfit_p(this%nwin))
    allocate(this%misfit_q(this%nwin))
    allocate(this%cc_max(this%nwin))
    allocate(this%select_meas(this%nwin))
    allocate(this%imeas(this%nwin))
    allocate(this%residuals(this%nwin))
    allocate(this%errors(this%nwin))
    this%select_meas = .true.
    this%misfits = 0.0_dp
    this%total_misfit = 0.0_dp

  end subroutine initialize
end module mt_tt_misfit