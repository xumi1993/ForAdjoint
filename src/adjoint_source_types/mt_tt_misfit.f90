module mt_tt_misfit
  use config
  use signal
  use adj_config

  implicit none

  type, extends(AdjointMeasurement) :: MTTTMisfit
    real(kind=dp), private, dimension(:), allocatable :: tshift, dlna, sigma_dt, sigma_dlna, &
                                            misfit_p, misfit_q, cc_max
    real(kind=dp) :: dt
    integer :: nlen_f
    logical :: is_mtm
  contains
    procedure :: calc_adjoint_source, check_time_series_acceptability
    procedure, private :: initialize
  end type MTTTMisfit

contains

  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:), allocatable :: s, d, adj_tw_p, adj_tw_q
    integer :: iwin, nlen, nb, ne, nlen_win
    logical :: is_mtm


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
    nlen = size(dat)
    if (allocated(this%adj_src)) deallocate(this%adj_src)
    allocate(this%adj_src(nlen))
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
    end do
  end subroutine calc_adjoint_source

  function check_time_series_acceptability(this, iwin, nlen_win) result(is_acceptable)
    class(MTTTMisfit), intent(in) :: this
    integer, intent(in) :: nlen_win, iwin
    logical :: is_acceptable

    ! Check if the time shift is within acceptable limits
    if (this%tshift(iwin) <= this%dt) then
      is_acceptable = .false.
    elseif (min_cycle_in_window * min_period > nlen_win) then
      is_acceptable = .false.
    else
      is_acceptable = .true.
    end if
  end function check_time_series_acceptability

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