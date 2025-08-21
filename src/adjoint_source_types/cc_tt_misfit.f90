module cc_tt_misfit
  use config
  use signal
  use adj_config
  use cross_correlate
  implicit none

  type, extends(AdjointMeasurement) :: CCTTMisfit
    real(kind=dp), private, dimension(:), allocatable :: tshift, dlna, sigma_dt, sigma_dlna, &
                                            misfit_p, misfit_q, cc_max
  contains
    procedure :: calc_adjoint_source, cc_measure_select
    procedure, private :: initialize
  end type CCTTMisfit

contains
  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(CCTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:), allocatable :: s, d, adj_tw_p, adj_tw_q
    integer :: iwin, nlen, nb, ne, nlen_win

    ! num of measurements
    if (size(windows, 2) /= 2) then
      write(*,*) 'Error: windows must have two columns (start and end times)'
      error stop
    end if
    this%nwin = size(windows, 1)
    ! allocate measurement arrays
    call this%initialize()
  
    ! number of time samples
    nlen = size(dat)
    if (allocated(this%adj_src)) deallocate(this%adj_src)
    allocate(this%adj_src(nlen))
    this%adj_src = 0.0_dp

    ! loop over windows
    do iwin = 1, this%nwin
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
      
      ! select measurements based on criteria
      this%cc_max(iwin) = cc_max_coef(d, s)
      call this%cc_measure_select(iwin)
      if (.not. this%select_meas(iwin)) cycle

      ! calculate adjoint source
      call calc_cc_adjsrc(s, this%tshift(iwin), this%dlna(iwin), dt, &
                          this%sigma_dt(iwin), this%sigma_dlna(iwin), &
                          this%misfit_p(iwin), this%misfit_q(iwin), &
                          adj_tw_p, adj_tw_q)

      ! add windwow to adjoint source
      call window_taper(adj_tw_p, taper_percentage, itaper_type)
      call window_taper(adj_tw_q, taper_percentage, itaper_type)

      select case (imeasure_type)
        case(2) ! CC-TT
          this%imeas(iwin) = IMEAS_CC_TT ! CC-TT
          this%total_misfit = this%total_misfit + this%misfit_p(iwin)
          this%misfits(iwin) = this%misfit_p(iwin)
          this%residuals(iwin) = this%tshift(iwin)
          this%errors(iwin) = this%sigma_dt(iwin)
          this%adj_src(nb:ne) = adj_tw_p
        case(3) ! CC-DLNA
          this%imeas(iwin) = IMEAS_CC_DLNA ! CC-DLNA
          this%total_misfit = this%total_misfit + this%misfit_q(iwin)
          this%misfits(iwin) = this%misfit_q(iwin)
          this%residuals(iwin) = this%dlna(iwin)
          this%errors(iwin) = this%sigma_dlna(iwin)
          this%adj_src(nb:ne) = adj_tw_q
      end select

    end do

  end subroutine calc_adjoint_source

  subroutine initialize(this)
    class(CCTTMisfit), intent(inout) :: this

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

  subroutine cc_measure_select(this, iwin)
    class(CCTTMisfit), intent(inout) :: this
    real(kind=dp) :: cc_max
    integer, intent(in) :: iwin

    ! cc_max = cc_max_coef(d, s)
    if ((this%cc_max(iwin) < CC_MIN) .or. (this%tshift(iwin) < TSHIFT_MIN) &
                                     .or. (this%tshift(iwin) > TSHIFT_MAX) &
                                     .or. (this%dlna(iwin) < DLNA_MIN) &
                                     .or. (this%dlna(iwin) > DLNA_MAX)) then
      this%select_meas(iwin) = .false.
      this%misfit_p(iwin) = 0.0_dp
      this%misfit_q(iwin) = 0.0_dp
      this%misfits(iwin) = 0.0_dp
    end if
  end subroutine cc_measure_select

end module cc_tt_misfit