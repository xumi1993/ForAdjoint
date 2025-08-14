module cc_tt_misfit
  use config
  use signal
  use adj_config
  use cross_correlate
  implicit none

  type, extends(AdjointMeasurement) :: CCTTMisfit
    real(kind=dp), dimension(:), allocatable :: tshift, dlna, sigma_dt, sigma_dlna, &
                                                misfit_p, misfit_q
  contains
    procedure :: calc_adjoint_source
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
    ! allocate windows
    ! allocate(this%window_chi(this%nwin, NCHI))
    ! this%window_chi = 0.0_dp
    allocate(this%misfits(this%nwin))
    allocate(this%tshift(this%nwin))
    allocate(this%dlna(this%nwin))
    allocate(this%sigma_dt(this%nwin))
    allocate(this%sigma_dlna(this%nwin))
    allocate(this%misfit_p(this%nwin))
    allocate(this%misfit_q(this%nwin))
    this%total_misfit = 0.0_dp

    ! number of time samples
    nlen = size(dat)
    allocate(this%adj_src(nlen))
    this%adj_src = 0.0_dp

    ! loop over windows
    do iwin = 1, this%nwin
      this%imeas = 2 ! CC-TT
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
          this%total_misfit = this%total_misfit + this%misfit_p(iwin)
          this%misfits(iwin) = this%misfit_p(iwin)
          this%adj_src(nb:ne) = adj_tw_p
        case(3) ! CC-DLNA
          this%total_misfit = this%total_misfit + this%misfit_q(iwin)
          this%misfits(iwin) = this%misfit_q(iwin)
          this%adj_src(nb:ne) = adj_tw_q
      end select

    end do

  end subroutine calc_adjoint_source

end module cc_tt_misfit