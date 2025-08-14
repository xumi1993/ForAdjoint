module waveform_misfit
  use config
  use signal
  implicit none

  type, extends(AdjointMeasurement) :: WaveformMisfit
  contains
    procedure :: calc_adjoint_source
  end type WaveformMisfit

contains
  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(WaveformMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:), allocatable :: s, d, adj_tw
    real(kind=dp) :: misfit
    integer :: iwin, nlen, nb, ne, nlen_win

    ! num of measurements
    if (size(windows, 2) /= 2) then
      write(*,*) 'Error: windows must have two columns (start and end times)'
      error stop
    end if
    this%nwin = size(windows, 1)
    allocate(this%misfits(this%nwin))
    this%misfits = 0.0_dp
    this%total_misfit = 0.0_dp

    ! number of time samples
    nlen = size(dat)
    allocate(this%adj_src(nlen))
    this%adj_src = 0.0_dp

    ! loop over windows
    do iwin = 1, this%nwin
      this%imeas = 1 ! Waveform diff
      call get_window_info(windows(iwin,:), dt, nb, ne, nlen_win)
      s = syn(nb:ne)
      d = dat(nb:ne)

      ! taper the windows
      call window_taper(s, taper_percentage, itaper_type)
      call window_taper(d, taper_percentage, itaper_type)

      adj_tw = s - d

      call window_taper(adj_tw, taper_percentage, itaper_type)

      ! calculate waveform misfit and adjoint source
      misfit = 0.5_dp * simpson((s - d)**2) * dt
      this%misfits(iwin) = misfit
      this%total_misfit = this%total_misfit + misfit
      this%adj_src(nb:ne) = adj_tw

    end do

  end subroutine calc_adjoint_source

end module waveform_misfit