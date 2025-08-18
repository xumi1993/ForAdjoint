module adj_config
  use config
  use signal
  implicit none

  type :: AdjointMeasurement
    real(kind=dp), dimension(:), allocatable :: adj_src
    real(kind=dp), dimension(:,:), allocatable :: window_chi
    real(kind=dp), dimension(:), allocatable :: misfits, residuals, errors
    real(kind=dp) :: total_misfit
    integer :: nwin
    integer, dimension(:), allocatable :: imeas
    logical, dimension(:), allocatable :: select_meas

  end type AdjointMeasurement

  real(kind=dp) :: min_period, max_period, taper_percentage=0.3_dp, &
                    dt_sigma_min=1.0_dp, dlna_sigma_min=0.5_dp, transfunc_waterlevel=1e-10_dp, &
                    water_threshold=0.02_dp, mt_nw=4.0_dp, phase_step=1.5_dp, &
                    dt_fac=2.0_dp, err_fac=2.5_dp, dt_max_scale=0.5_dp, TSHIFT_MIN=-5.0_dp, & 
                    TSHIFT_MAX=5.0_dp, DLNA_MIN=-1.5_dp, DLNA_MAX=1.5_dp, CC_MIN=0.8_dp
  integer :: itaper_type=4, imeasure_type=1, lnpt=15, min_cycle_in_window=3, &
              num_taper=5
  ! itaper_type: 1=Hanning, 2=Hamming, 3=cos, 4=cos^10
  ! imeasure_type: 1=Wavefrom-diff, 2=CC-TT, 3=CC-DLNA, 4=CC-TT-MT, 5=CC-DLNA-MT, 11=RF

end module adj_config