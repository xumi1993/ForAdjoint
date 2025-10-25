module vpmodel
  use config
  use signal, only: interp2
  implicit none

  real(kind=dp), private, parameter :: EARTH_RADIUS = 6371.0_dp ! Earth's radius in km
  real(kind=dp), private, parameter :: max_dep = 6360.0_dp, min_dep = -10.0_dp, dr = 5.0_dp, &
                                       min_gcarc = -1.0_dp, max_gcarc = 180.0_dp
  integer, private, parameter :: n_gcarc = 1200

  type :: velmodel_type
    integer :: n_layers, n_gcarc
    real(kind=dp), allocatable :: depths(:) ! depth at bottom of each layer (km)
    real(kind=dp), allocatable :: radius(:) ! radius of each layer (km)
    real(kind=dp), allocatable :: gcarc(:)  ! great circle arc (degrees)
    real(kind=dp), allocatable :: vp(:)     ! P-wave velocity in each layer (km/s)
    real(kind=dp), allocatable :: slowness(:) ! slowness in each layer (s/km)
    real(kind=dp), dimension(:,:), allocatable :: a, b, c, slowness_2d ! 2D slowness grid for ray tracing
  contains
    procedure :: init => initialize, get_travel_time
    procedure, private :: read_from_file, get_ak135
  end type velmodel_type
  
contains
  subroutine initialize(this, model_name)
    class(velmodel_type), intent(inout) :: this
    character(len=*), intent(in) :: model_name
    real(kind=dp) :: min_radius, max_radius, dgcarc
    integer :: i

    ! initialize radius and gcarc arrays
    min_radius = EARTH_RADIUS - max_dep
    max_radius = EARTH_RADIUS - min_dep
    this%n_layers = int((max_radius - min_radius)/dr) + 1
    this%radius = [(min_radius + (i - 1) * dr, i=1, this%n_layers)]
    this%depths = EARTH_RADIUS - this%radius

    ! Initialize great circle arc array (min_gcarc to max_gcarc)
    this%n_gcarc = n_gcarc
    dgcarc = (max_gcarc - min_gcarc)/real(this%n_gcarc - 1, kind=dp)
    this%gcarc = [(min_gcarc + (i - 1) * dgcarc, i = 1, this%n_gcarc)]
    this%gcarc = this%gcarc * 3.141592653/180.0

    if (trim(model_name) == 'AK135') then
      call this%get_ak135()
    else
      call this%read_from_file(trim(model_name))
    end if

    allocate(this%slowness_2d(this%n_layers, this%n_gcarc))
    allocate(this%a(this%n_layers, this%n_gcarc))
    allocate(this%b(this%n_layers, this%n_gcarc))
    allocate(this%c(this%n_layers, this%n_gcarc))
    this%a = 1.0_dp
    this%b = 1.0_dp
    this%c = 0.0_dp

    do i = 1, this%n_gcarc
      this%slowness_2d(:, i) = this%slowness
    end do

  end subroutine initialize
  
  function get_travel_time(this, dsrc, drec, gcarc_rec) result(ttime)
    class(velmodel_type), intent(in) :: this
    real(kind=dp), intent(in) :: dsrc ! source radius (km)
    real(kind=dp), dimension(:), intent(in) :: drec ! receiver radii (km)
    real(kind=dp), dimension(:), intent(in) :: gcarc_rec ! receiver great circle arc (degrees)
    real(kind=dp), dimension(:), allocatable :: ttime ! travel times (s)

    real(kind=dp), dimension(:), allocatable :: rrec
    real(kind=dp), dimension(:,:), allocatable :: tt_field
    real(kind=dp) :: rsrc, tsrc
    integer :: nrec, i

    allocate(tt_field(this%n_layers, this%n_gcarc))
    tt_field = 0.0_dp

    tsrc = 0.0_dp ! assuming source at 0 degrees
    rsrc = EARTH_RADIUS - dsrc
    call FSM_UW_PS_sphe_2d(this%radius, this%gcarc, & 
                          this%n_layers, this%n_gcarc, this%a, this%b, this%c, &
                          tt_field, this%slowness_2d, rsrc, tsrc)

    nrec = size(drec)
    rrec = EARTH_RADIUS - drec
    allocate(ttime(nrec))
    do i = 1, nrec
      ttime(i) = interp2(this%radius, this%gcarc, tt_field, rrec(i), gcarc_rec(i)*3.141592653/180.0)
    end do

  end function get_travel_time

  subroutine read_from_file(this, filename)
    class(velmodel_type), intent(inout) :: this
    character(len=*), intent(in) :: filename
    character(len=MAX_STRING_LEN) :: line
    integer :: i, j, ios
    real(kind=dp) :: vs, rho
    real(kind=dp), dimension(:), allocatable :: depths, vp
    integer :: unit, ndep

    ! Open the file for reading
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, "Error opening velocity model file: ", filename
      stop
    end if

    ! First pass: count the number of layers
    ndep = 0
    do
      read(unit, *, iostat=ios) line
      if (ios /= 0) exit
      ndep = ndep + 1
    end do

    ! Allocate arrays based on the number of layers
    allocate(depths(ndep))
    allocate(vp(ndep))

    ! Rewind the file to read the data into arrays
    rewind(unit)
    
    ! Second pass: read the data into arrays
    do i = 1, ndep
      read(unit, *, iostat=ios) line
      read(line, *, iostat=ios) depths(i), vp(i), vs, rho
      if (ios /= 0) then
        print *, "Error reading velocity model in ", trim(filename), " at line ", i
        stop
      end if
    end do
    close(unit)

    ! Set up the model grid
    allocate(this%vp(this%n_layers))
    allocate(this%slowness(this%n_layers))

    ! Interpolate vp to the desired depth grid
    do i = 1, this%n_layers
      if (this%depths(i) <= depths(1)) then
        this%vp(i) = vp(1)
      elseif (this%depths(i) >= depths(ndep)) then
        this%vp(i) = vp(ndep)
      else
        ! Linear interpolation
        do j = 1, ndep - 1
          if (this%depths(i) >= depths(j) .and. this%depths(i) < depths(j + 1)) then
            this%vp(i) = vp(j) + &
                         (vp(j + 1) - vp(j)) * (this%depths(i) - depths(j)) / (depths(j + 1) - depths(j))
            exit
          end if
        end do
      end if
    end do
    this%slowness = 1.0_dp / this%vp
  end subroutine read_from_file

  subroutine AK135(dep, slowness)
    real(kind=dp), intent(in) :: dep
    real(kind=dp), intent(out) :: slowness
    real(kind=dp), parameter :: ak135_depths(136) = [ &
      0.0, 20.0, 20.0, 35.0, 35.0, 77.5, 120.0, 165.0, 210.0, 210.0, 260.0, 310.0, 360.0, 410.0, 410.0, 460.0, &
      510.0, 560.0, 610.0, 660.0, 660.0, 710.0, 760.0, 809.5, 859.0, 908.5, 958.0, 1007.5, 1057.0, 1106.5, &
      1156.0, 1205.5, 1255.0, 1304.5, 1354.0, 1403.5, 1453.0, 1502.5, 1552.0, 1601.5, 1651.0, 1700.5, 1750.0, &
      1799.5, 1849.0, 1898.5, 1948.0, 1997.5, 2047.0, 2096.5, 2146.0, 2195.5, 2245.0, 2294.5, 2344.0, 2393.5, &
      2443.0, 2492.5, 2542.0, 2591.5, 2640.0, 2690.0, 2740.0, 2740.0, 2789.67, 2839.33, 2891.5, 2891.5, 2939.33, &
      2989.66, 3039.99, 3090.32, 3140.66, 3190.99, 3241.32, 3291.65, 3341.98, 3392.31, 3442.64, 3492.97, 3543.3, &
      3593.64, 3643.97, 3694.3, 3744.63, 3794.96, 3845.29, 3895.62, 3945.95, 3996.28, 4046.62, 4096.95, 4147.28, &
      4197.61, 4247.94, 4298.27, 4348.6, 4398.93, 4449.26, 4499.6, 4549.93, 4600.26, 4650.59, 4700.92, 4801.58, &
      4851.91, 4902.24, 4952.58, 5002.91, 5053.24, 5103.57, 5153.5, 5153.5, 5204.61, 5255.32, 5306.04, 5356.75, &
      5407.46, 5458.17, 5508.89, 5559.6, 5610.31, 5661.02, 5711.74, 5813.16, 5863.87, 5914.59, 5965.3, 6016.01, &
      6066.72, 6117.44, 6168.15, 6218.86, 6269.57, 6320.29, 6371.0 ]
    real(kind=dp), parameter :: ak135_vp(136) = [ &
      5.8, 5.8, 6.5, 6.5, 8.04, 8.045, 8.05, 8.175, 8.3, 8.3, 8.4825, 8.665, 8.8475, 9.03, 9.36, 9.528, &
      9.696, 9.864, 10.032, 10.2, 10.79, 10.9229, 11.0558, 11.1353, 11.2221, 11.3068, 11.3896, 11.4705, 11.5495, 11.6269, &
      11.7026, 11.7766, 11.8491, 11.92, 11.9895, 12.0577, 12.1245, 12.1912, 12.255, 12.3185, 12.3819, 12.4426, 12.5031, &
      12.5631, 12.6221, 12.6804, 12.7382, 12.7956, 12.8526, 12.9096, 12.9668, 13.0222, 13.0783, 13.1336, 13.1894, 13.2465, &
      13.3018, 13.3585, 13.4156, 13.4741, 13.5312, 13.59, 13.6494, 13.6494, 13.653, 13.6566, 13.6602, 8.0, 8.0382, &
      8.1283, 8.2213, 8.3122, 8.4001, 8.4861, 8.5692, 8.6496, 8.7283, 8.8036, 8.8761, 8.9461, 9.0138, &
      9.0792, 9.1426, 9.2042, 9.2634, 9.3205, 9.376, 9.4297, 9.4814, 9.5306, 9.5777, 9.6232, 9.6673, &
      9.71, 9.7513, 9.7914, 9.8304, 9.8682, 9.9051, 9.941, 9.9761, 10.0103, 10.0439, 10.0768, 10.1415, &
      10.1739, 10.2049, 10.2329, 10.2565, 10.2745, 10.2854, 10.289, 11.0427, 11.0585, 11.0718, 11.085, 11.0983, &
      11.1166, 11.1316, 11.1457, 11.159, 11.1715, 11.1832, 11.1941, 11.2134, 11.2219, 11.2295, 11.2364, 11.2424, &
      11.2477, 11.2521, 11.2557, 11.2586, 11.2606, 11.2618, 11.2622 ]
    integer :: i
    real(kind=dp) :: vel
    logical :: found

    vel = ak135_vp(1) ! default value
    found = .false.
    
    if (dep <= ak135_depths(1)) then
      vel = ak135_vp(1)
      found = .true.
    elseif (dep >= ak135_depths(size(ak135_depths))) then
      vel = ak135_vp(size(ak135_vp))
      found = .true.
    else
      do i = 1, size(ak135_depths) - 1
        if (dep >= ak135_depths(i) .and. dep < ak135_depths(i+1)) then
          vel = ak135_vp(i) + (ak135_vp(i+1) - ak135_vp(i)) * (dep - ak135_depths(i)) / (ak135_depths(i+1) - ak135_depths(i))
          found = .true.
          exit
        end if
      end do
    end if
    
    slowness = 1.0_dp / vel
  end subroutine AK135

  subroutine get_ak135(this)
    class(velmodel_type), intent(inout) :: this
    integer :: i
    real(kind=dp) :: slowness_val

    allocate(this%vp(this%n_layers))
    allocate(this%slowness(this%n_layers))

    do i = 1, this%n_layers
      call AK135(this%depths(i), slowness_val)
      this%slowness(i) = slowness_val
      this%vp(i) = 1.0_dp / slowness_val
    end do
  end subroutine get_ak135


end module vpmodel