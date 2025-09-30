module vpmodel
  use config
  use signal, only: interp2
  implicit none

  real(kind=dp), private, parameter :: EARTH_RADIUS = 6371.0_dp ! Earth's radius in km
  integer, private, parameter :: nt = 1000
  real(kind=dp), private, parameter :: max_dep = 670.0_dp, min_dep = -10.0_dp, dr = 1.0_dp

  type :: velmodel_type
    integer :: n_layers, n_gcarc = nt
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
    integer :: i

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
      ttime(i) = interp2(this%radius, this%gcarc, tt_field, rrec(i), gcarc_rec(i))
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
    this%n_layers = int((max_dep - min_dep)/dr) + 1
    this%depths = [(min_dep + (i - 1) * dr, i=1, this%n_layers)]
    this%radius = EARTH_RADIUS - this%depths
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

  subroutine AK135(dep,slowness)
    real(kind=dp), intent(in) :: dep
    real(kind=dp), intent(out) :: slowness
    real(kind=dp) :: velocity(54),layer(54),r1,vel
    integer :: ii,idx0

    velocity    = [ 5.800,5.800,6.500,6.500,8.040,8.045,8.050,8.175,8.300, &
                  & 8.300,8.483,8.665,8.847,9.030,9.360,9.528,9.696,9.864, &
                  & 10.032,10.200,10.790,10.923,11.056,11.135,11.222,11.307, &
                  & 11.390,11.470,11.550,11.627,11.703,11.777,11.849,11.920, &
                  & 11.990,12.058,12.125,12.191,12.255,12.318,12.382,12.443, &
                  & 12.503,12.563,12.622,12.680,12.738,12.796,12.853,12.910, &
                  & 12.967,13.022,13.078,13.134]
    layer       = [ 0.0,  20.0 ,20.0 ,35.0, 35.0, 77.5, 120.0,165.0,210.0, &
                  & 210.0,260.0,310.0,360.0,410.0,410.0,460.0,510.0,560.0, &
                  & 610.0, 660.0, 660.0, 710.0, 760.0, 809.5, 859.5, 908.5, &
                  &  958.0,1007.5,1057.0,1106.5,1156.0,1205.5,1255.0,1304.5, &
                  & 1354.0,1403.5,1453.0,1502.5,1552.0,1601.5,1651.0,1700.5, &
                  & 1750.0,1799.5,1849.0,1898.5,1948.0,1997.5,2047.0,2096.5, &
                  & 2146.0,2195.5,2245.0,2294.5]

    if (dep <= layer(1)) then
      vel = velocity(1)
    elseif (dep > layer(28)) then
      vel = velocity(28)
    else
      do ii=1,27
        if (dep>layer(ii) .and. dep<= layer(ii+1)) then
          idx0=ii
          exit
        end if
      end do
      r1 = min(1.0d0, (dep-layer(ii))/(layer(ii+1)-layer(ii)))
      vel = r1*velocity(idx0+1) + (1-r1)*velocity(idx0)
    end if

    slowness = 1.0/vel
  end subroutine AK135

  subroutine get_ak135(this)
    class(velmodel_type), intent(inout) :: this
    integer :: i
    real(kind=dp) :: slowness_val

    this%n_layers = int((max_dep - min_dep)/dr) + 1
    this%depths = [(min_dep + (i - 1) * dr, i=1, this%n_layers)]
    this%radius = EARTH_RADIUS - this%depths

    allocate(this%vp(this%n_layers))
    allocate(this%slowness(this%n_layers))

    do i = 1, this%n_layers
      call AK135(this%depths(i), slowness_val)
      this%slowness(i) = slowness_val
      this%vp(i) = 1.0_dp / slowness_val
    end do
    
  end subroutine get_ak135


end module vpmodel