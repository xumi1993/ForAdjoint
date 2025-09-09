module multitaper
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  public :: staper
contains

  ! Slepian - Thomson multi-taper procedure
  subroutine staper(nt, nev, v, ndim, a, w)
    integer, intent(in) :: nt, nev, ndim
    real(kind=dp), dimension(:,:), intent(out), allocatable :: v
    real(kind=dp), dimension(:), intent(out), allocatable :: a, w

    real(kind=dp), parameter :: PI = 3.141592653589793d0
    real(kind=dp) :: r2, om, com, hn, asav, rbd, dc, sm, s, sn, vmax, fw
    integer :: i, j, k, m, nxi, lh, lp1, neven, nodd, ntot, kk, kmax, nlow, nup

    r2 = sqrt(2.0d0)

    if (nt < 2) return
    allocate(v(ndim, nev))
    allocate(a(ndim), w(ndim))
    v = 0.0d0
    a = 0.0d0
    w = 0.0d0

    nxi = mod(nt, 2)
    lh = (nt/2) + nxi
    lp1 = nt + 1
    fw = dble(nev) / 2.0d0
    om = 2.0d0 * PI * fw / nt
    com = cos(om)
    hn = 0.5d0 * dble(lp1)

    do i = 1, lh
      a(i) = com * (i-hn)**2
      w(i) = 0.5d0 * dble(i * (nt-i))
    end do

    if (nxi == 0) then
      asav = a(lh) - w(lh)
      a(lh) = a(lh) + w(lh)
      rbd = 1.0d0 / (a(lh) + w(lh-1))
    else
      asav = w(lh-1)
      rbd = 1.0d0 / (w(lh) + w(lh-1))
      w(lh-1) = r2 * w(lh-1)
    end if

    do i = 1, lh
      a(i+lh) = w(i) * rbd
      w(i) = a(i+lh)**2
      a(i) = a(i) * rbd
    end do

    neven = max((nev+1)/2, 1)
    nodd = nev - neven

    ! Even tapers
    call tsturm(nt, lh, a, a(lh+1:), w, neven, v, ndim, w(lh+1:), 0)
    do i = 1, neven
      k = 2*i - 1
      if (nxi == 1) v(lh, k) = r2 * v(lh, k)
      do j = 1, lh
        v(lp1-j, k) = v(j, k)
      end do
    end do
    if (nodd <= 0) goto 34

    ! Odd tapers
    if (nxi == 0) then
      a(lh) = asav * rbd
    else
      a(nt) = asav * rbd
      w(lh-1) = asav * asav
    end if
    call tsturm(nt, lh-nxi, a, a(lh+1:), w, nodd, v, ndim, w(lh+1:), 1)
    do i = 1, nodd
      k = 2*i
      if (nxi == 1) v(lh, k) = 0.0d0
      do j = 1, lh
        v(lp1-j, k) = -v(j, k)
      end do
    end do
34  ntot = neven + nodd

    ! Bandwidth retention
    dc = 2.0d0 * com
    sm = 0.0d0
    s = sin(om)
    w(1) = om / PI
    w(2) = s / PI
    do j = 3, nt
      sn = dc * s - sm
      sm = s
      s = sn
      w(j) = s / (PI * dble(j-1))
    end do

    do m = 1, ntot
      vmax = abs(v(1, m))
      kmax = 1
      do kk = 2, lh
        if (abs(v(kk, m)) <= vmax) exit
        kmax = kk
        vmax = abs(v(kk, m))
      end do
      a(m) = 0.0d0
      nlow = kmax - 1
      do j = 1, nlow
        a(m) = a(m) + w(j+1) * v(nlow+1-j, m)
      end do
      nup = nt - nlow
      do j = 1, nup
        a(m) = a(m) + w(j) * v(nlow+j, m)
      end do
      a(m) = a(m) / v(kmax, m)
    end do

    return
  end subroutine staper

  ! Sturm sequence for tridiagonal matrix eigenvalues/eigenvectors
  subroutine tsturm(nt, n, a, b, w, nev, r, ndim, ev, ipar)
    integer, intent(in) :: nt, n, nev, ndim, ipar
    real(kind=dp), intent(in) :: a(:), b(:), w(:)
    real(kind=dp), intent(inout) :: r(ndim, *)
    real(kind=dp), intent(inout) :: ev(:)

    real(kind=dp), parameter :: epsi = 1.0d-15, epsi1 = 5.0d-15
    real(kind=dp) :: q, el, elam, u, umeps, x, ddot, rnorm
    integer :: i, j, ik, iag, m, jk, jm1

    if (n <= 0 .or. nev <= 0) return
    umeps = 1.0d0 - epsi
    ev(:nev) = -1.0d0
    u = 1.0d0
    do ik = 1, nev
      if (ik > 1) u = ev(ik-1) * umeps
      el = min(ev(ik), u)
      do
        elam = 0.5d0 * (u + el)
        if (abs(u - el) <= epsi1) exit
        iag = 0
        q = a(1) - elam
        if (q >= 0.0d0) iag = iag + 1
        do i = 2, n
          if (q == 0.0d0) then
            x = abs(b(i-1)) / epsi
          else
            x = w(i-1) / q
          end if
          q = a(i) - elam - x
          if (q >= 0.0d0) iag = iag + 1
          if (iag > nev) exit
        end do
        if (iag >= ik) then
          u = elam
        else
          el = elam
        end if
      end do
      call root(u, el, elam, a, b, w, n, ik)
      ev(ik) = elam
      jk = 2*ik + ipar - 1
      r(1, jk) = 1.0d0
      r(2, jk) = -(a(1) - ev(ik)) / b(1)
      ddot = 1.0d0 + r(2, jk)**2
      jm1 = 2
      do j = 3, n
        r(j, jk) = -((a(jm1) - ev(ik)) * r(jm1, jk) + b(j-2) * r(j-2, jk)) / b(jm1)
        ddot = ddot + r(j, jk)**2
        jm1 = j
      end do
      rnorm = sqrt(dble(nt) / (2.0d0 * ddot))
      do j = 1, n
        r(j, jk) = r(j, jk) * rnorm
      end do
    end do
    return
  end subroutine tsturm

  ! Newton's method for refining eigenvalues
  subroutine root(u, el, elam, a, bb, w, n, ik)
    real(kind=dp), intent(inout) :: u, el, elam
    real(kind=dp), intent(in) :: a(:), bb(:), w(:)
    integer, intent(in) :: n, ik

    real(kind=dp), parameter :: epsi = 1.0d-15, epsi1 = 5.0d-15
    real(kind=dp) :: an, b, bm, bn, del, x
    integer :: i, iag

    do
      elam = 0.5d0 * (u + el)
      if (abs(u - el) <= 1.5d0 * epsi1) return
      an = a(1) - elam
      b = 0.0d0
      bn = -1.0d0 / an
      iag = 0
      if (an >= 0.0d0) iag = iag + 1
      do i = 2, n
        if (an == 0.0d0) then
          x = abs(bb(i-1)) / epsi
        else
          x = w(i-1) / an
        end if
        an = a(i) - elam - x
        if (an == 0.0d0) an = epsi
        bm = b
        b = bn
        bn = ((a(i) - elam) * b - bm * x - 1.0d0) / an
        if (an >= 0.0d0) iag = iag + 1
      end do
      if (iag == ik) then
        el = elam
      else
        u = elam
      end if
      del = 1.0d0 / bn
      if (abs(del) <= epsi1) del = sign(epsi1, del)
      elam = elam - del
      if (elam >= u .or. elam <= el) cycle
    end do
    return
  end subroutine root

end module multitaper