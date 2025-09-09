module fftpack
    implicit none

    type, public :: fft_cls
        contains
            procedure :: fft
            procedure :: ifft
            procedure :: fftfreq
    end type fft_cls

contains

    function fft(this, x, nft) result(fpx)
        class(fft_cls) :: this
        integer :: nft
        double precision, dimension(nft) :: x
        complex, dimension(nft) :: fpx
        fpx = cmplx(x, 0)
        call fft_raw(nft, fpx, -1)
    end function fft

    function ifft(this, fpx, nft) result(x)
        class(fft_cls) :: this
        integer :: nft
        double precision, dimension(nft) :: x
        complex, dimension(nft) :: fpx, fp_norm
        fp_norm = fpx/nft
        call fft_raw(nft, fp_norm, +1)
        x = dble(real(fp_norm))
    end function ifft

    subroutine fft_raw(n, x, ind)
        integer :: n, ind, j, i, kmax, istep,m, k
        complex :: temp, theta
        complex, intent(inout) :: x(n)

        j=1
        do i=1, n
            if (i < j) then 
                temp=x(j)
                x(j)=x(i)
                x(i)=temp
            endif
            m = n/2
            do while (.true.)
                if (j > m) then
                    j=j-m
                    m=m/2
                    if (m < 2) exit
                else
                    exit
                endif
            enddo
            j=j+m
        enddo
        kmax=1
        do while(kmax < n)
            istep=kmax*2
            do k=1,kmax
            theta=cmplx(0.0,3.141592653*ind*(k-1)/kmax)
                do i=k, n, istep
                    j=i+kmax
                    temp=x(j)*cexp(theta)
                    x(j)=x(i)-temp
                    x(i)=x(i)+temp
                enddo
            enddo
            kmax=istep
        enddo
    end subroutine fft_raw

    function fftfreq(this, n, d) result(freqs)
        class(fft_cls) :: this
        integer, intent(in) :: n         
        double precision, intent(in) :: d        
        double precision, allocatable :: freqs(:) 

        integer :: i, N_half

        allocate(freqs(n))
        N_half = n / 2

        do i = 1, N_half
            freqs(i) = (i - 1) / (n * d)
        end do

        do i = N_half + 1, n
            freqs(i) = - (n - i + 1) / (n * d)
        end do

    end function fftfreq
end module