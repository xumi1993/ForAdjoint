module fftpack
    implicit none
    integer, private, parameter :: dp = kind(1.0d0)

    type, public :: fft_cls
        contains
            procedure :: fft
            procedure :: ifft
            procedure :: ifft_complex
            procedure :: fftfreq
            procedure :: hilbert
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

    function ifft_complex(this, fpx, nft) result(x)
        class(fft_cls) :: this
        integer :: nft
        complex(kind=dp), dimension(nft) :: x
        complex, dimension(nft) :: fpx, fp_norm
        fp_norm = fpx/nft
        call fft_raw(nft, fp_norm, +1)
        x = cmplx(fp_norm, kind=dp)
    end function ifft_complex

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

    function hilbert(this, x, n_opt) result(xa)
        ! FFT-based computation of the analytic signal
        ! The analytic signal is calculated by filtering out the negative frequencies and
        ! doubling the amplitudes of the positive frequencies in the FFT domain
        
        class(fft_cls) :: this
        real(kind=dp), dimension(:), intent(in) :: x     ! Input real signal
        integer, intent(in), optional :: n_opt           ! Number of Fourier components
        complex(kind=dp), dimension(:), allocatable :: xa ! Analytic signal output
        
        ! Local variables
        integer :: n, n_input, i
        real(kind=dp), dimension(:), allocatable :: x_padded
        complex(kind=dp), dimension(:), allocatable :: xf, h_weights
        
        n_input = size(x)
        
        ! Set the FFT size
        if (present(n_opt)) then
            n = n_opt
        else
            n = n_input
        end if
        
        if (n <= 0) then
            write(*,*) 'Error: N must be positive'
            return
        end if
        
        ! Prepare input data with zero padding if necessary
        allocate(x_padded(n))
        x_padded = 0.0_dp
        if (n_input <= n) then
            x_padded(1:n_input) = x(1:n_input)
        else
            x_padded(1:n) = x(1:n)
        end if
        
        ! Compute FFT
        allocate(xf(n))
        xf = this%fft(x_padded, n)
        
        ! Create Hilbert filter weights
        allocate(h_weights(n))
        h_weights = (0.0_dp, 0.0_dp)
        
        if (mod(n, 2) == 0) then
            ! Even length
            h_weights(1) = (1.0_dp, 0.0_dp)          ! DC component
            h_weights(n/2 + 1) = (1.0_dp, 0.0_dp)    ! Nyquist frequency
            do i = 2, n/2
                h_weights(i) = (2.0_dp, 0.0_dp)      ! Positive frequencies
            end do
            ! Negative frequencies remain zero
        else
            ! Odd length
            h_weights(1) = (1.0_dp, 0.0_dp)          ! DC component
            do i = 2, (n + 1)/2
                h_weights(i) = (2.0_dp, 0.0_dp)      ! Positive frequencies
            end do
            ! Negative frequencies remain zero
        end if
        
        ! Apply Hilbert filter
        do i = 1, n
            xf(i) = xf(i) * h_weights(i)
        end do
        
        ! Compute inverse FFT to get analytic signal
        allocate(xa(n))
        xa = this%ifft_complex(cmplx(xf, kind=4), n)
        
        ! Clean up
        deallocate(x_padded, xf, h_weights)
        
    end function hilbert
end module