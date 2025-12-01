program main
    implicit none
    real(kind=8) :: x, dx, S=0
    integer :: i, N

    N = int(1e8)

    dx = 1.d0/real(N-1)

    do i=1,N
        x = real(i-1)*dx
        S = S + F(x)*dx
    enddo

    print *,S

    contains
        function F(x)
            real(kind=8) :: F
            real(kind=8), intent(in) :: x

            F = 4.d0/(1.d0+x*x)

        end function F

end program