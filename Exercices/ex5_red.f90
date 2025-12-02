program main
    implicit none
    real(kind=8) :: x, dx, S
    integer(kind=8) :: i, N,OMP_GET_NUM_THREADS

    N = 5e11

    dx = 1.d0/real(N-1)

    !$OMP PARALLEL PRIVATE(x,i) SHARED(S)
    S = 0
    !$OMP DO SCHEDULE(STATIC, N/OMP_GET_NUM_THREADS()) REDUCTION(+ :S)
    do i=1,N
        x = real(i-1)*dx
        S = S + F(x)*dx
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! S = S*dx

    print *,S

    contains
        function F(x)
            real(kind=8) :: F
            real(kind=8), intent(in) :: x

            F = 4.d0/(1.d0+x*x)

        end function F

end program