program main
    implicit none
    real(kind=8) :: x, dx, S_th, S=0
    integer :: i, N, OMP_GET_NUM_THREADS

    N = int(1e8)

    dx = 1.d0/real(N-1)

    !$OMP PARALLEL PRIVATE(x,S_th,i) SHARED(S)
    S_th = 0
    !$OMP DO SCHEDULE(STATIC, N/OMP_GET_NUM_THREADS())
    do i=1,N
        x = real(i-1)*dx
        S_th = S_th + F(x)*dx
    enddo
    !$OMP END DO
    !$OMP ATOMIC
    S = S + S_th
    !$OMP END ATOMIC
    !$OMP END PARALLEL

    print *,S

    contains
        function F(x)
            real(kind=8) :: F
            real(kind=8), intent(in) :: x

            F = 4.d0/(1.d0+x*x)

        end function F

end program