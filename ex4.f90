program main
    implicit none
    real, dimension(:,:,:), allocatable :: A,B
    real :: sum
    integer :: i,j,k,dx,dy,dz, N=128

    ALLOCATE(A(N,N,N))
    ALLOCATE(B(N,N,N))

    B = 0.

    
    do k=1,N
        do j=1,N
            do i=1,N
                A(i,j,k) = real(i+j+k)
            enddo
        enddo
    enddo

    !$OMP PARALLEL SHARED(B,N,A) PRIVATE(i,j,k,dx,dy,dz,sum)
    !$OMP DO
    do k=3,N-2
        do j=3,N-2
            do i=3,N-2
                sum = 0
                do dz=k-2, k+2
                    do dy=j-2, j+2
                        do dx=i-2,i+2
                            sum = sum + A(dx,dy,dz)
                        enddo
                    enddo
                enddo
                B(i,j,k) = sum
            enddo
        enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    print *,B(N/2, N/2, N/2)

    DEALLOCATE(A)
    DEALLOCATE(B)

end program