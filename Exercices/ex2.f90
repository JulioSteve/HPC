program main
    implicit none
    real, dimension(:,:,:), allocatable :: A
    integer :: i,j,k,N=1000

    ALLOCATE(A(N,N,N))
    A = 1.
    
    open(unit=10, file="data.bin", form="unformatted")
    write(10) A
    ! do i=1,N
    !     do j=1,N
    !         do k=1,N
    !             write(10) A(i,j,k)
    !         enddo
    !     enddo
    ! enddo
    close(10)
    DEALLOCATE(A)
end program