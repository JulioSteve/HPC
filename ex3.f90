program main
    implicit none
    integer :: i

    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1,10
        write (*,*) i
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
end program