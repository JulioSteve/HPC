program main
    implicit none
    integer :: N,i,j
    integer,parameter :: Tsteps=int(3)
    real :: m,mtot,x,y,z,dist0, omega, dt, t, G, vcm
    real,dimension(:,:),allocatable :: pos,vel
    integer,dimension(8) :: seed
    
    seed = [324549107,1517762999,2043478039,169521159,-1935932375,513922432,1075905891,-1147996143]
    call random_seed(put=seed)

    print *,"Please, enter a number of bodies:"
    read(*,*) N
    m = 1.
    omega = 1.
    G = 1.
    mtot = m*N
    
    allocate(pos(N,3))
    allocate(vel(N,3))
    do i=1,N
        dist0 = 1.
        do while (dist0>=1)
            call random_number(x)
            call random_number(y)
            call random_number(z)
            x = 2*x-1
            y = 2*y-1
            z = 2*z-1
            
            dist0 = sqrt(x*x+y*y+z*z)
        end do
        pos(i,1) = x
        pos(i,2) = y
        pos(i,3) = z
        vel(i,1) = -y*omega
        vel(i,2) = x*omega
        vel(i,3) = 0
    end do
    
    vcm = 0
    do i=1,N
        vcm = vcm + vel(i,1)*vel(i,1) + vel(i,2)*vel(i,2) + vel(i,3)*vel(i,3)
    enddo
    vcm = sqrt(vcm)/N
    do i=1,N
        vel(i,:) = vel(i,:)-vcm
    enddo

    t = 0.
    dt = 5e-4
    open(unit=10,file="nbody.dat")
    write(10,*) "Time steps:",Tsteps
    do j=0,Tsteps
        if (mod(j,1)==0) then
            do i=1,N
                write(10,'(F8.5,3F10.2,3F10.2)') t,pos(i,:),vel(i,:)
            enddo
        endif
        call leapfrog(dt,pos,vel)
        t = t+dt
    enddo
    close(10)
    

    contains
        subroutine force(f,pos1,pos2)
            real,dimension(3),intent(in) :: pos1,pos2
            real,dimension(3),intent(out) :: f
            real :: fscal, eps
            eps=0.
            
            fscal = (sqrt((pos1(1)-pos2(1))**2+&
            (pos1(2)-pos2(2))**2+(pos1(3)-pos2(3))**2+eps**2)**3)

            fscal = G*m*m/fscal

            f = fscal*(pos2-pos1)
        end subroutine force
        
        subroutine leapfrog(dt,r,v)
            real,intent(in) :: dt
            real,dimension(N,3),intent(inout) :: r,v
            real,dimension(3) :: f
            integer :: i,j

            do i=1,N-1
                f = 0
                do j=i+1,N
                    call force(f,r(i,:),r(j,:))

                    v(i,:) = v(i,:) + dt/m*f
                    v(j,:) = v(j,:) - dt/m*f
                enddo
                
                r(i,:) = r(i,:) + v(i,:)*dt
            enddo
            r(N,:) = r(N,:) + v(N,:)*dt
        end subroutine leapfrog
end program main