program main
    implicit none
    integer :: i,j

    integer,parameter :: Tsteps=int(1e3), N=int(1e3)
    real, parameter :: dt=5e-3, m=1./N, mtot=m*N, omega=1., G=1., R=1.

    real :: x,y,z,dist0, t
    real,dimension(N,3) :: pos,vel
    real, dimension(3) :: vcm
    integer,dimension(8) :: seed
    
    seed = [324549107,1517762999,2043478039,169521159,-1935932375,513922432,1075905891,-1147996143]
    call random_seed(put=seed)
    
    do i=1,N
        dist0 = R
        do while (dist0>=R)
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
    
    do j=1,3
        vcm(j) = sum(vel(:,j))/N
        vel(:,j) = vel(:,j)-vcm(j)
    enddo

    t = 0.
    open(unit=10,file="settings.dat")
    write(10,*) "Time steps, number of bodies, dt"
    write(10,*) Tsteps, N, dt
    close(10)

    open(unit=11,file="pos.bin",form="unformatted")
    open(unit=12,file="vel.bin",form="unformatted")
    do j=0,Tsteps
        if (mod(j,10)==0) then
            write(11) pos(:,:)
            write(12) vel(:,:)
        endif
        call leapfrog(dt,pos,vel)
        t = t+dt
    enddo
    close(11)
    close(12)

    contains
        subroutine force(f,pos1,pos2)
            real,dimension(3),intent(in) :: pos1,pos2
            real,dimension(3),intent(out) :: f
            real :: fscal, eps
            eps=0.99
            
            fscal = (sqrt((pos1(1)-pos2(1))**2+&
            (pos1(2)-pos2(2))**2+(pos1(3)-pos2(3))**2+eps*eps)**3)

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