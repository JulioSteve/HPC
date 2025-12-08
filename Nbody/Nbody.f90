    program main
        implicit none
        integer :: i,j, modu

        integer,parameter :: Tsteps=int(5e3), N=int(5e2)
        real, parameter :: dt=2e-3, m=1./N, mtot=m*N, omega=1., G=1., R=1., eps=0.1

        real :: x,y,z,dist0, t, Ep, Ek

        real,dimension(N,3) :: pos,vel
        real, dimension(3) :: vcm
        integer,dimension(8) :: seed
        
        seed = [324549107,1517762999,2043478039,169521159,-1935932375,513922432,1075905891,-1147996143]
        call random_seed(put=seed)

        modu = int(1./(dt*60))
        if (modu==0) then
            modu=1
        endif

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
        open(unit=10,file="settings.dat", status="replace")
        write(10,*) "Time steps, number of bodies, dt, modulo (every simulation results are printed at each modulo)"
        write(10,*) Tsteps, N, dt, modu
        close(10)

        open(unit=11,file="pos.bin",form="unformatted", access="stream", status="replace")
        open(unit=12,file="vel.bin",form="unformatted", access="stream", status="replace")
        open(unit=13, file="energies.bin", form="unformatted", access="stream", status="replace")

        Ek = 0.5*m*sum(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2)+vel(:,3)*vel(:,3)) !initial kinetic energy
        call initialpot(Ep, pos) !we do it just to calculate the initial Potential Energy, the function does not modify anything but Ep
        ! we are forced to do so to automatize the next loops, doing the leapfrog first then writing
        write(11) pos(:,:)
        write(12) vel(:,:)
        write(13) Ek, Ep

        do j=1,Tsteps
            call leapfrog(pos,vel,Ep)
            if (mod(j,modu)==0) then
                Ek = 0.5*m*sum(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2)+vel(:,3)*vel(:,3))

                write(11) pos(:,:)
                write(12) vel(:,:)
                write(13) Ek, Ep
            endif
            t = t+dt
        enddo
        close(11)
        close(12)
        close(13)

        contains
            subroutine initialpot(potential,r)
                real,dimension(N,3),intent(in) :: r
                real,intent(out) :: potential
                real :: dx,dy,dz
                integer :: i,j

                potential=0
                
                do i=1,N-1
                    do j=i+1,N
                        dx = r(j,1)-r(i,1)
                        dy = r(j,2)-r(i,2)
                        dz = r(j,3)-r(i,3)

                        potential = potential + 1./sqrt(dx*dx+dy*dy+dz*dz+eps*eps)
                    enddo
                enddo

                potential = -G*m*m*potential
            
            end subroutine initialpot


            subroutine accelerations(A,potential,r)
                real,dimension(N,3),intent(in) :: r
                real,dimension(N,3),intent(out) :: A
                real,intent(out) :: potential
                real :: dx,dy,dz, r3, Ax, Ay, Az
                integer :: i,j

                A = 0
                potential=0
                
                do i=1,N-1
                    do j=i+1,N
                        dx = r(j,1)-r(i,1)
                        dy = r(j,2)-r(i,2)
                        dz = r(j,3)-r(i,3)
                        
                        r3 = (dx*dx+dy*dy+dz*dz+eps*eps)**1.5

                        Ax = G*m*dx/r3
                        Ay = G*m*dy/r3
                        Az = G*m*dz/r3

                        A(i,1) = A(i,1) + Ax
                        A(i,2) = A(i,2) + Ay
                        A(i,3) = A(i,3) + Az
                        
                        ! We calculated f_ij which appears in every next terms (with opposite sign by Newton 3rd law)
                        ! So we just fill every next term so we don't calculate these terms again
                        A(j,1) = A(j,1) - Ax
                        A(j,2) = A(j,2) - Ay
                        A(j,3) = A(j,3) - Az

                        potential = potential + 1./sqrt(dx*dx+dy*dy+dz*dz+eps*eps) ! We calculate the potential energy there to optimize
                    enddo
                enddo

                potential = -G*m*m*potential

            end subroutine accelerations


            subroutine leapfrog(r,v,pot)
                real,dimension(N,3),intent(inout) :: r,v
                real,dimension(N,3) :: A, r_half
                real,intent(out) :: pot
                
                r_half = r + v*dt/2.             !semi-drift

                call accelerations(A,pot,r_half)

                v = v + A*dt                     !kick

                r = r_half + v*dt/2.             !semi-drift

            end subroutine leapfrog


    end program main