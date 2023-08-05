subroutine distance(x1,y1,z1,x2,y2,z2,d) 
    implicit none
    real, intent(in) :: x1, y1, x2, y2, z1, z2
    real             :: d
    
    d = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
end subroutine distance

subroutine f_grav(m1,m2,d, Force)
    implicit none
    real, intent(in) :: m1, m2, d
    real             :: G = 6.67*10e-11
    real             :: Force
    
    Force = G * (m1 * m2 / d **2)
end subroutine f_grav

program n_body
    use ISO_FORTRAN_ENV, only : REAL128
    implicit none
    character(len=*), parameter :: OUT_FILE = 'n_body.txt' ! Output file.
    integer,          parameter :: N = 10, Cycle = 10    ! Number of particles & Cycles         
    integer                     :: i, j, fu 

    real                        :: u1(N), u2(N), mu_mass, sigma_mass, z_BM(N), mass(N)
    real                        :: pi = 4 * atan(1.0)

    real                        :: x(N), y(N), z(N)
    real                        :: distance_values(N)
    real, dimension(N)          :: vx, vy, vz, x1, y1, z1, x2, y2, z2

    ! Creating Mass Array for Objects - Standard Distribition 
    call random_number(u1)
    call random_number(u2)

    z_BM       = sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    mu_mass    = log(10e25) - 0.5 * log((10e17 / 10e16)**2)
    sigma_mass = sqrt(log((10e30 / 10e20)**2))
    mass       = exp(mu_mass + sigma_mass * z_BM)

    ! Creating initial XYZ Positions for Objects
    
    call random_number(x)
    call random_number(y)
    call random_number(z)

    x = x * 1000 + 1
    y = y * 1000 + 1
    z = z * 1000 + 1 

    ! Creating initial velocities
    vx = abs(vx * 0)
    vy = abs(vy * 0)
    vz = abs(vz * 0)

    
    open (newunit=fu, file=OUT_FILE, status='replace', action='write', form='formatted')       
        do i = 1, N
            write(fu, *) mass(i), x(i), y(i), z(i), vx(i), vy(i), vz(i)
            x1 = x(i)
            y1 = y(i)
            z1 = z(i)

            do j = i, N
                x2 = x(j)
                y2 = y(j)
                z2 = z(j)

                call distance(x1,y1,z1,x2,y2,z2,distance_values)
            end do 

        end do

    
        
    
    close(fu)

end program n_body




    