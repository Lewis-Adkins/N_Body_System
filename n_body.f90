
program n_body
    use ISO_FORTRAN_ENV, only: REAL128
    implicit none
    character(len=*), parameter :: OUT_FILE = 'n_body.txt' ! Output file.
    integer, parameter          :: N = 100, Cycle = 1000           ! Number of particles & Cycles         
    integer                     :: i, j, Frame, idx, fu 

    real                        :: mu_mass, sigma_mass 
    real                        :: pi = 4 * atan(1.0), G = 6.674e-11
    real                        :: t = 0.5 ! time interval 

    real, dimension(N*(N-1))    :: x_distances, y_distances, z_distances, r, combined_masses, Force
    real, dimension(N)          :: vx, vy, vz, Final_x_Force, Final_y_Force, Final_z_Force, x ,y , z, mass, z_BM, u1, u2, ax, ay, az
    real, dimension(N)          :: new_x_position, new_y_position, new_z_position, new_vx, new_vy, new_vz

    real , dimension(N*(N-1))   :: unit_x, unit_y, unit_z, Force_x, Force_y, Force_z
    real , dimension(N-1,N)     :: Force_x_Mat, Force_y_Mat, Force_z_Mat
  
    ! Creating Mass Array for Objects - Standard Distribution 
    call random_number(u1)
    call random_number(u2)

    z_BM       = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2)
    mu_mass    = log(10.0e25) - 0.5 * log((10.0e17 / 10.0e16)**2)
    sigma_mass = sqrt(log((10.0e30 / 10.0e20)**2))
    mass       = exp(mu_mass + sigma_mass * z_BM) / 10e10

     ! Creating initial XYZ Positions for Objects
    
    call random_number(x)
    call random_number(y)
    call random_number(z)

    x = x * 10e3 + 1
    y = y * 10e3 + 1
    z = z * 10e3 + 1 

    ! Creating initial velocities
    vx = 0.0
    vy = 0.0
    vz = 0.0



    x_distances = 0.0
    y_distances = 0.0
    z_distances = 0.0

    ! Open output file
    open(newunit=fu, file=OUT_FILE, status='replace', action='write', form='formatted')       
    !write(fu, *) 'Mass', 'x', 'y', 'z', 'Fx', 'Fy', 'Fz', 'ax', 'ay', 'az'

    do i = 1, N   
        write(fu, *) mass(i), ',', x(i), ',', y(i), ',', z(i), ',', &
        Final_x_Force(i), ',', Final_y_Force(i), ',', Final_z_Force(i), ',', &
        ax(i), ',', ay(i), ',', az(i)                 
    end do  

    do Frame = 1, Cycle    
        
        idx = 1   
        do i = 1, N   
            
            do j = 1, N     
                
                if (j /= i) then
                     
                    x_distances(idx) = x(i) - x(j)        
                    y_distances(idx) = y(i) - y(j)      
                    z_distances(idx) = z(i) - z(j)
                       
                    r = sqrt(x_distances ** 2 + y_distances ** 2 + z_distances ** 2)
       
                    combined_masses(idx) = G * mass(i) * mass(j)
           
                    Force_x = -x_distances * combined_masses/(r **3)          
                    Force_y = -y_distances * combined_masses/(r **3)                  
                    Force_z = -z_distances * combined_masses/(r **3)          
                                                
                    idx = idx + 1
                          
                end if      
            end do   
        end do

 
        Force_x_Mat = RESHAPE(Force_x, [N-1,N])  
        Force_y_Mat = RESHAPE(Force_y, [N-1,N]) 
        Force_z_Mat = RESHAPE(Force_z, [N-1,N])
        
        !print*, 'x_distances', x_distances
        !print*, 'r = ', r
        !print*, 'G*m1*m2 = ',combined_masses
        !print*, 'x = ',x_distances
        !print*,'Force_x = ', Force_x
         
        idx = 1
        
        do i = 1, N
   
            Final_x_Force(i) = SUM(Force_x_Mat(:,i))  
            Final_y_Force(i) = SUM(Force_y_Mat(:,i)) 
            Final_z_Force(i) = SUM(Force_z_Mat(:,i))
        
      
            ax(i)            = Final_x_Force(i) / mass(i)   
            ay(i)            = Final_y_Force(i) / mass(i)   
            az(i)            = Final_z_Force(i) / mass(i) 
            
            vx(i) = vx(i) + ax(i) * t
            vy(i) = vy(i) + ay(i) * t
            vz(i) = vz(i) + az(i) * t
       
            new_x_position(idx) = x(i) + vx(i) * t 
            new_y_position(idx) = y(i) + vy(i) * t 
            new_z_position(idx) = z(i) + vz(i) * t 

            idx = idx + 1
         
        end do
    
        do i = 1, N
            x(i) = new_x_position(i)
            y(i) = new_y_position(i)
            z(i) = new_z_position(i)
        end do

        !print*, 'ax = ',ax
        !print*, 'vx = ', vx
        

    do i = 1, N   
        write(fu, *) mass(i), ',', x(i), ',', y(i), ',', z(i), ',', &
        Final_x_Force(i), ',', Final_y_Force(i), ',', Final_z_Force(i), ',', &
        ax(i), ',', ay(i), ',', az(i)                 
    end do  
      
    end do
        
    close(fu)

end program n_body