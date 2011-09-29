program Laplace2D_all
! Compares Jacobi and G-S methods with SOR

use tecplot_write
implicit none

    real(rprec), parameter :: pi = 3.1415926535
    real(rprec), parameter :: xmin = 0.0, xmax = 2.0*pi
    real(rprec), parameter :: ymin = 0.0, ymax = 2.0*pi
    real(rprec), parameter :: dx = (xmax-xmin)/(64-1), dy = (ymax-ymin)/(64-1)
    real(rprec), parameter :: errmax = 1.0e-3
    real(rprec), parameter :: wj = 1.0, wgs = 1.0


    integer, dimension(2) :: kmax
    real(rprec), pointer, dimension(:,:) :: u0, u, unew
    real(rprec), pointer, dimension(:) :: x, y
    real(rprec) :: error, du, ran3
    integer :: Nx, Ny, i, j, k, m, n, p
    real clock_start, clock_end


! Call timer
    call cpu_time(clock_start)

! Determine number of gridpoints (0:N) where BCs are applied at 0 and N
    Nx = (xmax-xmin)/dx
    Ny = (ymax-ymin)/dy
    write(*,*) 'Nx,Ny',Nx,Ny

! Allocate arrays, initialize
    allocate(u(0:Nx,0:Ny),unew(0:Nx,0:Ny),u0(1:Nx-1,1:Ny-1))
    allocate(x(0:Nx),y(0:Ny))

! Create x,y arrays
    do i=0,Nx
        x(i) = xmin + i*dx
    enddo
    do j=0,Ny
        y(j) = ymin + j*dy
    enddo  

! Set u0 using random initial condition
    do j=1,Ny-1
    do i=1,Nx-1 
        u0(i,j) = ran3(i*j*clock_start)
    enddo
    enddo

! Apply BCs
    u(0,:) = 0.0
    u(Nx,:) = 0.0
    u(:,Ny) = 0.0
    do i=1,Nx-1
        u(i,0) = sin(10.*x(i)) + sin(20.*x(i)) + sin(30.*x(i)) + sin(40.*x(i))
    enddo
    unew = u


! Jacobi Iterate
    ! Initialize interior of u
        u(1:(Nx-1),1:(Ny-1)) = u0
    ! Iterate
        k = 0
        error = 100.
        do while (error .gt. errmax)
            error = 0.
            do i=1,Nx-1
                do j=1,Ny-1
                    du = wj * (0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)) - u(i,j))
                    unew(i,j) = u(i,j) + du
                    error = error + abs(du)
                enddo
            enddo
            error = error / (Nx*Ny)
            k = k + 1
            u = unew
        enddo
    ! Store value of kmax
        if (error .lt. 100) then
            kmax(1) = k
        else
            kmax(1) = 0.      ! solution did not converge
        endif

! Call timer
    call cpu_time(clock_end)
    write ( *, * ) 'Jacobi Elapsed CPU time = ', clock_end - clock_start

! Call timer
    call cpu_time(clock_start)

! Gauss-Siedel Iterate
    ! Initialize interior of u
        u(1:(Nx-1),1:(Ny-1)) = u0
    ! Iterate
        k = 0
        error = 100.
        do while (error .gt. errmax)
            error = 0.
            do i=1,Nx-1
                do j=1,Ny-1
                    du = wgs * (0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)) - u(i,j))
                    u(i,j) = u(i,j) + du
                    error = error + abs(du)
                enddo
            enddo
            error = error/(Nx*Ny)
            k = k + 1
        enddo
    ! Store value of kmax
        if (error .lt. 100) then
            kmax(2) = k
        else
            kmax(2) = 0.      ! solution did not converge
        endif

! Call timer
    call cpu_time(clock_end)
    write ( *, * ) 'GS Elapsed CPU time = ', clock_end - clock_start

! Write kmax to screen

    777 format (a5,f6.4,a7,es14.7,a6,f7.4,a7,i4) 
    write(*,777) '  w= ',wj,'err= ',errmax,'  u0= ',u0,'  k_j= ',kmax(1)

    write(*,777) '  w= ',wgs,'err= ',errmax,'  u0= ',u0,' k_gs= ',kmax(2)


end program Laplace2D_all
