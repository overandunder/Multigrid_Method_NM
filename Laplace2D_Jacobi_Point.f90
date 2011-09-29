program Laplace2D_Jacobi_Point
! Solves 2D Laplace eqn using Jacobi point method with SOR

use tecplot_write
implicit none

    real(rprec), parameter :: pi = 3.1415926535
    real(rprec), parameter :: xmin = 0.0, xmax = 2.0*pi
    real(rprec), parameter :: ymin = 0.0, ymax = 2.0*pi
    real(rprec), parameter :: dx = (xmax-xmin)/128, dy = (ymax-ymin)/128
    real(rprec), parameter :: errmax = 1.0e-5
    real(rprec), parameter :: wj = 1.0

    real(rprec), pointer, dimension(:,:) :: u, unew
    real(rprec), pointer, dimension(:) :: x, y
    real(rprec) :: error, du, ran3
    integer :: Nx, Ny, i, j, k
    real clock_start, clock_end

! Empty output directory
    call system("rm ./J_output/*")

! Call timer
    call cpu_time(clock_start)

! Determine number of gridpoints (0:N) where BCs are applied at 0 and N
    Nx = (xmax-xmin)/dx
    Ny = (ymax-ymin)/dy
    write(*,*) 'Nx,Ny',Nx,Ny

! Allocate arrays, initialize
    allocate(u(0:Nx,0:Ny),unew(0:Nx,0:Ny))
    allocate(x(0:Nx),y(0:Ny))

! Create x,y arrays
    do i=0,Nx
        x(i) = xmin + i*dx
    enddo
    do j=0,Ny
        y(j) = ymin + j*dy
    enddo  

! Set u using random initial condition
    do j=1,Ny-1
    do i=1,Nx-1 
        u(i,j) = ran3(i*j*clock_start)
    enddo
    enddo
    call write_u(0)

! Apply BCs
    u(0,:) = 0.0
    u(Nx,:) = 0.0
    u(:,Ny) = 0.0
    do i=1,Nx-1
        u(i,0) = sin(4.*x(i))
    enddo
    unew = u

! Jacobi Iterate
    k = 0
    open(unit=1,file='J_output/J_residual.dat',action='write',position='rewind')
    write(1,*) 'iteration number, residual'
    close(1)

    ! Compute error
        call error_u(k)
        open(unit=1,file='J_output/J_residual.dat',action='write',position='append')
        write(1,*) k,error
        close(1)

    do while (error .gt. errmax)
        ! Iterate once
            do i=1,Nx-1
                do j=1,Ny-1
                    du = wj * (0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)) - u(i,j))
                    unew(i,j) = u(i,j) + du
                enddo
            enddo
            k = k + 1
            u = unew

        !Write every 100 u's to file
            if (mod(k,100).eq.0) then
                call write_u(k)
            endif

        ! Compute error
            call error_u(k)
            open(unit=1,file='J_output/J_residual.dat',action='write',position='append')
            write(1,*) k,error
            close(1)
    enddo

    ! Check that solution converged
        if (error .gt. 100) k = 0    

! Call timer
    call cpu_time(clock_end)
    write ( *, * ) 'Jacobi Elapsed CPU time = ', clock_end - clock_start

! Write kmax to screen
    777 format (a5,f6.4,a7,es14.7,a7,i4) 
    write(*,777) '  w= ',wj,'err= ',errmax,'  k_j= ',k

contains

!*************************************************************
subroutine error_u(loop)
!*************************************************************
! 
integer, intent(in) :: loop

    error = 0.
    do j=1,Ny-1
    do i=1,Nx-1
        error = error + ( u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4.*u(i,j) )**2
    enddo
    enddo
    error = sqrt(error)
    write(*,*) 'Iter=',loop,'    ERROR=', error

end subroutine error_u

!*************************************************************
subroutine write_u(loop)
!*************************************************************
! 
integer, intent(in) :: loop
real(rprec) :: loop2
 character(64) :: fname, temp

    fname = 'J_output/J_u'

    if (loop.le.9) then
        write(temp,'(i1)') loop
    elseif (loop.le.999) then
        write(temp,'(i3)') loop
    elseif (loop.le.9999) then
        write(temp,'(i4)') loop                
    elseif (loop.le.99999) then
        write(temp,'(i5)') loop
    endif                

    fname = trim(fname) // temp
    fname = trim(fname) // '.dat'

    loop2 = 1.*real(loop)
    
    call write_tecplot_header_ND(fname,'rewind',3,(/Nx-1,Ny-1/),'"x","y","J_u"',1,1,loop2)
    call write_real_data_2D(fname,'append','formatted',1,Nx-1,Ny-1, &
        (/u(1:Nx-1,1:Ny-1)/),0,x(1:Nx-1),y(1:Ny-1))    

end subroutine write_u

end program Laplace2D_Jacobi_Point
