program Laplace2D_multigrid
! 2D Laplace, cell-centered
! Nx=Ny=N (must be a power of 2) and domain must be square, too.
! Solved with multigrid (GS as smoother) V-cycle

use types
use tecplot_write

implicit none

    real(rprec), parameter :: pi = 3.1415926535
    real(rprec), parameter :: xmin = 0.0, xmax = 2.0*pi
    real(rprec), parameter :: ymin = 0.0
    integer, parameter :: Ntot=128, nlev = 1
    real(rprec), parameter :: errmax = 1.0e-3
    
    integer :: i, j, k, Ntmp, loop
    character(64) :: fname, temp
    real :: clock_start, clock_end
    real(rprec) :: ymax, error

! Empty output directory
    call system("rm ./output/*")

! Start the clock
    call cpu_time(clock_start)

! Allocate array for number of levels specified  
    if (nlev .gt. floor(1+log(real(Ntot))/log(2.))) then
        write(*,*) 'Number of levels is too large.'
        stop
    endif
    nullify(all_levels_t%level_t)
    allocate(all_levels_t%level_t(nlev)) 

! Set basic variables for each level
    ymax = ymin + (xmax - xmin)     !to ensure domain is square

    nullify(all_levels_t%u)
    allocate(all_levels_t%u(0:Ntot+1,0:Ntot+1))

    do k=1,nlev
        Ntmp = Ntot / 2.0**(k-1)
        write(*,*) 'level, N',k,Ntmp
        all_levels_t%level_t(k)%N = Ntmp
        all_levels_t%level_t(k)%dx2 = ((xmax-xmin)/Ntmp)**2       

        call alloc(k,Ntmp)
    enddo

! Determine initial condition for u (randomly)
    call random_ic()

! Cycle
    loop = 0
    open(unit=1,file='output/residual.dat',action='write',position='rewind')
    write(1,*) 'iteration number, residual'
    close(1)

    ! Compute error
        call error_u(loop)
        open(unit=1,file='output/residual.dat',action='write',position='append')
        write(1,*) loop,error
        close(1)

    do while (error .gt. errmax)
        !Cycle once
            call fine_to_coarse(1,nlev)
            call coarse_to_fine(nlev,1)
            loop = loop+1

        !Write every 100 u's to file
            if (mod(loop,100).eq.0) then
                call write_u(loop)
            endif

        ! Compute error
            call error_u(loop)
            open(unit=1,file='output/residual.dat',action='write',position='append')
            write(1,*) loop,error
            close(1)
    enddo

! Stop the clock
    call cpu_time(clock_end)
    write ( *, * ) 'Elapsed CPU time = ', clock_end - clock_start

contains

!*************************************************************
subroutine alloc(k,Ntmp)
!*************************************************************
! 
!  This subroutine nullifies and allocates the arrays for a 
!  single level, 'k'.  'Ntmp' represents the number of gridpoints
!  in each direction.  The variables below are defined in types.f90
!
integer, intent(in) :: k, Ntmp

    nullify(all_levels_t%level_t(k)%x)
    nullify(all_levels_t%level_t(k)%y)
    nullify(all_levels_t%level_t(k)%BCtop)
    nullify(all_levels_t%level_t(k)%BCbot)
    nullify(all_levels_t%level_t(k)%BCleft)
    nullify(all_levels_t%level_t(k)%BCright)
    nullify(all_levels_t%level_t(k)%RHS)
    nullify(all_levels_t%level_t(k)%RHSlow)
    nullify(all_levels_t%level_t(k)%e)
    nullify(all_levels_t%level_t(k)%elow)

    allocate(all_levels_t%level_t(k)%x(Ntmp))
    allocate(all_levels_t%level_t(k)%y(Ntmp))
    allocate(all_levels_t%level_t(k)%BCtop(Ntmp))
    allocate(all_levels_t%level_t(k)%BCbot(Ntmp))
    allocate(all_levels_t%level_t(k)%BCleft(Ntmp))
    allocate(all_levels_t%level_t(k)%BCright(Ntmp))
    allocate(all_levels_t%level_t(k)%RHS(Ntmp*Ntmp))
    allocate(all_levels_t%level_t(k)%RHSlow(Ntmp*Ntmp))
    allocate(all_levels_t%level_t(k)%e(0:Ntmp+1,0:Ntmp+1))
    allocate(all_levels_t%level_t(k)%elow(0:Ntmp+1,0:Ntmp+1))

    ! Set x and y arrays
    do i=1,Ntmp
        all_levels_t%level_t(k)%x(i) = xmin + (i-0.5)*(xmax-xmin)/Ntmp
    enddo

    do j=1,Ntmp
        all_levels_t%level_t(k)%y(j) = ymin + (j-0.5)*(xmax-xmin)/Ntmp
    enddo

    ! Set BC arrays
    do i=1,Ntmp
        all_levels_t%level_t(k)%BCtop(i) = 0.
        all_levels_t%level_t(k)%BCbot(i) = sin(4.*all_levels_t%level_t(k)%x(i))
    enddo

    do j=1,Ntmp
        all_levels_t%level_t(k)%BCleft(j) = 0.
        all_levels_t%level_t(k)%BCright(j) = 0.
    enddo 

    ! Set RHS array
    all_levels_t%level_t(k)%RHS = 0.
    all_levels_t%level_t(k)%RHSlow = 0.

    !write(*,*) 'level',all_levels_t%level_t(k)%x  

end subroutine alloc


!*************************************************************
subroutine random_ic()
!*************************************************************
! 
!  This subroutine creates a random initial condition for the
!  array 'u' at level 1
!
integer :: i,j
real(rprec) :: ran3
real(rprec), pointer, dimension(:,:) :: p_u => null()
    
    p_u => all_levels_t%u
    
    ! Calculate
        do i=1,Ntot
        do j=1,Ntot 
            all_levels_t%u(i,j) = ran3(i*j*clock_start)
        enddo
        enddo

    ! Set values of u for i,j=0,Ntmp+1 using effective BCs (u0=2*uBC-u1)
        do i=1,Ntot   !(bottom and top)
            p_u(i,0) = 2.*all_levels_t%level_t(1)%BCbot(i) - p_u(i,1)
            p_u(i,Ntot+1) = 2.*all_levels_t%level_t(1)%BCtop(i) - p_u(i,Ntot)
        enddo
        do j=1,Ntot   !(left and right)
            p_u(0,j) = 2.*all_levels_t%level_t(1)%BCleft(j) - p_u(1,j)
            p_u(Ntot+1,j) = 2.*all_levels_t%level_t(1)%BCright(j) - p_u(Ntot,j)
        enddo

    ! Write to file
        call write_u(0)

end subroutine random_ic


!*************************************************************
subroutine fine_to_coarse(n_i,n_f)
!*************************************************************
! 
!  This subroutine performs the left half of the V-cycle from 
!  level 'n_i' down to level 'n_f' with initial guess of current
!  u at level n_i and current RHS from level n_i 
!
integer, intent(in) :: n_i,n_f
integer :: m
real(rprec), pointer, dimension(:,:) :: p_u => null()
    
    ! Set RHS for first level (RHS=Au)
        p_u => all_levels_t%u  
        do j=1,Ntot
        do i=1,Ntot
            all_levels_t%level_t(1)%RHS(i+(j-1)*Ntot) = p_u(i-1,j) + p_u(i+1,j) + p_u(i,j-1) + p_u(i,j+1) - 4.*p_u(i,j)
        enddo
        enddo

    ! Loop for each level
    do m=n_i,n_f
    
        ! Initial guess for error
        all_levels_t%level_t(m)%e = 0.

        ! Iterate with this initial guess
        call itsolve(m,all_levels_t%level_t(m)%N)
        !!write(*,*) 'Iterated at level ',m

        if (m .ne. n_f) then
            ! Calculate RHSlow which will be used by next-lowest level
            call residual(m,all_levels_t%level_t(m)%N)
            !!write(*,*) 'Resid calc at level ',m    

            ! Coarsen RHSlow @m and set it to RHS @(m-1)
            call restriction(m)
            !!write(*,*) 'Restricted at level ',m
        endif
    enddo

end subroutine fine_to_coarse


!*************************************************************
subroutine coarse_to_fine(n_f,n_i)
!*************************************************************
! 
!  This subroutine performs the right half of the V-cycle from 
!  level 'n_f' back up to level 'n_i' with initial guess of 
!  current u at level n_f and current RHS from level n_f 
!
integer, intent(in) :: n_i,n_f
integer :: m

    ! Iterate again at the bottom
        call itsolve(n_f,all_levels_t%level_t(n_f)%N)
        !!write(*,*) 'Iterated at level ',n_f

    ! Loop for each level
    do m=n_f,n_i+1,-1
        ! Use e @m to update e @(m-1)
        call prolongation(m)
        !!write(*,*) 'Prolonged at level ',m

        ! Iterate with corrected value of e to get an
        !   even better approximation which will be used
        !   to correct e at the next-finest level.
        !   RHS is unchanged from last time.
        call itsolve(m-1,all_levels_t%level_t(m-1)%N)
        !!write(*,*) 'Iterated at level ',m-1
    enddo

    ! Final update on u
        all_levels_t%u = all_levels_t%u - all_levels_t%level_t(1)%e

end subroutine coarse_to_fine


!*************************************************************
subroutine itsolve(p,Ntmp)
!*************************************************************
! 
!  This subroutine performs one iteration of Gauss-Seidel 
!  on the array 'e' for level 'p' with size 'Ntmp'^2 and source
!  term given by 'RHS'
!
integer, intent(in) :: p, Ntmp
integer :: i, j
real(rprec), pointer, dimension(:,:) :: p_e => null()
real(rprec), pointer, dimension(:) :: p_s => null()
real(rprec) :: dxsq
    
    p_e => all_levels_t%level_t(p)%e     
    p_s => all_levels_t%level_t(p)%RHS 
    dxsq = all_levels_t%level_t(p)%dx2  

    ! Perform iteration
        do j=1,Ntmp
        do i=1,Ntmp
            p_e(i,j) = 0.25*( p_e(i-1,j)+p_e(i+1,j)+p_e(i,j-1)+p_e(i,j+1) - dxsq*p_s(i+(j-1)*Ntmp) )
        enddo  
        enddo
   
    ! Update values of e for i,j=0,Ntmp+1 using effective BCs (u0=-u1)
        do i=1,Ntmp   !(bottom and top)
            p_e(i,0) = -1.*p_e(i,1)
            p_e(i,Ntmp+1) = -1.*p_e(i,Ntmp)
        enddo
        do j=1,Ntmp   !(left and right)
            p_e(0,j) = -1.*p_e(1,j)
            p_e(Ntmp+1,j) = -1.*p_e(Ntmp,j)
        enddo

end subroutine itsolve


!*************************************************************
subroutine residual(p,Ntmp)
!*************************************************************
! 
!  This subroutine computes the residual at level 'p' which 
!  will be the RHS of level 'p-1' once it is coarsened.  The
!  residual will be stored in 'RHSlow' at level 'p'.
!
!  RHSlow(p) = A*e(p) - RHS(p)  where A*u=0 is the Laplace eqn
!   then
!  RHS(p-1) = coarsened{ RHSlow(p) } in subroutine 'restriction'
!
integer, intent(in) :: p, Ntmp
integer :: i, j
real(rprec), pointer, dimension(:,:) :: p_e => null()
real(rprec), pointer, dimension(:) :: p_rhs => null(), p_rhslow => null()

    p_e => all_levels_t%level_t(p)%e
    p_rhs => all_levels_t%level_t(p)%RHS
    p_rhslow => all_levels_t%level_t(p)%RHSlow
  
    do j=1,Ntmp
    do i=1,Ntmp
        p_rhslow(i+(j-1)*Ntmp) = p_e(i-1,j) + p_e(i+1,j) + p_e(i,j-1) + p_e(i,j+1) - 4.*p_e(i,j)  &
                                - p_rhs(i+(j-1)*Ntmp)
    enddo
    enddo

end subroutine residual

!*************************************************************
subroutine restriction(p)
!*************************************************************
! 
!  This subroutine restricts the residual from level 'p' to 
!  level 'p+1'.  AKA: it takes RHSlow at level 'p' and coarsens
!  it to find RHS at level 'p+1'.  Here 1=p and 2=(p+1).
!
integer, intent(in) :: p
integer :: i, j, i2, j2, Ntmp1, Ntmp2
real(rprec), pointer, dimension(:) :: p_res1 => null(), p_res2 => null()

    p_res1 => all_levels_t%level_t(p)%RHSlow
    p_res2 => all_levels_t%level_t(p+1)%RHS
    Ntmp1 = all_levels_t%level_t(p)%N 
    Ntmp2 = all_levels_t%level_t(p+1)%N

! Add four nearest neighbors (good for conservation)
    do j2=1,Ntmp2
    do i2=1,Ntmp2
        i = 1+2*(i2-1)
        j = 1+2*(j2-1)
        p_res2(i2+(j2-1)*Ntmp2) =  (p_res1(i+(j-1)*Ntmp1)+p_res1(i+1+(j-1)*Ntmp1)+ &
                                    p_res1(i+j*Ntmp1)+p_res1(i+1+j*Ntmp1))
    enddo
    enddo

end subroutine restriction

!*************************************************************
subroutine prolongation(p)
!*************************************************************
! 
!  This subroutine interpolates the solution 'e' from level 'p'
!  to level 'p-1'.  Here 1=p and 0=(p-1)
!
!  The value of 'e' at level 'p' is first refined and stored 
!  in the array 'elow' at level 'p-1'.  This 'elow' is then
!  added to 'e' at level 'p-1' as the update.
!
integer, intent(in) :: p
integer :: i0, j0, i, j, Ntmp0, Ntmp1
real(rprec), pointer, dimension(:,:) :: p_e0 => null(), p_e1 => null()
real(rprec), pointer, dimension(:,:) :: p_elow0 => null()

    p_e0 => all_levels_t%level_t(p-1)%e
    p_elow0 => all_levels_t%level_t(p-1)%elow
    p_e1 => all_levels_t%level_t(p)%e
    Ntmp0 = all_levels_t%level_t(p-1)%N 
    Ntmp1 = all_levels_t%level_t(p)%N

! Take 1/4 value of closest neighbor (to whom it contributed during restriction)
    do j0=2,Ntmp0-1
    do i0=2,Ntmp0-1
        i = i0/2.   ! nearest gridpoint at left
        j = j0/2.   ! nearest gridpoint below
        i = i + mod(i0,2)
        j = j + mod(j0,2)
        p_elow0(i0,j0) = 0.25*p_e1(i,j)
    enddo
    enddo

! Finalize the correction; new value for 'e' at level 'p-1'
    p_e0 = p_e0 - p_elow0

end subroutine prolongation



!*************************************************************
subroutine error_u(loop)
!*************************************************************
! 
integer, intent(in) :: loop
real(rprec), pointer, dimension(:,:) :: p_u => null()

    p_u => all_levels_t%u  
    error = 0.
    do j=1,Ntot
    do i=1,Ntot
        error = error + ( p_u(i-1,j) + p_u(i+1,j) + p_u(i,j-1) + p_u(i,j+1) - 4.*p_u(i,j) )**2
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

    fname = 'output/u'

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
    
    call write_tecplot_header_ND(fname,'rewind',3,(/Ntot,Ntot/),'"x","y","u"',1,1,loop2)
    call write_real_data_2D(fname,'append','formatted',1,Ntot,Ntot, &
        (/all_levels_t%u(1:Ntot,1:Ntot)/),0,all_levels_t%level_t(1)%x,all_levels_t%level_t(1)%y)    


end subroutine write_u

end program Laplace2D_multigrid
