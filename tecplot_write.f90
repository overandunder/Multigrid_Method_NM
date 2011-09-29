module tecplot_write

use types, only:rprec

contains

!*************************************************************
subroutine write_real_data_3D(fname, write_pos, write_fmt, nvars, &
  imax, jmax, kmax, vars, ibuff, x,y,z)
!*************************************************************
! 
!  This subroutine variables the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. An example in which vars 
!  is to be specified as:
!    (/ u, v, w /)
!  
!  Inputs:
!  fname (char) - file to write to
!  write_pos (char) - postition to write in file : 'append' or 'rewind'
!  write_fmt (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  jmax (int) - size of 2nd dimension of variables in vars  
!  kmax (int)- size of 3rd dimension of variables in vars
!  vars (real, vector) - vector contaning variables to write
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction
!     2 - buffer on j direction
!     3 - buffer on k direction
!     4 - buffer on i,j direction
!     5 - buffer on i,k direction
!     6 - buffer on j,k direction
!     7 - buffer on i,j,k directions
!  x,y,z (real, vector, optional) - vectors containing x,y,z coordinates 
!
implicit none

    character(*), intent(in) :: fname, write_pos, write_fmt
    integer, intent(in) :: nvars, imax, jmax, kmax
    real(rprec), intent(in), dimension(nvars*imax*jmax*kmax) :: vars
    integer, intent(in) :: ibuff
    real(rprec), intent(in), dimension(:), optional :: x,y,z

    logical :: coord_pres

    integer :: i,j,k,n
    integer :: i0, j0, k0, imax_buff, jmax_buff, kmax_buff

    !integer, allocatable, dimension(:) :: ikey_x, ikey_y, ikey_z
    integer, allocatable, dimension(:,:,:,:) :: ikey_vars
    !real(rprec), allocatable, dimension(:,:,:,:) :: vars_dim

!  Check if spatial coordinates are specified
coord_pres=.false.
if(present(x) .and. present(y) .and. present(z)) coord_pres = .true.

if( ibuff == 0 ) then

  imax_buff = imax
  jmax_buff = jmax
  kmax_buff = kmax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
  jmax_buff = jmax
  kmax_buff = kmax

elseif( ibuff == 2 ) then

  imax_buff = imax
  jmax_buff = jmax + 1
  kmax_buff = kmax
  
elseif( ibuff == 3 ) then

  imax_buff = imax
  jmax_buff = jmax
  kmax_buff = kmax + 1
  
elseif( ibuff == 4 ) then

  imax_buff = imax + 1
  jmax_buff = jmax + 1
  kmax_buff = kmax  
  
elseif( ibuff == 5 ) then

  imax_buff = imax + 1
  jmax_buff = jmax
  kmax_buff = kmax + 1  
  
elseif( ibuff == 6 ) then

  imax_buff = imax
  jmax_buff = jmax + 1
  kmax_buff = kmax + 1  
  
elseif( ibuff == 7 ) then

  imax_buff = imax + 1
  jmax_buff = jmax + 1
  kmax_buff = kmax + 1  
  
endif

allocate(ikey_vars(nvars,imax_buff,jmax_buff,kmax_buff)) 


do n=1,nvars

  do k=1,kmax_buff
  
    k0 = buff_indx(k,kmax)
    
    do j = 1, jmax_buff
  
      j0 = buff_indx(j,jmax)
    
      do i = 1, imax_buff

        i0 = buff_indx(i,imax)
      
        ikey_vars(n,i,j,k) = (n-1)*imax*jmax*kmax + (k0-1)*imax*jmax + (j0-1)*imax + i0
        
      enddo
    
    enddo
    
  enddo
  
enddo

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)
  
!  Write the data
select case(write_fmt)

  case('formatted')
  
    !  Specify output format; may want to use a global setting
    	
    if (coord_pres) then
	    
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2,'(1e)') x(i)
	        enddo 
	      enddo
	    enddo
      
      do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2,'(1e)') y(j)
	        enddo 
	      enddo
	    enddo
      
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2,'(1e)') z(k)
	        enddo 
	      enddo
	    enddo      
      
    endif      
    
    do n=1, nvars
	    do k=1, kmax_buff
        do j=1, jmax_buff
          do i=1, imax_buff
            write(2,'(1e)') vars(ikey_vars(n,i,j,k))
	        enddo 
	      enddo
	    enddo 
    enddo

  case('unformatted')
  
    if (coord_pres) then
	  
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2) x(i)
	        enddo 
	      enddo
	    enddo
      
      do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2) y(j)
	        enddo 
	      enddo
	    enddo
      
	    do k=1, kmax_buff
  	    do j=1, jmax_buff
    	    do i=1, imax_buff
            write(2) z(k)
	        enddo 
	      enddo
	    enddo  
	  
    endif
    
    do n=1, nvars
	    do k=1, kmax_buff
        do j=1, jmax_buff
          do i=1, imax_buff
            write(2) vars(ikey_vars(n,i,j,k))
	        enddo 
	      enddo
	    enddo 
    enddo
	
end select

close(2)
  
deallocate(ikey_vars)

return
end subroutine write_real_data_3D

!*************************************************************
subroutine write_real_data_2D(fname, write_pos, write_fmt, nvars, &
  imax, jmax, vars, ibuff, x, y)
!*************************************************************
! 
!  This subroutine variables the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. An example in which vars 
!  is to be specified as:
!    (/ u, v, w /)
!  
!  Inputs:
!  fname (char) - file to write to
!  write_pos (char) - postition to write in file : 'append' or 'rewind'
!  write_fmt (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  jmax (int) - size of 2nd dimension of variables in vars 
!  vars (real, vector) - vector contaning variables to write
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction
!     2 - buffer on j direction
!     3 - buffer on i and j directions
!  x,y (real, vector, optional) - vectors containing x,y coordinates 
!
implicit none

    character(*), intent(in) :: fname, write_pos, write_fmt
    integer, intent(in) :: nvars, imax, jmax
    real(rprec), intent(in), dimension(nvars*imax*jmax) :: vars
    integer, intent(in) :: ibuff
    real(rprec), intent(in), dimension(:), optional :: x, y

    logical :: coord_pres

    integer :: i,j,n
    integer :: i0, j0, imax_buff, jmax_buff
    !integer, allocatable, dimension(:) :: ikey_x, ikey_y
    integer, allocatable, dimension(:,:,:) :: ikey_vars

!  Check if spatial coordinates specified
coord_pres=.false.
if(present(x) .and. present(y)) coord_pres = .true.

if( ibuff == 0 ) then

  imax_buff = imax
  jmax_buff = jmax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
  jmax_buff = jmax

elseif( ibuff == 2 ) then

  imax_buff = imax
  jmax_buff = jmax + 1
  
elseif( ibuff == 3 ) then

  imax_buff = imax + 1
  jmax_buff = jmax + 1
  
endif

allocate(ikey_vars(nvars,imax_buff,jmax_buff)) 
  
do n=1,nvars

  do j = 1, jmax_buff
  
    j0 = buff_indx(j,jmax)
    
    do i = 1, imax_buff

      i0 = buff_indx(i,imax)
      
      ikey_vars(n,i,j) = (n-1)*imax*jmax + (j0-1)*imax + i0
    
    enddo
    
  enddo
  
enddo

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)

!  Write the data
select case(write_fmt)
  case('formatted')
  
    !  Specify output format; may want to use a global setting
    	
    if (coord_pres) then

	    do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2,'(1e)') x(i)
	      enddo 
	    enddo
    	  
      do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2,'(1e)') y(j)
	      enddo 
	    enddo
 
    endif
	  
    do n=1, nvars
      do j=1, jmax_buff
	      do i=1,imax_buff
          write(2,'(1e)') vars(ikey_vars(n,i,j))
	      enddo 
	    enddo
    enddo

  case('unformatted')
  
    if (coord_pres) then
	  
	  
	    do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2) x(i)
	      enddo 
	    enddo
    	  
      do j=1, jmax_buff
  	    do i=1,imax_buff
          write(2) y(j)
	      enddo 
	    enddo
 
    endif

    do n=1, nvars
      do j=1, jmax_buff
	      do i=1,imax_buff
          write(2) vars(ikey_vars(n,i,j))
	      enddo 
	    enddo
    enddo    
	
end select

close(2)

deallocate(ikey_vars)

return
end subroutine write_real_data_2D

!*************************************************************
subroutine write_real_data_1D(fname, write_pos, write_fmt, nvars, &
  imax, vars, ibuff, x)
!*************************************************************
! 
!  This subroutine writes the variables given by vars to the
!  specified file, fname. The number of variables can be arbitrary
!  but must be specified by nvars. An example in which vars 
!  is to be specified as:
!    (/ u, v, w /)
!  
!  Inputs:
!  fname (char) - file to write to
!  write_pos (char) - postition to write in file : 'append' or 'rewind'
!  write_fmt (char) - Fotran format flag : 'formatted' or 'unformatted'
!  nvars (int) - number of variables contained in vars
!  imax (int) - size of 1st dimension of variables in vars
!  vars (real, vector) - vector contaning variables to write
!  ibuff (int) - flag for adding buffer region due to periodicity
!     0 - no buffer region
!     1 - buffer on i direction (i = 1 = imax + 1)
!  x (real,vector,optional) - vector containing x coordinates 
!

implicit none

    character(*), intent(in) :: fname, write_pos, write_fmt
    integer, intent(in) :: nvars, imax
    real(rprec), intent(in), dimension(nvars*imax) :: vars
    integer, intent(in) :: ibuff
    real(rprec), intent(in), dimension(:), optional :: x

    logical :: coord_pres

    integer :: i,n
    integer :: i0, imax_buff
    !integer, allocatable, dimension(:) :: ikey_x
    integer, allocatable, dimension(:,:) :: ikey_vars


!  Check if spatial coordinates specified
if(present(x)) coord_pres = .true.

if( ibuff == 0 ) then

  imax_buff = imax
  
elseif( ibuff == 1 ) then

  imax_buff = imax + 1
  
endif

allocate(ikey_vars(nvars,imax_buff)) 
  
do n=1,nvars
 
  do i = 1, imax_buff

    i0 = buff_indx(i,imax)
      
    ikey_vars(n,i) = (n-1)*imax + i0
    
  enddo
  
enddo

open (unit = 2,file = fname, status='unknown',form=write_fmt, &
  action='write',position=write_pos)
   
!  Write the data
select case(write_fmt)
  case('formatted')
  
    if (coord_pres) then
	  
      do i=1,imax_buff
        write(2,'(1e)') x(i)
      enddo 
	  
    endif
    
    do n=1, nvars
      do i=1,imax_buff
        write(2,'(1e)') vars(ikey_vars(n,i))
      enddo
	  enddo 

  case('unformatted')
  
    if (coord_pres) then
	  
      do i=1,imax_buff
        write(2) x(i)
      enddo 
	  
    endif
    
    do n=1, nvars
      do i=1,imax_buff
        write(2) vars(ikey_vars(n,i))
      enddo
	  enddo 

end select

close(2)

deallocate(ikey_vars)

return
end subroutine write_real_data_1D

!*************************************************************
subroutine write_tecplot_header_xyline(fname, write_pos, var_list)
!*************************************************************
!  The purpose of this routine is to write Tecplot header information
!  for xy-line formatted files
! 
!  Inputs:
!  fname (char)     - file name to write to
!  write_pos (char) - position in file to write data: append or rewind
!  var_list	(char)  - string containing variable names: Ex. "x", "u"
!

implicit none

    character(*), intent(in) :: fname, write_pos, var_list

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position=write_pos)

write(2,'(1a)') 'variables = ' // var_list

close(2)

return
end subroutine write_tecplot_header_xyline


!*************************************************************
subroutine write_tecplot_header_ND(fname, write_pos, nvars, &
  domain_size, var_list, zone, data_type, soln_time)
!*************************************************************
!  The purpose of this routine is to write Tecplot header information
!  for 1D, 2D, and 3D data files.
!
!  NOTE: domain_size needs to be specified as a vector even for 1D data.
!  Ex: 
!    1D : (/ Nx /)
!    2D : (/ Nx, Ny /)
!    3D : (/ Nx, Ny, Nz /)
! 
!  Inputs:
!  fname (char)     - file name to write to
!  write_pos (char) - position in file to write data: append or rewind
!  nvars (int)      - number of variables
!  domain_size (int, vector) - vector containing the diminsions of the data.
!  var_list	(char)  - string containing variable names: Ex. '"x", "u"'
!  zone (int)       - zone number
!  date_type (int) 	- specify Tecplot data type (precision): 1 - single, 2 - double
!  soln_time (real, optional) - time stamp
!
implicit none

    character(*), intent(in) :: fname, write_pos
    integer, intent(in) :: nvars
    integer, dimension(:), intent(in) :: domain_size
    character(*), intent(in) :: var_list
    integer, intent(in) :: zone, data_type
    real(rprec), optional, intent(in) :: soln_time

    character(120) :: tec_dt_str, tec_dat_str
    integer :: ndims

!  Get number of dimensions for data
ndims = size(domain_size,1)

if(ndims == 1) then
  write(tec_dat_str,"(1a,i9,1a,i9,1a,i9)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1)
elseif(ndims == 2) then
  write(tec_dat_str,"(1a,i9,1a,i6,1a,i6,1a,i6)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1),', j=', domain_size(2)
elseif(ndims == 3) then
  write(tec_dat_str,"(1a,i9,1a,i3,1a,i3,1a,i3,1a,i3)") 'ZONE T="', &
    zone,'", DATAPACKING=BLOCK, i=', domain_size(1),', j=', domain_size(2),', k=', domain_size(3)
endif

!  Create Tecplot DT string
call tecplot_data_type_str(nvars, data_type, tec_dt_str)

open (unit = 2,file = fname, status='unknown',form='formatted', &
  action='write',position=write_pos)

!  Write variable list
write(2,'(1a)') 'variables = ' // var_list
!  Write data layout size information
write(2,'(1a)') tec_dat_str
!  Write Tecplot data type for each variable
write(2,'(1a)') tec_dt_str

if (present (soln_time)) then
write(2,'(1a,e)') 'solutiontime=', soln_time
endif

close(2)

  
return
end subroutine write_tecplot_header_ND

!*************************************************************
subroutine tecplot_data_type_str(nvars, data_type, tec_dt_str)
!*************************************************************
implicit none

    integer, intent(in) ::  nvars, data_type
    character(120), intent(OUT) :: tec_dt_str
    character(7) :: tec_dt

    integer :: n

!  Specify single or double precision
if(data_type == 1) then
  tec_dt = ' SINGLE'
elseif(data_type == 2) then
  tec_dt = ' DOUBLE'
endif

!  Create DT string
tec_dt_str = 'DT=('
do n=1, nvars
  tec_dt_str = trim(adjustl(tec_dt_str)) // tec_dt
enddo
  tec_dt_str = trim(adjustl(tec_dt_str)) // ')'

return
end subroutine tecplot_data_type_str

!**********************************************************************
integer function buff_indx(i,imax)
!**********************************************************************
!  This function returns the physical index associated with the buffer 
!  region for the specified i and imax. 
!  For i = imax + 1 -> 1 is returned otherwise i is returned
implicit none

integer, intent(in) :: i,imax

if(i == imax + 1) then
  buff_indx = 1
else
  buff_indx = i
endif
  
return
end function buff_indx

!**********************************************************************

end module tecplot_write
