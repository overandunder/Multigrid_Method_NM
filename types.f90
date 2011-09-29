module types

public :: rprec
integer,parameter :: rprec = kind (1.d0) 

type sublevel
    integer :: N
    real(rprec) :: dx2
    real(rprec), pointer, dimension(:) :: x, y
    real(rprec), pointer, dimension(:) :: BCtop, BCbot, BCleft, BCright
    real(rprec), pointer, dimension(:) :: RHS, RHSlow
    real(rprec), pointer, dimension(:,:) :: e, elow
end type sublevel

type all_levels
    real(rprec), pointer, dimension(:,:) :: u
    type(sublevel), pointer, dimension(:) :: level_t
end type all_levels

type(all_levels) :: all_levels_t

end module types

