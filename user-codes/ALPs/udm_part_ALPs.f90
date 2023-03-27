!=======================================================================
module udm_part_ALPs
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none

!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "ALPs"
!-----------------------------------------------------------------------
double precision, allocatable,save:: Parameters(:)
integer, save :: kf_this
!-----------------------------------------------------------------------
double precision ga

contains

!=======================================================================
subroutine initialize
!=======================================================================
ga=Parameters(2) ! a gamma gamma coupling [1/MeV]
end subroutine initialize

!=======================================================================
function mass()
!        Unit: MeV
!=======================================================================
double precision mass
mass = Parameters(1)
return
end

!=======================================================================
function lifetime()
!        Unit: Second
!=======================================================================
double precision lifetime
if(mass()==0d0 .or. ga==0d0) then
    lifetime=1e+30
else
    lifetime=1.32E-19/(mass()**3 * ga**2) ! second. mass[MeV], coupling[1/MeV].
endif
return
end

!=======================================================================
subroutine decay
!=======================================================================
! udm_Kin:  Kinetic Energy
!-----------------------------------------------------------------------
call two_body_decay_uniform(22,22)
!-----------------------------------------------------------------------
end subroutine decay



















!=======================================================================
!     DO NOT CHANGE BELOW
!=======================================================================
subroutine caller(action,index)
!=======================================================================
integer action,index
if( action .eq.  2 ) call check_match_name(index)
if(udm_part_name(index) .ne. Name) return
if( action .eq.  1 ) call fill_Parameters_and_kf(index)
if( action .eq.  3 ) call initialize
if( action .eq. 11 .and. udm_kf_for_11 .eq. kf_this ) call fill_lifetime
if( action .eq. 12 .and. udm_kf_for_12 .eq. kf_this ) call fill_mass
if( action .eq. 21 .and. udm_kf_for_21 .eq. kf_this ) call decay
end subroutine caller

!=======================================================================
subroutine check_match_name(index)
!=======================================================================
integer index
if(udm_part_name(index) .eq. Name) udm_integer = udm_integer+1
end subroutine check_match_name

!=======================================================================
subroutine fill_Parameters_and_kf(index)
!=======================================================================
integer index
integer i
allocate( Parameters( udm_part_param_nMax ) )
do i=1,udm_part_param_nMax
  Parameters(i)=udm_part_param(index,i)
enddo
kf_this=udm_part_kf(index)
write(*,"(a,I7,a,a)") ," Set [ User defined particle ]: kf=",kf_this," ",trim(Name)
end subroutine fill_Parameters_and_kf

!=======================================================================
subroutine fill_lifetime
!=======================================================================
udm_lifetime = lifetime()
end subroutine fill_lifetime

!=======================================================================
subroutine fill_mass
!=======================================================================
udm_mass = mass()
end subroutine fill_mass

!=======================================================================
end module udm_part_ALPs
!=======================================================================



