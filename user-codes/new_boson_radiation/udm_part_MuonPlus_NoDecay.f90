! Particle Template Version = 1.0

!=======================================================================
module udm_part_MuonPlus_NoDecay
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be private.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "MuonPlus_NoDecay" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
integer kf_this ! The kf-code (particle ID) of this particle entered in input files.
contains

!=======================================================================
double precision function lifetime() ! Unit: Second
!=======================================================================
! The value should be greater than 0.
! udm_Kin: Kinetic energy of decaying particle (initial particle)
! get_cell(1): Cell ID where the particle decays
!=======================================================================
lifetime=1d+9
return
end
























!=======================================================================
!     DO NOT CHANGE BELOW
!=======================================================================
subroutine caller(action,index)
!=======================================================================
integer action,index
if( action .eq.  2 ) call check_match_name(index)
if(udm_part_name(index) .ne. Name) return
if( action .eq.  1 ) call fill_Parameters_and_kf(index)
if( action .eq. 11 .and. udm_kf_for_11 .eq. kf_this ) call fill_lifetime
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
end module udm_part_MuonPlus_NoDecay
!=======================================================================



