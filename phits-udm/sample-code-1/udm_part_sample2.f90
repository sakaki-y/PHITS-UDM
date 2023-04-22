! Particle Template Version = 1.0

!=======================================================================
module udm_part_sample_2
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be public.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "my_particle_2" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
integer kf_this ! The kf-code (particle ID) of this particle entered in input files.
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
! integer i,j
! double precision x,y
contains



!=======================================================================
subroutine print_comments
!=======================================================================
! print*,"****************************************"
! print*,"[Reference] for ", trim(Name)
! print*,"    ............"
! print*,"****************************************"
end subroutine print_comments


!=======================================================================
subroutine initialize
! This subroutine is called only once at the beginning of a calculation.
!=======================================================================
end subroutine initialize



!=======================================================================
double precision function mass() ! Unit: MeV
!=======================================================================
mass=Parameters(1) ! MeV
return
end

!=======================================================================
double precision function lifetime() ! Unit: Second
! The value should be greater than 0.
!=======================================================================
lifetime=Parameters(2) ! second
return
end

!=======================================================================
subroutine decay
!=======================================================================
! [Variables available in this subroutine]
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------
if(0.5 > get_random_0to1()) then
  call two_body_decay_uniform(12,12)      ! X -> 2 neutrinos (50%)
else
  call three_body_decay_uniform(12,12,12) ! X -> 3 neutrinos (50%)
endif
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
end module udm_part_sample_2
!=======================================================================



