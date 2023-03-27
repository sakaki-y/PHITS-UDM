! Particle Template Version = 1.0

!=======================================================================
module udm_part_MuonPlus
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be private.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "MuonPlus" ! This 'Name' is used in input files.
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
! print*,"[Reference] for ", Name
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
return
end

!=======================================================================
double precision function lifetime() ! Unit: Second
! The value should be greater than 0.
!=======================================================================
integer iblz1,iblz2
common /tlgeom/ iblz1,iblz2
if(iblz1==1) then
  lifetime=2.2e-10
else
  lifetime=2.2e-6
endif
return
end

!=======================================================================
double precision function get_x(kf)
!=======================================================================
integer kf,i
double precision x,P

x=get_random_0to1()
if    (kf==12 .or. kf==-12) then
  get_x=4d0*x**3 - 3d0*x**4
  return
elseif(kf==11 .or. kf==-11) then
  get_x=4d0*x**3 - 3d0*x**4
  return
elseif(kf==14 .or. kf==-14) then
  get_x=2d0*x**3 - 1d0*x**4
  return
else
  print*,"error: get_x"
  stop
endif

end

!=======================================================================
subroutine decay
!=======================================================================
! [Variables available in this subroutine]
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------
integer kf,i
double precision x,costh,sinth,phi,gamma,beta,Enu,mmu

call initialize_udm_event_info
set_final_state_number = 2

mmu=get_mass(13)
gamma=1d0+udm_Kin/mmu
beta=sqrt(dabs(1d0-1d0/gamma**2))

do i=1,2

  if(i==1) kf=-14 ! nu_mu_bar
  ! if(i==2) kf= 12 ! nu_e
  if(i==2) kf=-11 ! e+

  x=get_x(kf)
  costh=get_random(-1d0,1d0)
  sinth=sqrt(dabs(1d0-costh**2))
  phi=get_random(0d0,3.14159d0*2d0)

  Enu=x*mmu/2d0

  set_kf                 (i) = kf
  set_Total_Energy_in_MeV(i) = Enu*(gamma + gamma*beta*costh)
  set_Px_in_MeV          (i) = Enu*sinth*cos(phi)
  set_Py_in_MeV          (i) = Enu*sinth*sin(phi)
  set_Pz_in_MeV          (i) = Enu*(gamma*beta + gamma*costh)

enddo

call fill_final_state_for_decay

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
end module udm_part_MuonPlus
!=======================================================================



