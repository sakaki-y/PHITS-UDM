! Particle Template Version = 1.0

!=======================================================================
module udm_part_kill
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be private.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "kill" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
! integer kf_this ! The kf-code (particle ID) of this particle entered in input files.
integer n_kill ! Number of kill set
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
integer , allocatable, save :: kill_kf(:)
integer , allocatable, save :: kill_cell(:)
double precision , allocatable, save :: kill_Emin(:)
double precision , allocatable, save :: kill_Emax(:)
integer warning_proton
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
integer i
integer kf,cell
double precision Emin, Emax

warning_proton=0

if(n_kill > 0) then
  allocate( kill_kf  (n_kill) )
  allocate( kill_cell(n_kill) )
  allocate( kill_Emin(n_kill) )
  allocate( kill_Emax(n_kill) )
  do i=1,n_kill
    kf=nint(Parameters(4*(i-1)+1))
    kill_kf  (i)=kf
    kill_cell(i)=nint(Parameters(4*(i-1)+2))
    kill_Emin(i)=     Parameters(4*(i-1)+3)
    kill_Emax(i)=     Parameters(4*(i-1)+4)
  enddo
else
  open(1,file="kill.inp",status="old",err=999)
  n_kill=0
  do
888 continue
    read(1,*,err=888,end=777) kf,cell,Emin,Emax
    n_kill=n_kill+1
  enddo
777 continue
  allocate( kill_kf  (n_kill) )
  allocate( kill_cell(n_kill) )
  allocate( kill_Emin(n_kill) )
  allocate( kill_Emax(n_kill) )
  rewind(1)
  i=0
  do
666 continue
    read(1,*,err=666) kf,cell,Emin,Emax
    i=i+1
    kill_kf  (i)=kf
    kill_cell(i)=cell
    kill_Emin(i)=Emin
    kill_Emax(i)=Emax
    if(i==n_kill) exit
  enddo
  close(1)
endif

print*,"'kill proton' is not available now."
print "(A17,A10,A11,A11)","kf","cell","Emin","Emax"
do i=1,n_kill
  print "(A,I10,I10,E11.3,E11.3)", " kill: ",kill_kf(i),kill_cell(i),kill_Emin(i),kill_Emax(i)
enddo

return

999 continue
print*,"Error: kill.  Make kill.inp"
stop
return

end subroutine initialize



! !=======================================================================
! double precision function mass() ! Unit: MeV
! !=======================================================================
! return
! end

!=======================================================================
double precision function lifetime() ! Unit: Second
! The value should be greater than 0.
!=======================================================================
integer i

lifetime=0d+0

if(udm_kf_for_11==2212) return

do i=1,n_kill
  if(kill_kf  (i)==udm_kf_for_11 .or. kill_kf  (i)==0) then 
  if(kill_cell(i)==get_cell(1)) then 
  if(kill_Emin(i)> udm_Kin) then 
  if(kill_Emax(i)< udm_Kin) then 
    lifetime=1d-30
    return
  endif
  endif
  endif
  endif
enddo



return
end

!=======================================================================
subroutine decay
!=======================================================================
! [Variables available in this subroutine]
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------
call initialize_udm_event_info
set_final_state_number = 1
set_kf                 (1) = 12
set_Total_Energy_in_MeV(1) = 1d-6
set_Px_in_MeV          (1) = 1d-6
set_Py_in_MeV          (1) = 0d0
set_Pz_in_MeV          (1) = 0d0
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
! if( action .eq. 11 .and. udm_kf_for_11 .eq. kf_this ) call fill_lifetime
  if( action .eq. 11                                  ) call fill_lifetime
! if( action .eq. 12 .and. udm_kf_for_12 .eq. kf_this ) call fill_mass
! if( action .eq. 21 .and. udm_kf_for_21 .eq. kf_this ) call decay
  if( action .eq. 21                                  ) call decay
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
! kf_this=udm_part_kf(index)
! write(*,"(a,I7,a,a)") ," Set [ User defined particle ]: kf=",kf_this," ",trim(Name)
n_kill=udm_part_kf(index)
write(*,"(a,a)") ," Set [ User defined particle ]: ",trim(Name)
end subroutine fill_Parameters_and_kf

!=======================================================================
subroutine fill_lifetime
!=======================================================================
udm_lifetime = lifetime()
end subroutine fill_lifetime

! !=======================================================================
! subroutine fill_mass
! !=======================================================================
! udm_mass = mass()
! end subroutine fill_mass

!=======================================================================
end module udm_part_kill
!=======================================================================



