! Interaction Template Version = 1.0

!=======================================================================
module udm_int_neutrino
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be private.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "neutrino" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
integer, parameter :: num_initial = 2 ! The number of incident particles causing this interaction.
integer, save :: kf_initial(num_initial) = (/ 12, -14 /) ! The kf-codes (particle IDs) of the incident particles causing this interaction.
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
double precision, parameter :: KinDataMin = 100d3 ! 100 GeV

integer, parameter :: nKinMax = 100
integer         , save :: nKin_array(num_initial)=0
double precision, save ::  Kin_array(num_initial,nKinMax)=-1d0
integer         , save ::    nEvents(num_initial,nKinMax)=0

double precision Bias
integer nEventsPerFile
contains


!=======================================================================
subroutine print_comments
!=======================================================================
print*,"****************************************"
print*,"[Reference] for ", trim(Name)
print*,"The GENIE Neutrino Monte Carlo Generator"
print*,"Nucl.Instrum.Meth.A 614 (2010) 87-104"
print*,"https://inspirehep.net/literature/820590"
print*,"****************************************"
end subroutine print_comments


!=======================================================================
subroutine initialize
! This subroutine is called only once at the beginning of the calculation.
!=======================================================================
integer iNeu,nE,iKin
double precision Kin

Bias=Parameters(1)

do iNeu=1, num_initial
  open(1,file=trim(udm_phits_path)//"/src-udm/udm_int_neutrino/data_udm/"&
    &//trim(udm_int2char(kf_initial(iNeu)))//"/n_events_list.dat",status="old")
  read(1,*) nKin_array(iNeu)
  do iKin=1, nKin_array(iNeu)
    read(1,*) Kin_array(iNeu,iKin), nEvents(iNeu,iKin)
  enddo
  close(1)
enddo

open(1,file=trim(udm_phits_path)//"/src-udm/udm_int_neutrino/data_udm/"&
  &//"/nEventsPerFile.dat",status="old")
read(1,*) nEventsPerFile
close(1)

end subroutine initialize



!=======================================================================
double precision function Xsec_per_atom(Kin,Z,A)
! Integrated cross section per an atom.
! Unit: barn (10^-24 cm2)
implicit none
double precision Kin ! Kinetic energy of incident particle [MeV]
           integer Z ! Atomic number of target atom
           integer A ! Mass   number of target atom
!-----------------------------------------------------------------------
! [Variables available in this function]
! udm_kf_incident : The kf-codes (particle IDs) of the incident particles.
!=======================================================================
double precision fac

if(Kin < KinDataMin) then
  Xsec_per_atom=0.0
  return
endif

fac=Bias*(A**0.92)*Kin

if     (udm_kf_incident ==  12) then
  Xsec_per_atom=fac*1e-17
else if(udm_kf_incident == -12) then
  Xsec_per_atom=fac*1e-17*0.5d0
else if(udm_kf_incident ==  14) then
  Xsec_per_atom=fac*1e-17
else if(udm_kf_incident == -14) then
  Xsec_per_atom=fac*1e-17*0.5d0
else
  print*,"error"
  stop
endif

return
end




!=======================================================================
subroutine generate_final_state
!=======================================================================
! Subroutine to determine final state information (4-momenta, etc.).
! In this example, the X and e+- momenta of the final state are determined.
!-----------------------------------------------------------------------
! [Variables available in this subroutine]
! udm_kf_incident : The kf-codes (particle IDs) of the incident particles.
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------
integer i,j,nFS,iEvent,iKinData,iNeu,kf
double precision KinData,rKout,Kout,the,phi,mass,Etot,absp
double precision dKinMin,dKin
integer iKin_,iFile,iev

! get iNeu
iNeu=0
do i = 1, num_initial
  if( udm_kf_incident==kf_initial(i) ) then
    iNeu=i
  endif
enddo
if(iNeu==0) then
  print*,"ERROR in generate_final_state"
  stop
endif


! get iKin
dKinMin=1d+10
do iKin_ = 1, nKin_array(iNeu)
  dKin=dabs(udm_Kin - Kin_array(iNeu,iKin_))
  if(dKin<dKinMin) then
    dKinMin=dKin
    iKinData=iKin_
     KinData=Kin_array(iNeu,iKin_)
  endif
enddo


iEvent=get_random_int(1,nEvents(iNeu,iKinData))
iFile=1+iEvent/nEventsPerFile
iev  =1+mod(iEvent,nEventsPerFile)
print*,"iEvent__",iEvent,iFile,iev

open(1,file=trim(udm_phits_path)//"/src-udm/udm_int_neutrino/data_udm/"&
  &//trim(udm_int2char(udm_kf_incident))//"/"&
  &//trim(udm_int2char(nint(KinData)))//"/"&
  &//trim(udm_int2char(iFile))//".dat",status="old",err=999)

do i=1,iev-1
  read(1,*) nFS
  do j=1,nFS
    read(1,'()')
  enddo
enddo

call initialize_udm_event_info
read(1,*) set_final_state_number
do i=1,set_final_state_number
  read(1,*) rKout,the,phi,kf
  Kout=rKout*udm_Kin


  ! calculate kinematic variables
  phi=phi+get_random(0.0d0,2.0*3.1415d0)
  mass=get_mass(kf)
  Etot=Kout+mass
  absp=sqrt(dabs(Etot**2-mass**2))

  ! fill
  print*,kf,Etot,absp,mass,the,phi
  set_kf                 (i) = kf
  set_Total_Energy_in_MeV(i) = Etot
  set_Px_in_MeV          (i) = absp*sin(the)*cos(phi)
  set_Py_in_MeV          (i) = absp*sin(the)*sin(phi)
  set_Pz_in_MeV          (i) = absp*cos(the)
enddo
close(1)

call multiply_weight_by(1d0/Bias)
call fill_final_state
return

999 continue
print*,"ERROR: no file, udm_int_neutrino"
stop

end subroutine generate_final_state

























!=======================================================================
!     DO NOT CHANGE BELOW
!=======================================================================
subroutine caller(action,index)
!=======================================================================
integer action,index
if( action .eq.  2 ) call check_match_name(index)
if(udm_int_name(index) .ne. Name) return
if( action .eq.  1 ) call fill_Parameters(index)
if( action .eq.  3 ) call initialize
if( action .eq.  4 ) call print_comments
if(.not. match_initial()) return
if( action .eq. 11 ) call calc_averaged_Xsec
if( action .eq. 21 ) call generate_final_state
end subroutine caller

!=======================================================================
function match_initial()
!=======================================================================
logical match_initial
integer i
match_initial = .false.
do i = 1, num_initial
  if(udm_kf_incident .eq. kf_initial(i)) then
    match_initial = .true.
    exit
  endif
enddo
return
end

!=======================================================================
subroutine check_match_name(index)
!=======================================================================
integer index
udm_logical = udm_logical .or. (udm_int_name(index) .eq. Name)
end subroutine check_match_name

!=======================================================================
subroutine fill_Parameters(index)
!=======================================================================
integer index
integer i
allocate( Parameters( udm_int_param_nMax ) )
do i=1,udm_int_param_nMax
  Parameters(i)=udm_int_param(index,i)
enddo
print*,"Set [ User defined interaction ]: ",trim(Name)
end subroutine fill_Parameters

!=======================================================================
subroutine calc_averaged_Xsec
!     Averaged cross section for mixed atoms
!=======================================================================
integer i
do i = 1, num_nuclide
  udm_sigt = udm_sigt + mat_ratio(i)*Xsec_per_atom( udm_Kin, mat_Z(i), mat_A(i) )
enddo
end subroutine calc_averaged_Xsec

!=======================================================================
function get_hit_nuclide_Z_A(Kin)
!=======================================================================
integer i, get_hit_nuclide_Z_A(2)
double precision Kin, randTMP, XsecSum, tmp
!     ------------------------------------------------------------------
XsecSum=0d0
do i = 1, num_nuclide
  XsecSum = XsecSum + mat_ratio(i)*Xsec_per_atom( Kin, mat_Z(i), mat_A(i) )
enddo
!     ------------------------------------------------------------------
randTMP=get_random_0to1()
tmp=0d0
do i = 1, num_nuclide
  tmp = tmp + ( mat_ratio(i)*Xsec_per_atom( Kin, mat_Z(i), mat_A(i) ) ) / XsecSum
  if(randTMP .lt. tmp) then ! nuclide=i is accepted.
    get_hit_nuclide_Z_A(1)=mat_Z(i) ! Z of hit nuclide
    get_hit_nuclide_Z_A(2)=mat_A(i) ! A of hit nuclide
    return
  endif
enddo
!     ------------------------------------------------------------------
print*,"[CAUTION] Something wrong. hit_nuclide was not choosed."
get_hit_nuclide_Z_A(1)=mat_Z(1)
get_hit_nuclide_Z_A(2)=mat_A(1)
return
end

!=======================================================================
end module udm_int_neutrino
!=======================================================================






