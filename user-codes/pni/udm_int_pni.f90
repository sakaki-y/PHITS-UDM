! Interaction Template Version = 1.0

!=======================================================================
module udm_int_pni
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be private.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "pni" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
integer, parameter :: num_initial = 1 ! The number of incident particles causing this interaction.
integer, save :: kf_initial(num_initial) = (/ 22 /) ! The kf-codes (particle IDs) of the incident particles causing this interaction.
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
integer, parameter :: nZ = 6
integer, parameter :: ZminUseData = 26
integer, parameter :: Z_array(nZ) = (/ 26,41,50,60,74,82 /)
double precision, parameter :: KinDataMin = 140d0
double precision, parameter :: KinDataMax = 1d+6

integer, parameter :: nkfMax = 2500
integer, parameter :: nKinMax = 200
integer         , save :: nKin_array(nZ)=0
double precision, save ::  Kin_array(nZ,nKinMax)=-1d0
integer         , save ::   kf_array(nZ,nKinMax,nkfMax)=0
double precision, save ::  mul_array(nZ,nKinMax,nkfMax)=-1d0
contains


!=======================================================================
subroutine print_comments
!=======================================================================
! print*,"****************************************"
print*,"Note:",trim(Name)
print*,"Z containing data: 26,41,50,60,74,82"
print*,"Energy range: 140 MeV - 1 TeV"
! print*,"****************************************"
end subroutine print_comments


!=======================================================================
subroutine initialize
! This subroutine is called only once at the beginning of the calculation.
!=======================================================================
integer i,iZ,Z,iKin,nkf,ikf,kf
double precision Kin,mul

! Read:
!   nKin_array
!    Kin_array
!     kf_array
!    mul_array
do iZ=1, nZ
  Z=Z_array(iZ)
  open(1,file=trim(udm_phits_path)//"/src-udm/udm_int_pni_data/data_for_phits/"//trim(udm_int2char(Z))//"/multi.dat",status="old")
  read(1,*) nKin_array(iZ)
  do iKin=1, nKin_array(iZ)
    read(1,*) Kin, nkf
    Kin_array(iZ,iKin)=Kin
    do ikf=1, nkf
      read(1,*) kf_array(iZ,iKin,ikf), mul_array(iZ,iKin,ikf)
    enddo
  enddo
  close(1)
enddo

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
double precision phxs
Xsec_per_atom=phxs(Z,A-Z,Kin)

! Comment out:
! if( action .eq. 11 ) call calc_averaged_Xsec

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
integer Z,dZ,i,itmp,iZ,iZ_,iKin,iKin_,ikf,kf,kfZ,kfA,mul1,nx,ny,ix,iy,imul,isamp
double precision dKinMin,dKin,mul,mul2,Kin,xmin,xmax,ymin,ymax,x,y,dx,dy
double precision, allocatable :: P(:)
! double precision ::  P(20000)=0d0
double precision Kout,theta,phi,mass,Etot,absp
!-----------------------------------------------------------------------

! Add:
! if( action .eq. -22 ) call generate_final_state

! return to normal pni routine
! mat_Z(1) is true Z
if(udm_Kin < KinDataMin) return
if(udm_Kin > KinDataMax) return
if(mat_Z(1) < ZminUseData) return

! get Z and iZ
dZ=1000
do iZ_ = 1, nZ
  if(  abs(mat_Z(1)-Z_array(iZ_)) < dZ ) then
    dZ=abs(mat_Z(1)-Z_array(iZ_))
    iZ=iZ_
     Z=Z_array(iZ_) ! Z holding data
  endif
enddo

! get iKin
dKinMin=1d+10
do iKin_ = 1, nKin_array(iZ)
  dKin=dabs(udm_Kin - Kin_array(iZ,iKin_))
  if(dKin<dKinMin) then
    dKinMin=dKin
    iKin=iKin_
     Kin=Kin_array(iZ,iKin_)
  endif
enddo

! initialize
call initialize_udm_event_info
set_final_state_number = 0

! sampling for each kf
do ikf=1,nkfMax

  ! get kf, ikf
  kf = kf_array(iZ,iKin,ikf)
  mul=mul_array(iZ,iKin,ikf)
  if(kf==0) exit ! finish sampling

  ! get multiplicity
  mul1=int(mul)
  mul2=mul-mul1
  if(mul2 > get_random_0to1()) mul1=mul1+1

  ! "mul1" is the multiplicity of this "kf" 
  if(mul1==0) cycle

  ! search table
  open(1,file=trim(udm_phits_path)//"/src-udm/udm_int_pni_data/data_for_phits/"&
    &//trim(udm_int2char(Z))//"/"&
    &//trim(udm_int2char(nint(Kin)))//"/"&
    &//trim(udm_int2char(kf))//".dat",status="old",err=999)

  ! get grid size
  read(1,*) nx,ny
  if(mod(nx*ny,10) .ne. 0) then
    print*,"*** error: udm_int_pni: mod(nx*ny,10) .ne. 0"
    stop
  endif

  ! read table for this Kin
  allocate( P(nx*ny) )
  read(1,*) P
  close(1)

  ! sampling mul1 times
  do imul=1, mul1

    do isamp=1,1000000
      ix=get_random_int(1,nx)
      iy=get_random_int(1,ny)
      itmp=(ix-1)*ny+iy
      if( P(itmp) > get_random_0to1() ) exit
    enddo

    ! calculate x,y
    ! x = -sqrt(-log(Kout/Kin))
    ! y =  sqrt(theta)
    ! <=>
    ! Kin*exp(-(-x)**2) = Kout
    ! y**2              = theta
    xmin=-sqrt(-log(0.1/Kin))
    xmax=0.0
    ymin=0.0
    ymax=sqrt(3.1415)
    dx=(xmax-xmin)/nx
    dy=(ymax-ymin)/ny
    x=get_random(xmin+dx*(ix-1), xmin+dx*ix)
    y=get_random(ymin+dy*(iy-1), ymin+dy*iy)

    ! calculate Kout,theta
    Kout=udm_Kin*exp(-(-x)**2)
    theta=y**2
    phi=get_random(0.0d0,2.0*3.1415d0)
    mass=get_mass(kf)
    Etot=Kout+mass
    absp=sqrt(dabs(Etot**2-mass**2))

    ! fill
    set_final_state_number=set_final_state_number+1
    set_kf                 (set_final_state_number) = kf
    set_Total_Energy_in_MeV(set_final_state_number) = Etot
    set_Px_in_MeV          (set_final_state_number) = absp*sin(theta)*cos(phi)
    set_Py_in_MeV          (set_final_state_number) = absp*sin(theta)*sin(phi)
    set_Pz_in_MeV          (set_final_state_number) = absp*cos(theta)
  enddo

  ! finish x,y sampling
  deallocate( P )
  cycle

999 continue
  ! not x,y sampling (Heavy Nucleus)

  kfZ=kf/1000000
  kfA=mod(kf,1000000)
  if(kfZ > mat_Z(1)) cycle
  if(kfA > mat_A(1)) cycle

  do imul=1, mul1
    ! calculate Kout,theta
    Kout=0.1d0
    theta=get_random(0.0d0,3.1415d0)
    phi=get_random(0.0d0,2.0*3.1415d0)
    mass=get_mass(kf)
    Etot=Kout+mass
    absp=sqrt(dabs(Etot**2-mass**2))

    ! go to fill
    set_final_state_number=set_final_state_number+1
    set_kf                 (set_final_state_number) = kf
    set_Total_Energy_in_MeV(set_final_state_number) = Etot
    set_Px_in_MeV          (set_final_state_number) = absp*sin(theta)*cos(phi)
    set_Py_in_MeV          (set_final_state_number) = absp*sin(theta)*sin(phi)
    set_Pz_in_MeV          (set_final_state_number) = absp*cos(theta)
  enddo
enddo
! finish sampling for all kf


if(set_final_state_number > 0) call fill_final_state
udm_logical=.true. ! final state is generated


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
! if( action .eq. 11 ) call calc_averaged_Xsec
if( action .eq. 21 ) call generate_final_state
if( action .eq. -22 ) call generate_final_state
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
end module udm_int_pni
!=======================================================================






