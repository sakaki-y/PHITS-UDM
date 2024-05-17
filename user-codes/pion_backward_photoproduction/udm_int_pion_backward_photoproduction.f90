!=======================================================================
module udm_int_pion_backward_photoproduction
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
!-----------------------------------------------------------------------
integer, parameter :: num_initial = 1
integer, save :: kf_initial(num_initial) = (/22/)
character(len=99), parameter :: Name = "pion_backward_photoproduction"
!-----------------------------------------------------------------------
double precision, allocatable,save:: Parameters(:)
!-----------------------------------------------------------------------
integer, parameter :: ndata=9
double precision, save :: kGeVL(ndata)
double precision, save :: dsdOL(ndata)
double precision, parameter :: anglePionCM_min = 140 ! c.m.系におけるpionの角度の最小値 [degree]。
double precision, parameter :: mypi = 3.1415926535d0
double precision, save :: cut
!-----------------------------------------------------------------------
contains


!=======================================================================
subroutine print_comments
!=======================================================================
! print*,"****************************************"
! print*,"[Reference] for ", trim(Name)
! print*,"****************************************"
end subroutine print_comments


!=======================================================================
subroutine initialize
!=======================================================================
kGeVL=(/0.2, 0.25, 0.3, 0.4, 0.6, 0.79, 1. , 2.2, 7.  /) ! 入射光子エネルギー(Lab frame)
dsdOL=(/9.5, 13. , 15., 4.7, 2.5, 2.5 , 2.2, 0.1, 0.01/) ! dsigma/dOmega, dOmega = dcosTheta dPhi [micro-barn/sr]
cut=cos(anglePionCM_min/180d0*mypi)
end subroutine initialize


!=======================================================================
function Xsec_per_atom(Ein,Z_atom,A_atom)
!     Return value: Total cross section per an atom.
!     Unit: barn
!     Arguments:
!       Ein: Kinetic energy of inital particle [MeV]
!       Z_atom : Atomic number of the target atom
!       A_atom : Mass   number of the target atom
!=======================================================================
implicit none
double precision Xsec_per_atom, Ein
integer Z_atom, A_atom
! -----
double precision kGeV, shadowing, dsdO
! double precision phxs

kGeV=Ein/1000d0 ! MeV -> GeV
shadowing=A_atom*(1d0-0.072d0*log(dble(A_atom)))
dsdO=shadowing*1.0d-6*udm_interp(kGeV,ndata,kGeVL,dsdOL,"lin-log") ! micro-barn -> barn
Xsec_per_atom=2d0*mypi*(1d0+cut)*dsdO

! do kGeV=0.1d0, 8.0d0, 0.05d0
!   dsdO=shadowing*1.0d-6*udm_interp(kGeV,ndata,kGeVL,dsdOL,"lin-log") ! micro-barn -> barn
!   Xsec_per_atom=2d0*mypi*(1d0+cut)*dsdO
!   print*,kGeV,Xsec_per_atom/phxs(Z_atom,A_atom-Z_atom,Ein)
! enddo
! stop

return
end

!=======================================================================
subroutine generate_final_state
! udm_Kin       : mother's Kinetic Energy [MeV]
! udm_kf_incident : mother's kf-code
!=======================================================================
double precision costh_cm, sinth_cm, s, EN1_cm, pN1_cm, gamma, gamma_beta, phi
double precision EN2   , pxN2   , pyN2   , pzN2
double precision EN2_cm, pxN2_cm, pyN2_cm, pzN2_cm, pN2_cm
integer i
double precision, parameter :: mN1 = 938.272081 ! MeV
double precision, parameter :: mN2 = 939.565413 ! MeV
double precision, parameter :: mM  = 139.57061  ! MeV
! ----------------------------------------
! photon N1 --> M N2
! ----------------------------------------

! ----------------------------------------
! sampling neutron angle in c.m. frame
costh_cm = -get_random(-1d0,cut)
sinth_cm = sqrt(dabs(1d0-costh_cm**2))

! ----------------------------------------
! kinematics
s=mN1**2+2*udm_Kin*mN1
EN1_cm=sqrt(s)/2*(1+mN1**2/s)
EN2_cm=sqrt(s)/2*(1+mN2**2/s-mM**2/s)
pN1_cm=sqrt(dabs(EN1_cm**2-mN1**2))
pN2_cm=sqrt(dabs(EN2_cm**2-mN2**2))
gamma     =EN1_cm/mN1
gamma_beta=pN1_cm/mN1

phi=get_random(0d0,2*mypi)
pxN2_cm=pN2_cm*sinth_cm*cos(phi)
pyN2_cm=pN2_cm*sinth_cm*sin(phi)
pzN2_cm=pN2_cm*costh_cm

EN2 =gamma*EN2_cm + gamma_beta*pzN2_cm
pxN2=                          pxN2_cm
pyN2=                          pyN2_cm
pzN2=gamma_beta*EN2_cm + gamma*pzN2_cm

! ----------------------------------------
! fill only neutron
call initialize_udm_event_info
set_final_state_number = 1
set_kf                 (1) = 2112
set_Total_Energy_in_MeV(1) = EN2
set_Px_in_MeV          (1) = pxN2
set_Py_in_MeV          (1) = pyN2
set_Pz_in_MeV          (1) = pzN2
call fill_final_state
return
! --------------------------------------------------------------

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
function choose_hit_nuclide_Z_A(Ein)
!=======================================================================
integer i, choose_hit_nuclide_Z_A(2)
double precision Ein, randTMP, XsecSum, tmp
!     ------------------------------------------------------------------
XsecSum=0d0
do i = 1, num_nuclide
  XsecSum = XsecSum + mat_ratio(i)*Xsec_per_atom( Ein, mat_Z(i), mat_A(i) )
enddo
!     ------------------------------------------------------------------
randTMP=get_random_0to1()
tmp=0d0
do i = 1, num_nuclide
  tmp = tmp + ( mat_ratio(i)*Xsec_per_atom( Ein, mat_Z(i), mat_A(i) ) ) / XsecSum
  if(randTMP .lt. tmp) then ! nuclide=i is accepted.
    choose_hit_nuclide_Z_A(1)=mat_Z(i) ! Z of hit nuclide
    choose_hit_nuclide_Z_A(2)=mat_A(i) ! A of hit nuclide
    return
  endif
enddo
!     ------------------------------------------------------------------
print*,"[CAUTION] Something wrong. hit_nuclide was not choosed."
choose_hit_nuclide_Z_A(1)=mat_Z(1)
choose_hit_nuclide_Z_A(2)=mat_A(1)
return
end

!=======================================================================
end module udm_int_pion_backward_photoproduction
!=======================================================================






