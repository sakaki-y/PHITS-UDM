!=======================================================================
module udm_int_K0_prod
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
!-----------------------------------------------------------------------
integer, parameter :: num_initial = 1
integer, save :: kf_initial(num_initial) = (/22/)
character(len=99), parameter :: Name = "K0_prod"
!-----------------------------------------------------------------------
double precision, allocatable,save:: Parameters(:)
!-----------------------------------------------------------------------
integer, parameter :: ndata=15
double precision, save :: WL(ndata)
double precision, save :: sigmaL(ndata)
double precision, save :: xL(11)
double precision, save :: yL(11)
!-----------------------------------------------------------------------
contains


!=======================================================================
subroutine print_comments
!=======================================================================
print*,"****************************************"
print*,"[Reference] for ", trim(Name)
print*,"A. M. Boyarski, et.al.,"
print*,"Phys. Rev. Lett. 22, 1131 – Published 26 May 1969"
print*,"https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.22.1131"
print*,"****************************************"
end subroutine print_comments


!=======================================================================
subroutine initialize
!=======================================================================
WL    =(/1.690,1.749,1.821,1.862,1.893,1.940,1.986,2.036,2.134,2.280,2.430,2.582,2.743,2.887,3.000/)
sigmaL=(/0.000,0.417,0.756,0.928,0.970,0.876,0.798,0.788,0.882,0.956,0.956,0.842,0.645,0.440,0.320/)
xL=(/0.000, 0.031, 0.085, 0.143, 0.203, 0.257, 0.327, 0.394, 0.467, 0.537, 0.600/)
yL=(/0.519, 0.639, 0.742, 0.798, 0.797, 0.761, 0.695, 0.606, 0.491, 0.410, 0.339/)
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

double precision eGeV,W
double precision, parameter :: mpGeV = 0.938272081

eGeV=(Ein+0d0)/1000d0
W=dsqrt(mpGeV**2+2.*mpGeV*eGeV)
if(W > WL(ndata)) then
  Xsec_per_atom=6.0d-6/eGeV**2*Z_atom
  return
endif

Xsec_per_atom=udm_interp(W,ndata,WL,sigmaL,"lin-lin")*1.0d-6*Z_atom
return
end

!=======================================================================
function distribution(T)
!     Normalzation is not important
!     dP/dT
!=======================================================================
double precision distribution,T

if( T >= 0.600 ) then
  distribution=0.339*exp(-3.*T)/exp(-3.*0.600)
  return
endif

distribution=udm_interp(T,11,xL,yL,"lin-lin")
return
end

!=======================================================================
function distribution_max()
!=======================================================================
double precision distribution_max
distribution_max=0.798
return
end

!=======================================================================
subroutine generate_final_state
! udm_Kin       : mother's Kinetic Energy [MeV]
! udm_kf_incident : mother's kf-code
!=======================================================================
double precision eGeV,sqrts,Eg_cm, EN_cm,Eh_cm, pN_cm,ph_cm, gamma,gamma_beta, Tmin,Tmax,T
double precision y,ymax
double precision costh_cm, Eh_lab,ph_lab,costh_lab,theta,phi
double precision pxGeV,pyGeV,pzGeV
integer Z_A_hit(2), Z_hit, A_hit
integer i, n_sampling_max
double precision, parameter :: mN = 0.938272081
double precision, parameter :: mh = 0.497611
double precision, parameter :: mX = 1.18937

! ----------------------------------------
! initial parameters
eGeV=(udm_Kin+0d0)/1000d0
sqrts=sqrt(mN**2+2.*eGeV*mN)
Eg_cm=sqrts/2.*(1.       -mN**2 /sqrts**2)
EN_cm=sqrts/2.*(1.       +mN**2 /sqrts**2)
Eh_cm=sqrts/2.*(1.+(mh**2-mX**2)/sqrts**2)
pN_cm=sqrt(dabs(EN_cm**2-mN**2))
ph_cm=sqrt(dabs(Eh_cm**2-mh**2))
gamma     =EN_cm/mN
gamma_beta=pN_cm/mN

! ----------------------------------------
! MC parameters range and Area
Tmin=-mh**2+2.*Eg_cm*(Eh_cm-ph_cm)
Tmax=-mh**2+2.*Eg_cm*(Eh_cm+ph_cm)

! ----------------------------------------------------------------
Z_A_hit = choose_hit_nuclide_Z_A(udm_Kin)
Z_hit = Z_A_hit(1)
A_hit = Z_A_hit(2)

! ----------------------------------------------------------------
! rejection sampling
n_sampling_max = 100000
do i = 1, n_sampling_max

  ! --------------------------------------------------------------
  ! generate kinetic variables
  T=get_random(Tmin,Tmax)

  ! --------------------------------------------------------------
  ! calculating differential cross section
  y   =distribution(T)
  ymax=distribution_max()

  ! --------------------------------------------------------------
  ! Rejection sampling (as example)
  if(y/ymax > get_random_0to1()) then

    call initialize_udm_event_info

    ! K0 4-momentum
    costh_cm=(Eh_cm-(T+mh**2)/(2.*Eg_cm))/ph_cm
    Eh_lab=gamma*Eh_cm+gamma_beta*ph_cm*costh_cm
    if(Eh_lab**2-mh**2 < 0d0) then
      print*,"Caution: gen_ALP. Minus in Sqrt. Skip."
      cycle
    endif
    ph_lab=dsqrt(dabs(Eh_lab**2-mh**2))
    costh_lab=(gamma_beta*Eh_cm+gamma*ph_cm*costh_cm)/ph_lab
    theta=acos(costh_lab)
    phi=get_random(0d0,2.*3.14159d0)

    pxGeV=ph_lab*sin(theta)*cos(phi)
    pyGeV=ph_lab*sin(theta)*sin(phi)
    pzGeV=ph_lab*costh_lab

    set_final_state_number = 1
    if(get_random_0to1() > 0.5) then
      set_kf(1) =  311
    else
      set_kf(1) = -311
    endif
    set_Total_Energy_in_MeV(1) = Eh_lab*1000d0
    set_Px_in_MeV          (1) = pxGeV *1000d0
    set_Py_in_MeV          (1) = pyGeV *1000d0
    set_Pz_in_MeV          (1) = pzGeV *1000d0

!   [Additional Quantities]
!   set_isomer_level            (?) = ?? ! (Default=0)   0: Not nuclear isomer, 1: 1st nuclear isomer, 2: 2nd nuclear isomer 
!   set_excitation_energy_in_MeV(?) = ?? ! (Default=0.0) Excitation energy [MeV] 

    call fill_final_state
    return
  endif
  ! --------------------------------------------------------------
enddo
! ----------------------------------------------------------------

print*, "[CAUTION] sampling failed. Increase 'n_sampling_max'.", Name

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
end module udm_int_K0_prod
!=======================================================================






