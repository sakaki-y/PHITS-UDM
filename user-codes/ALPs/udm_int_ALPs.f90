!=======================================================================
module udm_int_ALPs
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none

!-----------------------------------------------------------------------
integer, parameter :: num_initial = 1
integer, save :: kf_initial(num_initial) = (/22/)
character(len=99), parameter :: Name = "ALPs"
!-----------------------------------------------------------------------
double precision, allocatable,save:: Parameters(:)
!-----------------------------------------------------------------------
integer, save :: kfX
double precision, save :: mX_MeV, gX_MeV
!-----------------------------------------------------------------------

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
!=======================================================================
kfX=int(Parameters(1))
mX_MeV=get_mass(kfX)
gX_MeV=Parameters(2)
end subroutine initialize


!=======================================================================
function Xsec_per_atom(Ein,Z_atom,A_atom)
!     Return value: Total cross section per an atom.
!     Unit: barn
!     Arguments:
!       Ein: Kinetic energy of initA_atoml particle
!       Z_atom : Atomic number of the target atom
!       A_atom : Mass   number of the target atom
!=======================================================================
implicit none
double precision Xsec_per_atom, Ein
integer Z_atom, A_atom
double precision eGeV,m,g,a,d,tmin,tmax,x,ymin,ymax,zmin,zmax
double precision, parameter :: pi = 3.14159d0
double precision, parameter :: alpha = 0.007299d0
double precision, parameter :: hbarc = 0.0197d0 ! GeV*sqrt(barn)
double precision, parameter :: me = 0.511d-3

eGeV=Ein/1000d0
m=mX_MeV/1000d0 ! mass [GeV]
g=gX_MeV*1000d0 ! coupling [1/GeV]
a=111d0*Z_atom**(-1./3.)/me
d=0.164d0*A_atom**(-2./3.)
tmin=m**4/(4.0*eGeV**2)
tmax=tmin+eGeV**2*pi**2
x=a**2*d
ymin=1+a**2*tmin
ymax=1+a**2*tmax
zmin=tmin+d
zmax=tmax+d

Xsec_per_atom= &
& g**2 *alpha *Z_atom**2 *hbarc**2 /8. &
& *(x/(x-1.))**2 &
& *((1.+2.*ymin/(x-1.))*log(zmin/zmax*ymax/ymin) &
&  -(1.-zmin/zmax)*((x-1)/ymax+2.))

return
end

!=======================================================================
function distribution(Ein,Z_atom,A_atom,lnt)
!     NormalZ_atomation is not important
!=======================================================================
double precision distribution,Ein,lnt
integer Z_atom,A_atom
double precision eGeV,m,a,t,a2t,d,tmin
double precision, parameter :: me = 0.511d-3

eGeV=Ein/1000d0
m=mX_MeV/1000d0 ! mass [GeV]

a=111d0*Z_atom**(-1./3.)/me
t=exp(lnt)
a2t=a**2*t
d=0.164d0*A_atom**(-2./3.)
tmin=m**4/(4.0*eGeV**2)
distribution=a2t/(1.+a2t)/(1+t/d)*(1-tmin/t)

return
end

!=======================================================================
function distribution_max(Ein,Z_atom,A_atom)
!=======================================================================
double precision distribution_max, Ein
integer Z_atom,A_atom
distribution_max=1d+0
return
end

!=======================================================================
subroutine generate_final_state
!=======================================================================
double precision Ein, xmin, xmax, x, y, ymax
integer i, n_sampling_max, kf
integer Z_A_hit(2), Z_hit, A_hit
double precision Etot,pabs

double precision, parameter :: pi = 3.14159d0
double precision, parameter :: alpha = 0.007299d0
double precision, parameter :: hbarc = 0.0197d0 ! GeV*sqrt(barn)
double precision, parameter :: me = 0.511d-3

double precision m,g,a,d,tmin,tmax,eGeV,lntmin,lntmax,lnt,t,theta,rMN,pE,absp,rphi,px,py,pz

! udm_Kin       : mother's Kinetic Energy [MeV]
! udm_kf_incident : mother's kf-code

Ein = udm_Kin ! Initial Kinetic Energy [MeV]
Z_A_hit = choose_hit_nuclide_Z_A(Ein)
Z_hit = Z_A_hit(1)
A_hit = Z_A_hit(2)

! ----------------------------------------
! MC parameters range and Area
eGeV=Ein/1000d0
m=mX_MeV/1000d0 ! mass [GeV]
g=gX_MeV*1000d0 ! coupling [1/GeV]
a=111d0*Z_hit**(-1./3.)/me
d=0.164d0*A_hit**(-2./3.)
tmin=m**4/(4.*eGeV**2)
tmax=tmin+eGeV**2*pi**2
lntmin=log(tmin)
lntmax=log(tmax)

! ----------------------------------------------------------------
! sampling
n_sampling_max = 100000
do i = 1, n_sampling_max

  ! --------------------------------------------------------------
  ! generate kinetic variables
  lnt=get_random(lntmin,lntmax)
  ! --------------------------------------------------------------
  ! calculating differential cross section
  y   =distribution    (Ein,Z_hit,A_hit,lnt)
  ymax=distribution_max(Ein,Z_hit,A_hit)

  ! --------------------------------------------------------------
  ! Rejection sampling (as example)
  if(y/ymax .gt. get_random_0to1()) then

    call initialize_udm_event_info

    t=exp(lnt)
    theta=dsqrt((t-tmin)/eGeV**2)
    rMN=0.94*A_hit ! GeV
    pE=eGeV-eGeV**2*theta**2/2./rMN-m**4/8./rMN/eGeV**2
    absp=dsqrt(pE**2-m**2)
    rphi = get_random(0d0,2.*pi)
    px=absp*sin(theta)*cos(rphi)
    py=absp*sin(theta)*sin(rphi)
    pz=absp*cos(theta)

    set_final_state_number = 1
    set_kf                 (1) = kfX
    set_Total_Energy_in_MeV(1) = pE*1000d0
    set_Px_in_MeV          (1) = px*1000d0
    set_Py_in_MeV          (1) = py*1000d0
    set_Pz_in_MeV          (1) = pz*1000d0

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
end module udm_int_ALPs
!=======================================================================






