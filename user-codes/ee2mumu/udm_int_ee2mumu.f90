!=======================================================================
module udm_int_ee2mumu
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none

!-----------------------------------------------------------------------
integer, parameter :: num_initial = 1
integer, save :: kf_initial(num_initial) = (/-11/)
character(len=99), parameter :: Name = "ee2mumu"
!-----------------------------------------------------------------------
double precision, allocatable,save:: Parameters(:)
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
double precision me, mmu, Elab, plab, s, pcm, Ecm, r
double precision, parameter :: pi = 3.14159d0
double precision, parameter :: alpha = 0.007299d0
double precision, parameter :: hbarc = 0.0197d0 ! GeV*sqrt(barn)

if(Ein <= 43727d0) then
  Xsec_per_atom=0d0
  return
endif

me =get_mass(11)/1000d0 ! e,  GeV
mmu=get_mass(13)/1000d0 ! mu, GeV
Elab=Ein/1000d0+me ! GeV, total energy
plab=sqrt(abs(Elab**2-me**2)) ! GeV
s=me**2+me**2+2d0*Elab*me ! GeV
Ecm=sqrt(s)/2d0
r=mmu**2/Ecm**2

if(1-r <= 0d0) then
  Xsec_per_atom=0d0
  return
endif

Xsec_per_atom= &
& Z_atom*hbarc**2 &
& *4d0*pi*alpha**2/3d0/s &
& *sqrt(1d0-r) &
& *(1d0+r/2d0) 

! print*,"Xsec_per_atom__",Xsec_per_atom

return
end

!=======================================================================
function distribution(Ein,Z_atom,A_atom,costhcm)
!     Normaliation is not important
!=======================================================================
double precision distribution,Ein,costhcm
integer Z_atom,A_atom
double precision me, mmu, Elab, plab, s, pcm, Ecm, r

me =get_mass(11)/1000d0 ! e,  GeV
mmu=get_mass(13)/1000d0 ! mu, GeV
Elab=Ein/1000d0+me ! GeV, total energy
plab=sqrt(abs(Elab**2-me**2)) ! GeV
s=me**2+me**2+2d0*Elab*me ! GeV
Ecm=sqrt(s)/2d0
pcm=sqrt(abs(Ecm**2-me**2))
r=mmu**2/Ecm**2

distribution=(1d0+r)+(1d0-r)*costhcm**2

return
end

!=======================================================================
function distribution_max(Ein,Z_atom,A_atom)
!=======================================================================
double precision distribution_max, Ein
integer Z_atom,A_atom
distribution_max=2d+0
return
end

!=======================================================================
subroutine generate_final_state
!=======================================================================
double precision Ein, xmin, xmax, x, y, ymax
integer i, n_sampling_max, kf
integer Z_A_hit(2), Z_hit, A_hit

double precision me, mmu, Elab, plab, s, pcm, Ecm, r
double precision gamma, gammabeta
double precision costhcm_min,costhcm_max,costhcm,thetacm,phicm
double precision Emucm , pxmucm ,pymucm ,pzmucm , pmucm
double precision Emulab, pxmulab,pymulab,pzmulab
double precision, parameter :: pi = 3.14159d0

! udm_Kin         : incident particles's Kinetic Energy [MeV]
! udm_kf_incident : incident particles's kf-code

Ein = udm_Kin ! Initial Kinetic Energy [MeV]
Z_A_hit = choose_hit_nuclide_Z_A(Ein)
Z_hit = Z_A_hit(1)
A_hit = Z_A_hit(2)

me =get_mass(11)/1000d0 ! e,  GeV
mmu=get_mass(13)/1000d0 ! mu, GeV
Elab=Ein/1000d0+me ! GeV, total energy
plab=sqrt(abs(Elab**2-me**2)) ! GeV
s=me**2+me**2+2d0*Elab*me ! GeV
Ecm=sqrt(s)/2d0
pcm=sqrt(abs(Ecm**2-me**2))
gamma    =Ecm/me
gammabeta=pcm/me
r=mmu**2/Ecm**2

! ----------------------------------------
! MC parameters range and Area
costhcm_min=-1d0
costhcm_max= 1d0

! ----------------------------------------------------------------
! sampling
n_sampling_max = 100000
do i = 1, n_sampling_max

  ! --------------------------------------------------------------
  ! generate kinetic variables
  costhcm=get_random(costhcm_min,costhcm_max)
  ! --------------------------------------------------------------
  ! calculating differential cross section
  y   =distribution    (Ein,Z_hit,A_hit,costhcm)
  ymax=distribution_max(Ein,Z_hit,A_hit)

  ! --------------------------------------------------------------
  ! Rejection sampling (as example)
  if(y/ymax .gt. get_random_0to1()) then

    call initialize_udm_event_info

    thetacm=acos(costhcm)
    ! print*,"thetacm__",thetacm/pi
    phicm=get_random(0d0,2d0*pi)

    Emucm =Ecm
    pmucm =sqrt(abs(Emucm**2-mmu**2))
    pxmucm=pmucm*sin(thetacm)*cos(phicm)
    pymucm=pmucm*sin(thetacm)*sin(phicm)
    pzmucm=pmucm*cos(thetacm)

    ! ------------------------------------
    set_final_state_number = 2
    ! ------------------------------------
    Emulab =gamma*Emucm+gammabeta*pzmucm
    pxmulab=                      pxmucm
    pymulab=                      pymucm
    pzmulab=gammabeta*Emucm+gamma*pzmucm

    set_kf                 (1) = 13
    set_Total_Energy_in_MeV(1) = Emulab *1000d0
    set_Px_in_MeV          (1) = pxmulab*1000d0
    set_Py_in_MeV          (1) = pymulab*1000d0
    set_Pz_in_MeV          (1) = pzmulab*1000d0
    ! ------------------------------------
    pxmucm=-pxmucm
    pymucm=-pymucm
    pzmucm=-pzmucm

    Emulab =gamma*Emucm+gammabeta*pzmucm
    pxmulab=                      pxmucm
    pymulab=                      pymucm
    pzmulab=gammabeta*Emucm+gamma*pzmucm

    set_kf                 (2) = -13
    set_Total_Energy_in_MeV(2) = Emulab *1000d0
    set_Px_in_MeV          (2) = pxmulab*1000d0
    set_Py_in_MeV          (2) = pymulab*1000d0
    set_Pz_in_MeV          (2) = pzmulab*1000d0
    ! ------------------------------------

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
end module udm_int_ee2mumu
!=======================================================================






