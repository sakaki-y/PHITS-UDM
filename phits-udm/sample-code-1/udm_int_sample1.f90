! Interaction Template Version = 1.0

!=======================================================================
module udm_int_sample_1
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be public.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "my_interaction_1" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
integer, parameter :: num_initial = 2 ! The number of incident particles causing this interaction.
integer, save :: kf_initial(num_initial) = (/ 11, -11 /) ! The kf-codes (particle IDs) of the incident particles causing this interaction.
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
! integer i,j
! double precision x,y
integer kf_X
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
! This subroutine is called only once at the beginning of the calculation.
!=======================================================================
kf_X=Parameters(1) ! kf-codes (particle ID) of X (outgoing particle).
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

if(Kin < 100.0) then
  Xsec_per_atom=0.0
  return
endif

if     (udm_kf_incident ==  11) then
  Xsec_per_atom=1e-6*Z
else if(udm_kf_incident == -11) then
  Xsec_per_atom=2e-6*Z
else
  print*,"error"
  stop
endif

return
end




!=======================================================================
!=======================================================================
! An example of sampling final state variables
!=======================================================================
!=======================================================================
double precision function distribution(Kin,Z,A, Kout)
! ----------------------------------------------------------------------
! Kout distribution for each inital state variables (Kin, Z, A).
! Normalization is not important, and this value should be less than 1 for any (Kin, Z, A).
implicit none
! ----------------------------------------------------------------------
! [Initial state variables]
double precision Kin ! Kinetic energy of incident particle [MeV]
           integer Z ! Atomic number of target atom
           integer A ! Mass   number of target atom
! ----------------------------------------------------------------------
! [Final state variables]
double precision Kout ! Kinetic energy of outgoing particle X.
!-----------------------------------------------------------------------
! [Variables available in this function]
! udm_kf_incident : The kf-codes (particle IDs) of the incident particles.
!=======================================================================
if(udm_kf_incident == 11) then
  distribution=Kout/Kin ! Kout distribution when the incident particle is 'electron'.
else
  distribution=1.0      ! Kout distribution when the incident particle is 'positron'.
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
integer Z_A_hit(2), Z, A
integer i, n_sampling_max
double precision Kout, Kout_min, Kout_max
double precision P, R
double precision m_X, E_X, p_X, theta_X
double precision m_e, E_e, p_e, theta_e
!-----------------------------------------------------------------------
! You may use Z and A for final state variables.
! Z and A of a target material are automatically obtained using the function 'get_hit_nuclide_Z_A'.
Z_A_hit=get_hit_nuclide_Z_A(udm_Kin)
Z=Z_A_hit(1)
A=Z_A_hit(2)
!-----------------------------------------------------------------------


! ================================== 
! = Your original sampling section =
! ==================================
! Range of sampling variable
Kout_min=0.0d0
Kout_max=udm_Kin-get_mass(kf_X)
! ----------------------------------------------------------------
! Start sampling
n_sampling_max = 100000
do i = 1, n_sampling_max
  ! --------------------------------------------------------------
  ! Randomly generate 'Kout' in the range.
  Kout=get_random(Kout_min, Kout_max)
  ! --------------------------------------------------------------
  ! Rejection sampling method.
  P=distribution(udm_Kin, Z, A, Kout) ! The likelihood of this 'Kout' happening.
  R=get_random_0to1()                 ! Random number from 0 to 1.
  ! If P > R, this Kout is accepted, otherwise rejected.
  if(P > R) then
    ! ----------------------------------------------------------------
    ! Accepted!!!
    ! Initialization is required before filling the final state information.
    call initialize_udm_event_info
    ! ----------------------------------------------------------------

    ! ------------------------------------
    ! Set number of final states to be recorded in event history.
    set_final_state_number = 2
    ! ------------------------------------

    ! ------------------------------------
    ! [Important Notice]
    ! Here, you set the final state momenta to the following "set_*" array, 
    ! assuming the direction of the momentum of the incident particle to be 
    ! positive along the "Z-axis". The final state momenta are automatically 
    ! rotated based on the actual direction outside of this subroutine.
    ! ------------------------------------

    m_X=get_mass(kf_X)
    E_X=Kout+m_X
    p_X=sqrt(E_X**2-m_X**2)
    theta_X=get_random(0.0d0,0.1d0)
    ! ------------------------------------
    ! Set 4-momentum of X
    set_kf                 (1) = kf_X
    set_Total_Energy_in_MeV(1) = E_X
    set_Px_in_MeV          (1) = p_X*sin(theta_X)
    set_Py_in_MeV          (1) = 0.0d0
    set_Pz_in_MeV          (1) = p_X*cos(theta_X)

    m_e=get_mass(udm_kf_incident)
    E_e=udm_Kin+m_e-E_X
    p_e=sqrt(E_e**2-m_e**2)
    theta_e=get_random(0.0d0,0.1d0)
    ! ------------------------------------
    ! Set 4-momentum of e+-
    set_kf                 (2) = udm_kf_incident
    set_Total_Energy_in_MeV(2) = E_e
    set_Px_in_MeV          (2) = p_e*sin(theta_e)
    set_Py_in_MeV          (2) = 0.0d0
    set_Pz_in_MeV          (2) = p_e*cos(theta_e)
    ! ------------------------------------

    ! [Additional Quantities]
    ! set_isomer_level            (?) = ?? ! (Default=0)   0: Not nuclear isomer, 1: 1st nuclear isomer, 2: 2nd nuclear isomer 
    ! set_excitation_energy_in_MeV(?) = ?? ! (Default=0.0) Excitation energy [MeV] 

    ! ------------------------------------
    ! The final states are recorded.
    call fill_final_state
    ! ------------------------------------
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
end module udm_int_sample_1
!=======================================================================






