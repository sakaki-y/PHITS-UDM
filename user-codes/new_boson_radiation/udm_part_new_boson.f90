! Particle Template Version = 1.0

!=======================================================================
module udm_part_new_boson
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be private.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "new_boson" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
integer kf_this ! The kf-code (particle ID) of this particle entered in input files.
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
integer boson_type
integer :: N_caution = 0
double precision m_phi, g_e, g_mu

double precision, parameter :: alpha=0.00729735d0
double precision, parameter :: pi=3.14159265d0
contains



!=======================================================================
subroutine print_comments
!=======================================================================
! print*,"****************************************"
! print*,"[Reference] for ", trim(Name)
! print*,"Yu-Sheng Liu, David McKeen, Gerald A. Miller"
! print*,"Phys.Rev.D 95 (2017) 3, 036010."
! print*,"https://inspirehep.net/literature/1487721"
! print*,"Yu-Sheng Liu, Gerald A. Miller"
! print*,"Phys.Rev.D 96 (2017) 1, 016004."
! print*,"https://inspirehep.net/literature/1598132"
! print*,"****************************************"
end subroutine print_comments


!=======================================================================
subroutine initialize
! This subroutine is called only once at the beginning of a calculation.
!=======================================================================
! read parameters
boson_type=Parameters(1)
m_phi     =Parameters(2)
g_e       =Parameters(3)
g_mu      =Parameters(4)
end subroutine initialize



!=======================================================================
double precision function mass() ! Unit: MeV
!=======================================================================
mass=m_phi
return
end

!=======================================================================
function asin_comp(x)
complex(kind(0d0)) asin_comp
double precision x
complex(kind(0d0)) :: i = (0d0,1d0)
asin_comp=-i*cdlog(i*x+cdsqrt((1d0,0d0)*(1d0-x**2)))
return
end

!=======================================================================
double precision function fS(t)
double precision t
fS = 1d0/(4d0*t**2) * cdabs( 1d0 + (1d0-1d0/t)*asin_comp(dsqrt(t))**2 )**2
return
end
!=======================================================================
double precision function fP(t)
double precision t
complex(kind(0d0)) braket
braket = t + cdsqrt( (t**2-t)*(1d0,0d0) )
fP = 1d0/(64d0*t**2) * cdabs( cdlog( (1d0-2d0*braket) )**2 )**2
return
end

!=======================================================================
double precision function Gamma1(m)
! phi -> n photon (n>=2)
double precision m, epsilon

! assign 'epsilon'
if    ( dabs(1d0-m/get_mass(11)) < 0.01 ) then
  epsilon = g_e  / dsqrt(4d0*pi*alpha) 
elseif( dabs(1d0-m/get_mass(13)) < 0.01 ) then
  epsilon = g_mu / dsqrt(4d0*pi*alpha) 
else
  print*,"Error: udm_int_new_boson_radiation, mass = ",m
  stop
endif

! assign 'Gamma1'
if    (boson_type==1) then
  Gamma1 = epsilon**2 * alpha**3/(4d0*pi**2) * m_phi**3/m**2 * fS(m_phi**2/4d0/m**2)
elseif(boson_type==2) then
  Gamma1 = epsilon**2 * alpha**3/(4d0*pi**2) * m_phi**3/m**2 * fP(m_phi**2/4d0/m**2)
elseif(boson_type==3) then
  Gamma1 = epsilon**2 * alpha**4/(2**7*3**6*5**2*pi**3) * m_phi**9/m**8 &
         & * (17d0/5d0 + 67d0/42d0*m_phi**2/m**2 + 128941d0/246960d0*m_phi**4/m**4)
elseif(boson_type==4) then
  Gamma1 = epsilon**2 * 127d0*alpha**5/(2**11*3**8*5**4*7**2*pi**4) * m_phi**13/m**12
else
  print*,"error: udm_part_new_boson, Gamma1"
  stop
endif

return
end

!=======================================================================
double precision function Gamma2(m)
! phi -> l l
double precision m, epsilon

! assign 'epsilon'
if    ( dabs(1d0-m/get_mass(11)) < 0.01 ) then
  epsilon = g_e  / dsqrt(4d0*pi*alpha) 
elseif( dabs(1d0-m/get_mass(13)) < 0.01 ) then
  epsilon = g_mu / dsqrt(4d0*pi*alpha) 
else
  print*,"Error: udm_int_new_boson_radiation, mass = ",m
  stop
endif

! assign 'Gamma1'
if    (boson_type==1) then
  Gamma2 = epsilon**2 * alpha/2d0 * m_phi * (1d0-4d0*m**2/m_phi**2)**(3d0/2d0)
elseif(boson_type==2) then
  Gamma2 = epsilon**2 * alpha/2d0 * m_phi * (1d0-4d0*m**2/m_phi**2)**(1d0/2d0)
elseif(boson_type==3) then
  Gamma2 = epsilon**2 * alpha/3d0 * m_phi * (1d0-4d0*m**2/m_phi**2)**(1d0/2d0) * (1d0+2d0*m**2/m_phi**2)
elseif(boson_type==4) then
  Gamma2 = epsilon**2 * alpha/3d0 * m_phi * (1d0-4d0*m**2/m_phi**2)**(3d0/2d0)
else
  print*,"error: udm_part_new_boson, Gamma2"
  stop
endif

return
end

!=======================================================================
double precision function lifetime() ! Unit: Second
! The value should be greater than 0.
!=======================================================================
double precision Gamma
Gamma = Gamma1(get_mass(11)) + Gamma1(get_mass(13))
if(m_phi > 2d0*get_mass(11)) Gamma = Gamma + Gamma2(get_mass(11))
if(m_phi > 2d0*get_mass(13)) Gamma = Gamma + Gamma2(get_mass(13))

! MeVInv2Sec = hbarc/c = 6.58209d-22 MeV * second
lifetime= 6.58209d-22 / Gamma
return
end

!=======================================================================
subroutine decay
!=======================================================================
! [Variables available in this subroutine]
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------
double precision tot, Br_photon, Br_e, Br_mu, P

Br_photon = Gamma1(get_mass(11)) + Gamma1(get_mass(13))

Br_e=0d0; Br_mu=0d0;
if(m_phi > 2d0*get_mass(11)) Br_e  = Gamma2(get_mass(11))
if(m_phi > 2d0*get_mass(13)) Br_mu = Gamma2(get_mass(13))

tot = Br_photon + Br_e + Br_mu
Br_photon = Br_photon/tot
Br_e      = Br_e     /tot
Br_mu     = Br_mu    /tot

P=get_random_0to1()
if (Br_photon      > P) then
  if    (boson_type==1) then
    call two_body_decay_uniform(22,22)
  elseif(boson_type==2) then
    call two_body_decay_uniform(22,22)
  elseif(boson_type==3) then
    call three_body_decay_uniform(22,22,22)
  elseif(boson_type==4) then
    if(N_caution==0) then
      print*,"phi -> 4 gammas has not been implemented yet."
      print*,"Tentatively changed to 3-body decay."
      N_caution=1
    endif
    call three_body_decay_uniform(22,22,22)
  else
    print*,"error: udm_part_new_boson, decay"
  endif
else if(Br_photon+Br_e > P) then
  call two_body_decay_uniform(11,-11)
else
  call two_body_decay_uniform(13,-13)
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
end module udm_part_new_boson
!=======================================================================



