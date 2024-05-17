! Interaction Template Version = 1.0

!=======================================================================
module udm_int_new_boson_radiation
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be public.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "new_boson_radiation" ! This 'Name' is used in input files.
double precision, allocatable, save :: Parameters(:) ! Parameters entered in input files.
integer, parameter :: num_initial = 4 ! The number of incident particles causing this interaction.
integer, save :: kf_initial(num_initial) = (/ 11, -11, 13, -13 /) ! The kf-codes (particle IDs) of the incident particles causing this interaction.
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
integer kf_phi
integer boson_type, nZ
integer, allocatable, save :: Z_array(:)
double precision m_phi, g_e, g_mu, Kmax
double precision dXSmax_this, dvolume
double precision PiVEC(4),PfVEC(4),PifVEC(4),pVEC(4),pprimeVEC(4),kVEC(4),vVEC(4),qVEC(4)

integer, parameter :: nK=100
double precision K_array(nK)
double precision, allocatable, save :: XS    (:,:,:) ! kf_in, Z, K
double precision, allocatable, save :: dXSmax(:,:,:) ! kf_in, Z, K

double precision, parameter :: alpha=0.00729735d0
double precision, parameter :: pi=3.14159265d0
double precision, parameter :: dalton=931.4941d0 ! MeV

contains


!=======================================================================
subroutine print_comments
!=======================================================================
print*,"****************************************"
print*,"[Reference] for ", trim(Name)
print*,"Yu-Sheng Liu, David McKeen, Gerald A. Miller"
print*,"Phys.Rev.D 95 (2017) 3, 036010."
print*,"https://inspirehep.net/literature/1487721"
print*,"Yu-Sheng Liu, Gerald A. Miller"
print*,"Phys.Rev.D 96 (2017) 1, 016004."
print*,"https://inspirehep.net/literature/1598132"
print*,"****************************************"
end subroutine print_comments



!=======================================================================
subroutine initialize
! This subroutine is called only once at the beginning of the calculation.
!=======================================================================
double precision m ! mass of injection particle
double precision E ! total energy of injection particle
integer          Z ! atomic number of target particle
integer in,iZ,iK
integer inANTI,itmp
double precision dK,Kmin
logical output_XS
character(len=20) boson_type_c
character(len=199) OutputName
!-----------------------------------------------------------------------
! boson_type = 1(Scalar), 2(Pseudoscalar), 3(Vector), 4(Axial-vector)
!-----------------------------------------------------------------------

! ----------------------------------------------------
! read parameters
kf_phi    =Parameters(1) ! boson's kf-codes (particle ID).
boson_type=Parameters(2)
m_phi     =Parameters(3)
g_e       =Parameters(4)
g_mu      =Parameters(5)
Kmax      =Parameters(6)*1.001d0
nZ        =Parameters(7)
allocate(Z_array(nZ))
do iZ=1, nZ
  Z_array(iZ)=Parameters(7+iZ)
enddo
if(           Parameters(7+nZ+1) < -0.5d0) then
  output_XS=.true.
else
  output_XS=.false.
endif

if(.true.) then
print*,"-------------------------------"
print*,"Parameters: new_boson_radiation"
print*,"-------------------------------"
print*,"kf-code of boson =", kf_phi
if    (boson_type==1) then; boson_type_c="Scalar"
elseif(boson_type==2) then; boson_type_c="Pseudoscalar"
elseif(boson_type==3) then; boson_type_c="Vector"
elseif(boson_type==4) then; boson_type_c="Axial-vector"
else
  print*,"[error]: boson_type should be 1, 2, 3, or 4."
endif
print*,"boson type =", boson_type_c
print*,"mass [MeV] =", m_phi
print*,"g_e  =", g_e
print*,"g_mu =", g_mu
print*,"Kmax =", Kmax/1.001d0
write(*,"(a)",advance='no') " Target atoms ="
do iZ=1, nZ
  write(*,"(I4)",advance='no') Z_array(iZ)
enddo
print*
print*,"-------------------------------"
endif

! ----------------------------------------------------
! others
Kmin=max(m_phi-m,0d0)
if(Kmax <= Kmin) then
  print*,"Error: E(max) < (boson mass). Set larger E(max)."
  stop
endif
dK=(Kmax-Kmin)/(nK-1)
do iK=1,nK
  K_array(iK)=Kmin+(iK-1)*dK
enddo
! ----------------------------------------------------
! calculate cross section
allocate(XS    (num_initial,nZ,nK))
allocate(dXSmax(num_initial,nZ,nK))
XS    =-1d0
dXSmax=-1d0
do in=1,num_initial
  ! Find the index of the antiparticle for which the cross section has already been calculated.
  inANTI=0
  do itmp=1,in-1
    if( kf_initial(itmp)==-kf_initial(in) ) inANTI=itmp
  enddo
do iZ=1,nZ
do iK=1,nK
  if(inANTI > 0) then ! "XS for the anti-particle has already been calculated."
    XS    (in,iZ,iK)=XS    (inANTI,iZ,iK)
    dXSmax(in,iZ,iK)=dXSmax(inANTI,iZ,iK)
  else
    m=get_mass(kf_initial(in)) ! mass of injection particle
    E=K_array(iK)+m            ! total energy of injection particle
    Z=Z_array(iZ)              ! atomic number of target particle
    XS    (in,iZ,iK)=sigma(m,E,Z) ! When 'sigma(...)' is calculated, 'dXSmax_this' is assigned.
    dXSmax(in,iZ,iK)=dXSmax_this
  endif
enddo
enddo
enddo

if(output_XS) then
  print*,"Output:"
  do in=1,num_initial
    do iZ=1,nZ
    OutputName="XS_"//trim(boson_type_c)&
      &//"_kf"//trim(udm_int2char(kf_initial(in)))&
      &//"_Z"//trim(udm_int2char(Z_array(iZ)))&
      &//"_mass"//trim(udm_double2char("(ES9.3)",m_phi))&
      &//"_ge"//trim(udm_double2char("(ES9.3)",g_e))&
      &//"_gmu"//trim(udm_double2char("(ES9.3)",g_mu))&
      &//".dat"
    print*,"--> ",trim(OutputName)
    open(1,file=trim(OutputName),status='replace')
      write(1,"(A)") "# (kinetic energy of incident)  (cross section [barn])"
      do iK=1,nK
        write(1,"(ES11.3,ES11.3)") K_array(iK), XS(in,iZ,iK)
      enddo
    close(1)
    enddo
  enddo
  stop
endif

end subroutine initialize







!=======================================================================
double precision function dot(VEC1,VEC2)
!!! mostly plus metric !!! (-,+,+,+)
double precision VEC1(4), VEC2(4)
dot = -VEC1(1)*VEC2(1) + VEC1(2)*VEC2(2) + VEC1(3)*VEC2(3) + VEC1(4)*VEC2(4)
return
end
!=======================================================================
double precision function FormFactor(t,Z)
double precision t, aFac, dFac
integer Z
double precision A
character(len=99) tmp
call udm_get_atomic_information(Z,A,tmp,tmp)
aFac=112d0*Z**(-1d0/3d0)/get_mass(11)
dFac=0.164d0 * 1d+6 * A**(-2d0/3d0)
FormFactor=Z*(aFac**2*t)/(1d0+aFac**2*t)/(1d0+t/dFac)
return
end
!=======================================================================
double precision function dsigma(m,E,Z)
double precision m ! mass of injection particle
double precision E ! total energy of injection particle
integer          Z ! atomic number of target particle
integer i
double precision Aatom
double precision xmin, xmax, tmin, tmax, x, theta, Mnuc, Qminus, Qplus
double precision Vabs, Epr_Ef, Ek, Pk, rootin, q0, Qabs, c0, c1, c2
double precision xmod, thetamod, phiq, tmod
double precision stil, util, t2, t
double precision A23, P2, A1, A2, A2m, A3, A3m, A4, A5, A5m, Jacobian, epsilon
double precision costhq, costhq1, costhq2, sinthq
double precision :: MeVInv2toBarn = 388.09d0 ! [barn * MeV^2]. hbarc^2
logical ok_costhq1, ok_costhq2
character(len=99) tmp

! -------------------------------------------------------
! initialization
dsigma=0d0
xmin=m_phi/E
xmax=1d0-m/E
call udm_get_atomic_information(Z,Aatom,tmp,tmp)
Mnuc=dalton*Aatom
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PiVEC(1)=Mnuc;  pVEC(1)=E;
PiVEC(2)=0d0 ;  pVEC(2)=0d0;
PiVEC(3)=0d0 ;  pVEC(3)=0d0;
PiVEC(4)=0d0 ;  pVEC(4)=dsqrt(dabs(E**2-m**2));
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -------------------------------------------------------
! sampling
xmod    =get_random(log(xmin/(1d0-xmin)), log(xmax/(1d0-xmax)))
thetamod=get_random(0d0, dsqrt(pi))
phiq    =get_random(0d0, 2d0*pi)
! ... then get tmod
x=1d0/(1d0+exp(-xmod))
theta=thetamod**2
Ek=x*E
Pk=dsqrt(dabs(Ek**2-m_phi**2))
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kVEC(1)=Ek           ;  vVEC(1)=kVEC(1)-pVEC(1);
kVEC(2)=Pk*sin(theta);  vVEC(2)=kVEC(2)-pVEC(2);
kVEC(3)=0d0          ;  vVEC(3)=kVEC(3)-pVEC(3);
kVEC(4)=Pk*cos(theta);  vVEC(4)=kVEC(4)-pVEC(4);
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Vabs=dsqrt(vVEC(2)**2 + vVEC(3)**2 + vVEC(4)**2)
Epr_Ef=Mnuc+E-Ek
util=2d0*dot(pVEC,kVEC) + m_phi**2
rootin = util**2 + 4d0*Mnuc*util*Epr_Ef + 4d0*Mnuc**2*Vabs**2
! if(rootin > 0d0) print*,">___",rootin
! if(rootin < 0d0) print*,"<___",rootin
if(rootin < 0d0) then
  dvolume=0d0
  return ! dsigma = 0d0
endif

Qminus=(Vabs*(util+2d0*Mnuc*Epr_Ef)-Epr_Ef*dsqrt(rootin)) / (2d0*Epr_Ef**2-2d0*Vabs**2)
 Qplus=(Vabs*(util+2d0*Mnuc*Epr_Ef)+Epr_Ef*dsqrt(rootin)) / (2d0*Epr_Ef**2-2d0*Vabs**2)
if(Qminus<0.2d0*Mnuc) then ! To avoid numerical error
  tmin=Qminus**2
else
  tmin=2d0*Mnuc*(dsqrt(Mnuc**2+Qminus**2)-Mnuc)
endif
if(Qplus<0.2d0*Mnuc) then ! To avoid numerical error
  tmax=Qplus**2
else
  tmax=2d0*Mnuc*(dsqrt(Mnuc**2+Qplus**2)-Mnuc)
endif
tmod=get_random(log(tmin), log(tmax))
t=exp(tmod)
! -------------------------------------------------------
dvolume=(log(xmax/(1d0-xmax))-log(xmin/(1d0-xmin))) * dsqrt(pi) * 2d0*pi * (log(tmax)-log(tmin))
! -------------------------------------------------------

! -------------------------------------------------------
! check if the sampled variables are physical.
q0=-t/(2d0*Mnuc)
Qabs=dsqrt(q0**2+t)
c0=(m**2+t+dot(vVEC,vVEC)+2d0*q0*vVEC(1))/(2d0*Qabs)
c1=cos(phiq)*vVEC(2)/c0
c2=          vVEC(4)/c0
costhq1=(c2+dabs(c1)*dsqrt(c1**2+c2**2-1d0))/(c1**2+c2**2)
costhq2=(c2-dabs(c1)*dsqrt(c1**2+c2**2-1d0))/(c1**2+c2**2)
ok_costhq1 = -1d0<costhq1 .and. costhq1<1d0
ok_costhq2 = -1d0<costhq2 .and. costhq2<1d0
if    ( ok_costhq1 .and. (.not. ok_costhq2)) then
  ! print*,1,costhq1, costhq2
  costhq=costhq1
elseif( ok_costhq2 .and. (.not. ok_costhq1)) then
  ! print*,2,costhq1, costhq2
  costhq=costhq2
elseif( ok_costhq1 .and.        ok_costhq2 ) then
  ! print*,3,costhq1, costhq2
  if(get_random_0to1() < 0.5d0) then
    costhq=costhq1
  else
    costhq=costhq2
  endif
else
  return
  ! dsigma is 0d0
  ! dvolume has a non-zero value
endif
sinthq=dsqrt(1d0-costhq**2)
! -------------------------------------------------------

! -------------------------------------------------------
! reconstruct 4-momenta from the sampled variables
qVEC(1) = q0
qVEC(2) = Qabs * sinthq * cos(phiq)
qVEC(3) = Qabs * sinthq * sin(phiq)
qVEC(4) = Qabs * costhq

pprimeVEC(1) = pVEC(1) + qVEC(1) - kVEC(1)
pprimeVEC(2) = pVEC(2) + qVEC(2) - kVEC(2)
pprimeVEC(3) = pVEC(3) + qVEC(3) - kVEC(3)
pprimeVEC(4) = pVEC(4) + qVEC(4) - kVEC(4)

PifVEC(1) = 2d0*PiVEC(1) - qVEC(1)
PifVEC(2) = 2d0*PiVEC(2) - qVEC(2)
PifVEC(3) = 2d0*PiVEC(3) - qVEC(3)
PifVEC(4) = 2d0*PiVEC(4) - qVEC(4)
! -------------------------------------------------------

! -------------------------------------------------------
! amplitude
P2  = dot(PifVEC,PifVEC)
stil=-2d0*dot(pprimeVEC,kVEC) + m_phi**2
t2  = 2d0*dot(pprimeVEC,pVEC) + 2d0*m_phi**2

A1 = -(stil+util)**2/(stil*util)*P2 -4d0*t/(stil*util)*dot(PifVEC,kVEC)**2
A2 = (stil+util)**2/(stil*util)**2
A2m= (stil-util)**2/(stil*util)**2
A3 = P2*t +4d0*( (util*dot(PifVEC,pVEC)+stil*dot(PifVEC,pprimeVEC))/(stil+util) )**2
A3m= P2*t +4d0*( (util*dot(PifVEC,pVEC)+stil*dot(PifVEC,pprimeVEC))/(stil-util) )**2
A4 = -2d0*(stil**2+util**2)/(stil*util)*P2
A5 = 8d0*t/(stil*util) * (dot(PifVEC,pVEC)**2 + dot(PifVEC,pprimeVEC)**2 - (t2+m_phi**2)/2d0*P2)
A5m= 8d0*t/(stil*util) * (dot(PifVEC,pVEC)**2 + dot(PifVEC,pprimeVEC)**2 - (t2-m_phi**2)/2d0*P2)

if    (boson_type==1) then
  A23 = A1-A2*(m_phi**2-4d0*m**2)*A3
elseif(boson_type==2) then
  A23 = A1-A2*m_phi**2*A3
elseif(boson_type==3) then
  A23 = A4-A5-2d0*A2*(m_phi**2+2d0*m**2)*A3
elseif(boson_type==4) then
  A23 = A4-A5m+4d0*m**2/m_phi**2*A1-2d0*A2m*(m_phi**2-4d0*m**2)*A3m
else
  print*,"error: udm_int_new_boson_radiation"
  stop
endif
! -------------------------------------------------------

! -------------------------------------------------------
! Finally, calculate dsigma
if    ( dabs(1d0-m/get_mass(11)) < 0.01 ) then
  epsilon = g_e  / dsqrt(4d0*pi*alpha) 
elseif( dabs(1d0-m/get_mass(13)) < 0.01 ) then
  epsilon = g_mu / dsqrt(4d0*pi*alpha) 
else
  print*,"Error: udm_int_new_boson_radiation, mass = ",m
  stop
endif

Jacobian = x*(1d0-x) * 2d0*dsqrt(theta) * t

dsigma = Jacobian * epsilon**2 * alpha**3 * Pk * E / dsqrt(dabs(E**2-m**2)) / Vabs
dsigma = dsigma * FormFactor(t,Z)**2/t**2 / (8d0*Mnuc**2) / (2d0*pi) *A23
dsigma = dsigma * MeVInv2toBarn
! -------------------------------------------------------

return
end


!=======================================================================
double precision function sigma(m,E,Z)
double precision m ! mass of injection particle
double precision E ! total energy of injection particle
integer          Z ! atomic number of target particle
integer i
integer :: Nsample = 100000
double precision integrant, volume, dsigma_this

integrant=0d0
volume   =0d0
dXSmax_this=0d0

! Monte-Carlo integration
do i=1, Nsample
  dsigma_this = dsigma(m,E,Z)
  if(dsigma_this > dXSmax_this) dXSmax_this=dsigma_this
  integrant = integrant + dsigma_this
  volume    = volume    + dvolume
enddo

! averaged
integrant=integrant/Nsample
volume   =volume   /Nsample

sigma=volume*integrant

return
end


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
integer in ! index of injection particle
integer iZ ! index of Z array

do in=1, num_initial
  if(kf_initial(in)==udm_kf_incident) exit
enddo

do iZ=1, nZ
  if(Z_array(iZ)==Z) exit
enddo

Xsec_per_atom=udm_interp(Kin, nK, K_array, XS(in,iZ,:), "lin-lin")

return
end




!=======================================================================
subroutine generate_final_state
!=======================================================================
! Subroutine to determine final state information (4-momenta, etc.).
!-----------------------------------------------------------------------
! [Variables available in this subroutine]
! udm_kf_incident : The kf-codes (particle IDs) of the incident particles.
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------
double precision m ! mass of injection particle
double precision E ! total energy of injection particle
integer          Z ! atomic number of target particle
integer Z_A_hit(2), A
integer i, n_sampling_max, in, iZ
double precision tmp, dXSmax_interp, azimuth
!-----------------------------------------------------------------------
! You may use Z and A for final state variables.
! Z and A of a target material are automatically obtained using the function 'get_hit_nuclide_Z_A'.
Z_A_hit=get_hit_nuclide_Z_A(udm_Kin)
Z=Z_A_hit(1)
A=Z_A_hit(2)
!-----------------------------------------------------------------------
m=get_mass(udm_kf_incident)
E=udm_Kin+m
!-----------------------------------------------------------------------

! ----------------------------------------------------------------
! Start sampling
n_sampling_max = 100000
do i = 1, n_sampling_max
  ! --------------------------------------------------------------
  ! Rejection sampling method.
  ! When dsigma() is called, final state's 4-momenta are assigned automatically.
  tmp=dsigma(m, E, Z)
  if(tmp <= 0d0) cycle
  dXSmax_interp=udm_interp(udm_Kin, nK, K_array, dXSmax(in,iZ,:), "lin-lin")
  if(tmp > (10d0)*dXSmax_interp*get_random_0to1()) then ! (10d0) is a safety factor to avoid dsigma() > dXSmax(...).
    ! ----------------------------------------------------------------
    ! Accepted!!!
    ! Initialization is required before filling the final state information.
    call initialize_udm_event_info
    ! ----------------------------------------------------------------

    ! ------------------------------------
    ! Set number of final states to be recorded in event history.
    set_final_state_number = 2
    ! ------------------------------------

    azimuth=2d0*pi*get_random_0to1()

    ! ------------------------------------
    ! Set 4-momentum of new boson
    set_kf                 (1) = kf_phi
    set_Total_Energy_in_MeV(1) = kVEC(1)
    set_Px_in_MeV          (1) = kVEC(2)*cos(azimuth) - kVEC(3)*sin(azimuth)
    set_Py_in_MeV          (1) = kVEC(2)*sin(azimuth) + kVEC(3)*cos(azimuth)
    set_Pz_in_MeV          (1) = kVEC(4)

    ! ------------------------------------
    ! Set 4-momentum of incident fermion
    set_kf                 (2) = udm_kf_incident
    set_Total_Energy_in_MeV(2) = pprimeVEC(1)
    set_Px_in_MeV          (2) = pprimeVEC(2)*cos(azimuth) - pprimeVEC(3)*sin(azimuth)
    set_Py_in_MeV          (2) = pprimeVEC(2)*sin(azimuth) + pprimeVEC(3)*cos(azimuth)
    set_Pz_in_MeV          (2) = pprimeVEC(4)

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
end module udm_int_new_boson_radiation
!=======================================================================






