import os, sys, glob


check_mode=False


# ======================================================================
filename=[]
before  =[]
after   =[]


# ======================================================================
def check(f,b,a):
    if not os.path.exists(f):
        sys.exit("not found: "+f)
    with open(f) as FILE:
        s=FILE.read()
    if s.count(b)!=1:
        print("@@@ N =",s.count(b))
        print(f)
        print(b)
        print(a)
        sys.exit("ERROR")


# # ======================================================================
# f="AAA.f"
# b="""\
# """
# a="""\
# """
# if check_mode: check(f,b,a)
# filename.append(f)
# before  .append(b)
# after   .append(a)


# ======================================================================
f="dataup.f"
b="""\
     &           kid .ne. 0 ) )
"""
a="""\
     &          (kid .ne. 0 .or. 
     &          (900000 .le. abs(kf) .and. abs(kf) .le. 999999 ) )
     &                      ) )
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="dklos.f"
b="""\
*     output :  in common from dklos2, dklos3 and c9decay              *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)
"""
a="""\
*     output :  in common from dklos2, dklos3 and c9decay              *
*                                                                      *
*                                                                      *
************************************************************************

      use udm_Parameter ! y.sakaki 2021/9
      use udm_Utility   ! y.sakaki 2021/9
      use udm_Manager   ! y.sakaki 2021/9

      implicit real*8(a-h,o-z)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="dklos.f"
b="""\
      common /othdcp/ nodcp(2000,6), nodcn
      common /dcayp/  adcayp(20), bdcayp(20)
"""
a="""\
      common /othdcp/ nodcp(2000,6), nodcn
      common /dcayp/  adcayp(20), bdcayp(20)
      common /udmodel/ iudmodel 
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="dklos.f"
b="""\
         kdcay = 0
"""
a="""\
         kdcay = 0

C y.sakaki 2021/10 -->
*-----------------------------------------------------------------------
*        user defined particle decay
*-----------------------------------------------------------------------
C        "user defined particle" function can modify any particle decay pattern.
C        So, this function is written at the beginning of dklos.

         if(iudmodel .gt. 0) then

           udm_kf_for_21=kprj
           udm_Kin=eein
           set_decay_success=.false.
           call user_defined_particle(21)
           if(set_decay_success) goto 1000
           if( 900000 .le. abs(kprj) .and. abs(kprj) .le. 999999) then
             goto 1000
           endif
         endif
C <-- y.sakaki 2021/10
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="dklos.f"
b="""\
            write(6,'(''*** no data in dklos.f kf = '',i7)') kprj
"""
a="""\
            if(udm_counter_dklos .eq. 0) ! y.sakaki 2021/10
     &      write(6,'(''*** no data in dklos.f kf = '',i7)') kprj
            if(iudmodel .gt. 0) udm_counter_dklos=1 ! y.sakaki 2021/10
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="dklos.f"
b="""\
*     sumary decay mode
"""
a="""\

 1000 continue

*-----------------------------------------------------------------------
*     sumary decay mode
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
      use neutrino_mod, only : neutrino_Xsec
      use ion_track_structure, only : pts_flt2
"""
a="""\
      use neutrino_mod, only : neutrino_Xsec
      use ion_track_structure, only : pts_flt2
      use udm_Parameter
      use udm_Utility
      use udm_Manager
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
cfrtati 2021/12/17
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)
"""
a="""\
cfrtati 2021/12/17
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)

      common /udmodel/ iudmodel ! y.sakaki 2021/09
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
*-----------------------------------------------------------------------
*        photon for nuclear data
*-----------------------------------------------------------------------
"""
a="""\
! y.sakaki 2021/09 -->
*-----------------------------------------------------------------------
*        user defined interaction for photon
*-----------------------------------------------------------------------
               if( iudmodel .ge. 1 ) then
                 totudm_other=totttl
                !-----------------------------
                 totudm    =0d0
                 totudm_mul=0d0
                !-----------------------------
                 if(udm_int_nMax .gt. 0) call fill_mat_info(mat)
                 do i = 1, udm_int_nMax
                   udm_Kin=ein
                   udm_kf_incident=ktyp
                   udm_sigt=0d0
                   call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                   totudm_this = udm_sigt * denm(mat)
                   totudm_(i) = totudm_this
                   totudm     = totudm_this               + totudm
                   totudm_mul = totudm_this * udm_bias(i) + totudm_mul
                 enddo
                 totttl = totttl + totudm
               endif
! <-- y.sakaki 2021/09

*-----------------------------------------------------------------------
*        photon for nuclear data
*-----------------------------------------------------------------------
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
                     fpl = tstep
"""
a="""\

! y.sakaki 2021/09 -->
*-----------------------------------------------------------------------
*        user defined interaction for e+ e-
*-----------------------------------------------------------------------
               if( iudmodel .ge. 1 ) then
                 totudm_other=totttl
                !-----------------------------
                 totudm    =0d0
                 totudm_mul=0d0
                !-----------------------------
                 if(udm_int_nMax .gt. 0) call fill_mat_info(mat)
                 do i = 1, udm_int_nMax
                   udm_Kin=ein
                   udm_kf_incident=ktyp
                   udm_sigt=0d0
                   call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                   totudm_this = udm_sigt * denm(mat)
                   totudm_(i) = totudm_this
                   totudm     = totudm_this               + totudm
                   totudm_mul = totudm_this * udm_bias(i) + totudm_mul
                 enddo
                 totttl = totttl + totudm
               endif
! <-- y.sakaki 2021/09

                     fpl = tstep
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
*-----------------------------------------------------------------------
*        for high energy part
*-----------------------------------------------------------------------
"""
a="""\

C y.sakaki 2021/9 -->
*-----------------------------------------------------------------------
*        Interaction for user defined particle
*-----------------------------------------------------------------------
C             else if( (900000 .le. abs(ktyp) .and. abs(ktyp) .le. 999999) 
C      &             .and. ein .le. dnmax(11) ) then
            else if( 900000 .le. abs(ktyp) .and. abs(ktyp) .le. 999999) 
     &        then
*-----------------------------------------------------------------------
*        user defined interaction for new particles
*-----------------------------------------------------------------------
               if( iudmodel .ge. 1 ) then
                 totudm_other=totttl
                !-----------------------------
                 totudm    =0d0
                 totudm_mul=0d0
                !-----------------------------
                 if(udm_int_nMax .gt. 0) call fill_mat_info(mat)
                 do i = 1, udm_int_nMax
                   udm_Kin=ein
                   udm_kf_incident=ktyp
                   udm_sigt=0d0
                   call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                   totudm_this = udm_sigt * denm(mat)
                   totudm_(i) = totudm_this
                   totudm     = totudm_this               + totudm
                   totudm_mul = totudm_this * udm_bias(i) + totudm_mul
                 enddo
                 totttl = totttl + totudm
               endif
C <-- y.sakaki 2021/9

*-----------------------------------------------------------------------
*        for high energy part
*-----------------------------------------------------------------------
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
                        totttl = totttl + sigvpi + sigbrm + sigppd

                     endif
"""
a="""\
                        totttl = totttl + sigvpi + sigbrm + sigppd

                     endif

! y.sakaki 2021/09 -->
*-----------------------------------------------------------------------
*        user defined interaction for muon- muon+
*-----------------------------------------------------------------------
               if( iudmodel .ge. 1 ) then
                 totudm_other=totttl
                !-----------------------------
                 totudm    =0d0
                 totudm_mul=0d0
                !-----------------------------
                 if(udm_int_nMax .gt. 0) call fill_mat_info(mat)
                 do i = 1, udm_int_nMax
                   udm_Kin=ein
                   udm_kf_incident=ktyp
                   udm_sigt=0d0
                   call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                   totudm_this = udm_sigt * denm(mat)
                   totudm_(i) = totudm_this
                   totudm     = totudm_this               + totudm
                   totudm_mul = totudm_this * udm_bias(i) + totudm_mul
                 enddo
                 totttl = totttl + totudm
               endif
! <-- y.sakaki 2021/09

"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
*       dcyfl : decay flight length (cm)                               *
*                                                                      *
************************************************************************
"""
a="""\
*       dcyfl : decay flight length (cm)                               *
*                                                                      *
************************************************************************

      use udm_Parameter
      use udm_Manager

"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
*-----------------------------------------------------------------------

               dcyfl = 0.0d0
               tlife = -1.0

*-----------------------------------------------------------------------
"""
a="""\
      common /udmodel/ iudmodel ! y.sakaki 2022/10

*-----------------------------------------------------------------------

               dcyfl = 0.0d0
               tlife = -1.0

*-----------------------------------------------------------------------

C y.sakaki 2022/10 -->
C        "user defined particle" function can modify any particle lifetime.
C        So, this function is written at the beginning of "dcyfl" function.
         if(iudmodel .gt. 0) then
           udm_kf_for_11=kf
           udm_Kin=ein
           udm_lifetime=0d0
           call user_defined_particle(11) ! get 'udm_lifetime' for 'kf'
           if(udm_lifetime > 0d0) then
             tlife = udm_lifetime
             goto 100
           elseif(900000 .le. abs(kf) .and. abs(kf) .le. 999999) then
             tlife = 1d+30
             goto 100
           endif ! otherwise go back to the normal mode
         endif
C <-- y.sakaki 2022/10

"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
         do k = 1, ndcy

            if( kf .eq. kfdcy(k) ) then

               idcyc = 1
               return

            end if

         end do
"""
a="""\
         do k = 1, ndcy

            if( kf .eq. kfdcy(k) ) then

               idcyc = 1
               return

            end if

         end do

C y.sakaki 2021/10 -->
C          if( 900000 .le. abs(kf) .and. abs(kf) .le. 999999  ) then
C              idcyc = 1
C          endif
C <-- y.sakaki 2021/10

"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
************************************************************************
*                                                                      *
      subroutine getfltmu
*                                                                      *
*       for muon, get total flight path length at ec(no)               *
*       last modified by S.Abe on 2016/07/21                           *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************

      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
"""
a="""\
************************************************************************
*                                                                      *
      subroutine getfltmu
*                                                                      *
*       for muon, get total flight path length at ec(no)               *
*       last modified by S.Abe on 2016/07/21                           *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************

      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use udm_Parameter
      use udm_Utility
      use udm_Manager
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
!<-20211126murofushi delete
!--      dimension dnel(1), denh(1), zz(1), a(1), den(1)
!--      equivalence ( das, dnel, denh, zz, a, den )
!--      dimension edns(1)
!--      equivalence ( das, edns )
!--->

*-----------------------------------------------------------------------
*        only for muon
*-----------------------------------------------------------------------
"""
a="""\
C y.sakaki 2021/09   For user defined model
      common /udmodel/ iudmodel 
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

!<-20211126murofushi delete
!--      dimension dnel(1), denh(1), zz(1), a(1), den(1)
!--      equivalence ( das, dnel, denh, zz, a, den )
!--      dimension edns(1)
!--      equivalence ( das, edns )
!--->

*-----------------------------------------------------------------------
*        only for muon
*-----------------------------------------------------------------------
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
         totttl = totttl + sigvpi + sigbrm + sigppd

      endif
"""
a="""\
         totttl = totttl + sigvpi + sigbrm + sigppd

      endif

! y.sakaki 2021/09 -->
*-----------------------------------------------------------------------
*        user defined interaction for muon- muon+
*-----------------------------------------------------------------------
               if( iudmodel .ge. 1 ) then
                 totudm_other=totttl
                !-----------------------------
                 totudm    =0d0
                 totudm_mul=0d0
                !-----------------------------
                 if(udm_int_nMax .gt. 0) call fill_mat_info(mat)
                 do i = 1, udm_int_nMax
                   udm_Kin=ein
                   udm_kf_incident=ktyp
                   udm_sigt=0d0
                   call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                   totudm_this = udm_sigt * denm(mat)
                   totudm_(i) = totudm_this
                   totudm     = totudm_this               + totudm
                   totudm_mul = totudm_this * udm_bias(i) + totudm_mul
                 enddo
                 totttl = totttl + totudm
               endif
! <-- y.sakaki 2021/09

"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
               else if( ityp .eq. 11 .and. abs(kf) .lt .100 ) then
*                       --- lepton (neutrino) ---

                call neutrino_Xsec(kf, ein, hydro, mat, lem, sigt, sigh)
                totttl = totttl + sigt
                tothyd = tothyd + sigh

*-----------------------------------------------------------------------
*                       --- pion ---
*-----------------------------------------------------------------------
"""
a="""\
               else if( ityp .eq. 11 .and. abs(kf) .lt .100 ) then
*                       --- lepton (neutrino) ---

                call neutrino_Xsec(kf, ein, hydro, mat, lem, sigt, sigh)
                totttl = totttl + sigt
                tothyd = tothyd + sigh

! y.sakaki -->
*-----------------------------------------------------------------------
*        user defined interaction for neutrino
*-----------------------------------------------------------------------
               if( iudmodel .ge. 1 ) then
                 totudm_other=totttl
                !-----------------------------
                 totudm    =0d0
                 totudm_mul=0d0
                !-----------------------------
                 if(udm_int_nMax .gt. 0) call fill_mat_info(mat)
                 do i = 1, udm_int_nMax
                   udm_Kin=ein
                   udm_kf_incident=ktyp
                   udm_sigt=0d0
                   call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                   totudm_this = udm_sigt * denm(mat)
                   totudm_(i) = totudm_this
                   totudm     = totudm_this               + totudm
                   totudm_mul = totudm_this * udm_bias(i) + totudm_mul
                 enddo
                 totttl = totttl + totudm
               endif
! <-- ! y.sakaki

*-----------------------------------------------------------------------
*                       --- pion ---
*-----------------------------------------------------------------------
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="getflt.f"
b="""\
*-----------------------------------------------------------------------
*        decay width
*-----------------------------------------------------------------------

         if( tlife .gt. 0.0 ) then

               dcyfl = 1.0 / tlife / rcc
     &               / dsqrt((( ein / rtyp ) + 2.0 )
     &               * ( ein / rtyp ) )

         end if
"""
a="""\
*-----------------------------------------------------------------------
*        decay width
*-----------------------------------------------------------------------

         if( tlife .gt. 0.0 ) then

               dcyfl = 1.0 / tlife / rcc
     &               / dsqrt((( ein / rtyp ) + 2.0 )
     &               * ( ein / rtyp ) )

         end if

      if(iudmodel .gt. 0) then
        if( 0.0d0 < tlife .and. tlife < 1d-19 ) then
          dcyfl=1e+20
        endif
      endif
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="ggm06.f"
b="""\
      subroutine sctpni(eein,wgti,ireg,imat,icmm)
*                                                                      *
*        calculate a collision of a photon with an atom.               *
*        last modified by K.Niita on 2012/12/27                        *
*                                                                      *
************************************************************************
C --- add NS 2020.04 ---
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,iaevt, ibevt, jaevt, jbevt
     &                      ,aevts,aevtr,bevts,bevtr
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
C --- end add NS 2020.04 ---
      use GGMBANKMOD !FURUTA
      use NGSDATAMOD, only : bindeg
      use GGMARRAYMOD !2020ASTOM
!<-20211126murofushi add
      use moddas_material
!--->

      implicit real*8 (a-h,o-z)
"""
a="""\
      subroutine sctpni(eein,wgti,ireg,imat,icmm)
*                                                                      *
*        calculate a collision of a photon with an atom.               *
*        last modified by K.Niita on 2012/12/27                        *
*                                                                      *
************************************************************************
C --- add NS 2020.04 ---
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,iaevt, ibevt, jaevt, jbevt
     &                      ,aevts,aevtr,bevts,bevtr
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
C --- end add NS 2020.04 ---
      use GGMBANKMOD !FURUTA
      use NGSDATAMOD, only : bindeg
      use GGMARRAYMOD !2020ASTOM
!<-20211126murofushi add
      use moddas_material
!--->
      use udm_Parameter
      use udm_Utility
      use udm_Manager

      implicit real*8 (a-h,o-z)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="ggm06.f"
b="""\
C --- change to mod_counter NS 2020.04 ---
C --- delete NS 2020.04 --- 
C for USE_MOD_COUNTER 
CCC      common /kcount/ rncnt(40), rnint(200), rnintr(200),
CCC     &                rnpnt(40), rnpntr(40)
C --- end add NS 2020.04 ---
cABE end

cABE bugfix @2014/11/11
!nais changed
!!      common /pnixs/ totmpni,pnixsm(1000)
      common /pnixs/ totmpni,pnixsm(pnlmax),ipnikind(pnlmax),egypni,
     &               dlibmax
!$OMP THREADPRIVATE(/pnixs/)
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)
      common /nrfdir/ usave(3)
!$OMP THREADPRIVATE(/nrfdir/)

cABE 2019/03/13
      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)

"""
a="""\
C --- change to mod_counter NS 2020.04 ---
C --- delete NS 2020.04 --- 
C for USE_MOD_COUNTER 
CCC      common /kcount/ rncnt(40), rnint(200), rnintr(200),
CCC     &                rnpnt(40), rnpntr(40)
C --- end add NS 2020.04 ---
cABE end

cABE bugfix @2014/11/11
!nais changed
!!      common /pnixs/ totmpni,pnixsm(1000)
      common /pnixs/ totmpni,pnixsm(pnlmax),ipnikind(pnlmax),egypni,
     &               dlibmax
!$OMP THREADPRIVATE(/pnixs/)
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)
      common /nrfdir/ usave(3)
!$OMP THREADPRIVATE(/nrfdir/)

cABE 2019/03/13
      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)

      common /udmodel/ iudmodel 

"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="ggm06.f"
b="""\
cABE add @2014/08/13
            call cputime(22)
cABE end
cABE 2016/03/09
            if( imuinthit .eq. 1 ) then
             call jgdrqd_vpi(izm,inm,erg,igdrqd) ! judging function, GDR or QD or PD?
            else
             igdrqd = jgdrqd(izm,inm,erg) ! judging function, GDR or QD or PD?
            endif

            ipim = 0    ! S.Abe 2015/08/10

*-----------------------------------------------------------------------
*     photonuclear GDR reaction or NRF
*-----------------------------------------------------------------------

            if ( igdrqd .ge. 1 ) then ! GDR or NRF reaction occurs  2014/8/25 ogawa modifie"""
a="""\
cABE add @2014/08/13
            call cputime(22)
cABE end
cABE 2016/03/09
            if( imuinthit .eq. 1 ) then
             call jgdrqd_vpi(izm,inm,erg,igdrqd) ! judging function, GDR or QD or PD?
            else
             igdrqd = jgdrqd(izm,inm,erg) ! judging function, GDR or QD or PD?
            endif

            ipim = 0    ! S.Abe 2015/08/10




            if(iudmodel .ge. 1) then
              mat_Z(1)=izm
              mat_A(1)=izm+inm
              do i = 1, udm_int_nMax
                ! generate final states ------------
                udm_Kin=eein
                udm_kf_incident=22
                udm_logical=.false.
                call user_defined_interaction(-22,i) ! special action: -22
                ! ----------------------------------
                if(udm_logical) goto 999 ! final state is generated
              enddo
            endif






*-----------------------------------------------------------------------
*     photonuclear GDR reaction or NRF
*-----------------------------------------------------------------------

            if ( igdrqd .ge. 1 ) then ! GDR or NRF reaction occurs  2014/8/25 ogawa modifie"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="jam.f"
b="""\
        if(kc.ne.0) chau=chaf(kc,(3-isign(1,kf1))/2)
        return
"""
a="""\
        if(kc.ne.0) chau=chaf(kc,(3-isign(1,kf1))/2)
        return
C...User defined particles. y.sakaki 2021/9 -->
      elseif( 900000 .le. abs(kf1) .and. abs(kf1) .le. 999999 ) then
        if(kf1 .gt. 0) write(chau,'(i6)') kf1
        if(kf1 .lt. 0) write(chau,'(i7)') kf1
        return
C <-- y.sakaki 2021/9
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="main.f"
b="""\
      use partmod ! frtati 2021/10/05
!<-20220513murofushi add
!--->

      implicit real*8 (a-h,o-z)
"""
a="""\
      use partmod ! frtati 2021/10/05
!<-20220513murofushi add
!--->
C y.sakaki 2021/09  For user defined model
      use udm_Parameter
      use udm_Utility
      use udm_Manager

      implicit real*8 (a-h,o-z)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="main.f"
b="""\
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto
      common /tetf1/ nlat3,itetcl,itetvol,itetauto
*-----------------------------------------------------------------------
"""
a="""\
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto
      common /tetf1/ nlat3,itetcl,itetvol,itetauto
*-----------------------------------------------------------------------
C y.sakaki 2021/09   For user defined model
      common /udmodel/ iudmodel 
*-----------------------------------------------------------------------
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="main.f"
b="""\
*-----------------------------------------------------------------------
*     ovly12, ovly13 or ovly14
*-----------------------------------------------------------------------
"""
a="""\

C y.sakaki 2021/09 -->
*-----------------------------------------------------------------------
*     initialization For user defined model
*-----------------------------------------------------------------------
      if(iudmodel .ge. 1) then
*       ----------------------------------------------------------------
        call udm_initialize
        call udm_initialize_Utility
C       -------------
        do i=1,udm_int_nMax
          call user_defined_interaction(1,i) ! Fill Parameters
        enddo
        call user_defined_particle(1) ! Fill Parameters
C       'Initialize' should be after 'Fill Parameters'.
        do i=1,udm_int_nMax
          call user_defined_interaction(3,i) ! Initialize
        enddo
        call user_defined_particle(3) ! Initialize
C       -------------
C       Check whether the entered modules exist
        call udm_check(1)
        call udm_check(2)
C       ----
        do i=1,udm_int_nMax
          udm_logical=.false.
          call user_defined_interaction(2,i) ! check whether name match
          if(.not. udm_logical) then
            print*,
     &"In [ User defined interaction ], there is no module corresponding
     & to the entered Name -> ",trim(udm_int_name(i))
            stop
          endif
        enddo
C       ----
        udm_integer=0
        call user_defined_particle(2) ! Count number of existing modules (for particle)
        if(udm_integer .ne. udm_part_nMax) then
          print*,
     &   "In [ User defined particle ], there is no module corresponding
     & to the entered Name"
          stop
        endif
C       -------------
        do i=1,udm_int_nMax
          call user_defined_interaction(4,i) ! Print Comments
        enddo
        call user_defined_particle(4) ! Print Comments
*       ----------------------------------------------------------------
      endif
C <-- y.sakaki 2021/09

*-----------------------------------------------------------------------
*     ovly12, ovly13 or ovly14
*-----------------------------------------------------------------------
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="makefile"
b="""\
 TARGET   = ../phits_$(ENVFLAGS)
"""
a="""\
 TARGET   = ../phits_$(ENVFLAGS)-udm
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="makefile"
b="""\
	egs5mod.f \\
"""
a="""\
	egs5mod.f \\
	udm_Parameter.f \\
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="makefile"
b="""\
#=======================================================================
# for machine dependent or user defined source and analysis
"""
a="""\
#=======================================================================
# modules Fortran 90
#=======================================================================
SRCS0_90 = \
	$(wildcard udm_part_*.f90) $(wildcard udm_int_*.f90)

#=======================================================================
# for machine dependent or user defined source and analysis
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="makefile"
b="""\
#=======================================================================
# for param.inc
"""
a="""\
#=======================================================================
# user defined model
#=======================================================================
SRCS1_90 = udm_Manager.f90

#=======================================================================
# for param.inc
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="makefile"
b="""\
SRCS = $(SRCS0) $(SRCS1) $(SRCS2) $(SRCS3) $(SRCS4) \\
       $(SRCS5) $(SRCS6) $(SRCS7) $(SRCS8) $(SRCS9) \\
       $(SRCS10) $(SRCS11) 

OBJS = $(SRCS:.f=.o) $(SRCSf90:.f90=.o)
"""
a="""\
OBJS = $(SRCS0:.f=.o) \\
       $(SRCS0_90:.f90=.o) \\
       $(SRCS1:.f=.o) \\
       $(SRCS1_90:.f90=.o) \\
       $(SRCS2:.f=.o) \\
       $(SRCS3:.f=.o) \\
       $(SRCS4:.f=.o) \\
       $(SRCS5:.f=.o) \\
       $(SRCS6:.f=.o) \\
       $(SRCS7:.f=.o) \\
       $(SRCS8:.f=.o) \\
       $(SRCS9:.f=.o) \\
       $(SRCS10:.f=.o) \\
       $(SRCS11:.f=.o) \\
       $(SRCSf90:.f90=.o)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
      use ion_track_structure, only : itsreac
!<-20211126murofushi add
      use moddas_material
!--->
      implicit real*8(a-h,o-z)
"""
a="""\
      use ion_track_structure, only : itsreac
!<-20211126murofushi add
      use moddas_material
!--->
      use udm_Parameter
      use udm_Utility
      use udm_Manager
      implicit real*8(a-h,o-z)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
*-----------------------------------------------------------------------

cABE change @2014/08/13
            call cputime(6)
"""
a="""\
C y.sakaki 2021/09   For user defined model
      common /udmodel/ iudmodel 

*-----------------------------------------------------------------------

cABE change @2014/08/13
            call cputime(6)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
! y.sakaki 2019/09/26 -->
c         if( ipnint .ge. 1 .and. pnimul .gt. 1.d0 ) then
         if( ( ipnint  .ge. 1 .and. pnimul .gt. 1.d0 ) .or.
     &       ( igmuppd .ge. 1 .and. gmumul .gt. 1.0d0 ) ) then
"""
a="""\
! y.sakaki 2021/09 -->
         if( ( ipnint  .ge. 1 .and. pnimul .gt. 1.d0 ) .or.
     &       ( igmuppd .ge. 1 .and. gmumul .gt. 1.0d0 ) .or.
     &         iudmodel .ge. 1  ) then
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
            totmul = totpni * pnimul + totgmm * gmumul + totgam + totegs
"""
a="""\
            totmul = totpni * pnimul + totgmm * gmumul + totgam + totegs
     &             + totudm_mul
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
            prbg3 = totgmm * gmumul / totmul
"""
a="""\
            prbg3 = totgmm * gmumul / totmul
            prbg_udp = totudm_mul / totmul
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
            elseif( r. lt. prbg3 + prbg1 ) then      ! y.sakaki 2021/05/06

                  jcoll = 19
"""
a="""\
            elseif( r. lt. prbg3 + prbg1 ) then      ! y.sakaki 2021/05/06

                  jcoll = 18
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
                  call sctgmm(eein,wgti,ireg,imat)
                  goto 1000
! <-- y.sakaki 2019/09/26
"""
a="""\
                  call sctgmm(eein,wgti,ireg,imat)
                  goto 1000
! <-- y.sakaki 2019/09/26

! y.sakaki 2021/9 -->
*-----------------------------------------------------------------------
            elseif( r .lt. prbg3 + prbg1 + prbg_udp ) then
            ! Choose a process from user-defined processes ( totudm_mul_(i) )
            ! Calculate weight for the choosed process 
              tmp=0d0
              r_udp=rn(0)
              call fill_mat_info(imat)
              do i = 1, udm_int_nMax
                tmp=tmp + (totudm_(i) * udm_bias(i)) / totudm_mul
                udm_Kin=eein
                udm_kf_incident=ktyp
                udm_sigt=0d0
                call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                if(udm_sigt==0d0) cycle
      ! If the charged particle is an incident particle, the energy may have
      ! changed compared to that used to calculate totudm_(i) in getflt.f. 
      ! Processes with a cross section of 0 for the changed energy are eliminated. 
      ! For example, this can be a problem for production processes in narrow 
      ! resonance processes. This could be addressed by
      ! (1) Shortening the following step lengths
      !     chard (Default value = 0.1. For electrons and positrons)
      !     deltm (Default value = 20.12345. For others)
      ! (2) Comment out 'if(udm_sigt==0d0) cycle' and address the problem in the udm_int* file.
                if(r_udp .lt. tmp) then ! index=i is accepted.
                  jcoll = 15
                  wtmg_udp = totmul / delsig / udm_bias(i)
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg_udp
                  oldwt = wt(ibkwt+no,ipomp+1)
                  wgti  = wt(ibkwt+no,ipomp+1)
                  ! generate final states ------------
                  udm_Kin=eein
                  udm_kf_incident=ktyp
                  call user_defined_interaction(21,i)
                  ! ----------------------------------
                  goto 1000
                endif
              enddo
              print*,"Skipped: User Defined Interaction"
! <-- ! y.sakaki 2021/9
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
*---------------------------------------------------------------------
*     electron for EGS5
*---------------------------------------------------------------------
cc H.Iwase

         if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         eein .le. dnmax(12) .and. iegsemi .ne. 0 ) then
"""
a="""\
*---------------------------------------------------------------------
*     electron for EGS5
*---------------------------------------------------------------------
cc H.Iwase

         if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         eein .le. dnmax(12) .and. iegsemi .ne. 0 ) then
     
! y.sakaki 2021/9 -->
            if(iudmodel .ge. 1) then
              totmul = totudm_other + totudm_mul
              prb_udm = totudm_mul / totmul
              if(rn(0) .lt. prb_udm) then
! ----------------------------------------------------------------------
            ! Choose a process from user-defined processes ( totudm_mul_(i) )
            ! Calculate weight for the choosed process 
              tmp=0d0
              r_udp=rn(0)
              call fill_mat_info(imat)
              do i = 1, udm_int_nMax
                tmp=tmp + (totudm_(i) * udm_bias(i)) / totudm_mul
                udm_Kin=eein
                udm_kf_incident=ktyp
                udm_sigt=0d0
                call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                if(udm_sigt==0d0) cycle
      ! If the charged particle is an incident particle, the energy may have
      ! changed compared to that used to calculate totudm_(i) in getflt.f. 
      ! Processes with a cross section of 0 for the changed energy are eliminated. 
      ! For example, this can be a problem for production processes in narrow 
      ! resonance processes. This could be addressed by
      ! (1) Shortening the following step lengths
      !     chard (Default value = 0.1. For electrons and positrons)
      !     deltm (Default value = 20.12345. For others)
      ! (2) Comment out 'if(udm_sigt==0d0) cycle' and address the problem in the udm_int* file.
                if(r_udp .lt. tmp) then ! index=i is accepted.
                  jcoll = 15
                  wtmg_udp = totmul / delsig / udm_bias(i)
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg_udp
                  oldwt = wt(ibkwt+no,ipomp+1)
                  wgti  = wt(ibkwt+no,ipomp+1)
                  ! generate final states ------------
                  udm_Kin=eein
                  udm_kf_incident=ktyp
                  call user_defined_interaction(21,i)
                  ! ----------------------------------
                  goto 1000
                endif
              enddo
              print*,"Skipped: User Defined Interaction"
! ----------------------------------------------------------------------
              end if
            endif
! <-- ! y.sakaki 2021/9
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
          goto 2000   ! S.Abe 2021/06/03, avoid to go to ncasc

         endif

*-----------------------------------------------------------------------
*     Reaction of Neutrino + Deuteron 
"""
a="""\

! y.sakaki 2021/9 -->
            if(iudmodel .ge. 1) then
              totmul = totudm_other + totudm_mul
              prb_udm = totudm_mul / totmul
              if(rn(0) .lt. prb_udm) then
! ----------------------------------------------------------------------
            ! Choose a process from user-defined processes ( totudm_mul_(i) )
            ! Calculate weight for the choosed process 
              tmp=0d0
              r_udp=rn(0)
              call fill_mat_info(imat)
              do i = 1, udm_int_nMax
                tmp=tmp + (totudm_(i) * udm_bias(i)) / totudm_mul
                udm_Kin=eein
                udm_kf_incident=ktyp
                udm_sigt=0d0
                call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                if(udm_sigt==0d0) cycle
      ! If the charged particle is an incident particle, the energy may have
      ! changed compared to that used to calculate totudm_(i) in getflt.f. 
      ! Processes with a cross section of 0 for the changed energy are eliminated. 
      ! For example, this can be a problem for production processes in narrow 
      ! resonance processes. This could be addressed by
      ! (1) Shortening the following step lengths
      !     chard (Default value = 0.1. For electrons and positrons)
      !     deltm (Default value = 20.12345. For others)
      ! (2) Comment out 'if(udm_sigt==0d0) cycle' and address the problem in the udm_int* file.
                if(r_udp .lt. tmp) then ! index=i is accepted.
                  jcoll = 15
                  wtmg_udp = totmul / delsig / udm_bias(i)
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg_udp
                  oldwt = wt(ibkwt+no,ipomp+1)
                  wgti  = wt(ibkwt+no,ipomp+1)
                  ! generate final states ------------
                  udm_Kin=eein
                  udm_kf_incident=ktyp
                  call user_defined_interaction(21,i)
                  ! ----------------------------------
                  goto 1000
                endif
              enddo
              print*,"Skipped: User Defined Interaction"
! ----------------------------------------------------------------------
              end if
            endif
! <-- ! y.sakaki 2021/9

          goto 2000   ! S.Abe 2021/06/03, avoid to go to ncasc

         endif

*-----------------------------------------------------------------------
*     Reaction of Neutrino + Deuteron 
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
*-----------------------------------------------------------------------
*     Non-Nucleons, Non-Elastic Collisions
*-----------------------------------------------------------------------
"""
a="""\
C sakaki -->
*-----------------------------------------------------------------------
*     Interaction of user-defined particles
*-----------------------------------------------------------------------

         if( 900000 .le. abs(ktyp) .and. abs(ktyp) .le. 999999 ) then

! y.sakaki 2021/9 -->
            if(iudmodel .ge. 1) then
              totmul = totudm_other + totudm_mul
              prb_udm = totudm_mul / totmul
              if(rn(0) .lt. prb_udm) then
! ----------------------------------------------------------------------
            ! Choose a process from user-defined processes ( totudm_mul_(i) )
            ! Calculate weight for the choosed process 
              tmp=0d0
              r_udp=rn(0)
              call fill_mat_info(imat)
              do i = 1, udm_int_nMax
                tmp=tmp + (totudm_(i) * udm_bias(i)) / totudm_mul
                udm_Kin=eein
                udm_kf_incident=ktyp
                udm_sigt=0d0
                call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                if(udm_sigt==0d0) cycle
      ! If the charged particle is an incident particle, the energy may have
      ! changed compared to that used to calculate totudm_(i) in getflt.f. 
      ! Processes with a cross section of 0 for the changed energy are eliminated. 
      ! For example, this can be a problem for production processes in narrow 
      ! resonance processes. This could be addressed by
      ! (1) Shortening the following step lengths
      !     chard (Default value = 0.1. For electrons and positrons)
      !     deltm (Default value = 20.12345. For others)
      ! (2) Comment out 'if(udm_sigt==0d0) cycle' and address the problem in the udm_int* file.
                if(r_udp .lt. tmp) then ! index=i is accepted.
                  jcoll = 15
                  wtmg_udp = totmul / delsig / udm_bias(i)
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg_udp
                  oldwt = wt(ibkwt+no,ipomp+1)
                  wgti  = wt(ibkwt+no,ipomp+1)
                  ! generate final states ------------
                  udm_Kin=eein
                  udm_kf_incident=ktyp
                  call user_defined_interaction(21,i)
                  ! ----------------------------------
                  goto 1000
                endif
              enddo
              print*,"Skipped: User Defined Interaction"
! ----------------------------------------------------------------------
              end if
            endif
! <-- ! y.sakaki 2021/9

           goto 1000

         endif
C <-- sakaki

*-----------------------------------------------------------------------
*     Non-Nucleons, Non-Elastic Collisions
*-----------------------------------------------------------------------
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="nreac.f"
b="""\
*-----------------------------------------------------------------------
*     Reaction of Neutrino + Deuteron 
*-----------------------------------------------------------------------

         if( abs(ktyp) .eq. 12 .or. abs(ktyp) .eq. 14 .or. 
     &    abs(ktyp) .eq. 16 ) then ! 2017/11/20 Ogawa Neutrino reaction
"""
a="""\
*-----------------------------------------------------------------------
*     Reaction of Neutrino + Deuteron 
*-----------------------------------------------------------------------

         if( abs(ktyp) .eq. 12 .or. abs(ktyp) .eq. 14 .or. 
     &    abs(ktyp) .eq. 16 ) then ! 2017/11/20 Ogawa Neutrino reaction






! y.sakaki 2021/9 -->
            if(iudmodel .ge. 1) then
              totmul = totudm_other + totudm_mul
              prb_udm = totudm_mul / totmul
              if(rn(0) .lt. prb_udm) then
! ----------------------------------------------------------------------
            ! Choose a process from user-defined processes ( totudm_mul_(i) )
            ! Calculate weight for the choosed process 
              tmp=0d0
              r_udp=rn(0)
              call fill_mat_info(imat)
              do i = 1, udm_int_nMax
                tmp=tmp + (totudm_(i) * udm_bias(i)) / totudm_mul
                udm_Kin=eein
                udm_kf_incident=ktyp
                udm_sigt=0d0
                call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
                if(udm_sigt==0d0) cycle
      ! If the charged particle is an incident particle, the energy may have
      ! changed compared to that used to calculate totudm_(i) in getflt.f. 
      ! Processes with a cross section of 0 for the changed energy are eliminated. 
      ! For example, this can be a problem for production processes in narrow 
      ! resonance processes. This could be addressed by
      ! (1) Shortening the following step lengths
      !     chard (Default value = 0.1. For electrons and positrons)
      !     deltm (Default value = 20.12345. For others)
      ! (2) Comment out 'if(udm_sigt==0d0) cycle' and address the problem in the udm_int* file.
                if(r_udp .lt. tmp) then ! index=i is accepted.
                  jcoll = 15
                  wtmg_udp = totmul / delsig / udm_bias(i)
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg_udp
                  oldwt = wt(ibkwt+no,ipomp+1)
                  wgti  = wt(ibkwt+no,ipomp+1)
                  ! generate final states ------------
                  udm_Kin=eein
                  udm_kf_incident=ktyp
                  call user_defined_interaction(21,i)
                  ! ----------------------------------
                  goto 1000
                endif
              enddo
              print*,"Skipped: User Defined Interaction"
! ----------------------------------------------------------------------
              end if
            endif
! <-- ! y.sakaki 2021/9







"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="read00.f"
b="""\
      use TDCHAINMOD, only: itdc,tdchain,talldcinit !FURUTA20200522
"""
a="""\
      use TDCHAINMOD, only: itdc,tdchain,talldcinit !FURUTA20200522
      use udm_Utility   ! S.ABe 2021/11/04
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="read00.f"
b="""\
cKN 2019/05/27
      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
"""
a="""\
C y.sakaki 2021/10   For user defined model
      common /udmodel/ iudmodel 

cKN 2019/05/27
      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="read00.f"
b="""\
            end if ! S.H. added (2017.8.11)
"""
a="""\
*-----------------------------------------------------------------------
*        [User Defined Interaction]
*-----------------------------------------------------------------------
cABE 2021/11/04

            else if( chcm(i1:i1+23) .eq.
     &                                 '[userdefinedinteraction]' ) then

               if( chcm(i1+24:i1+26) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call read_udinteract(jsn,jsi,dsin,idsi,ill,ilf,
     &                              jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [User Defined Particle]
*-----------------------------------------------------------------------
cABE 2021/11/04

            else if( chcm(i1:i1+20) .eq. '[userdefinedparticle]' ) then

               if( chcm(i1+21:i1+23) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call read_udpart(jsn,jsi,dsin,idsi,ill,ilf,
     &                          jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

            end if ! S.H. added (2017.8.11)
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="read00.f"
b="""\
      common /gmuppd/ igmuppd ! S.Abe 2019/11/07
"""
a="""\
      common /gmuppd/ igmuppd ! S.Abe 2019/11/07
      common /udmodel/ iudmodel ! y.sakaki 2021/09
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="read00.f"
b="""\
cAbe 2019/11/07
*-----------------------------------------------------------------------
*           multiplying factor for photon-induced muon pair production CS
"""
a="""\
C y.sakaki 2021/09
*-----------------------------------------------------------------------
*           user defined process
*-----------------------------------------------------------------------

                  iudmodel = mstz( 157 ) !!! <- iudmodel's "icdf" set in read02.f

                  if(igmuppd .gt. 0) then
                    nspred = 2
                    nspsgn = -1
                    if(mstz( 7 ) .ne. 0) nspred = iabs( mstz( 7 ) )
                  endif

                  if(iudmodel .gt. 0) then
                    i_udm_tmp = 0
                    do i=1,19
                      if(i .eq. 12 .or. i .eq. 13) cycle ! 12 and 13 are 0 by default.
                      if(wc1(i) .ne. 0.5d0) cycle ! Skip for set values.
                      if(i .eq. 2) then
                        ! If neutron's weight cutoff is set to 0, calculation time increase so much.
                        wc1(i) = 0.1d0/pnimul
                      else
                        wc1(i) = 0.0d0
                      endif                        
                      wc2(i) = 0.5d0 * wc1(i) * swtm(i)
                      i_udm_tmp = i_udm_tmp + 1
                    enddo
                    if(i_udm_tmp .ne. 0) 
     & write(*,*) "*** wc1(i) values are modified to alomost zero."
                  endif

cAbe 2019/11/07
*-----------------------------------------------------------------------
*           multiplying factor for photon-induced muon pair production CS
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="read02.f"
b="""\
      data icnu /359/
"""
a="""\
      data icnu /360/    ! increase by 1

*-----------------------------------------------------------------------
c y.sakaki 2021/09/02
      !!! Use the largest index (=icnu)
      !!! Set unused value to "icdf". This value is used as "mstz"'s index in read00.f.
      data chnm(360)/'iudmodel'/
      data icdf(360)/157/     ! 
      data icdl(360)/7/
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="utl02.f"
b="""\
      use NGSDATAMOD, only : bindeg
"""
a="""\
      use NGSDATAMOD, only : bindeg
      use udm_Parameter
      use udm_Manager
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)


# ======================================================================
f="utl02.f"
b="""\
*-----------------------------------------------------------------------

                  rmtyp = rmas(ityp)

            if( ityp .eq. 11 ) then

"""
a="""\
C y.sakaki 2021/09   For user defined model
      common /udmodel/ iudmodel 
*-----------------------------------------------------------------------

                  rmtyp = rmas(ityp)

            if( ityp .eq. 11 ) then

C y.sakaki 2021/9 -->
             if( iudmodel .gt. 0) then
               if(900000 .le. abs(ktyp) .and. abs(ktyp) .le. 999999)then
                 udm_kf_for_12=ktyp
                 udm_mass=-1d0
                 call user_defined_particle(12) ! get 'udm_mass' for 'ktyp'
                 if(udm_mass < 0d0) then
                   rmtyp = 0d0
                 else
                   rmtyp = udm_mass
                 endif
                 return
               endif
             endif
C <-- y.sakaki 2021/9
"""
if check_mode: check(f,b,a)
filename.append(f)
before  .append(b)
after   .append(a)














