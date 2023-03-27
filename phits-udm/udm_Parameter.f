************************************************************************
      module udm_Parameter
************************************************************************

      implicit none

      double precision, save :: udm_Kin
      integer         , save :: udm_kf_incident
      integer         , save :: udm_kf_for_11
      integer         , save :: udm_kf_for_12
      integer         , save :: udm_kf_for_21
      double precision, save :: udm_lifetime
      double precision, save :: udm_mass
      logical         , save :: udm_logical
      integer         , save :: udm_integer
      integer         , save :: udm_counter_dklos
      character(200)  , save :: udm_phits_path

C     From [ User defined interaction ]
      integer, save :: udm_int_nMax
      character(len=99), allocatable,save:: udm_int_name(:)
      double precision, allocatable,save:: udm_bias(:)
      double precision, allocatable,save:: udm_int_param(:,:)
      integer, parameter :: udm_int_param_nMax = 30

C     From [ User defined particle ]
      integer, save :: udm_part_nMax
      character(len=99), allocatable,save:: udm_part_name(:)
      integer, allocatable,save:: udm_part_kf(:)
      double precision, allocatable,save:: udm_part_param(:,:)
      integer, parameter :: udm_part_param_nMax = 30

C     For Material info
      integer, parameter :: num_nuclide_max = 200
      integer, save :: num_nuclide
      integer, save :: mat_Z(num_nuclide_max)
      integer, save :: mat_A(num_nuclide_max)
      double precision, save :: mat_ratio(num_nuclide_max)

C     For "tot" variables in getflt.f
      double precision, save :: udm_sigt, totudm, totudm_mul
      double precision, save :: totudm_other
      double precision, allocatable,save ::  totudm_(:)

C     For final states
      integer, parameter :: nFSmax = 10 ! maximal number of final states
      integer, save :: set_final_state_number
      integer, save :: set_kf(nFSmax)
      integer, save :: set_isomer_level(nFSmax) ! (0: Ground, 1,2: 1st, 2nd isomer)
      double precision, save :: set_Total_Energy_in_MeV(nFSmax)
      double precision, save :: set_Px_in_MeV(nFSmax)
      double precision, save :: set_Py_in_MeV(nFSmax)
      double precision, save :: set_Pz_in_MeV(nFSmax)
      double precision, save :: set_excitation_energy_in_MeV(nFSmax)
      logical         , save :: set_decay_success

      contains

************************************************************************
      subroutine udm_check(action)
************************************************************************
      integer action,i,j
      logical correct
      correct=.true.
C     ------------------------------------------------------------------
      if(action .eq. 1) then
C       [interaction]
        do i=1,udm_int_nMax
          do j=i+1,udm_int_nMax
            if(udm_int_name(i) .eq. udm_int_name(j)) then
              correct=.false.
              exit
            endif
          enddo
        enddo
        if(.not. correct) then
          print*,"In [ user defined interaction ], Name is duplicated."
          call abort
        endif
C       [particle]
        do i=1,udm_part_nMax
          do j=i+1,udm_part_nMax
            if(udm_part_name(i) .eq. udm_part_name(j)) then
              correct=.false.
              exit
            endif
          enddo
        enddo
        if(.not. correct) then
          print*,"In [ user defined particle ], Name is duplicated."
          call abort
        endif
C     ------------------------------------------------------------------
      elseif(action .eq. 2) then
        do i=1,udm_part_nMax
          do j=i+1,udm_part_nMax
            if(udm_part_kf(i) .eq. udm_part_kf(j)) then
              correct=.false.
              exit
            endif
          enddo
        enddo
        if(.not. correct) then
          print*,"In [ user defined particle ], kf is duplicated."
          call abort
        endif
C     ------------------------------------------------------------------
      endif
      end subroutine udm_check

************************************************************************
      subroutine udm_initialize
************************************************************************
      allocate( totudm_( udm_int_nMax ) )
      udm_counter_dklos=0
      end subroutine udm_initialize

************************************************************************
      subroutine initialize_udm_event_info
************************************************************************
      integer i
      do i=1,nFSmax
        set_final_state_number=0
        set_kf(i)=0
        set_isomer_level(i)=0
        set_Total_Energy_in_MeV(i)=0d0
        set_Px_in_MeV(i)=0d0
        set_Py_in_MeV(i)=0d0
        set_Pz_in_MeV(i)=0d0
        set_excitation_energy_in_MeV(i)=0d0
      enddo
      end subroutine initialize_udm_event_info

************************************************************************
      function get_ityp(kf)
************************************************************************
      integer get_ityp, kf
      get_ityp=11
      if(kf.eq.2212) get_ityp=1
      if(kf.eq.2112) get_ityp=2
      if(kf.eq. 211) get_ityp=3
      if(kf.eq. 111) get_ityp=4
      if(kf.eq.-211) get_ityp=5
      if(kf.eq. -13) get_ityp=6
      if(kf.eq.  13) get_ityp=7
      if(kf.eq. 321) get_ityp=8
      if(kf.eq. 311) get_ityp=9
      if(kf.eq.-321) get_ityp=10
C                    get_ityp=11
      if(kf.eq.  11) get_ityp=12
      if(kf.eq. -11) get_ityp=13
      if(kf.eq.  22) get_ityp=14
      if(kf.eq.1000002) get_ityp=15
      if(kf.eq.1000003) get_ityp=16
      if(kf.eq.2000003) get_ityp=17
      if(kf.eq.2000004) get_ityp=18
      if(kf.gt.2000004) get_ityp=19
      return
      end

************************************************************************
      function get_proton_neutron_number(kf)
************************************************************************
      integer get_proton_neutron_number(2)
      integer kf,Z,N
      get_proton_neutron_number(1)=0
      get_proton_neutron_number(2)=0
      if(kf.eq.2212) get_proton_neutron_number(1)=1
      if(kf.eq.2112) get_proton_neutron_number(2)=1
      if(kf.ge.1000002) then
        Z = kf / 1000000
        N = kf - kf / 1000000 * 1000000 - Z
        get_proton_neutron_number(1)=Z
        get_proton_neutron_number(2)=N
      endif
      return
      end

************************************************************************
      function udm_interp(X,ndata,Xdata,Ydata,type)
************************************************************************
      double precision udm_interp
      double precision X,Y
      integer ndata
      double precision,intent(in) :: Xdata(:)
      double precision,intent(in) :: Ydata(:)
      character(len=*) :: type
C     ------------------------------------------------------------------
      integer i,iX1,iX2
      double precision X1,X2,Y1,Y2
C     ------------------------------------------------------------------
      if(X .lt. Xdata(1)) then
        udm_interp=0d0
        return
      endif
C     -------------------
      if(X .gt. Xdata(ndata)) then
        udm_interp=0d0
        return
      endif
C     ------------------------------------------------------------------
      iX1=0
      do i = 2, ndata
        if(X .lt. Xdata(i)) then
          iX1=i-1
          exit
        endif
      enddo
C     ------------------------------------------------------------------
      if(0 .lt. iX1 .and. iX1 .lt. ndata) then
        iX2=iX1+1
        X1=Xdata(iX1); Y1=Ydata(iX1);
        X2=Xdata(iX2); Y2=Ydata(iX2);
C       --------------------------------
        if    (type .eq. "lin-lin") then
          Y=(Y2-Y1)/(X2-X1)*(X-X1)+Y1
C       --------------------------------
        elseif(type .eq. "lin-log") then
          if(Y1<=0d0 .or. Y2<=0d0) then
            Y=0d0
          else
            Y=exp((log(Y2)-log(Y1))/(X2-X1)*(X-X1)+log(Y1))
          endif
C       --------------------------------
        elseif(type .eq. "log-lin") then
          if(X1<=0d0 .or. X2<=0d0) then
            Y=0d0
          else
            Y=(Y2-Y1)/(log(X2)-log(X1))*(log(X)-log(X1))+Y1
          endif
C       --------------------------------
        elseif(type .eq. "log-log") then
          if(X1<=0d0 .or. X2<=0d0 .or. Y1<=0d0 .or. Y2<=0d0) then
            Y=0d0
          else
            Y=exp((log(Y2)-log(Y1))/(log(X2)-log(X1))*(log(X)-log(X1))
     &       +log(Y1))
          endif
C       --------------------------------
        else
          print*,"Caution: udm_interp in type"
          Y=(Y2-Y1)/(X2-X1)*(X-X1)+Y1
C       --------------------------------
        endif
        udm_interp=Y
        return
      else
        print*,"caution: udm_interp"
        udm_interp=0
        return
      endif
C     ------------------------------------------------------------------
      udm_interp=0
      return
      end
************************************************************************
      end module udm_Parameter
************************************************************************



















************************************************************************
      module udm_Utility
************************************************************************

      use udm_Parameter
      implicit real*8(a-h,o-z)
      contains

************************************************************************
      subroutine udm_initialize_Utility
************************************************************************
      common /paran/ icfn(100), ilfn(100), chfn(100) ! phits paths
      character chfn*100
      udm_phits_path=chfn(1)(1:ilfn(1))
      end subroutine udm_initialize_Utility















************************************************************************
*                                                                      *
      subroutine allocate_udinteract
*                                                                      *
*        Last Revised:     2021/11/04                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              Allocate udm_int_name, udm_bias and udm_int_param       *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      allocate(udm_int_name(udm_int_nMax))
      allocate(udm_bias(udm_int_nMax))
      allocate(udm_int_param(udm_int_nMax, udm_int_param_nMax))

      udm_int_name = " "
      udm_bias = 0.d0
      udm_int_param = 0.d0

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine read_udinteract(jsn,jsi,dsin,idsi,ill,ilf,
     &                           jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [User Defined Interaction] section of input files         *
*       modified by S.Abe on 2021/11/04                                *
*                                                                      *
************************************************************************

c      use UDINTERACTION

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      ierr  = 0

*-----------------------------------------------------------------------
*     read first line from jsi
*-----------------------------------------------------------------------

   40 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

      if( iskip .ne. 0 ) goto 40

*-----------------------------------------------------------------------
*     end of section
*-----------------------------------------------------------------------

      if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

         jpn = 1
         return

      end if

*-----------------------------------------------------------------------
*     read n_int
*-----------------------------------------------------------------------

      icl = i1

      if( chlw(icl:icl+4) .eq. 'n_int' ) then

         ic = inumc(chlw,il+1,i3,'=') + 1

            print*,"[n_int] ic,i3 = ",ic,i3
         if( ic .gt. i3 ) then
            ic=i3-1
            print*,"[n_int] ic,i3 = ",ic,i3,"(modified)"
         endif

         if( ic .gt. i3 ) goto 997

         icl = inumc(chlw,ic,i3,';') - 1

         ipm = i

         call onum(chlw,ic,icl,cvvv,ierr)

         if( ierr .ne. 0 ) goto 998

         udm_int_nMax = nint( cvvv )

      else

         goto 996

      endif

*-----------------------------------------------------------------------

      call allocate_udinteract

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

      i_udi = 0

  140 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

      if( iskip .ne. 0 ) goto 140

      if( chlw(i1:i1+3) .eq. 'name' ) goto 140      ! skip data sequence line

*-----------------------------------------------------------------------
*     end of section
*-----------------------------------------------------------------------

      if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

         jpn = 1
         return

      end if

*-----------------------------------------------------------------------

      i_udi = i_udi + 1

*-----------------------------------------------------------------------
*     read name
*-----------------------------------------------------------------------

      ic2  = i1
      ic = jnumc(chlw,ic2,i3)

      ntri = ic
      ntrf = 0
      ic = ic - 1

  151 ic = ic + 1

      if( chlw(ic:ic) .eq. ' ' ) then
         ntrf = ic - 1
         goto 161
      else if( ic .eq. i3 ) then
         goto 995
      else
         goto 151
      endif

  161 continue

      if( ntrf .lt. ntri ) goto 995
      if( ntrf - ntri + 1 .gt. 200 ) goto 994

      ii = 0
      do jj = ntri, ntrf
         ii = ii + 1
         udm_int_name(i_udi)(ii:ii) = chin(jj:jj)
      enddo

*-----------------------------------------------------------------------
*     read bias parameter
*-----------------------------------------------------------------------

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

      if( ierr .ne. 0 ) goto 999

      udm_bias(i_udi) = cvvv

      if( ic2 .ge. i3 ) goto 181

*-----------------------------------------------------------------------
*     read the other parameters
*-----------------------------------------------------------------------

      ipar = 0

  171 ipar = ipar + 1

      ic = jnumc(chlw,ic2,i3)

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

      if( ierr .ne. 0 ) goto 999

      udm_int_param(i_udi,ipar) = cvvv

      if( ic2 .ge. i3 ) goto 181
      goto 171

  181 continue

*-----------------------------------------------------------------------

      if( i_udi .lt. udm_int_nMax ) goto 140

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  994 continue

         m_err = 'Name should be less than 200 characters'//
     &           ' in [user defined interaction]'
         ErrCha = ''
         ErrID = 'L:0/R:read_udinteract/F:udinteract.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Only Name is defined in [user defined interaction]'
         ErrCha = ''
         ErrID = 'L:0/R:read_udinteract/F:udinteract.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'n_int is not defined in the first line'
         ErrCha = ''
         ErrID = 'L:0/R:read_udinteract/F:udinteract.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return


*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:0/R:read_udinteract/F:udinteract.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Value of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:0/R:read_udinteract/F:udinteract.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [user defined interaction] is wrong.'
         ErrCha = ''
         ErrID = 'L:0/R:read_udinteract/F:udinteract.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine

*-----------------------------------------------------------------------


************************************************************************
*                                                                      *
      subroutine allocate_udparticle
*                                                                      *
*        Last Revised:     2021/11/04                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              Allocate udm_part_name, udm_part_kf and udm_part_param                  *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      allocate(udm_part_name(udm_part_nMax))
      allocate(udm_part_kf(udm_part_nMax))
      allocate(udm_part_param(udm_part_nMax, udm_part_param_nMax))

      udm_part_name = " "
      udm_part_kf = 0
      udm_part_param = 0.d0

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine read_udpart(jsn,jsi,dsin,idsi,ill,ilf,
     &                           jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [User Defined Particle] section of input files            *
*       modified by S.Abe on 2021/11/04                                *
*                                                                      *
************************************************************************

c      use UDPARTICLE

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      ierr  = 0

*-----------------------------------------------------------------------
*     read first line from jsi
*-----------------------------------------------------------------------

   40 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

      if( iskip .ne. 0 ) goto 40

*-----------------------------------------------------------------------
*     end of section
*-----------------------------------------------------------------------

      if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

         jpn = 1
         return

      end if

*-----------------------------------------------------------------------
*     read n_part
*-----------------------------------------------------------------------

      icl = i1

      if( chlw(icl:icl+5) .eq. 'n_part' ) then

         ic = inumc(chlw,il+1,i3,'=') + 1

            print*,"[n_part] ic,i3 = ",ic,i3
         if( ic .gt. i3 ) then
            ic=i3-1
            print*,"[n_part] ic,i3 = ",ic,i3,"(modified)"
         endif

         if( ic .gt. i3 ) goto 997

         icl = inumc(chlw,ic,i3,';') - 1

         ipm = i

         call onum(chlw,ic,icl,cvvv,ierr)

         if( ierr .ne. 0 ) goto 998

         udm_part_nMax = nint( cvvv )

      else

         goto 996

      endif

*-----------------------------------------------------------------------

      call allocate_udparticle

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

      i_udp = 0

  140 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#%!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

      if( iskip .ne. 0 ) goto 140

      if( chlw(i1:i1+3) .eq. 'name' ) goto 140      ! skip data sequence line

*-----------------------------------------------------------------------
*     end of section
*-----------------------------------------------------------------------

      if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

         jpn = 1
         return

      end if

*-----------------------------------------------------------------------

      i_udp = i_udp + 1

*-----------------------------------------------------------------------
*     read name
*-----------------------------------------------------------------------

      ic2  = i1
      ic = jnumc(chlw,ic2,i3)

      ntri = ic
      ntrf = 0
      ic = ic - 1

  151 ic = ic + 1

      if( chlw(ic:ic) .eq. ' ' ) then
         ntrf = ic - 1
         goto 161
      else if( ic .eq. i3 ) then
         goto 995
      else
         goto 151
      endif

  161 continue

      if( ntrf .lt. ntri ) goto 995
      if( ntrf - ntri + 1 .gt. 200 ) goto 994

      ii = 0
      do jj = ntri, ntrf
         ii = ii + 1
         udm_part_name(i_udp)(ii:ii) = chin(jj:jj)
      enddo

*-----------------------------------------------------------------------
*     read bias parameter
*-----------------------------------------------------------------------

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

      if( ierr .ne. 0 ) goto 999

      udm_part_kf(i_udp) = nint( cvvv )

      if( ic2 .ge. i3 ) goto 181

*-----------------------------------------------------------------------
*     read the other parameters
*-----------------------------------------------------------------------

      ipar = 0

  171 ipar = ipar + 1

      ic = jnumc(chlw,ic2,i3)

      call snum(chlw,ic,i3,ic2,cvvv,ierr)

      if( ierr .ne. 0 ) goto 999

      udm_part_param(i_udp,ipar) = cvvv

      if( ic2 .ge. i3 ) goto 181
      goto 171

  181 continue

*-----------------------------------------------------------------------

      if( i_udp .lt. udm_part_nMax ) goto 140

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  994 continue

         m_err = 'Name should be less than 200 characters'//
     &           ' in [user defined particle]'
         ErrCha = ''
         ErrID = 'L:0/R:read_udpart/F:udparticle.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Only Name is defined in [user defined particle]'
         ErrCha = ''
         ErrID = 'L:0/R:read_udpart/F:udparticle.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'n_part is not defined in the first line'
         ErrCha = ''
         ErrID = 'L:0/R:read_udpart/F:udparticle.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return


*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:0/R:read_udpart/F:udparticle.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Value of parameter is wrong'
         ErrCha = ''
         ErrID = 'L:0/R:read_udpart/F:udparticle.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [user defined particle] is wrong.'
         ErrCha = ''
         ErrID = 'L:0/R:read_udpart/F:udparticle.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine

*-----------------------------------------------------------------------











************************************************************************
      subroutine fill_mat_info(mk)
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      dimension dnel(1), denh(1), zz(1), a(1), dens(1)
      equivalence ( das, dnel, denh, zz, a, dens )
      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      imat = mk

      num_nuclide=1

      do m = jmd(1+mk), jmd(1+mk+1) - 1
        ie = lme(2,m)
        izm = iza(m) / 1000       ! Z: proton number
        ims = iza(m) - 1000 * izm ! A: mass number
        icntm = m - jmd(1+mk) + 1
C         write(*,"(I5,I5,f12.6)") izm, ims, fme(m)

        ! -------------------------
        if( ims .eq. 0 ) then
C           print*,"@@@",izm,ims
          !--------------------------------------------------------------
          lemm = nint( dnel(kmat(imat)+1) )
          unnorm_AR_sum = 0.0d0 ! averaged cross section with the Unnormalized Abundance Ratio
          do lem = 1, lemm
            izt = zz(kmat(imat)+(lem-1)*3+3)
            if( izt .eq. izm ) then
              unnorm_AR = dens(kmat(imat)+(lem-1)*3+5) ! Unnormalized Abundance Ratio
              unnorm_AR_sum = unnorm_AR_sum + unnorm_AR
            end if
          end do
          !--------------------------------------------------------------
          do lem = 1, lemm
            izt = zz(kmat(imat)+(lem-1)*3+3)
            if( izt .eq. izm ) then
              unnorm_AR = dens(kmat(imat)+(lem-1)*3+5)
              sekt = sekt + unnorm_AR
              ims = nint( a(kmat(imat)+(lem-1)*3+4) ) ! mass number
C               write(*,"(A,I4,I5,f12.6)") "-",izt,ims,
C      &          unnorm_AR/unnorm_AR_sum*fme(m)
              mat_Z(num_nuclide)=izt
              mat_A(num_nuclide)=ims
              mat_ratio(num_nuclide)=unnorm_AR/unnorm_AR_sum*fme(m)
              num_nuclide=num_nuclide+1
              if(num_nuclide .gt. num_nuclide_max) then
                print*,
     &               "[CAUTION] Increase 'num_nuclide_max' 200 or more."
                return
              endif
            end if
          end do
        !--------------------------------------------------------------
        else
          mat_Z(num_nuclide)=izm
          mat_A(num_nuclide)=ims
          mat_ratio(num_nuclide)=fme(m)
          num_nuclide=num_nuclide+1
          if(num_nuclide .gt. num_nuclide_max) then
            print*,"[CAUTION] Increase 'num_nuclide_max' 200 or more."
            return
          endif
          endif
      end do
      num_nuclide=num_nuclide-1

      end subroutine fill_mat_info

************************************************************************
      function get_random(rmin,rmax)
************************************************************************
      implicit real*8 (a-h,o-z)
      get_random = rmin+(rmax-rmin)*unirn(dummy)
      return
      end

************************************************************************
      function get_random_0to1()
************************************************************************
      implicit real*8 (a-h,o-z)
      get_random_0to1 = unirn(dummy)
      return
      end

************************************************************************
      subroutine two_body_decay_uniform(kf1,kf2)
************************************************************************
      implicit real*8 (a-h,o-z)
      integer kf1, kf2
      double precision m0, m1, m2
      m0=get_mass(udm_kf_for_21)
      m1=get_mass(kf1)
      m2=get_mass(kf2)
C       if(udm_kf_for_21==999001)print*,udm_kf_for_21,"-->",kf1,kf2
C       if(udm_kf_for_21==999001)print*,m0,m1,m2
      if(m0 <= m1+m2) then
        print*,"*** [CAUTION] two_body_decay_uniform failed"
        print*,"*** ",udm_kf_for_21,"-->",kf1,kf2
        return
      endif
      call dklos2(
     &      get_ityp(udm_kf_for_21), 
     &               udm_kf_for_21, 
     &               udm_Kin, 
     &      get_ityp(kf1), kf1,
     &      get_ityp(kf2), kf2)
*-----------------------------------------------------------------------
      set_decay_success=.true.
*-----------------------------------------------------------------------
      end subroutine two_body_decay_uniform

************************************************************************
      subroutine three_body_decay_uniform(kf1,kf2,kf3)
************************************************************************
      implicit real*8 (a-h,o-z)
      integer kf1, kf2, kf3
      double precision m0, m1, m2, m3
      m0=get_mass(udm_kf_for_21)
      m1=get_mass(kf1)
      m2=get_mass(kf2)
      m3=get_mass(kf3)
C       if(udm_kf_for_21==999001)print*,udm_kf_for_21,"-->",kf1,kf2,kf3
C       if(udm_kf_for_21==999001)print*,m0,m1,m2,m3
      if(m0 <= m1+m2+m3) then
        print*,"*** [CAUTION] three_body_decay_uniform failed"
        print*,"*** ",udm_kf_for_21,"-->",kf1,kf2,kf3
        return
      endif
      call dklos3(
     &      get_ityp(udm_kf_for_21), 
     &               udm_kf_for_21, 
     &               udm_Kin, 
     &      get_ityp(kf1), kf1,
     &      get_ityp(kf2), kf2,
     &      get_ityp(kf3), kf3)
*-----------------------------------------------------------------------
      set_decay_success=.true.
*-----------------------------------------------------------------------
      end subroutine three_body_decay_uniform

************************************************************************
      function get_mass(kf)
************************************************************************
      double precision get_mass
      integer ityp, kf
      ityp=get_ityp(kf)
      get_mass=rmtyp(ityp,kf)
      return
      end

************************************************************************
      subroutine fill_final_state
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      parameter ( rpmass = 938.27d0, rnmass = 939.58d0 )
      real*8 event(20)

*-----------------------------------------------------------------------

      include 'param00.inc'
      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)

      integer iZN_i(2)
      double precision mass_i_MeV
      double precision absp_i
      double precision total_energy_MeV
      double precision kinetic_energy_MeV

*-----------------------------------------------------------------------

      do i=1,set_final_state_number

        nclsts = i

        kf_i  =set_kf(i)
        ityp_i=get_ityp(kf_i)
        iZN_i =get_proton_neutron_number(kf_i)

        iclusts(nclsts) = ipatf(ityp_i, kf_i)

        jclusts(0,nclsts) = 0
        jclusts(1,nclsts) = iZN_i(1)
        jclusts(2,nclsts) = iZN_i(2)
        jclusts(3,nclsts) = ityp_i
        jclusts(4,nclsts) = 0
        jclusts(5,nclsts) = ichgf(ityp_i, kf_i)
        jclusts(6,nclsts) = ibryf(ityp_i, kf_i)
        jclusts(7,nclsts) = kf_i
        jclusts(8,nclsts) = set_isomer_level(i)

        mass_i_MeV = rmtyp(ityp_i, kf_i)
        total_energy_MeV = set_Total_Energy_in_MeV(i)
        kinetic_energy_MeV = total_energy_MeV - mass_i_MeV
        absp_i = sqrt( set_Px_in_MeV(i)**2
     &               + set_Py_in_MeV(i)**2
     &               + set_Pz_in_MeV(i)**2 )
        qclusts(0,nclsts)  = 0.0
        qclusts(1,nclsts)  = set_Px_in_MeV(i) / absp_i
        qclusts(2,nclsts)  = set_Py_in_MeV(i) / absp_i
        qclusts(3,nclsts)  = set_Pz_in_MeV(i) / absp_i
        qclusts(4,nclsts)  = total_energy_MeV / 1000d0 ! GeV
        qclusts(5,nclsts)  = mass_i_MeV       / 1000d0 ! GeV
        qclusts(6,nclsts)  = set_excitation_energy_in_MeV(i)
        qclusts(7,nclsts)  = kinetic_energy_MeV
        qclusts(8,nclsts)  = 1.0 ! wgt / wg0
        qclusts(9,nclsts)  = 0.0
        qclusts(10,nclsts) = 0.0d0
        qclusts(11,nclsts) = 0.0d0
        qclusts(12,nclsts) = 0.0d0

      enddo
*-----------------------------------------------------------------------
      set_decay_success=.true.
*-----------------------------------------------------------------------
      end subroutine fill_final_state

************************************************************************
      subroutine multiply_weight_by(ratio)
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
*-----------------------------------------------------------------------
      wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * ratio
      end subroutine multiply_weight_by

************************************************************************
      end module udm_Utility
************************************************************************



