*=======================================================================
*                 S U B R O U T I N E  B 2 W A O S
*=======================================================================
*
*
*  1. PURPOSE
*
*     B2WAOS PRINTS OUT THE RESULTS OF A CALCULATION OF THE B2.5 EDGE
*     PLASMA CODE IN OLD-FASHIONED STYLE
*
*
*=======================================================================
*//DECLARATIONS//

      SUBROUTINE B2WAOS (NPRT, nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)

      use b2mod_types
      use b2mod_time
      use b2mod_constants
      use b2mod_diag
     . , only : label
      use b2mod_geo
      use b2mod_transport_namelist
!pb      use b2mod_plasma
      use b2mod_b2cmpa
      use b2mod_boundary_namelist
      use b2mod_residuals
      use b2mod_external
      use b2mod_rates
      implicit none

      integer, intent(in) :: nprt, nx, ny, ns
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny)

      integer ncall, model
      integer, save :: style = 1

*  -- LOCAL VARIABLES --
      INTEGER :: IX, IY, IS, IB, IXR, IYR
      REAL (KIND=R8) :: T0,RZO,CS,MACH,T1
      real (kind=R8) ::
     * fni(-1:nx,-1:ny,0:1), fne(-1:nx,-1:ny,0:1),
     2 ni(-1:nx,-1:ny,0:1),
     3 rz(-1:nx,-1:ny), pz(-1:nx,-1:ny)

*  -- PROCEDURES --
      REAL (KIND=R8) ::
     .      X, Y, XF, YF, RVX, RVY, RDX, RDY, RSX, RSY, PITX, PITY,
     .      SFEE, SFEI, SFNI(NBCD,0:NS)
      REAL (KIND=R8) :: sxhelp,syhelp
      real (kind=R8) :: XX(-1:NX,6), YY(-1:NY,6)
      real (kind=R8) :: b2tqca_phm0, b2treq_phm0
      real (kind=R8) ::
     * b2ssum, rfnax, rfnay, rfnix, rfniy, rfnex, rfney,
     &  rfhex, rfhey, rfhix, rfhiy, rfchx, rfchy
      EXTERNAL XERRAB
      X(IX)   = XX(IX,1)
      Y(IY)   = YY(IY,1)
      XF(IX)  = XX(IX,2)
      YF(IY)  = YY(IY,2)
      RVX(IX) = XX(IX,3)
      RVY(IY) = YY(IY,3)
      RDX(IX) = XX(IX,4)
      RDY(IY) = YY(IY,4)
      RSX(IX) = XX(IX,5)
      RSY(IY) = YY(IY,5)
      PITX(IX) = XX(IX,6)
      PITY(IY) = YY(IY,6)
      rfnax(ix,is) = b2ssum(ny+2,fna(ix,-1,0,is),nx+2)
      rfnay(iy,is) = b2ssum(nx+2,fna(-1,iy,1,is),1)
      rfhex(ix) = b2ssum(ny+2,fhe(ix,-1,0),nx+2)
      rfhey(iy) = b2ssum(nx+2,fhe(-1,iy,1),1)
      rfhix(ix) = b2ssum(ny+2,fhi(ix,-1,0),nx+2)
      rfhiy(iy) = b2ssum(nx+2,fhi(-1,iy,1),1)
      rfchx(ix) = b2ssum(ny+2,fch(ix,-1,0),nx+2)
      rfchy(iy) = b2ssum(nx+2,fch(-1,iy,1),1)
      rfnex(ix) = b2ssum(ny+2,fne(ix,-1,0),nx+2)
      rfney(iy) = b2ssum(nx+2,fne(-1,iy,1),1)
      rfnix(ix) = b2ssum(ny+2,fni(ix,-1,0),nx+2)
      rfniy(iy) = b2ssum(nx+2,fni(-1,iy,1),1)

      data ncall/0/, b2tqca_phm0/1.0_R8/, model/3/, b2treq_phm0/1.0_R8/


*=======================================================================
*//COMPUTATION//

*  ---------------------------------------------------------------------
*  PRELIMINARIES.
*  ---------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2waos')
*   ..test nget
      call xertst (1.le.nprt, 'faulty argument nprt')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

      if (ncall.eq.0) then
        call ipgeti ('b2tqca_model', model)
        call ipgetr ('b2tqca_phm0', b2tqca_phm0)
        call xertst (0.0_R8.le.b2tqca_phm0, 'faulty input phm0')
        call ipgeti ('b2tlmv_style', style)                          !srv 30.05.08
        call ipgetr ('b2treq_phm0', b2treq_phm0)
      end if

*  -- TEST EDGE VALUES OF SX AND SY --
      if (sum(abs(gs(-1,-1:ny,0))) > 0._R8)
     .   call xerrab ('b2waos--error, faulty edge gs(-1,:,0)')
      if (sum(abs(gs(-1:nx,-1,1))) > 0._R8)
     .   call xerrab ('b2waos--error, faulty edge gs(:,-1,1)')

*   ..compute auxiliaries
      call b2xpni (nx, ny, ns, na, ni)
      call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
      call b2xprz (nx, ny, ns, mp, am, na, rz)
      call b2xpfe (nx, ny, ns, ns_ext, qe,
     &  rza, za_ext, fna, fa_ext, fch, fne)
      call b2xpfi (nx, ny, ns, fna, fni)
*
      IXR = -1
      IYR = NY/2
      IF (NNCUT > 0) THEN
        IXR = rightcut(1)
        IYR = topcut(1) + 2
      END IF
      IYR = MIN(IYR,NY)

*  -- X-COORDINATE LABELS --

*     CELL CENTERS IN XX(0:NX+1,1), CELL FACES IN XX(0:NX,2);
*     INTEGRATE DISTANCE ACROSS CELL FACES AT IY=IYR.
      XX(IXR,1) = -0.5_R8*GS(IXR,IYR,1)/
     .             (PI*(CRX(IXR,IYR,2)+CRX(IXR,IYR,3)))
      XX(IXR,2) = 0
      DO IX = IXR-1, -1, -1
         XX(IX,2) = XX(IX+1,2) - GS(IX+1,IYR,1)
     .              /(PI*(CRX(IX+1,IYR,2)+CRX(IX+1,IYR,3)))
         XX(IX,1) = XX(IX,2) - 0.5_R8*GS(IX,IYR,1)
     .              /(PI*(CRX(IX,IYR,2)+CRX(IX,IYR,3)))
      END DO
      DO IX = IXR+1, NX
         XX(IX,1) = XX(IX-1,2) +
     .              0.5_R8*GS(IX,IYR,1)
     .              /(PI*(CRX(IX,IYR,2)+CRX(IX,IYR,3)))
         XX(IX,2) = XX(IX-1,2) + GS(IX,IYR,1)
     .              /(PI*(CRX(IX,IYR,2)+CRX(IX,IYR,3)))
      END DO
*     (XX(NX+1,2) IS ACTUALLY IMMATERIAL)

*  -- Y-COORDINATE LABELS --
*     CELL CENTERS IN YY(0:NY+1,1), CELL FACES IN YY(0:NY,2);
*     INTEGRATE DISTANCE ACROSS CELL FACES AT IX=IXR.
      YY(IYR,1) = -0.5_R8*GS(IXR+1,IYR,0)
     .            /(PI*(CRX(IXR,IYR,3)+CRX(IXR,IYR,1)))
      YY(IYR,2) = 0
      DO IY = IYR-1, -1, -1
         YY(IY,2) = YY(IY+1,2) - GS(IXR+1,IY+1,0)
     .              /(PI*(CRX(IXR,IY+1,3)+CRX(IXR,IY+1,1)))
         YY(IY,1) = YY(IY,2) - 0.5_R8*GS(IXR+1,IY,0)
     .              /(PI*(CRX(IXR,IY,3)+CRX(IXR,IY,1)))
      END DO
      DO IY = IYR+1, NY
         YY(IY,1) = YY(IY-1,2) +
     .              0.5_R8*GS(IXR+1,IY,0)
     .              /(PI*(CRX(IXR,IY,3)+CRX(IXR,IY,1)))
         YY(IY,2) = YY(IY-1,2) + GS(IXR+1,IY,0)
     .              /(PI*(CRX(IXR,IY,3)+CRX(IXR,IY,1)))
      END DO
*     (YY(NY+1,2) IS ACTUALLY IMMATERIAL)
*  ---------------------------------------------------------------------

*  -- TEST XX(*,1:2) AND YY(*,1:2) --
      DO IX = -1, NX-1
         IF (XX(IX,1) .GT. XX(IX,2) .OR.
     .         XX(IX,2) .GT. XX(IX+1,1)) THEN
            CALL XERRAB ('B2WAOS -- ERROR INVOLVING XX')
         END IF
      END DO
      DO IY = -1, NY-1
         IF (YY(IY,1) .GT. YY(IY,2) .OR.
     .         YY(IY,2) .GT. YY(IY+1,1)) THEN
            CALL XERRAB ('B2WAOS -- ERROR INVOLVING YY')
         END IF
      END DO


*  -- COMPUTE 1D X-REDUCTIONS --
      DO IX = -1, NX
*     -- 'VOL(IX)' --
         XX(IX,3) = SUM (VOL(IX,-1:NY))
*     -- 'DX(IX)' --
         T0 = 0._R8
         DO IY = -1, NY
            T0 = T0 + VOL(IX,IY)/HX(IX,IY)
         end do
         XX(IX,4) = XX(IX,3) / T0
*     -- 'SX(IX)' --
         XX(IX,5) = SUM (GS(IX,-1:NY,0))
*     -- 'PIT(IX)' --
         XX(IX,6) = SUM(VOL(IX,-1:NY)*BB(IX,-1:NY,0)/BB(IX,-1:NY,3))
     .              / XX(IX,3)
      end do

*  -- COMPUTE 1D Y-REDUCTIONS --
      DO IY = -1, NY
*     -- 'VOL(IY)' --
         YY(IY,3) = SUM(VOL(-1:NX,IY))
*     -- 'DY(IY)' --
         T0 = 0
         DO IX = -1, NX
            T0 = T0 + VOL(IX,IY)/HY(IX,IY)
         end do
         YY(IY,4) = YY(IY,3) / T0
*     -- 'SY(IY)' --
         YY(IY,5) = SUM(GS(-1:NX,IY,1))
*     -- 'PIT(IY)' --
         T0 = 0
         T1 = 0
         DO IX = -1, NX
            T0 = T0 + HX(IX,IY)
            T1 = T1 + HX(IX,IY) / (BB(IX,IY,0)/BB(IX,IY,3))
         end do
         YY(IY,6) = T0 / T1
      end do

*  ---------------------------------------------------------------------
*  PRINT THE STATE.
*  ---------------------------------------------------------------------
      WRITE (NPRT, 101) LABEL

      WRITE (NPRT, 1001) TIM

      WRITE (NPRT, 102) NX, NY, NS

      WRITE (NPRT, 103)
      DO IX = -1, NX
         T0 = SUM(HY(IX,-1:ny))
         IF (IX .LT. NX) THEN
            WRITE (NPRT, '(2X, I3, 1X, 1P, 7D11.3, 2X, 2D11.3)')
     .         IX, X(IX), RVX(IX), RDX(IX), T0, RVX(IX)/(RDX(IX)*T0),
     .         PITX(IX), RVX(IX)/RDX(IX), XF(IX), RSX(IX)
         ELSE IF (IX .EQ. NX) THEN
            WRITE (NPRT, '(2X, I3, 1X, 1P, 7D11.3)')
     .         NX, X(NX), RVX(NX), RDX(NX),
     .         T0, RVX(NX)/(RDX(NX)*T0),
     .         PITX(NX), RVX(NX)/RDX(NX)
         END IF
      END DO

      WRITE (NPRT, 104)
      DO IY = NY, -1, -1
         T0 = SUM(HX(-1:NX,IY))
         IF (IY .EQ. NY) THEN
            WRITE (NPRT, '(2X, I3, 1X, 1P, 7D11.3)')
     .         NY, Y(NY), RVY(NY), T0, RDY(NY),
     .         RVY(NY)/(T0*RDY(NY)),
     .         PITY(NY), RVY(NY)/RDY(NY)
         ELSE IF (IY .LT. NY) THEN
            WRITE (NPRT, '(2X, I3, 1X, 1P, 7D11.3, 2X, 2D11.3)')
     .         IY, Y(IY), RVY(IY), T0, RDY(IY), RVY(IY)/(T0*RDY(IY)),
     .         PITY(IY), RVY(IY)/RDY(IY), YF(IY), RSY(IY)
         END IF
      END DO

      WRITE (NPRT, '(/)')

!      CALL PHYPRT (NS, MI, ZI, NPRT)

      IF (MINVAL(AM(1:NS)) .LT. 0 .OR.
     .      MINVAL(ZN(1:NS)) .LT. 0) THEN
         CALL XERRAB ('B2WAOS -- FAULTY ARGUMENT MI, ZI')
      END IF

*  -- OUTPUT --
      do is = 0,ns-1
        if (is_neutral(is)) cycle
        write (nprt,
     .      '(//1x,A,I12,4(/1x,A,ES12.2),5(/1x,a,es12.2,a,i6))')
     .      ' SPECIES   ', is,
     .      ' MASS MU   ', am(is),
     .      ' CHARGE ZI ', zn(is),
     .      ' EL.TH.DIF ', PARM_HCE,
     .      ' IO.TH.DIF ', PARM_HCI(is),
     .      ' PAR.VIS   ', b2tqca_phm0, '  model flag ', model,
     .      ' TRA.VIS   ', PARM_VSA(is),'  model flag ', flag_vsa,
     .      ' DIF.NI    ', PARM_DNA(is),'  model flag ', flag_dna,
     .      ' DIF.PR    ', PARM_DPA(is),'  model flag ', flag_dpa,
     .      ' V.CONV    ', PARM_VLA(is),'  model flag ', flag_vla
      end do

      write (nprt,'(//1x,A,/1x,A,i12,6(/1x,a,es12.3))')
     .      ' PARAMETERS FOR ENERGY TRANSPORT:',
     .      ' FLUX LIMIT MODEL FLAG           ', style,
     .      ' FLUX LIMIT COEFFICIENT (FLIM)   ', cflme,
     .      ' FLUX LIMIT COEFFICIENT (FLIMI)  ', cflmi,
     .      ' FLUX LIMIT COEFFICIENT (FLIMV)  ', cflmv,
     .      ' EQUIPARTITION MULTIPLIER ..     ', b2treq_phm0


      WRITE (NPRT, 110) ARESHE, ARESHI, (ARESCO(IS), IS = 0, NS-1),
     .                  (ARESMO(IS), IS = 0, NS-1)

      WRITE (NPRT, 111)
      DO IX = -1, NX-1
         WRITE (NPRT, '(1PE12.3, E11.3, 2X, 2D11.3, (T50, 10D11.3))')
     .      XF(IX), RSX(IX), RFHEX(IX), RFHIX(IX), RFCHX(IX),
     .      RFNEX(IX), RFNIX(IX),
     .      (RFNAX(IX,IS), IS = 0, NS-1)
      ENDDO

      WRITE (NPRT, 112)
      DO IY = NY, 0, -1
         WRITE (NPRT, '(1PE12.3, E11.3, 2X, 2D11.3, (T50, 10D11.3))')
     .      YF(IY), RSY(IY), RFHEY(IY), RFHIY(IY), RFCHY(IY),
     .      RFNEY(IY), RFNIY(IY),
     .      (RFNAY(IY,IS), IS = 0, NS-1)
      END DO

      WRITE (NPRT,1112)

      WRITE (NPRT,*) ' FLUXES INTEGRATED ALONG BOUNDARY '
!      IF(NLCURR.EQ.0) THEN
        SFNI=0.D0
        DO IB=1,NBC
        SFEE=0.D0
        SFEI=0.D0
        IF (BCCHAR(IB) .EQ. 'S') THEN
        WRITE (NPRT, 113) IB, BCPOS(IB)
        DO IX = BCSTART(IB),BCEND(IB)
ctb       WRITE (NPRT, '(/1PE12.3, 4D11.3, 2X, (T59, 6D11.3))')
          syhelp=max(gs(ix,bcpos(ib)+1,1),1.D-20)
          WRITE (NPRT, '(/I4,8X, 4(1PE11.3), 2X, (T59, 6D11.3))')
     .      IX, TE(IX,BCPOS(IB))/EV, TI(IX,BCPOS(IB))/EV,
     .      NE(IX,BCPOS(IB)),
     .      PZ(IX,BCPOS(IB)), (NA(IX,BCPOS(IB),IS), IS = 0, NS-1)
          WRITE (NPRT, '(1PE12.3, 4D11.3, 2X, (T59, 6D11.3))')
     .      GS(IX,BCPOS(IB)+1,1), FHE(IX,BCPOS(IB)+1,1)/SYhelp,
     .      FHI(IX,BCPOS(IB)+1,1)/SYhelp,
     .      SUM(ZN(0:ns-1)*FNA(IX,BCPOS(IB)+1,1,0:NS-1))/syhelp,
     .         PO(IX,BCPOS(IB)),
     .      (FNA(IX,BCPOS(IB)+1,1,IS)/SYhelp, IS = 0, NS-1)
          WRITE (NPRT, '((1X, T13, 2ES11.3, (T59,6D11.3)))')
     .        FHE(IX,BCPOS(IB)+1,1)/(FNA(IX,BCPOS(IB)+1,1,1)+1.D-30)
     .                        /TE(IX,BCPOS(IB)),
     .        FHI(IX,BCPOS(IB)+1,1)/(FNA(IX,BCPOS(IB)+1,1,1)+1.D-30)
     .                        /TI(IX,BCPOS(IB)),
     .        (UA(IX,BCPOS(IB)+1,IS), IS = 0, NS-1),
     .        (UA(IX,BCPOS(IB),IS), IS = 0, NS-1)
          WRITE (NPRT, '(1PE12.3, 3D11.3, 2X, (T59, 6D11.3))')
     .      GS(IX,BCPOS(IB)+1,1),
     .      FHE(IX,BCPOS(IB)+1,1), FHI(IX,BCPOS(IB)+1,1),
     .      SUM(ZN(0:ns-1)*FNA(IX,BCPOS(IB)+1,1,0:NS-1)),
     .      (FNA(IX,BCPOS(IB)+1,1,IS), IS = 0, NS-1)
          SFEE=SFEE+FHE(IX,BCPOS(IB)+1,1)
          SFEI=SFEI+FHI(IX,BCPOS(IB)+1,1)
          DO IS = 0, NS-1
            if (.not.is_neutral(is))
     .       SFNI(IB,IS)=SFNI(IB,IS)+FNA(IX,BCPOS(IB)+1,1,IS)
          ENDDO
        END DO

        WRITE (NPRT,'(1X,3(A,1PD12.4))') ' SFEE = ',SFEE,' SFEI = ',
     .         SFEI,' SFNI = ',SUM(SFNI(IB,0:NS-1))
        ENDIF

        IF (BCCHAR(IB) .EQ. 'N') THEN
        WRITE (NPRT, 114) IB, BCPOS(IB)
        DO IX = BCSTART(IB),BCEND(IB)
          syhelp=max(gs(ix,bcpos(ib),1),1D-20)
          WRITE (NPRT, '(/I4,8X, 4(1PE11.3), 2X, (T59, 6D11.3))')
     .      IX, TE(IX,BCPOS(IB))/EV, TI(IX,BCPOS(IB))/EV,
     .      NE(IX,BCPOS(IB)),
     .      PZ(IX,BCPOS(IB)), (NA(IX,BCPOS(IB),IS), IS = 0, NS-1)
          WRITE (NPRT, '(1PE12.3, 4D11.3, 2X, (T59, 6D11.3))')
     .      GS(IX,BCPOS(IB),1), FHE(IX,BCPOS(IB),1)/SYhelp,
     .      FHI(IX,BCPOS(IB),1)/SYhelp,
     .      SUM(ZN(0:NS-1)*FNA(IX,BCPOS(IB),1,0:NS-1)) /
     .         SYhelp, Po(IX,BCPOS(IB)),
     .      (FNA(IX,BCPOS(IB),1,IS)/SYhelp, IS = 0, NS-1)
          WRITE (NPRT, '((1X, T13, 3ES11.3, (T59,6D11.3)))')
     .      FHE(IX,BCPOS(IB),1)/(FNA(IX,BCPOS(IB),1,1)+1.D-30)
     .                        /TE(IX,BCPOS(IB)),
     .      FHI(IX,BCPOS(IB),1)/(FNA(IX,BCPOS(IB),1,1)+1.D-30)
     .                        /TI(IX,BCPOS(IB)),
     .      (UA(IX,BCPOS(IB)-1,IS), IS = 0, NS-1),
     .      (UA(IX,BCPOS(IB),IS), IS = 0, NS-1)
          WRITE (NPRT, '(1PE12.3, 3D11.3, 2X, (T59, 6D11.3))')
     .      GS(IX,BCPOS(IB),1),
     .      FHE(IX,BCPOS(IB),1), FHI(IX,BCPOS(IB),1),
     .      SUM(ZN(0:NS-1)*FNA(IX,BCPOS(IB),1,0:NS-1)),
     .      (FNA(IX,BCPOS(IB),1,IS), IS = 0, NS-1)
          SFEE=SFEE-FHE(IX,BCPOS(IB),1)
          SFEI=SFEI-FHI(IX,BCPOS(IB),1)
          DO IS = 0, NS-1
            if (.not.is_neutral(is))
     .        SFNI(IB,IS)=SFNI(IB,IS)-FNA(IX,BCPOS(IB),1,IS)
          ENDDO
        END DO
        WRITE (NPRT,'(1X,3(A,1PD12.4))') ' SFEE = ',SFEE,' SFEI = ',
     .         SFEI,' SFNI = ',SUM(SFNI(IB,0:NS-1))
        ENDIF

        IF (BCCHAR(IB) .EQ. 'W') THEN
        WRITE (NPRT, 115) IB, BCPOS(IB)
        DO IY = BCEND(IB),BCSTART(IB), -1
          sxhelp=max(gs(bcpos(ib)+1,iy,0),1D-20)
          WRITE (NPRT, '(/I4,8X, 4(1PE11.3), 2X, (T59, 6D11.3))')
     .      IY, TE(BCPOS(IB),IY)/EV, TI(BCPOS(IB),IY)/EV,
     .      NE(BCPOS(IB),IY), PZ(BCPOS(IB),IY),
     .      (NA(BCPOS(IB),IY,IS), IS = 0, NS-1)
          WRITE (NPRT, '(1PE12.3, 4D11.3, 2X, (T59, 6D11.3))')
     .      GS(BCPOS(IB)+1,IY,0), FHE(BCPOS(IB)+1,IY,0)/SXhelp,
     .      FHI(BCPOS(IB)+1,IY,0)/SXhelp,
     .      SUM(ZN(0:NS-1)*FNA(BCPOS(IB)+1,IY,1,0:ns-1)) /
     .         SXhelp, PO(BCPOS(IB),IY),
     .      (FNA(BCPOS(IB)+1,IY,1,IS)/SXhelp, IS = 0, NS-1)
          RZO = 0.0_R8
          DO IS = 0, NS-1
            IF (IS_NEUTRAL(IS)) CYCLE
            RZO = RZO + AM(IS)*MP*NA(BCPOS(IB),IY,IS)
          END DO
          CS = SQRT(PZ(BCPOS(IB),IY) / RZO)
          MACH=-UA(BCPOS(IB),IY,1)/CS
          WRITE (NPRT, '((1X, T13, 4ES11.3, T70,6D11.3))')
     .      FHE(BCPOS(IB)+1,IY,0)/(FNA(BCPOS(IB)+1,IY,0,1)+1.E-30)
     .                        /TE(BCPOS(IB),IY),
     .      FHI(BCPOS(IB)+1,IY,0)/(FNA(BCPOS(IB)+1,IY,0,1)+1.E-30)
     .                        /TI(BCPOS(IB),IY),
     .      MACH,(UA(BCPOS(IB),IY,IS),IS=0,NS-1)
          WRITE (NPRT, '(1PE12.3, 3D11.3, 2X, (T59, 6D11.3))')
     .      GS(BCPOS(IB)+1,IY,0),
     .      FHE(BCPOS(IB)+1,IY,0),FHI(BCPOS(IB)+1,IY,0),
     .      SUM(ZN(0:NS-1)*FNA(BCPOS(IB)+1,IY,0,0:NS-1)),
     .      (FNA(BCPOS(IB)+1,IY,0,IS), IS = 0, NS-1)
          SFEE=SFEE+FHE(BCPOS(IB)+1,IY,0)
          SFEI=SFEI+FHI(BCPOS(IB)+1,IY,0)
          DO IS = 0, NS-1
            if (.not.is_neutral(is))
     .        SFNI(IB,IS)=SFNI(IB,IS)+FNA(BCPOS(IB)+1,IY,0,IS)
          ENDDO
        END DO
        WRITE (NPRT,'(1X,3(A,1PD12.4))') ' SFEE = ',SFEE,' SFEI = ',
     .         SFEI,' SFNI = ',SUM(SFNI(IB,0:NS-1))
        ENDIF

        IF (BCCHAR(IB) .EQ. 'E') THEN
        WRITE (NPRT, 116) IB, BCPOS(IB)
        DO IY = BCEND(IB),BCSTART(IB), -1
          ix = BCPOS(IB)
          sxhelp=max(gs(ix,iy,0),1D-20)
          WRITE (NPRT, '(/I4,8X, 4(1PE11.3), 2X, (T59, 6D11.3))')
     .      IY, TE(BCPOS(IB),IY)/EV, TI(BCPOS(IB),IY)/EV,
     .      NE(BCPOS(IB),IY),
     .      PZ(BCPOS(IB),IY), (NA(BCPOS(IB),IY,IS), IS = 0, NS-1)
          WRITE (NPRT, '(1PE12.3, 4D11.3, 2X, (T59, 6D11.3))')
     .      GS(IX,IY,0), FHE(IX,IY,0)/SXhelp,
     .      FHI(IX,IY,0)/SXhelp,
     .      SUM(ZN(0:NS-1)*FNA(IX,IY,0,0:NS-1)) /
     .         SXhelp, PO(BCPOS(IB),IY),
     .      (FNA(IX,IY,0,IS)/SXhelp, IS = 0, NS-1)
          RZO = 0.0_R8
          DO IS = 0, NS-1
            IF (IS_NEUTRAL(IS)) CYCLE
            RZO = RZO + AM(IS)*MP*NA(BCPOS(IB),IY,IS)
          END DO
          CS = SQRT(PZ(BCPOS(IB),IY) / RZO)
          MACH=UA(IX,IY,1)/CS
          WRITE (NPRT, '((1X, 5ES11.3, T70,6D11.3))')
     .    FHE(IX,IY,0)/(FNA(IX,IY,0,1)+1.E-30)/TE(BCPOS(IB),IY),
     .    FHI(IX,IY,0)/(FNA(IX,IY,0,1)+1.E-30)/TI(BCPOS(IB),IY),
     .    MACH,
     .    (UA(IX,IY,IS),IS=0,NS-1)
          WRITE (NPRT, '(1PE12.3, 3D11.3, 2X, (T59, 6D11.3))')
     .      GS(IX,IY,0),
     .      FHE(IX,IY,0), FHI(IX,IY,0),
     .      SUM(ZN(0:NS-1)*FNA(IX,IY,0,0:NS-1)),
     .      (FNA(IX,IY,0,IS), IS = 0, NS-1)
          SFEE=SFEE-FHE(IX,IY,0)
          SFEI=SFEI-FHI(IX,IY,0)
          DO IS = 0, NS-1
            if (.not.is_neutral(is))
     .        SFNI(IB,IS)=SFNI(IB,IS)-FNA(IX,IY,0,IS)
          ENDDO
        END DO
        WRITE (NPRT,'(1X,3(A,1PD12.4))') ' SFEE = ',SFEE,' SFEI = ',
     .         SFEI,' SFNI = ',SUM(SFNI(IB,0:NS-1))
        ENDIF
      ENDDO

      call b2prnti(NPRT, nx, ny, ns)

* ..return
      ncall = ncall+1
      call subend ()

      RETURN

*=======================================================================
*//FORMATS//

  101 FORMAT (/2X,A)

  102 FORMAT (/'  THE GRID CONTAINS NX * NY INTERIOR CELLS, ',
     .      'WHERE NX =', I3, ' AND NY =', I3, '.' /
     .      '  WITH BOUNDARY CELLS INCLUDED THE SIZE IS ',
     .      '(NX+2) * (NY+2).' /
     .      '  THE NUMBER OF ION SPECIES (NS) IS ', I2, '.')

  103 FORMAT (/'  THE GRID METRIC PARAMETERS REDUCED TO ONE ',
     .      'DIMENSION ARE AS FOLLOWS:' //
     .      '  A) FUNCTION OF IX; AVERAGED OR SUMMED OVER IY.' //
     .      T4, 'IX', T11, 'X(IX)', T21, 'VOL(IX)', T32, 'DX(IX)',
     .      T43, 'DY(IX)', T54, 'DZ(IX)', T65, 'PIT(IX)', T76,
     .      'SX(IX)', T88, 'X(IX+1/2)', T98, 'SX(IX+1/2)' /)

  104 FORMAT (/'  B) FUNCTION OF IY; AVERAGED OR SUMMED OVER IX.' //
     .      T4, 'IY', T11, 'Y(IY)', T21, 'VOL(IY)', T32, 'DX(IY)',
     .      T43, 'DY(IY)', T54, 'DZ(IY)', T65, 'PIT(IY)', T76,
     .      'SY(IY)', T88, 'Y(IY+1/2)', T98, 'SY(IY+1/2)' /)

  110 FORMAT (//'  RESEE/EI/CO :', 2X, 1P, 2D9.1, (T38, 10D9.1))

  111 FORMAT (//'  CROSS-SECTIONAL FLUXES FHE, FHI, FNI.' //
     .      '  1. INTEGRATED OVER Y AS FUNCTION OF X:' //
     .      T8, 'X+', T19, 'SX', T31, 'FHEX', T42, 'FHIX',
     .      T53, 'FCHX', T64, 'FNEX', T75, 'FNIX',
     .      T86, 'FNIX(1)', T97, 'FNIX(..)' /)

  112 FORMAT (/'  2. INTEGRATED OVER X AS FUNCTION OF Y:' //
     .      T8, 'Y+', T19, 'SY', T31, 'FHEY', T42, 'FHIY',
     .      T53, 'FCHY', T64, 'FNEY', T75, 'FNIY',
     .      T86, 'FNIY(1)', T97, 'FNIY(..)' /)

 1112 FORMAT (//'  STATE OF THE PLASMA ON THE BOUNDARIES,',
     .      ' AND PERPENDICULAR FLUX DENSITIES.' // )
  113 FORMAT (/' SOUTH EDGE:',I3, '   AT POSITION IY:', I3 //
     .      T8, 'IX', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'NI(1)', T74, 'NI(..)' /
     .      T8, 'SY', T16, 'FEEY/SY', T27, 'FEIY/SY', T38, 'FNEY/SY',
     .      T51, 'PHI', T62, 'FNIY/SY', T73, 'FNIY/SY'/
     .      T5, 'SINANGY',
     .      T16, 'DELTAE', T27, 'DELTAI', T39, 'VV', T52, 'VO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SY/SQ(1-SINYQ)',
     .      T16, 'FEEY', T27, 'FEIY', T38, 'FNEY',
     .      T62, 'FNIY', T73,  'FNIY')

  114 FORMAT (/' NORTH EDGE:',I3, '   AT POSITION IY:', I3 //
     .      T8, 'IX', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'NI(1)', T74, 'NI(..)' /
     .      T8, 'SY', T16, 'FEEY/SY', T27, 'FEIY/SY', T38, 'FNEY/SY',
     .      T51, 'PHI', T62, 'FNIY/SY', T73, 'FNIY/SY'/
     .      T5, 'SINANGY',
     .      T16, 'DELTAE', T27, 'DELTAI', T39, 'VV', T52, 'VO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SY/SQ(1-SINYQ)',
     .      T16, 'FEEY', T27, 'FEIY', T38, 'FNEY',
     .      T62, 'FNIY', T73,  'FNIY')

  115 FORMAT (/'  WEST EDGE:',I3, '   AT POSITION IX:', I3 //
     .      T8, 'IY', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'NI(1)', T74, 'NI(..)' /
     .      T8, 'SX', T16, 'FEEX/SX', T27, 'FEIX/SX', T38, 'FNEX/SX',
     .      T51, 'PHI', T62, 'FNIX/SX', T73, 'FNIX/SX'/
     .      T5, 'SINANGX',
     .      T16, 'DELTAE', T27, 'DELTAI', T38, 'MACH', T52 ,'UO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SX/SQ(1-SINXQ)',
     .      T16, 'FEEX', T27, 'FEIX', T38, 'FNEX',
     .      T62, 'FNIX', T73,  'FNIX')

  116 FORMAT (/'  EAST EDGE:',I3, '   AT POSITION IX:', I3 //
     .      T8, 'IY', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'NI(1)', T74, 'NI(..)' /
     .      T8, 'SX', T16, 'FEEX/SX', T27, 'FEIX/SX', T38, 'FNEX/SX',
     .      T51, 'PHI', T62, 'FNIX/SX', T73, 'FNIX/SX'/
     .      T5, 'SINANGX',
     .      T16, 'DELTAE', T27, 'DELTAI', T38, 'MACH', T52, 'UO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SX/SQ(1-SINXQ)',
     .      T16, 'FEEX', T27, 'FEIX', T38, 'FNEX',
     .      T62, 'FNIX', T73,  'FNIX')

 1001 FORMAT ('  TIME (SECONDS):', 2X, 1D11.4)

 1113 FORMAT (/' SOUTH EDGE:',I3, '   AT POSITION IY:', I3 //
     .      T8, 'IX', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'UP', T74, 'UD', T85, 'NI(1)', T96, 'NI(..)'/
     .      T8, 'SY', T16, 'FEEY/SY', T27, 'FEIY/SY', T38, 'FNEY/SY',
     .      T51, 'PHI', T63, 'DD', T74, 'CD',
     .      T82, 'FNIY/SY', T93, 'FNIY/SY'/
     .      T5, 'SINANGY',
     .      T16, 'DELTAE', T27, 'DELTAI', T39, 'VV', T52 ,'VO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SY/SQ(1-SINYQ)',
     .      T16, 'FEEY', T27, 'FEIY', T38, 'FNEY',
     .      T62, 'FNIY', T73,  'FNIY')

 1114 FORMAT (/' NORTH EDGE:',I3, '   AT POSITION IY:', I3 //
     .      T8, 'IX', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'UP', T74, 'UD', T85, 'NI(1)', T96, 'NI(..)' /
     .      T8, 'SY', T16, 'FEEY/SY', T27, 'FEIY/SY', T38, 'FNEY/SY',
     .      T51, 'PHI', T63, 'DD', T74, 'CD',
     .      T82, 'FNIY/SY', T93, 'FNIY/SY'/
     .      T5, 'SINANGY',
     .      T16, 'DELTAE', T27, 'DELTAI', T39, 'VV', T52,'VO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SY/SQ(1-SINYQ)',
     .      T16, 'FEEY', T27, 'FEIY', T38, 'FNEY',
     .      T62, 'FNIY', T73,  'FNIY')

 1115 FORMAT (/'  WEST EDGE:',I3, '   AT POSITION IX:', I3 //
     .      T8, 'IY', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'CC', T74, 'CP', T85, 'NI(1)', T96, 'NI(..)' /
     .      T8, 'SX', T16, 'FEEX/SX', T27, 'FEIX/SX', T38, 'FNEX/SX',
     .      T51, 'PHI', T63, 'UD', T74, 'CD',
     .      T82, 'FNIX/SX', T93, 'FNIX/SX'/
     .      T5, 'SINANGX',
     .      T16, 'DELTAE', T27, 'DELTAI', T38, 'MACH', T52,'UO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SX/SQ(1-SINXQ)',
     .      T16, 'FEEX', T27, 'FEIX', T38, 'FNEX',
     .      T62, 'FNIX', T73,  'FNIX')

 1116 FORMAT (/'  EAST EDGE:',I3, '   AT POSITION :', I3 //
     .      T8, 'IY', T19, 'TE', T30, 'TI', T41, 'NE', T52, 'PZ',
     .      T63, 'CC', T74, 'CP', T85, 'NI(1)', T96, 'NI(..)' /
     .      T8, 'SX', T16, 'FEEX/SX', T27, 'FEIX/SX', T38, 'FNEX/SX',
     .      T51, 'PHI', T63, 'UD', T74, 'CD',
     .      T82, 'FNIX/SX', T93, 'FNIX/SX'/
     .      T5, 'SINANGX',
     .      T16, 'DELTAE', T27, 'DELTAI', T38, 'MACH', T52,'UO',
     .      T63, 'UU(1)', T74,  'UU(..)' /
     .      T3, 'SX/SQ(1-SINXQ)',
     .      T16, 'FEEX', T27, 'FEIX', T38, 'FNEX',
     .      T62, 'FNIX', T73,  'FNIX')


*=======================================================================
*//END B2PRNT//

      END
*//B2PRNT//

*//B2PRNTI//
*=======================================================================
*                 S U B R O U T I N E  B 2 P R N T I
*=======================================================================
*
*
*  1. PURPOSE
*
*     B2PRNTI  PRINTS OUT THE RESIDUALS FORMATTED TO INTEGER FORMAT
*RHW-version
*
*=======================================================================
*//DECLARATIONS//

      SUBROUTINE B2PRNTI(NPRT, nx, ny, ns)

      use b2mod_types
      use b2mod_geo
      use b2mod_residuals
      implicit none

      integer, intent(in) :: NPRT, nx, ny, ns
      real (kind=R8) :: reshlp(-1:nx,-1:ny)
      integer :: iz (-1:nx,-1:ny)
      INTEGER :: IX,IY,IS
      EXTERNAL IPRI

      DO IS = 0, NS-1
       DO IX = -1, NX
        DO IY = -1, NY
         RESHLP(IX,IY)=RESCO(IX,IY,IS)
         RESHLP(-1,IY)=0.0
         RESHLP(NX,IY)=0.0
         RESHLP(IX,-1)=0.0
         RESHLP(IX,NY)=0.0
        END DO
       END DO
       CALL IPRI(NPRT,RESHLP(-1,-1),IZ(-1,-1),NX,NY,'RESCO   ')
      END DO

      DO IS = 0, NS-1
       DO IX = -1, NX
        DO IY = -1, NY
         RESHLP(IX,IY)=RESMO(IX,IY,IS)
         RESHLP(-1,IY)=0.0
         RESHLP(NX,IY)=0.0
         RESHLP(NX,IY)=0.0
         RESHLP(IX,-1)=0.0
         RESHLP(IX,NY)=0.0
        END DO
       END DO
       CALL IPRI(NPRT,RESHLP(-1,-1),IZ(-1,-1),NX,NY,'RESMO   ')
      END DO

      DO IX = -1, NX
       DO IY = -1, NY
        RESHLP(IX,IY)=RESPO(IX,IY)
        RESHLP(-1,IY)=0.0
        RESHLP(NX,IY)=0.0
        RESHLP(IX,-1)=0.0
        RESHLP(IX,NY)=0.0
       END DO
      END DO
      CALL IPRI(NPRT,RESHLP(-1,-1),IZ(-1,-1),NX,NY,'RESPO   ')

      DO IX = -1, NX
       DO IY = -1, NY
        RESHLP(IX,IY)=RESHE(IX,IY)
        RESHLP(-1,IY)=0.0
        RESHLP(NX,IY)=0.0
        RESHLP(IX,-1)=0.0
        RESHLP(IX,NY)=0.0
       END DO
      END DO
      CALL IPRI(NPRT,RESHLP(-1,-1),IZ(-1,-1),NX,NY,'RESHE   ')

      DO IX = -1, NX
       DO IY = -1, NY
        RESHLP(IX,IY)=RESHI(IX,IY)
        RESHLP(-1,IY)=0.0
        RESHLP(NX,IY)=0.0
        RESHLP(IX,-1)=0.0
        RESHLP(IX,NY)=0.0
       END DO
      END DO
      CALL IPRI(NPRT,RESHLP(-1,-1),IZ(-1,-1),NX,NY,'RESHI   ')

*-----------------------------------------------------------------------

      RETURN

*=======================================================================
*//END B2PRNTI//

      END
CRHW:BREIASLB.IPRI   AUSDRUCK VON 2-DIM. ARRAYS ALS INTEGER-ZAHLEN
C                    DIESE ROUTINE IST AN DIE BRAAMS-ROUTINEN ANGEPASST.
C                    DIE X-KOORDINATE (1. INDEX) LAEUFT SENKRECHT
C                    DIE Y-KOORDINATE (2. INDEX) LAEUFT WAAGRECHT
C
C
C    BEIM AUFRUF MUESSEN DEFINIERT SEIN:
C      Z    REAL*8 (0:NXD+1,0:NYD+1) ...WERTE DES ARRAYS
C      IZ   INT *4 (0:NXD+1,0:NYD+1) ...MUSS NUR DIMENSIONIERT SEIN
C      NXD                           ...DIMENSION DES 1. INDEX VON Z
C      NYD                           ...DIMENSION DES 2. INDEX VON Z
C      NX                            ...ANZAHL DER X-WERTE (1. INDEX)
C      NY                            ...ANZAHL DER Y-WERTE (2. INDEX)
C      TEXT CHAR*8                   ...TABELLENUEBERSCHRIFT
C
      SUBROUTINE IPRI (NPRT, Z, IZ, NX, NY, TEXT)
C
      USE B2MOD_TYPES
      IMPLICIT NONE
C
      INTEGER NX, NY, NPRT
      INTEGER IZ (-1:NX,-1:NY), ID (-1:39),
     .        IX, IY, NX1, NX2, NY1, NY2
      REAL(KIND=R8) Z(-1:NX,-1:NY), YMIN, YMAX, ZMAX, SIGF
      CHARACTER*8 TEXT, COTXT1, COTXT2
C
      DATA  COTXT1 /'TE      '/,  COTXT2 /'TI      '/
C
C
      YMAX = 0.0
      YMIN = 0.0
      ZMAX = 0.0
C
      DO IX=-1,NX
       DO IY=-1,NY
         YMAX = MAX (YMAX,Z(IX,IY))
         YMIN = MIN (YMIN,Z(IX,IY))
         ZMAX = MAX (ZMAX,ABS(Z(IX,IY)))
       END DO
      END DO
C
      DO 12 IX=-1,NX
CCCCC    JX(IX) = IX
C
        IF (ZMAX.GT.0.0)   THEN
           DO 13 IY=-1,NY
             SIGF = 1.0
             IF (Z(IX,IY).LT.0.0)   SIGF = -1.0
             IZ(IX,IY) = 100.0 * Z(IX,IY)/ZMAX + 0.5*SIGF
   13      CONTINUE
C
        ELSE
C
           DO 14 IY=-1,NY
             IZ(IX,IY) = 0
   14      CONTINUE
        ENDIF
C
   12 CONTINUE
C
      ZMAX = YMAX
      IF (ABS(YMIN).GT.YMAX)   ZMAX = YMIN
C
      IF (TEXT.EQ.COTXT1 .OR. TEXT.EQ.COTXT2)   ZMAX = ZMAX/1.6021D-19
C
      WRITE (NPRT,200) TEXT, ZMAX
  200 FORMAT (//1x,'   AUSDRUCK DES ARRAYS  >>>>>  ',A8,
     .        '  <<<<   ZMAX =',1PE13.3)
C
      NY2 = -2
C
   20 CONTINUE
      NX2 = -2
      NY1 = NY2+1
      NY2 = NY1+39
      NY2 = MIN0(NY2,NY)
C
      DO 21 IY=NY1,NY2
         ID(IY-NY1) = IY
   21 CONTINUE
C
      WRITE (NPRT,203) (ID(IY-NY1),IY=NY1,NY2)
      WRITE (NPRT,204)
C
   22 CONTINUE
      NX1 = NX2+1
      NX2 = NX1+48
      NX2 = MIN0(NX2,NX)
C
      DO 23 IX=NX1,NX2
      WRITE (NPRT,205) IX,(IZ(IX,IY),IY=NY1,NY2)
C==>> IF(IX.LT.NX2 .AND. MOD(IX,10).EQ.0)   WRITE(NPRT,202)
   23 CONTINUE
C
      IF(NX2.LT.NX)  GO TO 22
C
      IF(NY2.LT.NY)  GO TO 20
      RETURN
C
  202 FORMAT (6X,'|')
  203 FORMAT (' IY = | ',40I4)
  204 FORMAT ('  IX  |',/,1X,172('-'),/,'  IX  |')
  205 FORMAT (1X,I4,' | ',40I4)
CC201 FORMAT (/,' IX =',34I4,/,' IY',/)
C
      END

