
C =====================================================================
      SUBROUTINE PLOTCONTOUR(NCV,IS,
     . CR,CZ,FIELD,ZLOD,ZHID,PHYSICAL,LINLOG)
C =====================================================================
      use b2mod_types
      use b2mod_b2plot_values
      use b2mod_b2plot_debug
      use b2mod_b2plot_xy
      use b2mod_boundary_namelist
      use b2mod_subsys
      use b2us_data
      implicit none
#include "graphics_vars.h"
      INTEGER NCV,IS,LINLOG
cxpb  INTEGER NNCUT_L
      REAL (KIND=R8) :: CR(NCV),CZ(NCV),FIELD(NCV),ZLOD,ZHID
      REAL (KIND=R8), ALLOCATABLE :: FVX(:)
      CHARACTER*20 FMTZ
      CHARACTER*256 STRING
      CHARACTER*20 DEBUGSTRING
      LOGICAL PHYSICAL,DEBUG,LBLBR

      INTEGER LIWK,LRWK,LMAP,NOGRPS,NWRK
      PARAMETER (NWRK=50000)
      PARAMETER (LIWK=50000,LRWK=50000,LMAP=5000000,NOGRPS=6)
      REAL (KIND=R4) :: XWRK(NWRK),YWRK(NWRK),RWRK(LRWK)
      INTEGER IWRK(LIWK),MAP(LMAP),IAREA(NOGRPS),IGRP(NOGRPS)
      integer irold,nerro,nerrf
      INTEGER LOPN,LOEN,LOTN
      INTEGER NOPN,NOEN,NOTN
      PARAMETER (LOPN=4)
      PARAMETER (LOEN=5)
      PARAMETER (LOTN=4)
      INTEGER NPNT,NEDG,NTRI
      REAL (KIND=R4), ALLOCATABLE :: RPNT(:)
      INTEGER, ALLOCATABLE :: nodelist(:,:), edgelist(:,:),
     . IEDG(:),ITRI(:),IBEG(:),IBPT(:)
c$$$      INTEGER IPTP,IPTE,IPTT
      INTEGER LFIN(100),choose_color
      CHARACTER*20 LBLS(100)
      logical, allocatable :: done(:)

      integer ix,iy,ixstrt,ixend,iystrt,iyend
      integer nfmtz,noz,nozm,
     .        i,j,k,ib,ic,ie,it,je,jt,
     .        iCv,iCn,iFc,iVx,iVx1,iVx2,iVa,iVb
      real(kind=R4) :: rlo, rhi, zlo, zhi, zinc, zval, clv, special
      real(kind=R4) :: ctval(1000), rmax, rmin, zmax, zmin
      logical right_match_found, left_match_found
      integer no_available_colors
      real (kind=R4) :: ceil, sfloor
      external no_available_colors, ceil, sfloor, xertst
      intrinsic mod

      EXTERNAL FILL,CPDRPL

      logical first
      save first
      data first /.true./

      INTEGER LNGDEC
      PARAMETER(LNGDEC=200*200)
      REAL (KIND=R4) :: Z(LNGDEC)

      call subini ('plotcontour')
      if(first) then
        call alloc_b2mod_b2plot_xy(LNGDEC)
        first=.false.
      endif
      rwrk=0.0_R4   ! bug in NCAR: The variable (RWRK(209)) has an undefined value.
      z=0.0_R4

C ========== GLOBAL OPTIONS ==========================================
#ifdef DBG
      DEBUG=.true.
#else
      DEBUG=.false.
#endif
      LBLBR=.true.
C ====================================================================

      ZLO=real(ZLOD,R4)
      ZHI=real(ZHID,R4)
      rlo=real(field(1),R4)
      rhi=real(field(1),R4)
      DO iCv=2,nCv
        RLO=MIN(RLO,REAL(FIELD(iCv),R4))
        RHI=MAX(RHI,REAL(FIELD(iCv),R4))
      ENDDO
      IF(zlo.eq.zhi) then
        zlo=rlo
        zhi=rhi
      ENDIF

      if(zlo.eq.zhi) then
        zlo=min(zlo*0.9_R4,zlo*1.1_R4)
        zhi=max(zhi*0.9_R4,zhi*1.1_R4)
        if(zlo.eq.zhi) then
          zlo=-1.0_R4
          zhi=+1.0_R4
        endif
      endif

      WRITE(*,*) 'zlo =',ZLO,' zhi =',ZHI

C =========== Setting up Contour Levels and Labeling =================
      IF(mod(LINLOG,2).EQ.1) THEN   ! ==> do linear plot <====================
        IF(NOFCLS.EQ.0) THEN ! standard contouring
          write (*,*) 'STANDARD CONTOURING'
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ,.true.)
          NCONT=NOZ*NOZM
          if(b2plot_debug.ge.10)
     1     write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          do while(NCONT+3.GT.no_available_colors(0).AND.NOZM.GT.1)
            NOZM=NOZM/2
            NCONT=NOZ*NOZM
            if(b2plot_debug.ge.10)
     1       write (*,*) 'ncont =',ncont,' noz =',noz,' nozm =',nozm
          ENDDO
          ZINC=(ZHI-ZLO)/(NCONT)
          if(b2plot_debug.ge.10)
     1       write (*,*) 'zlo =',zlo,' zhi =',zhi,' zinc =',zinc
          DO I=1,NCONT+1
            ctval(i)=zlo+(i-1)*zinc
            write (*,*) 'contour level',i,
     1             '=',zlo+(i-1)*zinc,'in color',i
          ENDDO
        ELSE ! customised contouring
          write (*,*) 'CUSTOMISED CONTOURING'
          NCONT=NOFCLS-1
          if(ncont+3.gt.no_available_colors(0)) then
C do something to avoid this
          endif
          NOZM=NCONT/5
          if(b2plot_debug.ge.3)
     1     write (*,*) 'ncont =',ncont
          do i=1,ncont+1
            write (*,*) 'clval(',i,') =',clval(i)
          enddo
          CALL NICE(ZLO,ZHI,NOZ,NOZM,FMTZ,NFMTZ,.true.)
          ZINC=(ZHI-ZLO)/(NCONT)
          if(b2plot_debug.ge.10)
     1       write (*,*) 'zlo =',zlo,' zhi =',zhi,' zinc =',zinc
          DO I=1,NCONT+1
            ctval(i)=clval(i)
            write (*,*) 'contour level',i,'=',clval(i),'in color',i
          ENDDO
        ENDIF
        WRITE(STRING,'(a,3(1pg12.5,a))')
     1   'contoured from',ZLO,'by',ZINC,'to',ZHI
      ELSE  ! ==> do log plot <=======================================
        IF(NOFCLS.EQ.0) THEN ! standard contouring
          write (*,*) 'STANDARD CONTOURING'
          if(zhi.le.0.0_R4) then
            write(*,*) 'zhi <= 0 and logarithmic axes chosen'
            zhi=1.0_R4
          endif
          if(zlo.le.0.0_R4) zlo=10.0_R4**(-ndec)*zhi
          ZHI=10.0_R4**CEIL(LOG10(ZHI))
          ZLO=max(10.0_R4**SFLOOR(LOG10(ZLO)),10.0_R4**(-ndec)*zhi)
          write (*,*) 'zlo =',zlo,' zhi =',zhi
          NCONT=INT(CEIL(LOG10(ZHI))-SFLOOR(LOG10(ZLO)))
          FMTZ='(1PG7.1)'
          if(ncont*3+3.lt.no_available_colors(0)) then
            ncont=ncont*3
            do i=1,ncont+1
              zval=10.0_R4**((i-1)/3)*zlo
              if(mod(i-1,3).eq.1) zval=zval*2.0_R4
              if(mod(i-1,3).eq.2) zval=zval*5.0_R4
              ctval(i)=zval
            enddo
          else
            do i=1,ncont+1
              zval=10.0_R4**(i-1)*zlo
              ctval(i)=zval
            enddo
          endif
        ELSE ! customised contouring
          write (*,*) 'CUSTOMISED CONTOURING'
          ZLO=10.0_R4**SFLOOR(LOG10(ZLO))
          ZHI=10.0_R4**CEIL(LOG10(ZHI))
          NCONT=NOFCLS-1
          FMTZ='(1PG9.3)'
          DO I=1,NCONT+1
            ctval(i)=clval(i)
            write (*,*) 'contour level',i,'=',clval(i),' in color',i
          ENDDO
        ENDIF
        WRITE(STRING,'(a,2(1pg12.5,a))')
     1   'contoured from',ZLO,'to',ZHI
      ENDIF ! ==> end of lin/log plot <===============================

      IF(PHYSICAL) THEN  ! ========== plot physical domain

C Build triangle grid
C First, set the field data
        allocate(FVX(mpg%nVx))
        call intvertex(mpg%nCv,mpg%nVx,mpg,geo%vxVol,field,FVX)
        allocate(nodelist(1:mpg%nCi,1:maxval(mpg%cvVxP(:,2))))
        allocate(edgelist(1:mpg%nCi,1:maxval(mpg%cvFcP(:,2))))
        allocate(IBPT(1:mpg%nVx))
        allocate(RPNT(1:mpg%nVx*LOPN))
        i=0
        ib=0
        do iVx=1,mpg%nVx
          if(mpg%vxCvP(iVx,2).eq.0) cycle
          ib=ib+1
          IBPT(iVx)=ib
          i=i+1
          RPNT(i)=real(geo%vxX(iVx),R4)
          i=i+1
          RPNT(i)=real(geo%vxY(iVx),R4)
          i=i+1
          RPNT(i)=0.0_R4
          i=i+1
          RPNT(i)=real(FVX(iVx),R4)
        end do
        NOPN=ib
        NPNT=NOPN*LOPN
C Second, build the list of edges
        i=mpg%nFc
        j=mpg%nCi
        do iCv=1,mpg%nCi
          i=i+(mpg%cvFcP(iCv,2)-3)
          j=j+(mpg%cvVxP(iCv,2)-3)
        end do
        NOEN=i
        NEDG=NOEN*LOEN
        NOTN=j
        NTRI=NOTN*LOTN
        allocate(IEDG(NEDG))
        allocate(IBEG(NOEN))
        allocate(ITRI(NTRI))
        allocate(done(1:maxval(mpg%cvFcP(1:mpg%nCi,2))))
        do iCv=1,mpg%nCi
cxpb We assume convex control volumes
cxpb First we order the vertices and edges
          done=.false.
          iFc=mpg%cvFc(mpg%cvFcP(iCv,1))
          iVa=mpg%fcVx(iFc,1)
          iVb=mpg%fcVx(iFc,2)
          nodelist(iCv,1)=iVa
          nodelist(iCv,2)=iVb
          edgelist(iCv,1)=iFc
          ic=2
          ie=1
          done(1)=.true.
          do while (ic.lt.mpg%cvVxP(iCv,2))
            do k=2,mpg%cvFcP(iCv,2)
              if(ic.eq.mpg%cvVxP(iCv,2)) cycle
              if(done(k)) cycle
              right_match_found=.false.
              iFc=mpg%cvFc(mpg%cvFcP(iCv,1)+k-1)
              iVx1=mpg%fcVx(iFc,1)
              iVx2=mpg%fcVx(iFc,2)
              if (iVx1.eq.iVb) then
                right_match_found=.true.
                iVb=iVx2
              else if (iVx2.eq.iVb) then
                right_match_found=.true.
                iVb=iVx1
              end if
              if (right_match_found) then
                ic=ic+1
                nodelist(iCv,ic)=iVb
                ie=ie+1
                edgelist(iCv,ie)=iFc
                done(k)=.true.
              end if
            end do
          end do
          do k=2,mpg%cvFcP(iCv,2)
            if(done(k)) cycle
            iFc=mpg%cvFc(mpg%cvFcP(iCv,1)+k-1)
            ie=ie+1
            edgelist(iCv,ie)=iFc
          end do
        end do
        deallocate(done)
        allocate(done(3*j))
        done=.false.
        i=0
        ie=0
        je=mpg%nFc
cxpb Copy in the existing edges information
cxpb For now just assign the CV number as the triangle index
        do iCv=1,mpg%nCi
          do k=1,mpg%cvFcP(iCv,2)
            iFc=edgelist(iCv,k)
            if(done(iFc)) cycle
            ie=ie+1
            IBEG(iFc)=ie
            i=i+1
            IEDG(i)=(IBPT(nodelist(iCv,k))-1)*LOPN
            i=i+1
            if (k.lt.mpg%cvVxP(iCv,2)) then
              IEDG(i)=(IBPT(nodelist(iCv,k+1))-1)*LOPN
            else
              IEDG(i)=(IBPT(nodelist(iCv,1))-1)*LOPN
            end if
            if (k.lt.3) then
              j=k
            else if (k.eq.mpg%cvFcP(iCv,2)) then
              j=3
            else
              j=2
            end if
            if (mpg%fcVx(iFc,1).eq.nodelist(iCv,k)) then
              i=i+1
              if (mpg%fcCv(iFc,1).gt.mpg%nCi) then
                IEDG(i)=-1
              else
                IEDG(i)=(mpg%fcCv(iFc,1)-1)*LOTN+j
              end if
              i=i+1
              if (mpg%fcCv(iFc,2).gt.mpg%nCi) then
                IEDG(i)=-1
              else
                IEDG(i)=(mpg%fcCv(iFc,2)-1)*LOTN+j
              end if
            else
              i=i+1
              if (mpg%fcCv(iFc,2).gt.mpg%nCi) then
                IEDG(i)=-1
              else
                IEDG(i)=(mpg%fcCv(iFc,2)-1)*LOTN+j
              end if
              i=i+1
              if (mpg%fcCv(iFc,1).gt.mpg%nCi) then
                IEDG(i)=-1
              else
                IEDG(i)=(mpg%fcCv(iFc,1)-1)*LOTN+j
              end if
            end if
            i=i+1
            IEDG(i)=0
            done(iFc)=.true.
          end do
cxpb Now add any additional diagonals
          do k=4,mpg%cvFcP(iCv,2)
            je=je+1
            ie=ie+1
            IBEG(je)=ie
            i=i+1
            IEDG(i)=(IBPT(nodelist(iCv,1))-1)*LOPN
            i=i+1
            IEDG(i)=(IBPT(nodelist(iCv,k-1))-1)*LOPN
            i=i+1
            IEDG(i)=(iCv-1)*LOTN+3
            i=i+1
            IEDG(i)=(iCv-1)*LOTN+1
            i=i+1
            IEDG(i)=0
          end do
        end do
cxpb Now assign the triangles
        it=mpg%nCi
        ie=mpg%nFc
        do iCv=1,mpg%nCi
          do k=3,mpg%cvFcP(iCv,2)
            if(k.eq.3) then
              jt=iCv
            else
              it=it+1
              jt=it
              ie=ie+1
            end if
            ib=(jt-1)*LOTN
            if (k.eq.3) then
              ITRI(ib+1)=(IBEG(edgelist(iCv,1))-1)*LOEN
              ITRI(ib+2)=(IBEG(edgelist(iCv,k-1))-1)*LOEN
              if (k.eq.mpg%cvFcP(iCv,2)) then
                ITRI(ib+3)=(IBEG(edgelist(iCv,k))-1)*LOEN
              else
                ITRI(ib+3)=(IBEG(ie+1)-1)*LOEN
              end if
            else
              ITRI(ib+1)=(IBEG(ie)-1)*LOEN
              IEDG((IBEG(ie)-1)*LOEN+4)=ib+1
              ITRI(ib+2)=(IBEG(edgelist(iCv,k-1))-1)*LOEN
              if ((IBPT(nodelist(iCv,k-1))-1)*LOPN.eq.
     .            IEDG((IBEG(edgelist(iCv,k-1))-1)*LOEN+1)) then
                IEDG((IBEG(edgelist(iCv,k-1))-1)*LOEN+4)=ib+2
              else
                IEDG((IBEG(edgelist(iCv,k-1))-1)*LOEN+3)=ib+2
              end if
              if (k.lt.mpg%cvFcP(iCv,2)) then
                ITRI(ib+3)=(IBEG(ie+1)-1)*LOEN
                IEDG((IBEG(ie+1)-1)*LOEN+4)=ib+3
              else
                ITRI(ib+3)=(IBEG(edgelist(iCv,k))-1)*LOEN
                if ((IBPT(nodelist(iCv,k))-1)*LOPN.eq.
     .            IEDG((IBEG(edgelist(iCv,k))-1)*LOEN+1)) then
                  IEDG((IBEG(edgelist(iCv,k))-1)*LOEN+4)=ib+3
                else
                  IEDG((IBEG(edgelist(iCv,k))-1)*LOEN+3)=ib+3
                end if
              endif
            end if
            ITRI(ib+4)=0
            rmin = MIN(RPNT(IEDG(ITRI(ib+1)+1)+1),
     .                 RPNT(IEDG(ITRI(ib+1)+2)+1),
     .                 RPNT(IEDG(ITRI(ib+2)+1)+1),
     .                 RPNT(IEDG(ITRI(ib+2)+2)+1),
     .                 RPNT(IEDG(ITRI(ib+3)+1)+1),
     .                 RPNT(IEDG(ITRI(ib+3)+1)+1),
     .                 RPNT(IEDG(ITRI(ib+3)+2)+1))
            rmax = MAX(RPNT(IEDG(ITRI(ib+1)+1)+1),
     .                 RPNT(IEDG(ITRI(ib+1)+2)+1),
     .                 RPNT(IEDG(ITRI(ib+2)+1)+1),
     .                 RPNT(IEDG(ITRI(ib+2)+2)+1),
     .                 RPNT(IEDG(ITRI(ib+3)+1)+1),
     .                 RPNT(IEDG(ITRI(ib+3)+2)+1))
            zmin = MIN(RPNT(IEDG(ITRI(ib+1)+1)+2),
     .                 RPNT(IEDG(ITRI(ib+1)+2)+2),
     .                 RPNT(IEDG(ITRI(ib+2)+1)+2),
     .                 RPNT(IEDG(ITRI(ib+2)+2)+2),
     .                 RPNT(IEDG(ITRI(ib+3)+1)+2),
     .                 RPNT(IEDG(ITRI(ib+3)+2)+2))
            zmax = MAX(RPNT(IEDG(ITRI(ib+1)+1)+2),
     .                 RPNT(IEDG(ITRI(ib+1)+2)+2),
     .                 RPNT(IEDG(ITRI(ib+2)+1)+2),
     .                 RPNT(IEDG(ITRI(ib+2)+2)+2),
     .                 RPNT(IEDG(ITRI(ib+3)+1)+2),
     .                 RPNT(IEDG(ITRI(ib+3)+2)+2))
            IF (rmin.gt.rhi.or.rmax.gt.rlo.or.
     .          zmin.gt.zhi.or.zmax.le.zlo) ITRI(ib+4)=1
          end do
        end do

c$$$C Plot mesh for debugging
c$$$        CALL GSPLCI (11)
c$$$        DO 101 IT=0,NTRI-LOTN,LOTN
c$$$          CALL LINED (RPNT(IEDG(ITRI(IT+1)+1)+1),
c$$$     +                RPNT(IEDG(ITRI(IT+1)+1)+2),
c$$$     +                RPNT(IEDG(ITRI(IT+1)+2)+1),
c$$$     +                RPNT(IEDG(ITRI(IT+1)+2)+2))
c$$$          CALL LINED (RPNT(IEDG(ITRI(IT+2)+1)+1),
c$$$     +                RPNT(IEDG(ITRI(IT+2)+1)+2),
c$$$     +                RPNT(IEDG(ITRI(IT+2)+2)+1),
c$$$     +                RPNT(IEDG(ITRI(IT+2)+2)+2))
c$$$          CALL LINED (RPNT(IEDG(ITRI(IT+3)+1)+1),
c$$$     +                RPNT(IEDG(ITRI(IT+3)+1)+2),
c$$$     +                RPNT(IEDG(ITRI(IT+3)+2)+1),
c$$$     +                RPNT(IEDG(ITRI(IT+3)+2)+2))
c$$$  101   CONTINUE
c$$$C Label the points of the mesh, in aqua.
c$$$        CALL GSPLCI (7)
c$$$        CALL GSLWSC (2.)
c$$$C
c$$$        DO 114 IPTP=0,NPNT-LOPN,LOPN
c$$$          CALL LABLIT (IPTP/LOPN+1,.024,0.,RPNT(IPTP+1),RPNT(IPTP+2))
c$$$  114   CONTINUE
c$$$C Label the edges of the mesh, in deep sky blue.
c$$$        CALL GSPLCI (6)
c$$$        CALL GSLWSC (2.)
c$$$C
c$$$        DO 115 IPTE=0,NEDG-LOEN,LOEN
c$$$          CALL LABLIT (IPTE/LOEN+1,.012,
c$$$     +                 ATAN2(RPNT(IEDG(IPTE+2)+2)-RPNT(IEDG(IPTE+1)+2),
c$$$     +                       RPNT(IEDG(IPTE+2)+1)-RPNT(IEDG(IPTE+1)+1)),
c$$$     +                   .5*(RPNT(IEDG(IPTE+1)+1)+RPNT(IEDG(IPTE+2)+1)),
c$$$     +                   .5*(RPNT(IEDG(IPTE+1)+2)+RPNT(IEDG(IPTE+2)+2)))
c$$$C Label triangle number to the left with a negative value
c$$$          CALL LABLIT (-(IEDG(IPTE+3)/LOTN+1),.012,
c$$$     +                 ATAN2(RPNT(IEDG(IPTE+2)+2)-RPNT(IEDG(IPTE+1)+2),
c$$$     +                       RPNT(IEDG(IPTE+2)+1)-RPNT(IEDG(IPTE+1)+1)),
c$$$     +               .25*(3.*RPNT(IEDG(IPTE+1)+1)+RPNT(IEDG(IPTE+2)+1)),
c$$$     +               .25*(3.*RPNT(IEDG(IPTE+1)+2)+RPNT(IEDG(IPTE+2)+2)))
c$$$C Label triangle number to the right
c$$$          CALL LABLIT ( (IEDG(IPTE+4)/LOTN+1),.012,
c$$$     +                 ATAN2(RPNT(IEDG(IPTE+2)+2)-RPNT(IEDG(IPTE+1)+2),
c$$$     +                       RPNT(IEDG(IPTE+2)+1)-RPNT(IEDG(IPTE+1)+1)),
c$$$     +               .25*(RPNT(IEDG(IPTE+1)+1)+3.*RPNT(IEDG(IPTE+2)+1)),
c$$$     +               .25*(RPNT(IEDG(IPTE+1)+2)+3.*RPNT(IEDG(IPTE+2)+2)))
c$$$  115   CONTINUE
c$$$C Label the triangles of the mesh, in slate blue.
c$$$C Label any unused triangles in celeste.
c$$$C Write the edge numbers near them.
c$$$        CALL GSPLCI (4)
c$$$        CALL GSLWSC (2.)
c$$$C
c$$$        DO 116 IPTT=0,NTRI-LOTN,LOTN
c$$$          IF (ITRI(IPTT+4).EQ.0) THEN
c$$$            CALL GSPLCI (4)
c$$$          ELSE
c$$$            CALL GSPLCI (9)
c$$$          END IF
c$$$          CALL LABLIT (IPTT/LOTN+1,.018,0.,
c$$$     +                 (RPNT(IEDG(ITRI(IPTT+1)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+1)+2)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+2)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+2)+1))/6.,
c$$$     +                 (RPNT(IEDG(ITRI(IPTT+1)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+1)+2)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+2)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+2)+2))/6.)
c$$$          CALL LABLIT (ITRI(IPTT+1)/LOEN+1,.012,0.,
c$$$     +              (4.*RPNT(IEDG(ITRI(IPTT+1)+1)+1)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+1)+2)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+2)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+2)+1))/12.,
c$$$     +              (4.*RPNT(IEDG(ITRI(IPTT+1)+1)+2)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+1)+2)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+2)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+2)+2))/12.)
c$$$          CALL LABLIT (ITRI(IPTT+2)/LOEN+1,.012,0.,
c$$$     +                 (RPNT(IEDG(ITRI(IPTT+1)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+1)+2)+1)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+2)+1)+1)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+2)+2)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+2)+1))/12.,
c$$$     +                 (RPNT(IEDG(ITRI(IPTT+1)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+1)+2)+2)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+2)+1)+2)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+2)+2)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+3)+2)+2))/12.)
c$$$          CALL LABLIT (ITRI(IPTT+3)/LOEN+1,.012,0.,
c$$$     +                 (RPNT(IEDG(ITRI(IPTT+1)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+1)+2)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+1)+1)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+2)+1)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+3)+1)+1)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+3)+2)+1))/12.,
c$$$     +                 (RPNT(IEDG(ITRI(IPTT+1)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+1)+2)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+1)+2)+
c$$$     +                  RPNT(IEDG(ITRI(IPTT+2)+2)+2)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+3)+1)+2)+
c$$$     +               4.*RPNT(IEDG(ITRI(IPTT+3)+2)+2))/12.)
c$$$  116   CONTINUE
c$$$        CALL GSPLCI (1)
c$$$        CALL GSLWSC (1.)

C Initialize CONPAKT
        IF (CONTPLOT)
     >   CALL CTMESH (RPNT,NPNT,LOPN,
     ,                IEDG,NEDG,LOEN,
     ,                ITRI,NTRI,LOTN,
     ,                RWRK,LRWK,IWRK,LIWK)

C Initialise Areas
        CALL ARINAM(MAP, LMAP)
        CALL GSFAIS(1)

        IC=1
        X(IC)=XLO
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YLO
        if (b2plot_debug.ge.8) then
          write(*,*) 'Physical domain bounding box : '
          write(*,*) 'X : ', X(1:5)
          write(*,*) 'Y : ', Y(1:5)
        end if
        CALL AREDAM(MAP,X,Y,IC,2,-1,0)
cxpb Obtaining the contour outlines of the plasma region to be plotted
        deallocate(done)
        allocate(done(mpg%nCi+1:mpg%nCv))
        done=.false.
        do while (.not.all(done))
          do iCv = mpg%nCi+1, mpg%nCv
            if (done(iCv)) cycle
            iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            X(1)=real(geo%vxX(iVx1),R4)
            Y(1)=real(geo%vxY(iVx1),R4)
            X(2)=real(geo%vxX(iVx2),R4)
            Y(2)=real(geo%vxY(iVx2),R4)
            iVa = iVx1
            iVb = iVx2
            ic=2
            done(iCv)=.true.
            do while (.not.(iVa.eq.iVb.or.all(done)))
              do iCn = mpg%nCi+1,mpg%nCv
                if (done(iCn)) cycle
                left_match_found=.false.
                right_match_found=.false.
                iFc = mpg%cvFc(mpg%cvFcP(iCn,1))
                iVx1 = mpg%fcVx(iFc,1)
                iVx2 = mpg%fcVx(iFc,2)
                if (iVx1.eq.iVa) then
                  left_match_found=.true.
                  do i=ic,1,-1
                    X(i+1)=X(i)
                    Y(i+1)=Y(i)
                  end do
                  X(1)=real(geo%vxX(iVx2),R4)
                  Y(1)=real(geo%vxY(iVx2),R4)
                  iVa = iVx2
                  ic = ic+1
                else if (iVx1.eq.iVb) then
                  right_match_found=.true.
                  X(ic+1)=real(geo%vxX(iVx2),R4)
                  Y(ic+1)=real(geo%vxY(iVx2),R4)
                  iVb = iVx2
                  ic = ic+1
                else if (iVx2.eq.iVa) then
                  left_match_found=.true.
                  do i=ic,1,-1
                    X(i+1)=X(i)
                    Y(i+1)=Y(i)
                  end do
                  X(1)=real(geo%vxX(iVx1),R4)
                  Y(1)=real(geo%vxY(iVx1),R4)
                  iVa = iVx1
                  ic = ic+1
                else if (iVx2.eq.iVb) then
                  right_match_found=.true.
                  X(ic+1)=real(geo%vxX(iVx1),R4)
                  Y(ic+1)=real(geo%vxY(iVx1),R4)
                  iVb = iVx1
                  ic = ic+1
                end if
                done(iCn)=left_match_found.or.right_match_found
              end do
            end do
            CALL AREDAM(MAP,X,Y,IC,2,0,1)
            IF(.NOT.COLPRT.OR..NOT.CONTPLOT) CALL CURVED(X,Y,IC)
          end do
        end do
        deallocate(done)

        IF(LABELS) THEN
          CALL CTSETI('LLP - LINE LABEL POSITIONING',3)
          IF(LALONG) CALL CTSETI('LLO - LINE LABEL ORIENTATION',1)
        ELSE
          CALL CTSETI('LLP - LINE LABEL POSITIONING',0)
        ENDIF
        IF(INTPOL) THEN
          CALL CTSETR('T2D - TENSION OF THE 2D SPLINES',1.0)
          CALL CTSETI('PIC - POINT INTERPOLATION FLAG',-10)
          CALL IPOLATE(.FALSE.) ! for further plots
          WRITE (*,*) 'using spline interpolation'
        ELSE
          CALL CTSETR('T2D - TENSION OF THE 2D SPLINES',0.0)
          CALL CTSETI('PIC - POINT INTERPOLATION FLAG',-10)
        ENDIF
        CALL CTSETC('CFT - CONSTANT FIELD LABEL',' ')
        CALL CTSETI('CLS - CONTOUR LEVEL SELECTION FLAG',0)
        CALL CTSETI('NCL - NUMBER OF CONTOUR LEVELS',NCONT+1)
        DO I=1,NCONT+1
          CALL CTSETI('PAI - PARAMETER ARRAY INDEX',I)
          CALL CTSETR('CLV - CONTOUR LEVEL VALUE',CTVAL(I))
          CALL CTSETR('CLL - CONTOUR LINE LINE WIDTH',LWIDTH)
          CALL CTSETI('AIA - AREA ABOVE IDENTIFIER',I+1)
          CALL CTSETI('AIB - AREA BELOW IDENTIFIER',I)
          CALL CTSETI('CLU - CONTOUR LEVEL USE FLAG',3)
        END DO
        CALL CTSETC('ILT - INFORMATION LABEL TEXT STRING',
     1   TRIM(STRING))
        CALL CTSETR('ILX - INFORMATION LABEL X POSITION',0.5)
        CALL CTSETR('ILY - INFORMATION LABEL Y POSITION',0.95)
        CALL CTSETI('ILP - INFORMATION LABEL POSITION',-3)
        CALL CTSETI('GIC - GROUP IDENTIFIER FOR CONTOURLINES',1)
        CALL CTSETI('GIL - GROUP IDENTIFIER FOR LABELS',2)
        CALL CTSETI('GIS - GROUP IDENTIFIER FOR VSTRIPS',6)
        CALL CTSETI('NVS - NUMBER OF VSTRIPS',0)

        IF(CONTPLOT) THEN
C Put contour lines in the area map
          CALL CTCLAM(RPNT,IEDG,ITRI,RWRK,IWRK,MAP)
C Add label boxes to area map
          IF(LABELS) THEN
            CALL CTLBAM(RPNT,IEDG,ITRI,RWRK,IWRK,MAP)
          ENDIF
C Set fill style to solid, and fill contours if COLPRT==TRUE
          IF(COLPRT) THEN
            call entsr(irold,1)
            CALL ARPRAM(MAP,1,0,0)
            if(nerro(nerrf).ne.0) then
              call eprin
              call errof
              call retsr(irold)
              call subend ()
              return
            endif
            call retsr(irold)
            CALL ARSCAM(MAP,XWRK,YWRK,NWRK,IAREA,IGRP,NOGRPS,FILL)
          ENDIF
        ENDIF

        deallocate(RPNT)
        deallocate(IEDG,ITRI,IBEG,IBPT)
        deallocate(FVX)
        deallocate(nodelist,edgelist)

      ELSE IF (mpg%isClassicalGrid.eq.1) THEN   ! ======== plot computational domain
C Initialise Areas
        CALL ARINAM(MAP, LMAP)
        CALL GSFAIS(1)

        CALL CPSETI('MAP - MAPPING FLAG',4)
        IF(LABELS) THEN
          CALL CPSETI('LLP - LINE LABEL POSITIONING',3)
          IF(LALONG) CALL CPSETI('LLO - LINE LABEL ORIENTATION',1)
        ELSE
          CALL CPSETI('LLP - LINE LABEL POSITIONING',0)
        ENDIF
        IF(INTPOL) THEN
          CALL CPSETR('T2D - TENSION OF THE 2D SPLINES',1.0)
          CALL CPSETI('PIC - POINT INTERPOLATION FLAG',-10)
          CALL IPOLATE(.FALSE.) ! for further plots
          WRITE (*,*) 'using spline interpolation'
        ELSE
          CALL CPSETR('T2D - TENSION OF THE 2D SPLINES',0.0)
          CALL CPSETI('PIC - POINT INTERPOLATION FLAG',-10)
        ENDIF
        CALL CPSETC('CFT - CONSTANT FIELD LABEL',' ')
        CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',0)
        CALL CPSETI('NCL - NUMBER OF CONTOUR LEVELS',NCONT+1)
        DO I=1,NCONT+1
          CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
          CALL CPSETR('CLV - CONTOUR LEVEL VALUE',CTVAL(I))
          CALL CPSETR('CLL - CONTOUR LINE LINE WIDTH',LWIDTH)
          CALL CPSETI('AIA - AREA ABOVE IDENTIFIER',I+1)
          CALL CPSETI('AIB - AREA BELOW IDENTIFIER',I)
c          IF(MOD(I-1,3).EQ.0)
c     1     CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
          CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
        END DO
        CALL CPSETC('ILT - INFORMATION LABEL TEXT STRING',
     1   TRIM(STRING))
        CALL CPSETR('ILX - INFORMATION LABEL X POSITION',0.5)
        CALL CPSETR('ILY - INFORMATION LABEL Y POSITION',0.95)
        CALL CPSETI('ILP - INFORMATION LABEL POSITION',-3)

        SPECIAL=-1E30
        CALL CPSETR('SPV - SPECIAL VALUE',SPECIAL)

        IC=1
        X(IC)=XLO
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YHI
        IC=IC+1
        X(IC)=XHI
        Y(IC)=YLO
        IC=IC+1
        X(IC)=XLO
        Y(IC)=YLO
        if (b2plot_debug.ge.8) then
          write(*,*) 'Computational domain bounding box : '
          write(*,*) 'X : ', X(1:5)
          write(*,*) 'Y : ', Y(1:5)
        end if
        CALL AREDAM(MAP,X,Y,IC,2,-1,0)
cxpb Obtaining the computational domain contour
        IC=1
        X(IC)=real(max(-1.5_R8,values(1)),R4)
        Y(IC)=real(max(-1.5_R8,values(3)),R4)
        IC=IC+1
        X(IC)=real(max(-1.5_R8,values(1)),R4)
        Y(IC)=real(min(real(mpg%ny,R8)+0.5_R8,values(4)),R4)
        IC=IC+1
        X(IC)=real(min(real(mpg%nx,R8)+0.5_R8,values(2)),R4)
        Y(IC)=real(min(real(mpg%ny,R8)+0.5_R8,values(4)),R4)
        IC=IC+1
        X(IC)=real(min(real(mpg%nx,R8)+0.5_R8,values(2)),R4)
        Y(IC)=real(min(-1.5_R8,values(3)),R4)
        IC=IC+1
        X(IC)=real(min(-1.5_R8,values(1)),R4)
        Y(IC)=real(min(-1.5_R8,values(3)),R4)
        if (b2plot_debug.ge.8) then
          write(*,*) 'Computational domain data bounds : '
          write(*,*) 'X : ', X(1:5)
          write(*,*) 'Y : ', Y(1:5)
        end if
        CALL AREDAM(MAP,X,Y,IC,2,0,1)
        IF(.NOT.COLPRT.OR..NOT.CONTPLOT) CALL CURVED(X,Y,IC)

        CALL CPSETI('GIC - GROUP IDENTIFIER FOR CONTOURLINES',1)
        CALL CPSETI('GIL - GROUP IDENTIFIER FOR LABELS',2)
        CALL CPSETI('GIS - GROUP IDENTIFIER FOR VSTRIPS',6)
        CALL CPSETI('NVS - NUMBER OF VSTRIPS',0)

        IF(.NOT.CONTPLOT) GOTO 5152

        IC=0
        ixstrt=max(-1,min(mpg%nx,nint(values(1))))
        iystrt=max(-1,min(mpg%ny,nint(values(3))))
        ixend =max(-1,min(mpg%nx,nint(values(2))))
        iyend =max(-1,min(mpg%ny,nint(values(4))))
        DO IY=iystrt,iyend
          DO IX=ixstrt,ixend
            IC=IC+1
            X(IC)=real(ix,R4)    ! was real(CR(IX,0),R4)
            Y(IC)=real(iy,R4)    ! was real(CZ(IY,0),R4)
            if (mpg%imapCv(ix,iy).ne.0) then
              Z(IC)=real(FIELD(mpg%imapCv(IX,IY)),R4)
            else
              Z(IC)=0.0_R4
            end if
          ENDDO
        ENDDO
        NNX=ixend-ixstrt+1
        NNY=iyend-iystrt+1
        IF(IC.NE.NNX*NNY) WRITE(*,*) 'ic <> nnx*nny',IC,NNX,NNY
        IF (b2plot_debug.ge.8)
     >   WRITE(*,*) 'NX,NY,NNX,NNY',mpg%nx,mpg%ny,NNX,NNY
C         CALL CPCNRC(z,NNX,NNX,NNY,ZLO,ZHI,0.0,1,0,0)
        CALL CPSETI('SET - DO-SET-CALL FLAG',0)
        CALL CPSETI('MAP - MAPPING FLAG',1)
        CALL CPSETR('XC1 - X COORDINATE AT I=1', real(ixstrt,R4))
        CALL CPSETR('XCM - X COORDINATE AT I=M', real(ixend,R4))
        CALL CPSETR('YC1 - Y COORDINATE AT J=1', real(iystrt,R4))
        CALL CPSETR('YCN - Y COORDINATE AT J=N', real(iyend,R4))
        CALL CPRECT(z, nnx, nnx, nny, RWRK, LRWK, IWRK, LIWK)
C Add contours to area map
        CALL CPCLAM(z, RWRK, IWRK, MAP)
C Add label boxes to area map
        IF(LABELS) THEN
          CALL CPLBAM(z, RWRK, IWRK, MAP)
        ENDIF
C Set fill style to solid, and fill contours if COLPRT==TRUE
        IF(COLPRT) THEN
          call entsr(irold,1)
          CALL ARPRAM(MAP,1,0,0)
          if(nerro(nerrf).ne.0) then
            call eprin
            call errof
            call retsr(irold)
            call subend ()
            return
          endif
          call retsr(irold)
          CALL ARSCAM(MAP,XWRK,YWRK,NWRK,IAREA,IGRP,NOGRPS,FILL)
        ENDIF
C Draw Labels
        IF(LABELS) THEN
          CALL CPLBDR(z,RWRK,IWRK)
        ENDIF
C Draw Contours
        IF(CONTOURS) THEN
C            CALL CPCLDR(z,RWRK,IWRK)
          CALL CPCLDM(z,RWRK,IWRK,MAP,CPDRPL)
        ENDIF

 5152   CONTINUE

      ENDIF

      IF(LBLBR.AND.CONTPLOT.and.colprt) THEN
        DO I=1,NCONT+1
          LFIN(I)=choose_color(i-1,0)
          CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
          CALL CPGETR('CLV - CONTOUR LEVEL VALUES',CLV)
          if(nofcls.eq.0) then ! standard contouring
            WRITE (LBLS(I), FMTZ) CLV
          else
            if(abs(clv).ge.1e3.or.abs(clv).lt.1e-2) then
              write(LBLS(I),'(1p,e10.3)') CLV
            elseif(abs(clv).ge.1e2) then
              write(LBLS(I),'(f8.3)') CLV
            elseif(abs(clv).ge.1e1) then
              write(LBLS(I),'(f8.4)') CLV
            else
              write(LBLS(I),'(f8.5)') CLV
            endif
          endif
          CALL SHRINK_LABEL(LBLS(I))
        ENDDO
        i=ncont+2
        LFIN(I)=choose_color(i-1,0)
        LBLS(I:100)=' '
        CALL GSLWSC(1.0)
        if(noax.gt.0) then
          if(rlo.gt.clval(1).and.rhi.le.clval(ncont+1)) then
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),1,LBLS,
     +       NCONT+1,1)
          elseif(rlo.le.clval(1).and.rhi.le.clval(ncont+1)) then
            do i=ncont+2,2,-1
              lbls(i)=lbls(i-1)
            enddo
            lbls(1)=' '
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     +       NCONT+2,1)
          elseif(rlo.gt.clval(1).and.rhi.gt.clval(ncont+1)) then
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT+1,0.2/(0.2+0.8*lbsz),1.0,LFIN(2),0,LBLS,
     +       NCONT+2,1)
          else
            CALL LBLBAR(1,
     1       colorbar_left*XSIZE+XPOS,colorbar_right*XSIZE+XPOS,
     2       0.10*YSIZE+YPOS,0.90*YSIZE+YPOS,
     3       NCONT+2,0.2/(0.2+0.8*lbsz),1.0,LFIN,0,LBLS,
     +       NCONT+1,1)
          endif
        endif
      ENDIF

C Draw Mesh/Separatrix if wanted
      if(wall) call data_copy(is)
      call additional_elements(nCv,physical)

      IF(DEBUG) THEN
        DO I=1,2
          WRITE(DEBUGSTRING,'(I5)') I
          CALL ARDBPA(MAP,I,DEBUGSTRING)
        ENDDO
      ENDIF

      call subend ()
      RETURN
      END SUBROUTINE PLOTCONTOUR

c$$$      SUBROUTINE LABLIT (INUM,SIZE,ANGR,XPOS,YPOS)
c$$$      use b2mod_types
c$$$      implicit none
c$$$      INTEGER INUM
c$$$      REAL(KIND=R4) :: SIZE,ANGR,XPOS,YPOS
c$$$C
c$$$      INTEGER I
c$$$      REAL(KIND=R4) :: RTOD
c$$$      CHARACTER*8 CNUM
c$$$C
c$$$C Generate a label expressing the value of the integer INUM and place it
c$$$C on the frame at the position (XPOS,YPOS) in the current user system,
c$$$C using characters of size SIZE (in the fractional coordinate system),
c$$$C and an angle of ANGR (in radians).
c$$$C
c$$$C Conversion factor, radians to degrees.
c$$$C
c$$$      DATA RTOD / 57.2957795130823_R4 /
c$$$C
c$$$C Encode the integer into CNUM.
c$$$C
c$$$      WRITE (CNUM,'(I8)') INUM
c$$$C
c$$$C Write the label and return.  The loop finds the first non-blank in the
c$$$C label so that we can use sub-string notation in the call to PLCHMQ.
c$$$C
c$$$      DO 101 I=1,8
c$$$        IF (CNUM(I:I).NE.' ') THEN
c$$$          CALL PLCHMQ (XPOS,YPOS,CNUM(I:8),SIZE,RTOD*ANGR,0.)
c$$$          RETURN
c$$$        END IF
c$$$  101 CONTINUE
c$$$C
c$$$C If no non-blank character is found, return, having done nothing.
c$$$C
c$$$      RETURN
c$$$C
c$$$      END SUBROUTINE LABLIT
c$$$
