*-----------------------------------------------------------------------
*.specification

      subroutine write_traduitoutb2us(nfin,nCv,g,m,cvBp,cvBt,cvFpsi)
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nfin, nCv
      real (kind=R8) :: cvBp(nCv), cvBt(nCv), cvFpsi(nCv)
      type(geometry) :: g 
      type(mapping)  :: m  

*-----------------------------------------------------------------------
*.documentation
*
*     This routine writes a traduit.out.b2us-file from the geometry and
*     mapping structures.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, iVx, iFt, iFs, n2, idum(0:9)

*-----------------------------------------------------------------------
*.computation

*   ..write version header
      call cfverw (nfin, '03.002.000')

*   ..write basic mesh dimensions
      idum(0) = m%nCi
      idum(1) = m%nFc
      idum(2) = m%nVx
      idum(3) = m%nCg
      idum(4) = m%nFs
      idum(5) = m%nFt
      call cfwuin (nfin,6,idum,'nCi,nFc,nVx,nCg,nFs,nFt')
      idum(0) = m%nCmxVx
      idum(1) = m%nCmxFc
      idum(2) = m%nFmxCv
      idum(3) = m%nVmxCv
      idum(4) = m%nVmxFc
      call cfwuin (nfin,5,idum,'nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc')
      idum(0) = m%isClassicalGrid
      call cfwuin (nfin,1,idum,'isClassicalGrid')
      if (m%isClassicalGrid == 1) then
        idum(0) = m%nx
        idum(1) = m%ny
        idum(2) = m%nncut
        call cfwuin (nfin, 3, idum, 'nx,ny,nncut')
      endif

*   ..write vertex information
      write (nfin,'(a)') '*cf: Vx vxX vxY vxPsi vxBx vxBy vxFfbz'
      do iVx = 1, nnVx
        write(nfin,'(I11,6(1X,F15.8)') iv(iVx), g%vxX(iVx), g%vxY(iVx),
     &    g%vxFpsi(iVx), g%vxBb(iVx,0), g%vxBb(iVx,1), g%vxFfbz(iVx)
      enddo

*   ..write cell information
      write (nfin,'(a)') '*cf: cv cvVxP(:,1) cvVxP(:,2) cvX cvY psi bp bt cflags(:) cvReg cvFt'
      do iCv = 1, m%nCv
        write(nfin,'(2(I11,1X),I4,1X,5(F15.8,1X),2(I4,1X),I11') ic(iCv), m%cvVxP(iCv,1), 
     &    m%cvVxP(iCv,2), g%cvX(iCv), g%cvY(iCv), cvFpsi(iCv),
     &    cvBp(iCv), cvBt(iCv), m%cvLbl(iCv), m%cvReg(iCv), 
     &    m%cvFt(iCv)
      enddo
      call cfwuin(nfin,m%nCmxVx,m%cvVx,'cvVx')
      call cfwuin(nfin,m%nCmxFc,m%cvFc,'cvFc')

*   ..write face information
      write (nfin,'(a)') '*cf: fc fcVx(:,1) fcVx(:,2) fcLbl fcReg fcAligned'
      do iFc = 1, nnFc
        write(nfin,'(3(I11,1X),3(I4,1X)') if(iFc), m%fcVx(iFc,1), 
     &    m%fcVx(iFc,2), m%fcLbl(iFc), m%fcReg(iFc), 
     &    m%fcAligned(iFc)
      enddo

*   ..write flux tube information
      if (m%nFt .gt. 0) then
        write (nfin,'(a)') '*cf: ft ftCvP(:,1) ftCvP(:,2) ftFcP(:,1) ftFcP(:,2) ftReg'
        do iFt = 1, m%nFt
          write(nfin,'(5(I11,1X),I11)') ifcv(iFt), m%ftCvP(iFt,1), 
     &      m%ftCvP(iFt,2), m%ftFcP(iFt,1), m%ftFcP(iFt,2), 
     &      m%ftReg(iFt)
        enddo   
        call cfwuin(nfin,m%nCv,m%ftCv,'ftCv')
        call cfwuin(nfin,m%nFc,m%ftFc,'ftFc')
      endif

*   ..write flux surface information
      if (m%nFs .gt. 0) then
        write(nfin,'(a)') '*cf: fs fsFcP(:,1) fsFcP(:,2) fsPsi'
        do iFs = 1, m%nFs
          write(nfin,'(3(I11,1X),F15.8)') ifsv(iFs),  m%fsFcP(iFs,1),
     &      m%fsFcP(iFs,2), g%fsPsi(iFs)
        enddo
        call cfwuin(nfin,m%nFc,m%fsFc,'fsFc')
      endif

*   ..write information for back-mapping to structured grid
      if (m%isClassicalGrid.eq.1) then
        n2 = (m%nx+2)*(m%ny+2)
        call cwruin(nfin,n2,m%imapCv,'imapCv')
        call cwruin(nfin,n2,m%imapFcx,'imapFcx')
        call cwruin(nfin,n2,m%imapFcy,'imapFcy')
        call cwruin(nfin,n2,m%imapVx,'imapVx')
      endif

      return

*-----------------------------------------------------------------------
*.end write_traduitoutb2us

      end subroutine write_traduitoutb2us
