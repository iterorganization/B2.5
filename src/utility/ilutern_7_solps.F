      SUBROUTINE ILUTERN_7(N1, N2, N1D,                            !sv 11.01.00
     . ASO1, ASO, AWE, ACP, AEA, ANO, ANO1, RHS, SOL)              !sv 11.01.00
      use b2mod_types
      use b2mod_ma28_for_7diag
      use b2mod_indirect
      IMPLICIT NONE
      external bd1copy
C

*  -- INPUT ARGUMENTS --
      INTEGER N1, N2, N1D
      real (kind=R8) ::
     *   ASO1(N1D,N2), ASO(N1D,N2), AWE(N1D,N2), ACP(N1D,N2),      !sv 11.01.00
     .   AEA(N1D,N2), ANO(N1D,N2), ANO1(N1D,N2), RHS(N1D,N2)       !sv 11.01.00

*  -- OUTPUT ARGUMENTS --
      real (kind=R8) ::
     *   SOL(N1D,N2)

*  -- WORKSPACE ARGUMENTS --
      real (kind=R8) ::
     *   A1HN

*  -- LOCAL SCALARS --
      INTEGER N, I, J, K, IK, NZ
      real (kind=R8) ::
     *   x05baf,cpu_start,cpu_end
      integer iflag

      real (kind=R8) ::
     *   u
      parameter (u=1.0e0_R8)

c ma28 common blocks
      real (kind=R8) ::
     *   eps,rmin,resid
      integer irncp,icncp,minirn,minicn,irank
      logical abort1,abort2
      common /ma28fdc/ eps,rmin,resid,irncp,icncp,minirn,minicn,
     1 irank,abort1,abort2
      real (kind=R8) ::
     *   tol,themax,big,dxmax,errmax,dres,cgce
      integer ndrop,maxit,noiter,nsrch,istart
      logical lbig
      common /ma28hdc/ tol,themax,big,dxmax,errmax,dres,cgce,
     1 ndrop,maxit,noiter,nsrch,istart,lbig
c

*  -- PROCEDURES --
      real (kind=R8) ::
     *   b2sasum
      EXTERNAL XERRAB, b2sasum, x05baf

      logical lfirst
      integer ncall, nwrite
      save lfirst, ncall, nwrite
      data lfirst /.true./, ncall /0/, nwrite /0/

*=======================================================================
*//COMPUTATION//

*  ---------------------------------------------------------------------
*  TEST INPUTS AND SET AUXILIARIES.
*  ---------------------------------------------------------------------

      if(ncall.eq.0) then
         call ipgeti ('ma28_for_7d_nwrite', nwrite)
      endif
c     write(*,*) 'ncall,nwrite',ncall,nwrite

      cpu_start=x05baf()
      if(ncall.lt.1000) then
*  -- CHECK N1, N2, N1D, ITERS, ACPAR, ITPAR --
         IF (N1 .LT. 2 .OR. N2 .LT. 2 .OR. N1D .LT. N1) THEN
          CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT N1, N2, N1D') !sv 06.07.00
         END IF

*  -- CHECK EDGE VALUES OF ASO, AWE, AEA, ANO --
        do i = 1, N1
          do j = 1, N2
            IF (ASO(i,j) .GT. 0 .AND. bottomiy(i-2,j-2).eq.-2) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ASO')
            ENDIF
            IF (ASO1(i,j) .GT. 0 .AND. bottomiy(i-2,j-2).eq.-2) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ASO1')
            ENDIF
            IF (ASO1(i,j) .GT. 0 .AND. 
     &       bottomiy(bottomix(i-2,j-2),bottomiy(i-2,j-2)).eq.-2) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ASO1')
            ENDIF
            IF (ANO(i,j) .GT. 0 .AND. topiy(i-2,j-2).eq.N2-1) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ANO')
            ENDIF
            IF (ANO1(i,j) .GT. 0 .AND. topiy(i-2,j-2).eq.N2-1) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ANO1')
            ENDIF
            IF (ANO1(i,j) .GT. 0 .AND. 
     &       topiy(topix(i-2,j-2),topiy(i-2,j-2)).eq.N2-1) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ANO1')
            ENDIF
            IF (AWE(i,j) .GT. 0 .AND. leftix(i-2,j-2).eq.-2) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT AWE')
            ENDIF
            IF (AEA(i,j) .GT. 0 .AND. rightix(i-2,j-2).eq.N1-1) THEN
              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT AEA')
            ENDIF
          enddo
        enddo
      endif
      N = N1*N2
C    ZUWEISUNG
      NZ=0
      IK=0
      ICNH=0
      IRNH=0
      DO K=1,NC0
         IF (IJSH(K).EQ.-1) A1HN=ASO1(ISH(K),JSH(K))               !sv 26.01.00
         IF (IJSH(K).EQ.0) A1HN=ASO(ISH(K),JSH(K))
         IF (IJSH(K).EQ.1) A1HN=AWE(ISH(K),JSH(K))
         IF (IJSH(K).EQ.2) A1HN=ACP(ISH(K),JSH(K))
         IF (IJSH(K).EQ.3) A1HN=AEA(ISH(K),JSH(K))
         IF (IJSH(K).EQ.4) A1HN=ANO(ISH(K),JSH(K))
         IF (IJSH(K).EQ.5) A1HN=ANO1(ISH(K),JSH(K))                !sv 26.01.00
c       IF (A1HN.EQ.0.0e0_R8) THEN
         if(
     1   (IJSH(K).EQ.-1.and.jsh(k).eq.1).or.                       !sv 26.01.00
     1   (IJSH(K).EQ.-1.and.jsh(k).eq.2).or.                       !sv 26.01.00
     1   (IJSH(K).EQ.0.and.jsh(k).eq.1).or.
     1   (IJSH(K).EQ.1.and.ish(k).eq.1).or.
     1   (IJSH(K).EQ.3.and.ish(k).eq.n1).or.
     1   (IJSH(K).EQ.5.and.jsh(k).eq.n2-1).or.                     !sv 26.01.00
     1   (IJSH(K).EQ.5.and.jsh(k).eq.n2).or.                       !sv 26.01.00
     1   (IJSH(K).EQ.4.and.jsh(k).eq.n2)) then
         ELSE
          NZ=NZ+1
          A1H(NZ)=A1HN
          IRNH(NZ)=IRN(K)
          IK=IK+1
          ICNH(NZ)=ICN(K)
         ENDIF
      enddo
      DO I=1,N1
         DO J=1,N2
          R1(I+(J-1)*N1)=RHS(I,J)
          S1(I+(J-1)*N1)=0.0e0_R8
         enddo
      enddo
c-dpc
      if(nwrite.gt.0) then
        if(nwrite.gt.ncall) then
          if(ncall.eq.0) then
            write(*,*) 'Opening 78'                            !sv 11.01.00
            open(78,file='b2_7_diagonal_matrix')               !sv 11.01.00
          endif
          write(*,*) 'Writing 78'                              !sv 11.01.00
          write(78,*) n,' n'
          do i=1,n
            write(78,*) r1(i)                                  !sv 11.01.00
          enddo
          write(78,*) nz,' nz'                                 !sv 11.01.00
          do i=1,nz
            write(78,*) IRNH(i),ICNH(i),A1H(i)                 !sv 11.01.00
          enddo
        elseif(nwrite.eq.ncall) then
          write(*,*) 'Closing 78'                              !sv 11.01.00
          close(78)                                            !sv 11.01.00
        endif
      endif
c-dpc
      if(nltrsol.eq.0) then
      elseif(nltrsol.eq.1) then
         WRITE(6,*) 'ILUTERN'
      else
         WRITE(6,*) 'ILUTERN'
      endif
      if(lfirst) then
         write(*,*) 'ilutern_7: nz = ',nz                       !sv 11.01.00
         if(nz.gt.N_NONZEROS) stop 'Increase N_NONZEROS'
         iflag=0
      endif
      do i=1,nz
         irnh_save(i)=irnh(i)
         icnh_save(i)=icnh(i)
         a1h_save(i)=a1h(i)
      enddo
      if(.not.lfirst) then
         call ma28bdcopy
     1   (n,nz,a1h,solvdim,irnh,icnh,icnhkeep,ikeep,iw28,w,iflag)
         if(nltrsol.gt.0) write(*,*) 'iflag ',iflag
         if(iflag.lt.0) then
          do i=1,nz
            irnh(i)=irnh_save(i)
            icnh(i)=icnh_save(i)
            a1h(i)=a1h_save(i)
          enddo
         endif
      endif
      if(lfirst.or.iflag.lt.0) then
         call ma28adcopy
     1   (n,nz,a1h,solvdim,irnh,solvdim1,icnh,u,ikeep,iw28,w,iflag)
         if(nltrsol.gt.0)
     1   write(*,'(a,i6,i8,a,i8,i8,a,i8)') 
     .    'iflag,minirn/dim,minicn/dim ',
     2     iflag,minirn,' /',solvdim1,minicn,' /',solvdim
         do i=1,solvdim
          icnhkeep(i)=icnh(i)
         enddo
         lfirst=.false.
      endif
      call ma28idcopy
     1 (n,nz,a1h_save,irnh_save,icnh_save,solvdim,a1h,icnhkeep,ikeep,
     1 r1,s1,res,w,1,1e-10_R8,iflag)
      if(nltrsol.gt.0) write(*,'(a,i3,1p,2g15.8,i6)')
     1 'iflag,errmax,dres,noiter ',iflag,errmax,dres,noiter
      DO 50,K=1,N
        J=(K-1)/N1+1
        I=K-(J-1)*N1
        SOL(I,J)=S1(K)
  50  CONTINUE

*  ---------------------------------------------------------------------
      cpu_end=x05baf()
      ncall=ncall+1
      RETURN

      entry restart_ma28_for_7                              !sv 11.01.00
      lfirst=.true.
      write(*,*) 'Restarting ma28'
      return

*=======================================================================
*//END ILUTERN//

      END
