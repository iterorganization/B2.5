      SUBROUTINE ILUTERN(N1, N2, N1D, ASO, AWE, ACP, AEA, ANO,
     .      RHS, SOL, NAME)                                          !srv 18.03.10
      use b2mod_types
      use b2mod_ma28
      use b2mod_indirect
      IMPLICIT NONE
      external bd1
C

*  -- INPUT ARGUMENTS --
      INTEGER N1, N2, N1D
      character*(*) name                                             !srv 18.03.10
      REAL (KIND=R8) ::
     *   ASO(N1D,N2), AWE(N1D,N2), ACP(N1D,N2),
     .      AEA(N1D,N2), ANO(N1D,N2), RHS(N1D,N2)

*  -- OUTPUT ARGUMENTS --
      REAL (KIND=R8) ::
     *   SOL(N1D,N2)

*  -- WORKSPACE ARGUMENTS --
      REAL (KIND=R8) ::
     *   A1HN

*  -- LOCAL SCALARS --
      INTEGER N, I, J, K, IK, NZ
      integer ix, iy
      integer iflag

      real (kind=R8) ::
     *   u
      parameter (u=1.0_R8)

c ma28 common blocks
      real (kind=R8) ::
     *   eps,rmin,resid
      integer irncp,icncp,minirn,minicn,irank
      logical abort1,abort2
      common /ma28fd/ eps,rmin,resid,irncp,icncp,minirn,minicn,
     1 irank,abort1,abort2
      real (kind=R8) ::
     *   tol,themax,big,dxmax,errmax,dres,cgce
      integer ndrop,maxit,noiter,nsrch,istart
      logical lbig
      common /ma28hd/ tol,themax,big,dxmax,errmax,dres,cgce,
     1 ndrop,maxit,noiter,nsrch,istart,lbig
c

*  -- PROCEDURES --
      real (kind=R8) ::
     *   b2sasum
      EXTERNAL XERRAB, b2sasum

      logical lfirst
      integer ncall, nwrite
      integer, save :: kwrite_fmt = -1                                !srv 18.03.10
      integer, save :: kwrite_bnr = -1                                !srv 18.03.10
      integer, save :: b2swtest=0
      save lfirst, ncall, nwrite
      data lfirst /.true./, ncall /0/, nwrite /0/

*=======================================================================
*//COMPUTATION//

*  ---------------------------------------------------------------------
*  TEST INPUTS AND SET AUXILIARIES.
*  ---------------------------------------------------------------------

      if(ncall.eq.0) then
       call ipgeti ('ma28_nwrite', nwrite)
       call ipgeti ('ma28_ilutern_kwrite_fmt', kwrite_fmt)                    !srv 18.05.02
       call ipgeti ('ma28_ilutern_kwrite_bnr', kwrite_bnr)                    !srv 18.05.02
       call ipgeti ('b2swtest', b2swtest)
      endif
c     write(*,*) 'ncall,nwrite',ncall,nwrite

      if(ncall.lt.1000) then
*  -- CHECK N1, N2, N1D, ITERS, ACPAR, ITPAR --
        IF (N1 .LT. 2 .OR. N2 .LT. 2 .OR. N1D .LT. N1) THEN
          CALL XERRAB ('ILUTERN -- FAULTY ARGUMENT N1, N2, N1D')
        END IF

*  -- CHECK EDGE VALUES OF ASO, AWE, AEA, ANO --
        do i = 1, N1
          do j = 1, N2
            IF (ASO(i,j) .GT. 0 .AND. bottomiy(i-2,j-2).eq.-2) THEN
              CALL XERRAB ('ILUTERN -- FAULTY ARGUMENT ASO')
            ENDIF
            IF (ANO(i,j) .GT. 0 .AND. topiy(i-2,j-2).eq.N2-1) THEN
              CALL XERRAB ('ILUTERN -- FAULTY ARGUMENT ANO')
            ENDIF
            IF (AWE(i,j) .GT. 0 .AND. leftix(i-2,j-2).eq.-2) THEN
              CALL XERRAB ('ILUTERN -- FAULTY ARGUMENT AWE')
            ENDIF
            IF (AEA(i,j) .GT. 0 .AND. rightix(i-2,j-2).eq.N1-1) THEN
              CALL XERRAB ('ILUTERN -- FAULTY ARGUMENT AEA')
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
        IF (IJSH(K).EQ.0) A1HN=ASO(ISH(K),JSH(K))
        IF (IJSH(K).EQ.1) A1HN=AWE(ISH(K),JSH(K))
        IF (IJSH(K).EQ.2) A1HN=ACP(ISH(K),JSH(K))
        IF (IJSH(K).EQ.3) A1HN=AEA(ISH(K),JSH(K))
        IF (IJSH(K).EQ.4) A1HN=ANO(ISH(K),JSH(K))
c        IF (A1HN.EQ.0.0_R8) THEN
        if(
     1   (IJSH(K).EQ.0.and.jsh(k).eq.1).or.
     1   (IJSH(K).EQ.1.and.ish(k).eq.1).or.
     1   (IJSH(K).EQ.3.and.ish(k).eq.n1).or.
     1   (IJSH(K).EQ.4.and.jsh(k).eq.n2)) then
        ELSE
         IF(A1HN.NE.0.0_R8) THEN                                !srv 06.07.10 
          NZ=NZ+1
          A1H(NZ)=A1HN
          IRNH(NZ)=IRN(K)
          IK=IK+1
          ICNH(NZ)=ICN(K)
         ENDIF                                                  !srv 06.07.10
        ENDIF
      enddo
      DO I=1,N1
        DO J=1,N2
          R1(I+(J-1)*N1)=RHS(I,J)
          S1(I+(J-1)*N1)=0.0_R8
        enddo
      enddo

csw 12oct2012 in case of b2swtest=T --> reinit solver
c             (maybe necessary in case of multi-fluid?)
      if(b2swtest.eq.1) lfirst=.true.

      if(kwrite_bnr.eq.ncall) then                                   !srv 18.03.10 {
       write(*,*) 'Opening 76'
       open(76,file='b2_5_rhs_binary_'//trim(name)//'.dat',
     .         form='unformatted')                                   !srv 11.03.13
       write(76) n
       do i=1,n
        write(76) r1(i)
       enddo
       close(76)
       write(*,*) 'Opening 77'
       open(77,file='b2_5_mtr_binary_'//trim(name)//'.dat',
     .         form='unformatted')                                   !srv 11.03.13
       write(77) n, nz
       do i=1,nz
        write(77) IRNH(i),ICNH(i),A1H(i)
       enddo
       close(77)
      endif                             

      if(kwrite_fmt.eq.ncall) then      
       write(*,*) 'Opening 76'
       open(76,file='b2_5_rhs_'//trim(name)//'.dat')
       write(76,*) n
       do i=1,n
        write(76,*) r1(i)
       enddo
       close(76)
       open(77,file='b2_5_matrix_'//trim(name)//'.dat')
       write(77,*) n, nz
       do i=1,nz
        write(77,*) IRNH(i),ICNH(i),A1H(i)
       enddo
       close(77)
      endif                                                          !srv 18.03.10 }

c-dpc
      if(nwrite.gt.0) then
        if(nwrite.gt.ncall) then
          if(ncall.eq.0) then
            write(*,*) 'Opening 79'
            open(79,file='b2_matrix')
          endif
          write(*,*) 'Writing 79'
          write(79,*) 'Called from '//trim(name)
          write(79,*) n,' n'
          write(79,*) 'ASO AWE ACP AEA ANO RHS'
          do i=1,n
            iy = (i-1)/n1+1
            ix = i-(iy-1)*n1
            write(79,'(i6,2(i4,1x),1p,6(e15.8,1x))') i, ix-2, iy-2,
     .       aso(ix,iy), awe(ix,iy), acp(ix,iy), aea(ix,iy), ano(ix,iy),
     .       rhs(ix,iy)
          enddo
          write(79,*) nz,' nz'
          do i=1,nz
            write(79,'(2(i6,1x),1p,1e15.8)') IRNH(i),ICNH(i),A1H(i)
          enddo
        elseif(nwrite.eq.ncall) then
          write(*,*) 'Closing 79'
          close(79)
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
        write(*,*) 'ilutern: nz = ',nz
        call xertst (nz.le.N_NONZEROS, 'Increase N_NONZEROS')
        iflag=0
      endif
      do i=1,nz
        irnh_save(i)=irnh(i)
        icnh_save(i)=icnh(i)
        a1h_save(i)=a1h(i)
      enddo
      if(.not.lfirst) then
        call ma28bd
     1   (n,nz,a1h,solvdim,irnh,icnh,icnhkeep,ikeep,iw28,w,iflag)
        if(nltrsol.gt.0 .or. iflag.lt.0) write(*,*) 'iflag ',iflag   !srv 28.06.16
        if(iflag.lt.0.and.iflag.ne.-14) then   ! do not do this for duplicate elements
          do i=1,nz
            irnh(i)=irnh_save(i)
            icnh(i)=icnh_save(i)
            a1h(i)=a1h_save(i)
          enddo
        endif
      endif
      if(lfirst.or.(iflag.lt.0.and.iflag.ne.-14)) then   ! do not do this for duplicate elements
        call ma28ad
     1   (n,nz,a1h,solvdim,irnh,solvdim1,icnh,u,ikeep,iw28,w,iflag)
        if(nltrsol.gt.0 .or. iflag.lt.0)                             !srv 28.06.16
     1   write(*,'(a,i6,i8,a,i8,i8,a,i8)') 
     .    'iflag,minirn/dim,minicn/dim ',
     2     iflag,minirn,' /',solvdim1,minicn,' /',solvdim
        do i=1,solvdim
          icnhkeep(i)=icnh(i)
        enddo
        lfirst=.false.
      endif
      call ma28id
     1 (n,nz,a1h_save,irnh_save,icnh_save,solvdim,a1h,icnhkeep,ikeep,
     1 r1,s1,res,w,1,1e-10_R8,iflag)
      if(nltrsol.gt.0 .or. iflag.lt.0) write(*,'(a,i3,1p,2g15.8,i6)') !srv 28.06.16
     1 'iflag,errmax,dres,noiter ',iflag,errmax,dres,noiter

      if(iflag.lt.0) then                                            !srv 28.06.16 {

        do i=1,nz
          irnh(i) = irnh_save(i)
          icnh(i) = icnh_save(i)
          a1h(i)  = a1h_save(i)
        enddo

        call ma28ad
     1   (n,nz,a1h,solvdim,irnh,solvdim1,icnh,u,ikeep,iw28,w,iflag)
        if(nltrsol.gt.0 .or. iflag.lt.0) write(*,2)
     1     trim(name), iflag, minirn, solvdim1, minicn, solvdim
   2    format(' ma28ad for ',a13,': iflag =',i9/,
     1   ' minirn =',i9,' solvdim1 =',i9,' minicn =',i9,' solvdim =',i9)
        do i=1,solvdim
          icnhkeep(i)=icnh(i)
        enddo

        call ma28id
     1   (n,nz,a1h_save,irnh_save,icnh_save,solvdim,a1h,icnhkeep,ikeep,
     1    r1,s1,res,w,1,1e-10_R8,iflag)

        if(iflag.lt.0) then
          write(*,"(a)") 'Problem with solving '//trim(name))
        else
          write(*,"(3a)") ' Problem of ', trim(name), ' has been solved'
        endif
        if(nltrsol.gt.0 .or. iflag.lt.0)
     &   write(*,3)  trim(name), iflag, errmax, dres, noiter
   3    format(' ma28id for ',a13,' iflag =',i9,
     2         ' errmax =',es9.1,' dres =',es9.1,' noiter =',i2)

      endif                                                          !srv 28.06.16 }

      DO K=1,N
        J=(K-1)/N1+1
        I=K-(J-1)*N1
        SOL(I,J)=S1(K)
      ENDDO

*  ---------------------------------------------------------------------
      ncall=ncall+1
      RETURN

      entry restart_ma28
      lfirst=.true.
      write(*,*) 'Restarting ma28'
      return

*=======================================================================
*//END ILUTERN//

      END
