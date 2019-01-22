      SUBROUTINE ILUTERN_7(N1, N2, N1D,                            !srv 11.01.00 26.12.17
     . ASO1, ASO, AWE, ACP, AEA, ANO, ANO1, RHS, SOL, name)        !srv 11.01.00
      use b2mod_types
      use b2mod_ma28_for_7diag
      use b2mod_indirect
      use b2mod_openmp
      IMPLICIT NONE
#ifndef GFORTRAN
      external bd1copy3                                            !srv 26.12.17
#endif
C

*  -- INPUT ARGUMENTS --
      INTEGER N1, N2, N1D
      real (kind=R8) ::
     *   ASO1(N1D,N2), ASO(N1D,N2), AWE(N1D,N2), ACP(N1D,N2),      !srv 11.01.00
     .   AEA(N1D,N2), ANO(N1D,N2), ANO1(N1D,N2), RHS(N1D,N2)       !srv 11.01.00
      character*(*) name                                           !srv 18.05.02

*  -- OUTPUT ARGUMENTS --
      real (kind=R8) ::
     *   SOL(N1D,N2)

*  -- WORKSPACE ARGUMENTS --
      real (kind=R8) ::
     *   A1HN

*  -- LOCAL SCALARS --
      INTEGER N, I, J, K, IK, NZ
      real (kind=R8) :: cpu_start, cpu_end
      integer iflag

      real (kind=R8) ::
     *   u
      parameter (u=1.0e0_R8)
      real (kind=R8), save :: prec = 1.0e-10_R8                    !srv 16.01.18

c ma28 common blocks
      real (kind=R8) ::
     *   eps,rmin,resid
      integer irncp,icncp,minirn,minicn,irank
      logical abort1,abort2
      common /ma28fdc3/ eps,rmin,resid,irncp,icncp,minirn,minicn,  !srv 26.12.17
     1 irank,abort1,abort2
      real (kind=R8) ::
     *   tol,themax,big,dxmax,errmax,dres,cgce
      integer ndrop,maxit,noiter,nsrch,istart
      logical lbig
      common /ma28hdc3/ tol,themax,big,dxmax,errmax,dres,cgce,     !srv 26.12.17
     1 ndrop,maxit,noiter,nsrch,istart,lbig
!$OMP THREADPRIVATE(/ma28fdc3/,/ma28hdc3/)
c

*  -- PROCEDURES --
      EXTERNAL XERRAB

      logical lfirst
      integer ncall, nwrite
      integer, save :: kwrite_fmt = -1                             !srv 18.05.02
      integer, save :: kwrite_bnr = -1                             !srv 18.05.02
      save lfirst, ncall, nwrite
!$OMP THREADPRIVATE(lfirst, ncall, nwrite, kwrite_fmt, kwrite_bnr)
      data lfirst /.true./, ncall /0/, nwrite /0/

*=======================================================================
*//COMPUTATION//
*   ..subprogram start-up calls
      call subini ('ilutern_7')

*  ---------------------------------------------------------------------
*  TEST INPUTS AND SET AUXILIARIES.
*  ---------------------------------------------------------------------

      if(ncall.eq.0) then
        call ipgeti ('ma28_for_7d_nwrite', nwrite)
        call ipgeti ('ma28_for_7d_kwrite_fmt', kwrite_fmt)          !srv 18.05.02
        call ipgeti ('ma28_for_7d_kwrite_bnr', kwrite_bnr)          !srv 18.05.02
        call ipgetr ('ma28_for_7d_prec', prec)                      !srv 16.01.18
        if (in_parallel() .and. ( nwrite.ne.0 .or.
     &     kwrite_bnr.ne.-1 .or. kwrite_fmt.ne.-1)) then
!$OMP MASTER
         write(*,*) 'ilutern_7 OpenMP warning: no file IO ' ,
     &       'in parallel mode'
!$OMP END MASTER
         nwrite = 0
         kwrite_fmt = -1
         kwrite_bnr = -1
       endif 
#ifdef GFORTRAN
       call set_data_bd1copy3
       call set_data_bd2copy3
       call set_data_bd3copy3
#endif
      endif
c     write(*,*) 'ncall,nwrite',ncall,nwrite

      call cpu_time(cpu_start)
!      if(ncall.lt.1000) then                                       !srv 26.12.17 {
!*  -- CHECK N1, N2, N1D, ITERS, ACPAR, ITPAR --
!         IF (N1 .LT. 2 .OR. N2 .LT. 2 .OR. N1D .LT. N1) THEN
!           CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT N1, N2, N1D')  !srv 06.07.00
!         END IF

!*  -- CHECK EDGE VALUES OF ASO, AWE, AEA, ANO --
!        do i = 1, N1
!          do j = 1, N2
!            IF (ASO(i,j) .GT. 0 .AND. bottomiy(i-2,j-2).eq.-2) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ASO')
!            ENDIF
!            IF (ASO1(i,j) .GT. 0 .AND. bottomiy(i-2,j-2).eq.-2) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ASO1')
!            ENDIF
!            IF (ASO1(i,j) .GT. 0 .AND.
!     &       bottomiy(bottomix(i-2,j-2),bottomiy(i-2,j-2)).eq.-2) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ASO1')
!            ENDIF
!            IF (ANO(i,j) .GT. 0 .AND. topiy(i-2,j-2).eq.N2-1) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ANO')
!            ENDIF
!            IF (ANO1(i,j) .GT. 0 .AND. topiy(i-2,j-2).eq.N2-1) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ANO1')
!            ENDIF
!            IF (ANO1(i,j) .GT. 0 .AND.
!     &       topiy(topix(i-2,j-2),topiy(i-2,j-2)).eq.N2-1) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT ANO1')
!            ENDIF
!            IF (AWE(i,j) .GT. 0 .AND. leftix(i-2,j-2).eq.-2) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT AWE')
!            ENDIF
!            IF (AEA(i,j) .GT. 0 .AND. rightix(i-2,j-2).eq.N1-1) THEN
!              CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT AEA')
!            ENDIF
!          enddo
!        enddo
!      endif                                                      !srv 26.12.17 }
      N = N1*N2
C    ZUWEISUNG
      NZ=0
      IK=0
      ICNH=0
      IRNH=0
      DO K=1,NC0
         IF (IJSH(K).EQ.-1) A1HN=ASO1(ISH(K),JSH(K))               !srv 26.01.00
         IF (IJSH(K).EQ.0) A1HN=ASO(ISH(K),JSH(K))
         IF (IJSH(K).EQ.1) A1HN=AWE(ISH(K),JSH(K))
         IF (IJSH(K).EQ.2) A1HN=ACP(ISH(K),JSH(K))
         IF (IJSH(K).EQ.3) A1HN=AEA(ISH(K),JSH(K))
         IF (IJSH(K).EQ.4) A1HN=ANO(ISH(K),JSH(K))
         IF (IJSH(K).EQ.5) A1HN=ANO1(ISH(K),JSH(K))                !srv 26.01.00
c       IF (A1HN.EQ.0.0e0_R8) THEN
         if(
     1   (IJSH(K).EQ.-1.and.jsh(k).eq.1).or.                       !srv 26.01.00
     1   (IJSH(K).EQ.-1.and.jsh(k).eq.2).or.                       !srv 26.01.00
     1   (IJSH(K).EQ.0.and.jsh(k).eq.1).or.
     1   (IJSH(K).EQ.1.and.ish(k).eq.1).or.
     1   (IJSH(K).EQ.3.and.ish(k).eq.n1).or.
     1   (IJSH(K).EQ.5.and.jsh(k).eq.n2-1).or.                     !srv 26.01.00
     1   (IJSH(K).EQ.5.and.jsh(k).eq.n2).or.                       !srv 26.01.00
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
      if(kwrite_bnr.eq.ncall) then                                !srv 18.05.02 {
        write(*,*) 'Opening 83'
        open(83,file='b2_7_rhs_binary_'//trim(name),form='unformatted') !srv 11.03.13
        write(83) n
        do i=1,n
          write(83) r1(i)
        enddo
        close(83)
        write(*,*) 'Opening 82'
        open(82,file='b2_7_matrix_binary_'//trim(name),
     .       form='unformatted')                                  !srv 11.03.13
        write(82) n, nz
        do i=1,nz
          write(82) IRNH(i),ICNH(i),A1H(i)
        enddo
        close(82)
      endif                                                       !srv 18.05.02 }
      if(kwrite_fmt.eq.ncall) then                                !srv 18.05.02 {
        write(*,*) 'Opening 83'
        open(83,file='b2_7_rhs_'//trim(name))
        write(83,*) n
        do i=1,n
          write(83,*) r1(i)
        enddo
        close(83)
        open(82,file='b2_7_matrix_'//trim(name))
        write(82,*) n, nz
        do i=1,nz
          write(82,*) IRNH(i),ICNH(i),A1H(i)
        enddo
        close(82)
      endif                                                    !srv 18.05.02 }
cdpc
      if(nwrite.gt.0) then
        if(nwrite.gt.ncall) then
          if(ncall.eq.0) then
            write(*,*) 'Opening 81'                            !srv 11.01.00
            open(81,file='b2_7_diagonal_matrix')               !srv 11.01.00
          endif
          write(*,*) 'Writing 81'                              !srv 11.01.00
          write(81,*) n,' n'
          do i=1,n
            write(81,*) r1(i)                                  !srv 11.01.00
          enddo
          write(81,*) nz,' nz'                                 !srv 11.01.00
          do i=1,nz
            write(81,*) IRNH(i),ICNH(i),A1H(i)                 !srv 11.01.00
          enddo
        elseif(nwrite.eq.ncall) then
          write(*,*) 'Closing 81'                              !srv 11.01.00
          close(81)                                            !srv 11.01.00
        endif
      endif
cdpc
      if(nltrsol.eq.0) then
      elseif(nltrsol.eq.1) then
         WRITE(6,*) 'ILUTERN'
      else
         WRITE(6,*) 'ILUTERN'
      endif
      if(lfirst) then
         write(*,*) 'ilutern_7: nz = ',nz                      !srv 11.01.00
         call xertst(nz.le.N_NONZEROS, 'Increase N_NONZEROS')
         iflag=0
      endif
      do i=1,nz
         irnh_save(i)=irnh(i)
         icnh_save(i)=icnh(i)
         a1h_save(i)=a1h(i)
      enddo
      if(.not.lfirst) then
         call ma28bdcopy3                                          !srv 26.12.17
     1   (n,nz,a1h,solvdim,irnh,icnh,icnhkeep,ikeep,iw28,w,iflag)
         if(nltrsol.gt.0 .or. iflag.lt.0) write(*,*) 'iflag ',iflag
         if(iflag.lt.0) then
          do i=1,nz
            irnh(i)=irnh_save(i)
            icnh(i)=icnh_save(i)
            a1h(i)=a1h_save(i)
          enddo
         endif
      endif
      if(lfirst.or.iflag.lt.0) then
         call ma28adcopy3                                          !srv 26.12.17
     1   (n,nz,a1h,solvdim,irnh,solvdim1,icnh,u,ikeep,iw28,w,iflag)
         if(nltrsol.gt.0 .or. iflag.lt.0)
     1   write(*,'(a,i6,i8,a,i8,i8,a,i8)')
     .    'iflag,minirn/dim,minicn/dim ',
     2     iflag,minirn,' /',solvdim1,minicn,' /',solvdim
         do i=1,solvdim
          icnhkeep(i)=icnh(i)
         enddo
         lfirst=.false.
      endif
      call ma28idcopy3                                             !srv 26.12.17
     1 (n,nz,a1h_save,irnh_save,icnh_save,solvdim,a1h,icnhkeep,ikeep,
     1 r1,s1,res,w,1,prec,iflag)                                   !srv 16.01.18
      if(nltrsol.gt.0 .or. iflag.lt.0) write(*,'(1x,a,i3,1p,2g15.8,i6)')
     1 'iflag,errmax,dres,noiter ',iflag,errmax,dres,noiter

      if(iflag.lt.0) then                                          !srv 28.06.16 {

        do i=1,nz
          irnh(i) = irnh_save(i)
          icnh(i) = icnh_save(i)
          a1h(i)  = a1h_save(i)
        enddo

        call ma28adcopy3                                           !srv 26.12.17
     1   (n,nz,a1h,solvdim,irnh,solvdim1,icnh,u,ikeep,iw28,w,iflag)
        if(nltrsol.gt.0 .or. iflag.lt.0) write(*,2)
     1     trim(name), iflag, minirn, solvdim1, minicn, solvdim
    2     format(' ma28ad for ',a13,': iflag =',i9/,
     1     ' minirn =',i9,' solvdim1 =',i9,
     2     ' minicn =',i9,' solvdim =',i9)                         !srv 09.07.02
        do i=1,solvdim
          icnhkeep(i)=icnh(i)
        enddo

        call ma28idcopy3                                           !srv 26.12.17
     1   (n,nz,a1h_save,irnh_save,icnh_save,solvdim,a1h,icnhkeep,ikeep,
     1    r1,s1,res,w,1,prec,iflag)                                !srv 16.01.18

        if(iflag.lt.0) then
         write(*,"(a)") 'Problem with solving of '//trim(name)
        else
         write(*,"(3a)") ' Problem of ', trim(name), ' has been solved'
        endif
        if(nltrsol.gt.0 .or. iflag.lt.0)
     &   write(*,3)  trim(name), iflag, errmax, dres, noiter
    3    format(' ma28id for ',a13,': iflag =',i9,
     2          ' errmax =',es9.1,' dres =',es9.1,' noiter =',i2)  !srv 09.07.02 28.06.16 }

      endif                                                        !srv 28.06.16 }

      DO K=1,N
        J=(K-1)/N1+1
        I=K-(J-1)*N1
        SOL(I,J)=S1(K)
      ENDDO

*  ---------------------------------------------------------------------
      call cpu_time(cpu_end)
      ncall=ncall+1
      call subend ()
      RETURN

      entry restart_ma28_for_7                                     !srv 11.01.00
      lfirst=.true.
      write(*,*) 'Restarting ma28_for_7'                           !srv 02.04.18
      return

*=======================================================================
*//END ILUTERN_7//

      END
