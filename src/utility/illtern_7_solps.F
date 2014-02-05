      SUBROUTINE ILLTERN_7(N1, N2, N1D, nz, RHS, SOL)                  !sv 29.05.02 {
      use b2mod_types
	use b2mod_ysmp_sdrv                        
      use b2mod_indirect
      IMPLICIT NONE
C

*  -- INPUT ARGUMENTS --
      INTEGER N1, N2, N1D, NZ
      real (kind=R8) ::
     .       RHS(N1D,N2)

*  -- OUTPUT ARGUMENTS --
      real (kind=R8) ::
     *   SOL(N1D,N2)

*  -- LOCAL SCALARS --
      INTEGER N, I, J, K, l, m, IK
      real (kind=R8) ::
     *   x05baf,cpu_start,cpu_end
      integer iflag

c

*  -- PROCEDURES --
      INTRINSIC MOD
      real (kind=R8) ::
     *   sasum
      EXTERNAL XERRAB, sasum, x05baf

      logical lfirst
      integer ncall, nwrite_fmt, nwrite_bnr                           !sv 29.05.02
	integer       ::  ner, path                                     !sv 29.05.02
	integer, save ::  nltrsol                                       !sv 29.05.02
	real (kind=R8) :: time_begin, time_END                          !sv 29.05.02
      save lfirst, ncall, nwrite_fmt, nwrite_bnr                      !sv 29.05.02
      data lfirst /.true./, ncall /0/, nwrite_fmt /0/, nwrite_bnr /0/ !sv 29.05.02

*=======================================================================
*//COMPUTATION//

*  ---------------------------------------------------------------------
*  TEST INPUTS AND SET AUXILIARIES.
*  ---------------------------------------------------------------------
      if(ncall.eq.0) then
	 call ipgeti ('sdrv_for_7d_nwrite_fmt', nwrite_fmt)          !sv 29.05.02
	 call ipgeti ('sdrv_for_7d_nwrite_bnr', nwrite_bnr)          !sv 29.05.02
	 call ipgeti ('b2ux7p_nltrsol', nltrsol)                     !sv 29.05.02
      endif
c     write(*,*) 'ncall,nwrite',ncall,nwrite

      cpu_start=x05baf()
	iflag = 0                                                  !sv 16.04.02
      if(ncall.lt.10) then
*  -- CHECK N1, N2, N1D, ITERS, ACPAR, ITPAR --
	 IF (N1 .LT. 2 .OR. N2 .LT. 2 .OR. N1D .LT. N1) THEN
	  CALL XERRAB ('ILUTERN_7 -- FAULTY ARGUMENT N1, N2, N1D') !sv 06.07.00
	 END IF

      endif
      N = N1*N2
      DO I=1,N1
	 DO J=1,N2
	  R1(I+(J-1)*N1)=RHS(I,J)
	  S1(I+(J-1)*N1)=0.0e0_R8
	 enddo
      enddo
c-dpc
      if(nwrite_fmt.gt.0) then
	if(nwrite_fmt.gt.ncall) then
	  if(ncall.eq.0) then
	    write(*,*) 'Opening 78'                            !sv 11.01.00
	    write(*,*) 'Opening 79'                            !sv 11.01.00
	    open(78,file='b2_7_diag_symm_matrix')              !sv 11.01.00
	    open(79,file='b2_7_rhs')               !sv 11.01.00
	  endif
	  write(*,*) 'Writing 78'                              !sv 11.01.00
	  write(*,*) 'Writing 79'                              !sv 11.01.00
	  write(79,*) n,' n'
	  do i=1,n
	    write(79,*) r1(i)                                  !sv 11.01.00
	  enddo
	  write(78,*) ' nz =',nz
	  do k=1,n+1
	   write(78,*) k, ja(k), a(k), ia(k)
	  enddo
	  do k=n+2,nz
	   write(78,*) k, ja(k), a(k)
	  enddo
	elseif(nwrite_fmt.eq.ncall) then
	  write(*,*) 'Closing 78'                              !sv 11.01.00
	  close(78)                                            !sv 11.01.00
	  write(*,*) 'Closing 79'                              !sv 11.01.00
	  close(79)                                            !sv 11.01.00
	endif
      endif
      if(nwrite_bnr.gt.0) then
	if(nwrite_bnr.gt.ncall) then
	 if(ncall.eq.0) then
	  write(*,*) 'Opening 78'                            !sv 11.01.00
	  write(*,*) 'Opening 79'                            !sv 11.01.00
        open (78,file='b2_7_diag_symm_matrix_binary',form='unformatted') !srv 11.03.13
        open (79,file='b2_7_rhs_binary',form='unformatted')              !srv 11.03.13
	 endif
	  write(*,*) 'Writing 78'                              !sv 11.01.00
	  do i=1,n
	    write(79) r1(i)                                    !sv 11.01.00
	  enddo
	  write(78) ' kk =',nz
	  do k=1,n+1
	   write(78) k, ja(k), a(k), ia(k)
	  enddo
	  do k=n+2,nz
	   write(78) k, ja(k), a(k)
	  enddo
	elseif(nwrite_bnr.eq.ncall) then
	  write(*,*) 'Closing 78'                              !sv 11.01.00
	  write(*,*) 'Closing 79'                              !sv 11.01.00
	  close(78)                                            !sv 11.01.00
	  close(79)                                            !sv 11.01.00
	endif
      endif

      if(lfirst) then
C
        PATH = 2
!	  call cpu_time(time_begin)
        CALL  ODRV
     *     (N, IA,JA,A, P,IP, NSP,ISP, PATH, IFLAG)                      
C
!        call cpu_time(time_end)
!	  write (*,"(' Time for ODRV is',1pd9.2)") time_end - time_begin
	if(IFLAG.ne.0) then                                             !sv 13.02.03 {
       write(*,"(' ODRV ',': iflag =',i9,'  NSP =',I9)") iflag, NSP
	 stop
	endif                                                           !sv 13.02.03 }
C
	if(nltrsol.gt.0)                                                   !sv 09.07.02
     *    write(*,"(' ODRV ',': iflag =',i9,'  NSP =',I9)") iflag, NSP   !sv 09.07.02
C
        PATH = 4
!	  call cpu_time(time_begin)
C                             Perform symbolyc LDL' factorization (SSF)
        CALL SDRV
     *     (N, P,IP, IA,JA,A, r1, s1, NSP,ISP,RSP,ESP, PATH, IFLAG)
C
!        call cpu_time(time_end)

!	write (*,"(' Time for SDRV (SSF) is',1pd9.2)") time_end-time_begin
	if(IFLAG.ne.0) then                                             !sv 13.02.03 {
       write (*,"(' SDRV (SSF): ESP =',I9,' FLAG =',I9,'  NSP =',I9)")
     *       ESP, Iflag, NSP
	 stop
	endif                                                           !sv 13.02.03 }
	if(nltrsol.gt.0)                                                   !sv 09.07.02
     * write (*,"(' SDRV (SSF): ESP =',I9,' FLAG =',I9,'  NSP =',I9)")
     *       ESP, Iflag, NSP
	  lfirst=.false.
	 endif
C
       if(.not.lfirst) then
C
        PATH = 3
!	  call cpu_time(time_begin)
C       REORDER SYMMETRICALLY STORED MATRIX AS SPECIFIED BY INPUT PERMUTATION
        CALL  ODRV
     *     (N, IA,JA,A, P,IP, NSP,ISP, PATH, IFLAG)                      
C
!        call cpu_time(time_end)
C
        PATH = 6
!	  call cpu_time(time_begin)
C                             Perform numerical LDL' factorization
         CALL SDRV
     *     (N, P,IP, IA,JA,A, r1, s1, NSP,ISP,RSP,ESP, PATH, IFLAG)
C
!        call cpu_time(time_end)

!	write (6,"(' Time for SDRV (SNF) is',1pd9.2)") time_end-time_begin
	if(IFLAG.ne.0) then                                             !sv 13.02.03 {
       write (*,"(' SDRV (SNF): ESP =',I9,' FLAG =',I9,'  NSP =',I9)")
     *       ESP, Iflag, NSP
	 stop
	endif                                                           !sv 13.02.03 }
	if(nltrsol.gt.0)                                                !sv 09.07.02
     * write (*,"(' SDRV (SNF): ESP =',I9,' FLAG =',I9,'  NSP =',I9)")
     *       ESP, Iflag, NSP
C
        PATH = 3
!	  call cpu_time(time_begin)
C                                 Solve A * X = B
         CALL SDRV
     *     (N, P,IP, IA,JA,A, r1, s1, NSP,ISP,RSP,ESP, PATH, IFLAG)
C
!        call cpu_time(time_end)

!	write (6,"(' Time for SDRV (SNS) is',1pd9.2)") time_end-time_begin
	if(iflag.ne.0) then                                             !sv 13.02.03 {
       write (*,"(' SDRV (SNS): ESP =',I9,' FLAG =',I9,'  NSP =',I9)")
     *       ESP, Iflag, NSP
	 stop
	endif                                                           !sv 13.02.03 }
	if(nltrsol.gt.0)                                                   !sv 09.07.02
     * write (*,"(' SDRV (SNS): ESP =',I9,' FLAG =',I9,'  NSP =',I9)")
     *       ESP, Iflag, NSP
C
       endif
	  
	 
      DO 50,K=1,N
        J=(K-1)/N1+1
        I=K-(J-1)*N1
        SOL(I,J)=S1(K)
  50  CONTINUE

*  ---------------------------------------------------------------------
      cpu_end=x05baf()
      ncall=ncall+1
      RETURN

      entry restart_sdrv_for_7                             
      lfirst=.true.
      write(*,*) ' Restarting sdrv'
      return

*=======================================================================
*//END ILLTERN//

      END                                                          !sv 29.05.02 }
