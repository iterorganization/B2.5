      module b2mod_ilutern_9
      implicit none

      contains

      SUBROUTINE ILUTERN_9(N1, N2, N1D, NZ, RHS, SOL, name)      !srv 13.10.06 (
      use b2mod_types
      use b2mod_ma28_for_9diag
      use b2mod_indirect
      use b2mod_openmp
      IMPLICIT NONE
#ifndef GFORTRAN
      external bd1copy
#endif
C

*  -- INPUT ARGUMENTS --
      INTEGER N1, N2, N1D, NZ
      REAL (KIND=R8) :: RHS(N1D,N2)
      character*(*) name                                         !srv 18.05.02

*  -- OUTPUT ARGUMENTS --
      REAL (KIND=R8) ::
     *   SOL(N1D,N2)

*  -- LOCAL SCALARS --
      INTEGER N, I, J, K
#ifdef _OPENMP
      real (kind=R8) :: cpu_start, cpu_end
#else
      real (kind=R4) :: cpu_start, cpu_end
#endif
      integer iflag

      real (kind=R8) ::
     *   u
      parameter (u=1.0_R8)

c ma28 common blocks
      real (kind=R8) ::
     *   eps,rmin,resid
      integer irncp,icncp,minirn,minicn,irank
      logical abort1,abort2
      common /ma28fdc/ eps,rmin,resid,irncp,icncp,minirn,minicn,
     1 irank,abort1,abort2
      real (kind=R8) ::
     *   rtol,themax,big,dxmax,errmax,dres,cgce
      integer ndrop,maxit,noiter,nsrch,istart,ibig
      common /ma28hdc/ rtol,themax,big,dxmax,errmax,dres,cgce,
     1 ndrop,maxit,noiter,nsrch,istart,ibig
#ifdef _OPENMP
!$OMP THREADPRIVATE(/ma28fdc/, /ma28hdc/)
#endif
c

*  -- PROCEDURES --
      EXTERNAL XERRAB
#if ( defined(_OPENMP) && defined(CRAY) )
      real (kind=R8) :: omp_get_wtime
      external omp_get_wtime
#endif

      integer ncall, nwrite
      integer, save :: kwrite_fmt = -1                                !srv 18.05.02
      integer, save :: kwrite_bnr = -1                                !srv 18.05.02
      save ncall, nwrite
#ifdef _OPENMP
!$OMP THREADPRIVATE(ncall, nwrite)
!$OMP THREADPRIVATE(kwrite_fmt, kwrite_bnr)
#endif
      data ncall /0/, nwrite /0/

*=======================================================================
*//COMPUTATION//

*  ---------------------------------------------------------------------
*  TEST INPUTS AND SET AUXILIARIES.
*  ---------------------------------------------------------------------

      if(ncall.eq.0) then
        call ipgeti ('ma28_for_9d_nwrite', nwrite)
        call ipgeti ('ma28_for_9d_kwrite_fmt', kwrite_fmt)            !srv 18.05.02
        call ipgeti ('ma28_for_9d_kwrite_bnr', kwrite_bnr)            !srv 18.05.02
        if (in_parallel() .and. ( nwrite.ne.0 .or.
     &     kwrite_bnr.ne.-1 .or. kwrite_fmt.ne.-1)) then
#ifdef _OPENMP
!$OMP MASTER
#endif
         write(*,*) 'ilutern_9 OpenMP warning: no file IO ' ,
     &       'in parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
         nwrite = 0
         kwrite_fmt = -1
         kwrite_bnr = -1
       endif
#ifdef GFORTRAN
       call set_data_bd1copy
       call set_data_bd2copy
       call set_data_bd3copy
#endif
      endif
c     write(*,*) 'ncall,nwrite',ncall,nwrite

#ifdef _OPENMP
      cpu_start = omp_get_wtime()
#else
      call cpu_time(cpu_start)
#endif
      iflag=0
      if(ncall.lt.1000) then
*  -- CHECK N1, N2, N1D, ITERS, ACPAR, ITPAR --
        IF (N1 .LT. 2 .OR. N2 .LT. 2 .OR. N1D .LT. N1) THEN
          CALL XERRAB ('ILUTERN_9 -- FAULTY ARGUMENT N1, N2, N1D')
        END IF

      endif
      N = N1*N2
      DO I = 1, N1
        DO J = 1, N2
         R1_9p(I+(J-1)*N1) = RHS(I,J)
         S1_9p(I+(J-1)*N1) = 0.0_R8
        enddo
      enddo

      if(kwrite_bnr.eq.ncall) then                                  !srv 18.05.02 {
         write(*,*) 'Opening 76'
         open(76,file='b2_9_rhs_binary_'//trim(name),form='unformatted') !srv 11.03.13
         write(76) n
         do i=1,n
           write(76) r1_9p(i)
         enddo
         close(76)
         write(*,*) 'Opening 77'
         open(77,file='b2_9_matrix_binary_'//
     &       trim(name),form='unformatted')                         !srv 11.03.13
         write(77) n, nz
         do i=1,nz
           write(77) IRNH_9p(i),ICNH_9p(i),A1H_9p(i)
         enddo
         close(77)
      endif                                                         !srv 18.05.02 }
      if(kwrite_fmt.eq.ncall) then                                  !srv 18.05.02 {
        write(*,*) 'Opening 76'
        open(76,file='b2_9_rhs_'//trim(name))
        write(76,*) n
        do i=1,n
          write(76,*) r1_9p(i)
        enddo
        close(76)
        open(77,file='b2_9_matrix_'//trim(name))
        write(77,*) n, nz
        do i=1,nz
          write(77,*) IRNH_9p(i),ICNH_9p(i),A1H_9p(i)
        enddo
        close(77)
      endif                                                         !srv 18.05.02 }
cdpc
      if(nwrite.gt.0) then
        if(nwrite.gt.ncall) then
          if(ncall.eq.0) then
            write(*,*) 'Opening 79'
            open(79,file='b2_9_matrix')
          endif
          write(*,*) 'Writing 79'
          write(79,*) 'Called from '//trim(name)
          write(79,*) n,' n'
          do i=1,n
            write(79,*) r1_9p(i)
          enddo
          write(79,*) nz,' nz'
          do i=1,nz
            write(79,*) IRNH_9p(i),ICNH_9p(i),A1H_9p(i)
          enddo
        elseif(nwrite.eq.ncall) then
           write(*,*) 'Closing 79'
           close(79)
        endif
      endif
cdpc
      if(lfirst_ilu9) then
        write(*,*) 'ilutern_9: nz = ',nz
        call xertst (nz.le.N_NONZEROS_9p, 'Increase N_NONZEROS')
      endif

      do i=1,nz
        irnh_save_9p(i)=irnh_9p(i)
        icnh_save_9p(i)=icnh_9p(i)
        a1h_save_9p(i)=a1h_9p(i)
      enddo

      if(.not.lfirst_ilu9) then
        call ma28bdcopy
     1   (n,nz,a1h_9p,solvdim_9p,
     1    irnh_9p,icnh_9p,icnhkeep_9p,ikeep_9p,iw28_9p,w_9p,iflag)
        if(nltrsol_9p.gt.0 .or. iflag.lt.0)                          !srv 09.07.02 28.06.16
     *   write(*,"(' ma28bd for ',a13,': iflag =',i9)")
     *    trim(name), iflag                                          !srv 09.07.02 28.06.16
        if(iflag.lt.0) then
          do i=1,nz
            irnh_9p(i)=irnh_save_9p(i)
            icnh_9p(i)=icnh_save_9p(i)
            a1h_9p(i)=a1h_save_9p(i)
          enddo
        endif
      endif

      if(lfirst_ilu9.or.iflag.lt.0) then
        call ma28adcopy
     1   (n,nz,a1h_9p,solvdim_9p,irnh_9p,solvdim1_9p,icnh_9p,
     1    u,ikeep_9p,iw28_9p,w_9p,iflag)
        if(nltrsol_9p.gt.0 .or. iflag.lt.0) write(*,2)               !srv 09.07.02 28.06.16
     1     trim(name), iflag, minirn,
     1     solvdim1_9p, minicn, solvdim_9p                           !srv 09.07.02 28.06.16
    2  format(' ma28ad for ',a13,': iflag =',i9/,
     1 ' minirn =',i9,' solvdim1 =',i9,
     2 ' minicn =',i9,' solvdim =',i9)                               !srv 09.07.02
        do i=1,solvdim_9p
          icnhkeep_9p(i)=icnh_9p(i)
        enddo
        lfirst_ilu9=.false.
      endif

      call ma28idcopy
     1 (n,nz,a1h_save_9p,irnh_save_9p,icnh_save_9p,
     1  solvdim_9p,a1h_9p,icnhkeep_9p,ikeep_9p,
     1  r1_9p,s1_9p,res_9p,w_9p,1,1e-10_R8,iflag)
      if(nltrsol_9p.gt.0 .or. iflag.lt.0)                            !srv 09.07.02 28.06.16 {
     &   write(*,3)  trim(name), iflag, errmax, dres, noiter
    3    format(' ma28id for ',a13,
     1       ': iflag =',i9,' errmax =',es9.1,
     2       ' dres =',es9.1,' noiter =',i2)                         !srv 09.07.02 28.06.16 }

      if(iflag.lt.0) then                                            !srv 28.06.16 {

        do i=1,nz
          irnh_9p(i) = irnh_save_9p(i)
          icnh_9p(i) = icnh_save_9p(i)
          a1h_9p(i)  = a1h_save_9p(i)
        enddo

        call ma28adcopy
     1   (n,nz,a1h_9p,solvdim_9p,irnh_9p,solvdim1_9p,icnh_9p,
     1    u,ikeep_9p,iw28_9p,w_9p,iflag)
        if(nltrsol_9p.gt.0 .or. iflag.lt.0) write(*,2)
     1     trim(name), iflag, minirn, solvdim1_9p, minicn, solvdim_9p
        do i=1,solvdim_9p
          icnhkeep_9p(i)=icnh_9p(i)
        enddo

        call ma28idcopy
     1 (n,nz,a1h_save_9p,irnh_save_9p,icnh_save_9p,
     1  solvdim_9p,a1h_9p,icnhkeep_9p,ikeep_9p,
     1  r1_9p,s1_9p,res_9p,w_9p,1,1e-10_R8,iflag)

        if(iflag.lt.0) then
         write(*,"(a)") 'Problem with solving '//trim(name)
        else
         write(*,"(3a)") ' Problem of ', trim(name), ' has been solved'
        endif
        if(nltrsol_9p.gt.0 .or. iflag.lt.0)
     &   write(*,3)  trim(name), iflag, errmax, dres, noiter

      endif                                                          !srv 28.06.16 }

      do  k=1,n
        j=(k-1)/n1+1
        i=k-(j-1)*n1
        sol(i,j)=s1_9p(k)
      enddo

*  ---------------------------------------------------------------------
#ifdef _OPENMP
      cpu_end = omp_get_wtime()
#else
      call cpu_time(cpu_end)
#endif
      ncall=ncall+1
      RETURN
      END SUBROUTINE ILUTERN_9

*=======================================================================
*//END ILUTERN_9//                                                    !srv 13.10.06

      END MODULE B2MOD_ILUTERN_9
