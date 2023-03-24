      SUBROUTINE ILUTERN_US(NCV, NZ, RHS, SOL, name)      !srv 13.10.06 (
      use b2mod_types
      use b2mod_ma28_for_us
#ifdef _OPENMP
      use b2mod_openmp
#endif
      IMPLICIT NONE
C

*  -- INPUT ARGUMENTS --
      INTEGER NCV, NZ
      REAL (KIND=R8) :: RHS(NCV)
      character*(*) name                                         !srv 18.05.02

*  -- OUTPUT ARGUMENTS --
      REAL (KIND=R8) ::
     *   SOL(NCV)

*  -- LOCAL SCALARS --
      INTEGER N, I
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
c

*  -- PROCEDURES --
#ifndef _OPENMP
      INTRINSIC cpu_time
#else
#ifdef CRAY
      real (kind=R8) :: omp_get_wtime
      external omp_get_wtime
#endif
#endif
      EXTERNAL XERRAB

      integer ncall, nwrite
      integer, save :: kwrite_fmt = -1                                !srv 18.05.02
      integer, save :: kwrite_bnr = -1                                !srv 18.05.02
      save ncall, nwrite
      data ncall /0/, nwrite /0/

*=======================================================================
*//COMPUTATION//

*  ---------------------------------------------------------------------
*  TEST INPUTS AND SET AUXILIARIES.
*  ---------------------------------------------------------------------

      if(ncall.eq.0) then
        call ipgeti ('ma28_for_us_nwrite', nwrite)
        call ipgeti ('ma28_for_us_kwrite_fmt', kwrite_fmt)            !srv 18.05.02
        call ipgeti ('ma28_for_us_kwrite_bnr', kwrite_bnr)            !srv 18.05.02
      endif
c     write(*,*) 'ncall,nwrite',ncall,nwrite

#ifdef _OPENMP
      cpu_start = omp_get_wtime()
#else
      call cpu_time(cpu_start)
#endif
      iflag=0
      if(ncall.lt.1000) then
*  -- CHECK NCV, ITERS, ACPAR, ITPAR --
        IF (NCV .LT. 1) THEN
          CALL XERRAB ('ILUTERN_US -- FAULTY ARGUMENT NCV')
        END IF

      endif
      
      N  = NCV
      RR1 = RHS
      S1 = 0.0_R8

      if(kwrite_bnr.eq.ncall) then                                  !srv 18.05.02 {
         write(*,*) 'Opening 76'
         open(76,file='b2_us_rhs_binary_'//trim(name),
     &     form='unformatted') !srv 11.03.13
         write(76) n
         do i=1,nCv
           write(76) rr1(i)
         enddo
         close(76)
         write(*,*) 'Opening 77'
         open(77,file='b2_us_matrix_binary_'//
     &       trim(name),form='unformatted')                         !srv 11.03.13
         write(77) n, nz
         do i=1,nz
           write(77) IRNH(i),ICNH(i),A1H(i)
         enddo
         close(77)
      endif                                                         !srv 18.05.02 }
      if(kwrite_fmt.eq.ncall) then                                  !srv 18.05.02 {
        write(*,*) 'Opening 76'
        open(76,file='b2_us_rhs_'//trim(name))
        write(76,*) n
        do i=1,n
          write(76,*) rr1(i)
        enddo
        close(76)
        open(77,file='b2_us_matrix_'//trim(name))
        write(77,*) n, nz
        do i=1,nz
          write(77,*) IRNH(i),ICNH(i),A1H(i)
        enddo
        close(77)
      endif                                                         !srv 18.05.02 }
cdpc
      if(nwrite.gt.0) then
        if(nwrite.gt.ncall) then
          if(ncall.eq.0) then
            write(*,*) 'Opening 79'
            open(79,file='b2_us_matrix')
          endif
          write(*,*) 'Writing 79'
          write(79,*) 'Called from '//trim(name)
          write(79,*) n,' n'
          do i=1,n
            write(79,*) rr1(i)
          enddo
          write(79,*) nz,' nz'
          do i=1,nz
            write(79,*) IRNH(i),ICNH(i),A1H(i)
          enddo
        elseif(nwrite.eq.ncall) then
           write(*,*) 'Closing 79'
           close(79)
        endif
      endif
cdpc
      if(lfirst_ilu_us) then
        write(*,*) 'ilutern_us: nz = ',nz
        call xertst (nz.le.N_NONZEROS, 'Increase N_NONZEROS')
      endif

      do i=1,nz
        irnh_save(i)=irnh(i)
        icnh_save(i)=icnh(i)
        a1h_save(i)=a1h(i)
      enddo

      if(.not.lfirst_ilu_us) then
        call ma28bdcopy
     1   (n,nz,a1h,solvdim,irnh,icnh,icnhkeep,ikeep,iw28,w,iflag)
        if(nltrsol.gt.0 .or. iflag.lt.0)                             !srv 09.07.02 28.06.16
     *   write(*,"(' ma28bd for ',a13,': iflag =',i9)") 
     *    trim(name), iflag                                          !srv 09.07.02 28.06.16
        if(iflag.lt.0) then
          do i=1,nz
            irnh(i)=irnh_save(i)
            icnh(i)=icnh_save(i)
            a1h(i)=a1h_save(i)
          enddo
        endif
      endif

      if(lfirst_ilu_us.or.iflag.lt.0) then
        call ma28adcopy
     1   (n,nz,a1h,solvdim,irnh,solvdim1,icnh,u,ikeep,iw28,w,iflag)
        if(nltrsol.gt.0 .or. iflag.lt.0) write(*,2)                  !srv 09.07.02 28.06.16
     1     trim(name), iflag, minirn, solvdim1, minicn, solvdim      !srv 09.07.02 28.06.16
   2  format(' ma28ad for ',a13,': iflag =',i9/,
     1 ' minirn =',i9,' solvdim1 =',i9,
     2 ' minicn =',i9,' solvdim =',i9)                               !srv 09.07.02
        do i=1,solvdim
          icnhkeep(i)=icnh(i)
        enddo
        lfirst_ilu_us=.false. 
      endif

      call ma28idcopy
     1 (n,nz,a1h_save,irnh_save,icnh_save,solvdim,a1h,icnhkeep,ikeep,
     1 rr1,s1,res,w,1,1e-10_R8,iflag)
      if(nltrsol.gt.0 .or. iflag.lt.0)                               !srv 09.07.02 28.06.16 {
     &   write(*,3)  trim(name), iflag, errmax, dres, noiter
   3  format(' ma28id for ',a13,
     1       ': iflag =',i9,' errmax =',es9.1,
     2       ' dres =',es9.1,' noiter =',i2)                         !srv 09.07.02 28.06.16 }

      if(iflag.lt.0) then                                            !srv 28.06.16 {

        do i=1,nz
          irnh(i) = irnh_save(i)
          icnh(i) = icnh_save(i)
          a1h(i)  = a1h_save(i)
        enddo

        call ma28adcopy
     1   (n,nz,a1h,solvdim,irnh,solvdim1,icnh,u,ikeep,iw28,w,iflag)
        if(nltrsol.gt.0 .or. iflag.lt.0) write(*,2)
     1     trim(name), iflag, minirn, solvdim1, minicn, solvdim
        do i=1,solvdim
          icnhkeep(i)=icnh(i)
        enddo

        call ma28idcopy
     1 (n,nz,a1h_save,irnh_save,icnh_save,solvdim,a1h,icnhkeep,ikeep,
     1 rr1,s1,res,w,1,1e-10_R8,iflag)

        if(iflag.lt.0) then
         write(*,"(a)") 'Problem with solving '//trim(name)
        else
         write(*,"(3a)") ' Problem of ', trim(name), ' has been solved'
        endif
        if(nltrsol.gt.0 .or. iflag.lt.0)
     &   write(*,3)  trim(name), iflag, errmax, dres, noiter

      endif                                                          !srv 28.06.16 }

      sol = s1

*  ---------------------------------------------------------------------
#ifdef _OPENMP
      cpu_end = omp_get_wtime()
#else
      call cpu_time(cpu_end)
#endif
      ncall=ncall+1
      RETURN
*=======================================================================
*//END ILUTERN_US//                                                    !srv 13.10.06

      END SUBROUTINE ILUTERN_US

