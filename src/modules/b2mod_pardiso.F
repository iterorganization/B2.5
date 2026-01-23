      module b2mod_PARDISO
      use b2mod_types
      use b2mod_ma28_for_us !! the storage of the matrix in triplet format is in this module....perhaps move later on
      implicit none

      contains 


      subroutine SolvePARDISO(nCv, nz, b, sol, name, flag)

      use, intrinsic :: iso_c_binding, only: c_int
      implicit none

      !--------------------------------------
      ! TODO. At list of optional todos: 
      ! 1. The initialization phase can be moved to a separate routine and the parameters can be saved in the module. 
      ! 2. The solver settings can perhaps be further optimized for B2.5 matrices.
      ! 3. The shape of the stencils for each equation remains unchanged throughout the equations. It might be possible to save and store some 
      !    pre-computed info per equation to prevent re-computing it at every equation. 

      !--------------------------------------
      ! Arguments
      integer, intent(in) :: nCv, nz
      character*(*) :: name
      real(R8), intent(in) :: b(nCv)       !! right-hand side
      real(R8), intent(out) :: sol(nCv)    !! solution
      integer, intent(out) :: flag

      !--------------------------------------
      ! Internal MKL PARDISO variables
      integer :: ia_t(nz), ja_t(nz)        !! ia_t, ja_t, a_t = matrix in triplet format
      real(R8) :: a_t(nz)
      integer(R8) :: pt(64)  ! internal solver memory pointer
      integer :: iparm(64)
      integer :: maxfct, mnum, mtype, phase, nrhs
      integer :: error
      integer :: i
      integer :: ia(nCv+1), ja(nz)         !! ia, ja, a = matrix in CSR format, needed by PARDISO
      real(R8), dimension(:) :: a(nz)

#ifdef PARDISO
      !--------------------------------------
      ! Read in matrix from b2mod_ma28 (storage can perhaps be moved in the future)
      ia_t = irnh(1:nz)
      ja_t = icnh(1:nz)
      a_t = a1h(1:nz)


      !--------------------------------------
      ! Convert triplet to CSR format
      call TripletToCSR(nCv, nz, ia_t, ja_t, a_t, ia, ja, a)

      !--------------------------------------
      ! Initialization
      pt = 0_8
      iparm = 0
      maxfct = 1
      mnum = 1
      mtype = 11       ! Real unsymmetric
      phase = 0
      nrhs = 1
      error = 0

      !--------------------------------------
      ! Initialize default iparm
      iparm(1) = 1       ! use default parameters
      iparm(2) = 2       ! Fill-in reordering from METIS (optional)
      iparm(3) = 0
      iparm(4) = 0
      iparm(5) = 0
      iparm(6) = 0
      iparm(7) = 0
      iparm(8) = 0
      iparm(9) = 13      ! perturb the pivot elements
      iparm(10) = 1      ! scaling
      iparm(11) = 0
      iparm(12) = 0
      iparm(13) = 0
      iparm(14) = 0
      iparm(15) = 0
      iparm(16) = 0
      iparm(17) = -1     ! output: number of nonzeros in factor
      iparm(18) = -1     ! output: Mflops
      iparm(19) = 0
      iparm(20) = 0

      !--------------------------------------
      ! Override some default settings to match ma28 almost exactly
      iparm(10) = 13   ! pivoting strategy
      iparm(11) = 1    ! scaling
      iparm(13) = 1    ! matching-based ordering

      !--------------------------------------
      ! Step 1: Analysis
      phase = 11
      call pardiso(pt, maxfct, mnum, mtype, phase, nCv, a, ia, ja, 
     .          0, nrhs, iparm, 0, b, sol, error)

      if (error /= 0) then
        print *, 'SolvePARDISO: Analysis failed with error =', error
        flag = error
        return
      endif

      !--------------------------------------
      ! Step 2: Numerical factorization
      phase = 22
      call pardiso(pt, maxfct, mnum, mtype, phase, nCv, a, ia, ja, 
     .          0, nrhs, iparm, 0, b, sol, error)

      if (error /= 0) then
        print *, 'SolvePARDISO: Factorization failed with error =',error
        flag = error
        return
      endif

      !--------------------------------------
      ! Step 3: Solve
      phase = 33
      call pardiso(pt, maxfct, mnum, mtype, phase, nCv, a, ia, ja, 
     .          0, nrhs, iparm, 0, b, sol, error)

      if (error /= 0) then
        print *, 'SolvePARDISO: Solve failed with error =', error
        flag = error
        return
      endif

      flag = 0

      !--------------------------------------
      ! Step 4: Release internal memory
      phase = -1
      call pardiso(pt, maxfct, mnum, mtype, phase, nCv, a, ia, ja, 
     .          0, nrhs, iparm, 0, b, sol, error)
#else
      call xerrab('PARDISO not available! Link and re-compile with'//
     . ' PARDISO library')
#endif

      end subroutine SolvePARDISO


      subroutine TripletToCSR(n, nz, irn, jcn, a_triplet, ia, ja,
     .                         a_csr)

      implicit none
      integer, intent(in) :: n       ! Number of rows
      integer, intent(in) :: nz      ! Number of nonzeros
      integer, intent(in) :: irn(nz) ! Row indices (1-based)
      integer, intent(in) :: jcn(nz) ! Column indices (1-based)
      real(R8), intent(in) :: a_triplet(nz) ! Values

      integer, intent(out) :: ia(n+1) ! Row pointers (CSR)
      integer, intent(out) :: ja(nz)  ! Column indices (CSR)
      real(R8), intent(out) :: a_csr(nz) ! Values in CSR

      ! Local variables
      integer :: i, row
      integer, allocatable :: row_counts(:)
      integer :: ptr, pos

      !-----------------------------------
      ! Step 1: count number of nonzeros per row
      allocate(row_counts(n))
      row_counts = 0

      do i = 1, nz
        row = irn(i)
        row_counts(row) = row_counts(row) + 1
      end do

      !-----------------------------------
      ! Step 2: build row pointer array
      ia(1) = 1
      do i = 1, n
        ia(i+1) = ia(i) + row_counts(i)
      end do

      ! Step 3: temporary position array to track fill
      row_counts = 0

      !-----------------------------------
      ! Step 4: fill ja and a_csr arrays
      do i = 1, nz
       row = irn(i)
       ptr = ia(row) + row_counts(row)
       ja(ptr) = jcn(i)
       a_csr(ptr) = a_triplet(i)
       row_counts(row) = row_counts(row) + 1
      end do

      !-----------------------------------
      deallocate(row_counts)
      end subroutine TripletToCSR



      end module b2mod_PARDISO
