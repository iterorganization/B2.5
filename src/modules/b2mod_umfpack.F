      module b2mod_umfpack
      use b2mod_types
      use b2mod_ma28_for_us
      use, intrinsic :: iso_c_binding     !! binding with C
      implicit none

      ! UMFPACK variables - hard coded here...
      integer(c_int) :: umfpack_a = 0
      integer(c_int) :: umfpack_control = 20
      integer(c_int) :: umfpack_info = 90


      ! UMFPACK double integer routines
      ! ================================
      ! Note: the final authority of the implementation on the c-side of
      ! these routines lies with the manual provided by Dr. Tim Davis. 
      ! See the github repository 
      ! (https://github.com/DrTimothyAldenDavis/SuiteSparse, and the 
      ! UMFPACK package thereof) and the manual included there. Below,
      ! we provide a brief description in the interface of the expected
      ! arguments for each function. 
    
      interface 
        ! Defaults
        subroutine UmfpackDefaultsDI(control)  
     &       bind(c, name='umfpack_di_defaults')

            ! Description
            !============
            ! This routine simply loads in the default values for the 
            ! control variables used. Normally, there is no need to 
            ! change these. If this is the case, refer to the manual. 
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            real(c_double)      :: control(*)

        end subroutine

        ! Symbolic creator
        subroutine UmfpackSymbolicDI(nrow, ncol, Ap, Ai, Ax, symbolic, 
     &       control, info)  
     &       bind(c, name='umfpack_di_symbolic')

            ! Description
            !============
            ! This routine loads in the symbolic factorization, which is
            ! a sort of preprocessing step for the numerical 
            ! factorization. The arguments to be parsed are the number 
            ! of rows and columns of the matrix, the sparse matrix in 
            ! compressed column storage format (Ap is the column 
            ! pointer, Ai the row indices, and Ax the values), and the
            ! control and info arrays. The 'symbolic' argument here 
            ! should be a pointer, which, on exit, points to the 
            ! symbolic object used in the numerical counterpart of this
            ! routine. 
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Import
            import umfpack_control, umfpack_info

            ! Declare
            integer(c_int), value           :: nrow, ncol
            integer(c_int), intent(in)      :: Ap(*), Ai(*)
            real(c_double)                  :: 
     &           control(umfpack_control), info(umfpack_info)
            real(c_double), intent(in)      :: Ax(*)
            type(c_ptr)                     :: symbolic 

        end subroutine

        ! Symbolic destructor
        subroutine UmfpackFreeSymbolicDI(symbolic)  
     &       bind(c, name='umfpack_di_free_symbolic')

            ! Description
            !============
            ! This routine destroys the symbolic object
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            type(c_ptr)                     :: symbolic 

        end subroutine

        ! Numeric creator
        subroutine UmfpackNumericDI(Ap, Ai, Ax, symbolic, numeric, 
     &       control, info) 
     &       bind(c, name='umfpack_di_numeric')

            ! Description
            !============
            ! Constructs the object to factorize the matrix numerically,
            ! which is used in the solver. The input arguments are 
            ! largely the same as in the symbolic counterpart of this 
            ! routine, though the symbolic argument is now passed by
            ! value, whereas the numeric argument is passed as a 
            ! pointer. 

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int)                  :: Ap(*), Ai(*)
            real(c_double)                  :: control(*), info(*)
            real(c_double)                  :: Ax(*)
            type(c_ptr)                     :: numeric
            type(c_ptr), value              :: symbolic

        end subroutine

        ! Numeric destructor
        subroutine UmfpackFreeNumericDI(numeric)  
     &       bind(c, name='umfpack_di_free_numeric')

            ! Description
            !============
            ! Destroys the numeric object
            
            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            type(c_ptr)                     :: numeric 

        end subroutine

        ! Solver
        subroutine UmfpackSolveDI(sys, Ap, Ai, Ax, x, b, numeric, 
     &       control, info) 
     &       bind(c, name='umfpack_di_solve')

            ! Description
            !============
            ! General solver for a system ( Ax = b). Note that the 'A'
            ! in this system may be manipulated by adjusting the 'sys'
            ! argument. To solve a regular system, this should equal 
            ! UMFPACK_A (which is currently 0). This UMFPACK_A is a 
            ! variable that should be supplied by UMFPACK.

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int)                  :: Ap(*), Ai(*)
            integer(c_int), value           :: sys
            real(c_double)                  :: control(*), info(*)
            real(c_double)                  :: Ax(*), x(*), b(*)
            type(c_ptr), value              :: numeric

        end subroutine

        ! Info reporter - not used
        subroutine UmfpackReportInfoDI(control, info) 
     &       bind(c, name='umfpack_di_report_info')

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            real(c_double)                     :: info(*), control(*)
            
        end subroutine

        ! Coordinate to compressed column storage format
        subroutine UmfpackTriplet2ColDI(nrow, ncol, nz, Ti, Tj, Tx, 
     &       Ap, Ai, Ax, Map) 
     &       bind(c, name='umfpack_di_triplet_to_col')

            ! C module
            use, intrinsic :: iso_c_binding 

            ! Declare
            integer(c_int), value           :: nrow, ncol, nz
            integer(c_int), intent(in)      :: Ti(*), Tj(*)
            real(c_double), intent(in)      :: Tx(*)
            integer(c_int)                  :: Ap(*), Ai(*)
            real(c_double)                  :: Ax(*)
            integer(c_int)                  :: Map(*)

        end subroutine

      end interface   
    
      contains


    ! UMFPACK 
    !========
    ! The sparse solver 
      subroutine solve_umfpack(n, nnz, b, sol)

        ! Description
        !============
        ! This routine solves a sparse linear system by calling the 
        ! UMFPACK solvers. Here, it is assumed that the variables are
        ! in double precision. The input matrix A should be in a 
        ! MySparse type format and will be converted to compressed 
        ! column storage format (CCS) using the UMFPACK routines. The
        ! system is then solved by calling the required UMFPACK solver. 
        ! The rhs should be given in b, the solution will be returned in
        ! sol. Note that the solver here is a direct solver! 

        ! The result whether the solver is successful is parsed through
        ! 'info_solver'. It is zero when successful, and nonzero when
        ! not.

        ! Declare variables
        !==================
        ! Arguments
        integer, intent(in) 	:: n, nnz
        real(R8), intent(in)    :: b(n)
        real(R8), intent(out)   :: sol(n)

        ! Matrix in triplet format, as saved in b2mod_ma28
        integer			:: Arow(nnz), Acol(nnz)
        real(R8)                :: Aval(nnz)

        ! Auxiliary
        integer(c_int)     :: Ap(n+1), Ai(nnz), Map(nnz)
        integer(I8)        :: flag
        real(c_double)     :: Ax(nnz)

        real(c_double)     :: sol_c(n)
        real(c_double)     :: 
     &       control(umfpack_control), info(umfpack_info)

        type(c_ptr)        :: symbolic, numeric


        ! Read in matrix from b2mod_ma28 (storage can perhaps be moved in the future)
        Arow = irnh(1:nnz)-1
        Acol = icnh(1:nnz)-1
        Aval = a1h(1:nnz)

        ! Convert to CCS
        !===============
        ! First, convert the row and column indices to zero based 
        ! indexing

        ! Call converter
        call UmfpackTriplet2ColDI(n, n, nnz, Arow, 
     &       Acol, Aval, Ap, Ai, Ax, Map)

        ! Solve
        !======
        ! Solve using the C-routines of UMFPACK. The Double precision, 
        ! integer routines (_di_) are used 

        ! Set defaults
        call UmfpackDefaultsDI(control)

        ! Set symbolic factorization
        call UmfpackSymbolicDI(n, n, Ap, Ai, Ax, symbolic, 
     &       control, info)

        ! Compute numerical factorization
        call UmfpackNumericDI(Ap, Ai, Ax, symbolic, numeric, 
     &       control, info)

        ! Compute the solution
        call UmfpackSolveDI(umfpack_a, Ap, Ai, Ax, sol_c, b, numeric, 
     &       control, info)
        if (info(1) .ne. 0) then 
            print *, 'Linear solver may not have converged, info: ',
     &              info(1)
        end if

        ! Destroy symbolic
        call UmfpackFreeSymbolicDI(symbolic)

        ! Destroy numeric
        call UmfpackFreeNumericDI(numeric)

        ! Post-process
        !=============
        ! Cast solution
        sol = sol_c 

        ! Cast flag
        flag = nint(info(1))

      end subroutine

      
      end module b2mod_umfpack