*-----------------------------------------------------------------------
*.specification

      subroutine b2trca(nnucl, amti, N5, N6, N9,
     &  N10, N11, N12, nz2n, Ainv)                                      !som 05.02.2021
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl
      real (kind=R8) ::
     &  amti(0:nnucl-1), nz2n(0:nnucl-1),
     &  N5(0:nnucl-1,0:nnucl-1),
     &  N6(0:nnucl-1,0:nnucl-1),
     &  N9(0:nnucl-1,0:nnucl-1),
     &  N10(0:nnucl-1,0:nnucl-1), N11(0:nnucl-1,0:nnucl-1),
     &  N12(0:nnucl-1,0:nnucl-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: Ainv(0:2*nnucl-1,0:2*nnucl-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCA computes the inverted matrix for heat flux/R-moments.
*
*  3. description
*
*     This routine computes:
*     - Ainv: the inverted matrix
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i,j,dnnucl
      real (kind=R8) :: temp,
     &  G5(0:nnucl-1,0:nnucl-1),
     &  G6(0:nnucl-1,0:nnucl-1),
     &  G9(0:nnucl-1,0:nnucl-1),
     &  G10(0:nnucl-1,0:nnucl-1), G11(0:nnucl-1,0:nnucl-1),
     &  G12(0:nnucl-1,0:nnucl-1), A(0:2*nnucl-1,0:2*nnucl-1)
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trca')
*   ..main
      dnnucl=2*nnucl
*   ..G-objects for the martix
      do j = 0, nnucl-1
        G5(:,j)=nz2n(:)*nz2n(j)*N5(:,j)
        G6(:,j)=nz2n(:)*nz2n(j)*N6(:,j)
        G9(:,j)=nz2n(:)*nz2n(j)*N9(:,j)
        G10(:,j)=nz2n(:)*nz2n(j)*N10(:,j)
        G11(:,j)=nz2n(:)*nz2n(j)*N11(:,j)
        G12(:,j)=nz2n(:)*nz2n(j)*N12(:,j)
      enddo

*   ..make the matrix before the inversion
*   ..fill both diagonal and cross elements
      A(0:nnucl-1,0:nnucl-1)=G6(0:nnucl-1,0:nnucl-1)
      A(nnucl:dnnucl-1,0:nnucl-1)=7.0_R8*G10(0:nnucl-1,0:nnucl-1)
      do j = 0, nnucl-1
        A(0:nnucl-1,j+nnucl)=amti(0:nnucl-1)*G10(0:nnucl-1,j)
        A(nnucl:dnnucl-1,j+nnucl)=amti(j)*G12(0:nnucl-1,j)

*    ..addition for diagonal elements
        A(j,j)=A(j,j)+sum(G5(j,:))
        temp=sum(G9(j,:))
        A(j+nnucl,j)=A(j+nnucl,j)+7.0_R8*temp
        A(j,j+nnucl)=A(j,j+nnucl)+amti(j)*temp
        A(j+nnucl,j+nnucl)=A(j+nnucl,j+nnucl)+amti(j)*sum(G11(j,:))
      enddo

      call invert_matrix(A, Ainv, nnucl, nz2n)

      call subend ()
      return

      end subroutine b2trca
