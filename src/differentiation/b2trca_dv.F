*-----------------------------------------------------------------------
*.specification

      subroutine b2trca_dv(nnucl, amti, amtid, N5, N6, N9,
     &  N10, N11, N12, nz2n, nz2nd, Ainv, Ainvd, nbdirs)
      use b2mod_types
      use b2mod_subsys
      use b2mod_diffsizes
      implicit none
      external DGETRF
      external DGETRI
*   ..input arguments (unchanged on exit)
      integer nnucl, nbdirs
      real (kind=R8) ::
     &  amti(0:nnucl-1), nz2n(0:nnucl-1),
     &  N5(0:nnucl-1,0:nnucl-1),
     &  N6(0:nnucl-1,0:nnucl-1),
     &  N9(0:nnucl-1,0:nnucl-1),
     &  N10(0:nnucl-1,0:nnucl-1), N11(0:nnucl-1,0:nnucl-1),
     &  N12(0:nnucl-1,0:nnucl-1)
      real (kind=R8) ::
     &  amtid(nbdirs, 0:nnucl-1), nz2nd(nbdirs, 0:nnucl-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: Ainv(0:2*nnucl-1,0:2*nnucl-1)
      real (kind=R8) :: Ainvd(nbdirs,0:2*nnucl-1,0:2*nnucl-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCA_DV is the discretely differentiated version of b2trca
*     in tangent vector AD mode.
*
*  3. description
*
*     This routine computes:
*     - Ainv: the inverted matrix
*     - Ainvd: the directional derivative of Ainv, as Ainvd = -Ainv*Ad*Ainv
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
      real (kind=R8) :: tempd(nbdirsmax),
     &  G5d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G6d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G9d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G10d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G11d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G12d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  Ad(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
      real(kind=R8), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: info, nd
      real (kind=R8) :: dum(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trca_dv')
*   ..main
      dnnucl=2*nnucl
      G5d = 0.0_R8
      G6d = 0.0_R8
      G9d = 0.0_R8
      G10d = 0.0_R8
      G11d = 0.0_R8
      G12d = 0.0_R8
*   ..G-objects for the martix
      do j = 0, nnucl-1
        do nd = 1, nbdirs
          G5d(nd,:,j)=N5(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*nz2nd(nd,j))
          G6d(nd,:,j)=N6(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*nz2nd(nd,j))
          G9d(nd,:,j)=N9(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*nz2nd(nd,j))
          G10d(nd,:,j)=N10(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*
     &                 nz2nd(nd,j))
          G11d(nd,:,j)=N11(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*
     &                 nz2nd(nd,j))
          G12d(nd,:,j)=N12(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*
     &                 nz2nd(nd,j))
        end do
        G5(:,j)=nz2n(:)*nz2n(j)*N5(:,j)
        G6(:,j)=nz2n(:)*nz2n(j)*N6(:,j)
        G9(:,j)=nz2n(:)*nz2n(j)*N9(:,j)
        G10(:,j)=nz2n(:)*nz2n(j)*N10(:,j)
        G11(:,j)=nz2n(:)*nz2n(j)*N11(:,j)
        G12(:,j)=nz2n(:)*nz2n(j)*N12(:,j)
      enddo

*   ..make the matrix before the inversion
      Ad = 0.0_R8
*   ..fill both diagonal and cross elements
      do nd = 1, nbdirs
        Ad(nd,0:nnucl-1,0:nnucl-1)=G6d(nd,0:nnucl-1,0:nnucl-1)
        Ad(nd,nnucl:dnnucl-1,0:nnucl-1)=7.0_R8*
     &                                  G10d(nd,0:nnucl-1,0:nnucl-1)
      end do
      A(0:nnucl-1,0:nnucl-1)=G6(0:nnucl-1,0:nnucl-1)
      A(nnucl:dnnucl-1,0:nnucl-1)=7.0_R8*G10(0:nnucl-1,0:nnucl-1)
      do j = 0, nnucl-1
        temp=sum(G9(j,:))
        do nd = 1, nbdirs
          Ad(nd,0:nnucl-1,j+nnucl)=G10(0:nnucl-1,j)*amtid(nd,0:nnucl-1)+
     &                             amti(0:nnucl-1)*G10d(nd,0:nnucl-1,j)
          Ad(nd,nnucl:dnnucl-1,j+nnucl)=G12(0:nnucl-1,j)*amtid(nd,j)+
     &                                  amti(j)*G12d(nd,0:nnucl-1,j)

*        ..addition for diagonal elements
          Ad(nd,j,j)=Ad(nd,j,j)+SUM(G5d(nd,j,:))
          tempd(nd)=SUM(G9d(nd,j,:))
          Ad(nd,j+nnucl,j)=Ad(nd,j+nnucl,j)+7.0_R8*tempd(nd)
          Ad(nd,j,j+nnucl)=Ad(nd,j,j+nnucl)+temp*amtid(nd,j)+
     &                     amti(j)*tempd(nd)
          Ad(nd,j+nnucl,j+nnucl)=Ad(nd,j+nnucl,j+nnucl)+SUM(G11(j,:))*
     &                           amtid(nd,j)+amti(j)*SUM(G11d(nd,j,:))
        end do

        A(0:nnucl-1,j+nnucl)=amti(0:nnucl-1)*G10(0:nnucl-1,j)
        A(nnucl:dnnucl-1,j+nnucl)=amti(j)*G12(0:nnucl-1,j)

*    ..addition for diagonal elements
        A(j,j)=A(j,j)+sum(G5(j,:))
        A(j+nnucl,j)=A(j+nnucl,j)+7.0_R8*temp
        A(j,j+nnucl)=A(j,j+nnucl)+amti(j)*temp
        A(j+nnucl,j+nnucl)=A(j+nnucl,j+nnucl)+amti(j)*sum(G11(j,:))
      enddo

*    ..the matrix inversion
      Ainv=A
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(dnnucl, dnnucl, Ainv, dnnucl, ipiv, info)

      if (info /= 0) then
         print *, nz2n
         print *
         do i = 0, dnnucl-1
           print *, (A(i,j), j = 0, dnnucl-1 )
         enddo
         print *
#ifndef NAGFOR
         call sleep(3)
#endif
         flush(6)
         call xerrab ('Matrix is numerically singular!')
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(dnnucl, Ainv, dnnucl, ipiv, work, dnnucl, info)

      call xertst(info.eq.0, 'Matrix inversion failed!')

*    ..discrete version of d(inv(A)):
*      A*Ainv = I
*      d(A*Ainv) = Ad*Ainv + A*Ainvd = 0
*      Ainvd = -Ainv*Ad*Ainv
      Ainvd = 0.0_R8
      dum = 0.0_R8
      do nd = 1, nbdirs
*     ..first do dum = -Ainv*Ad
        do i = 0, nnucl-1
          do j = 0, nnucl-1
            dum(nd,i,j) = -sum(Ainv(i,:)*Ad(nd,:,j))
          end do
        end do
*     ..then do Ainvd = dum*Ainv
        do i = 0, nnucl-1
          do j = 0, nnucl-1
            Ainvd(nd,i,j) = sum(dum(nd,i,:)*Ainv(:,j))
          end do
        end do
      end do

      call subend ()
      return

      end subroutine b2trca_dv
