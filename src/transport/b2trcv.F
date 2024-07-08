*-----------------------------------------------------------------------
*.specification

      subroutine b2trcv(nnucl, amti,mdm, N3, N4, N13,
     & N14, N15, N16, nz2n, Cinv)                                       !som 05.02.2021
      use b2mod_types
      use b2mod_subsys
      implicit none
      external DGETRF
      external DGETRI
*   ..input arguments (unchanged on exit)
      integer nnucl
      real (kind=R8) :: amti(0:nnucl-1), mdm(0:nnucl-1,0:nnucl-1),
     &  nz2n(0:nnucl-1),
     &  N3(0:nnucl-1,0:nnucl-1),
     &  N4(0:nnucl-1,0:nnucl-1),
     &  N13(0:nnucl-1,0:nnucl-1), N14(0:nnucl-1,0:nnucl-1),
     &  N15(0:nnucl-1,0:nnucl-1), N16(0:nnucl-1,0:nnucl-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: Cinv(0:2*nnucl-1,0:2*nnucl-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCV computes the inverted matrix for the viscous stress tensor.
*
*  3. description
*
*     This routine computes:
*     - Cinv: the inverted matrix
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i,j,dnnucl
      real (kind=R8) :: temp,
     &  G3(0:nnucl-1,0:nnucl-1),
     &  G4(0:nnucl-1,0:nnucl-1),
     &  G13(0:nnucl-1,0:nnucl-1), G14(0:nnucl-1,0:nnucl-1),
     &  G15(0:nnucl-1,0:nnucl-1), G16(0:nnucl-1,0:nnucl-1),
     &  C(0:2*nnucl-1,0:2*nnucl-1)
      real(kind=R8), dimension(size(C,1)) :: work  ! work array for LAPACK
      integer, dimension(size(C,1)) :: ipiv   ! pivot indices
      integer :: info
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trcv')
*   ..main
      dnnucl=2*nnucl
*    ..G-objects for the matrix
      do j = 0, nnucl-1
        G3(:,j)=nz2n(:)*nz2n(j)*N3(:,j)
        G4(:,j)=nz2n(:)*nz2n(j)*N4(:,j)
        G13(:,j)=nz2n(:)*nz2n(j)*N13(:,j)
        G14(:,j)=nz2n(:)*nz2n(j)*N14(:,j)
        G15(:,j)=nz2n(:)*nz2n(j)*N15(:,j)
        G16(:,j)=nz2n(:)*nz2n(j)*N16(:,j)
      enddo

*    ..make the matrix before the inversion
*    ..fill both diagonal and cross elements
      C(0:nnucl-1,0:nnucl-1)=G4(0:nnucl-1,0:nnucl-1)
      C(nnucl:dnnucl-1,0:nnucl-1)=3.5_R8*mdm(0:nnucl-1,0:nnucl-1)*
     & G14(0:nnucl-1,0:nnucl-1)
      do j = 0, nnucl-1
        C(0:nnucl-1,j+nnucl)=amti(j)*G14(0:nnucl-1,j)
        C(nnucl:dnnucl-1,j+nnucl)=amti(0:nnucl-1)*G16(0:nnucl-1,j)

*    ..addition for diagonal elements
        C(j,j)=C(j,j)+sum(G3(j,:))
        temp=sum(G13(j,:))
        C(j+nnucl,j)=C(j+nnucl,j)+3.5_R8*temp
        C(j,j+nnucl)=C(j,j+nnucl)+amti(j)*temp
        C(j+nnucl,j+nnucl)=C(j+nnucl,j+nnucl)+amti(j)*sum(G15(j,:))
      enddo

*    ..the matrix inversion
      Cinv=C
      ! DGETRF computes an LU factorization of a general M-by-N matrix C
      ! using partial pivoting with row interchanges.
      call DGETRF(dnnucl, dnnucl, Cinv, dnnucl, ipiv, info)

      if (info /= 0) then
         print *, nz2n
         print *
         do i = 0, dnnucl-1
           print *, (C(i,j), j = 0, dnnucl-1 )
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
      call DGETRI(dnnucl, Cinv, dnnucl, ipiv, work, dnnucl, info)

      call xertst (info.eq.0, 'Matrix inversion failed!')

!       print *, 'C:'
!       do i = 0, dnnucl-1
!         print *, (C(i,j), j = 0, dnnucl-1 )
!       enddo
!       print *
!       print *, 'Cinv:'
!       do i = 0, dnnucl-1
!         print *, (Cinv(i,j), j = 0, dnnucl-1 )
!       enddo
!       print *
      call subend ()
      return

      end subroutine b2trcv
