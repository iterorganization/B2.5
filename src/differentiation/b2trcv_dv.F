*-----------------------------------------------------------------------
*.specification

      subroutine b2trcv_dv(nnucl, amti, amtid, mdm, N3, N4, N13,
     & N14, N15, N16, nz2n, nz2nd, Cinv, Cinvd, nbdirs)
      use b2mod_types
      use b2mod_subsys
      use b2mod_diffsizes
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
      real (kind=R8) :: amtid(nbdirsmax,0:nnucl-1),
     &  nz2nd(nbdirsmax,0:nnucl-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: Cinv(0:2*nnucl-1,0:2*nnucl-1),
     &  Cinvd(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCV_DV is the discretely differentiated version of b2trcv
*     in tangent vector AD mode.
*
*  3. description
*
*     This routine computes:
*     - Cinv: the inverted matrix
*     - Cinvd: the directional derivative of Cinv, as Cinvd = -Cinv*Cd*Cinv
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i,j,dnnucl,nbdirs
      real (kind=R8) :: temp,
     &  G3(0:nnucl-1,0:nnucl-1),
     &  G4(0:nnucl-1,0:nnucl-1),
     &  G13(0:nnucl-1,0:nnucl-1), G14(0:nnucl-1,0:nnucl-1),
     &  G15(0:nnucl-1,0:nnucl-1), G16(0:nnucl-1,0:nnucl-1),
     &  C(0:2*nnucl-1,0:2*nnucl-1)
      real (kind=R8) :: tempd(nbdirsmax),
     &  G3d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G4d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G13d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G14d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G15d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  G16d(nbdirsmax,0:nnucl-1,0:nnucl-1),
     &  Cd(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
      real(kind=R8), dimension(size(C,1)) :: work  ! work array for LAPACK
      integer, dimension(size(C,1)) :: ipiv   ! pivot indices
      integer :: info, nd
      real (kind=R8) :: dum(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trcv_dv')
*   ..main
      dnnucl=2*nnucl
      G3d = 0.0_R8
      G4d = 0.0_R8
      G13d = 0.0_R8
      G14d = 0.0_R8
      G15d = 0.0_R8
      G16d = 0.0_R8
*    ..G-objects for the matrix
      do j = 0, nnucl-1
        do nd = 1, nbdirs
          G3d(nd,:,j)=N3(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*nz2nd(nd,j))
          G4d(nd,:,j)=N4(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*nz2nd(nd,j))
          G13d(nd,:,j)=N13(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*
     &                 nz2nd(nd,j))
          G14d(nd,:,j)=N14(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*
     &                 nz2nd(nd,j))
          G15d(nd,:,j)=N15(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*
     &                 nz2nd(nd,j))
          G16d(nd,:,j)=N16(:,j)*(nz2n(j)*nz2nd(nd,:)+nz2n(:)*
     &                 nz2nd(nd,j))
        end do
        G3(:,j)=nz2n(:)*nz2n(j)*N3(:,j)
        G4(:,j)=nz2n(:)*nz2n(j)*N4(:,j)
        G13(:,j)=nz2n(:)*nz2n(j)*N13(:,j)
        G14(:,j)=nz2n(:)*nz2n(j)*N14(:,j)
        G15(:,j)=nz2n(:)*nz2n(j)*N15(:,j)
        G16(:,j)=nz2n(:)*nz2n(j)*N16(:,j)
      enddo

*    ..make the matrix before the inversion
      Cd = 0.0_R8
*    ..fill both diagonal and cross elements
      do nd = 1, nbdirs
        Cd(nd,0:nnucl-1,0:nnucl-1)=G4d(nd, 0:nnucl-1, 0:nnucl-1)
        Cd(nd,nnucl:dnnucl-1,0:nnucl-1)=mdm(0:nnucl-1,0:nnucl-1)*3.5_R8*
     &                                  G14d(nd,0:nnucl-1,0:nnucl-1)
      end do
      C(0:nnucl-1,0:nnucl-1)=G4(0:nnucl-1,0:nnucl-1)
      C(nnucl:dnnucl-1,0:nnucl-1)=3.5_R8*mdm(0:nnucl-1,0:nnucl-1)*
     & G14(0:nnucl-1,0:nnucl-1)
      do j = 0, nnucl-1
        temp = SUM(G13(j, :))
        do nd = 1, nbdirs
          Cd(nd,0:nnucl-1,j+nnucl)=G14(0:nnucl-1,j)*amtid(nd,j)+amti(j)*
     &                             G14d(nd,0:nnucl-1,j)
          Cd(nd,nnucl:dnnucl-1,j+nnucl)=
     &      G16(0:nnucl-1,j)*amtid(nd,0:nnucl-1)+
     &      amti(0:nnucl-1)*G16d(nd,0:nnucl-1,j)

*    ..addition for diagonal elements
          Cd(nd,j,j)=Cd(nd,j,j)+SUM(G3d(nd,j,:))
          tempd(nd)=SUM(G13d(nd,j,:))
          Cd(nd,j+nnucl,j)=Cd(nd,j+nnucl,j)+3.5_R8*tempd(nd)
          Cd(nd,j,j+nnucl)=Cd(nd,j,j+nnucl)+temp*amtid(nd,j)+
     &                     amti(j)*tempd(nd)
          Cd(nd,j+nnucl,j+nnucl)=Cd(nd,j+nnucl,j+nnucl)+
     &      SUM(G15(j,:))*amtid(nd,j)+amti(j)*SUM(G15d(nd,j,:))
        END DO
        C(0:nnucl-1,j+nnucl)=amti(j)*G14(0:nnucl-1,j)
        C(nnucl:dnnucl-1,j+nnucl)=amti(0:nnucl-1)*G16(0:nnucl-1,j)

*    ..addition for diagonal elements
        C(j,j)=C(j,j)+sum(G3(j,:))
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

*    ..discrete version of d(inv(C)):
*      C*Cinv = I
*      d(C*Cinv) = Cd*Cinv + C*Cinvd = 0
*      Cinvd = -Cinv*Cd*Cinv
      Cinvd = 0.0_R8
      dum = 0.0_R8
      do nd = 1, nbdirs
*     ..first do dum = -Cinv*Cd
        do i = 0, nnucl-1
          do j = 0, nnucl-1
            dum(nd,i,j) = -sum(Cinv(i,:)*Cd(nd,:,j))
          end do
        end do
*     ..then do Cinvd = dum*Cinv
        do i = 0, nnucl-1
          do j = 0, nnucl-1
            Cinvd(nd,i,j) = sum(dum(nd,i,:)*Cinv(:,j))
          end do
        end do
      end do

      call subend ()
      return

      end subroutine b2trcv_dv
