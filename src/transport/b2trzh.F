*-----------------------------------------------------------------------
*.specification

      subroutine b2trzh (nx, ny, ns, nnucl, znmax, na, ti)          !som 05.02.2021  
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_constants     
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nnucl, znmax
      real (kind=R8) ::
     &  na(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRZH computes the transport coefficient using 21N-moment approach, which  
*     is discussed in [V. M. Zhdanov. Transport Processes in Multicomponent Plasma,
*     vol. 44. 10 2002.]
*
*  3. description 
*
*     This routine computes:
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer ncall
      integer ix, iy, is, inucl, iz, j, 
     & nucl2s(0:nnucl-1,0:znmax), znnucl_int(0:nnucl-1) 
      real (kind=R8) :: zn_temp, am_temp, amti(0:nnucl-1),
     & amnucl(0:nnucl-1), znnucl(0:nnucl-1), z2n(0:nnucl-1),
     & zeffne, nz2n(0:nnucl-1),  
     & N2(0:nnucl-1,0:nnucl-1), N5(0:nnucl-1,0:nnucl-1), 
     & N6(0:nnucl-1,0:nnucl-1), N8(0:nnucl-1,0:nnucl-1),
     & N9(0:nnucl-1,0:nnucl-1),
     & N10(0:nnucl-1,0:nnucl-1), N11(0:nnucl-1,0:nnucl-1),
     & N12(0:nnucl-1,0:nnucl-1),
     & G2(0:nnucl-1,0:nnucl-1), G8(0:nnucl-1,0:nnucl-1),
     & mu(0:nnucl-1,0:nnucl-1), mum(0:nnucl-1,0:nnucl-1),
     & rmu(0:nnucl-1,0:nnucl-1),
     & mdm(0:nnucl-1,0:nnucl-1), Ainv(0:2*nnucl-1,0:2*nnucl-1),
     & Ainv_mdf(0:nnucl-1,0:2*nnucl-1)

      save ncall
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trzh')
*   ..set internal parameters on first call
!     if (ncall.eq.0) then

!       print *, 'checkpoint_31'
!       print *, znmax  
*   ..calculate the atom number and the atomic mass of each type of ions
        inucl=1
        znnucl(0)=zn(0)
        amnucl(0)=am(0)         
        zn_temp=zn(0)
        am_temp=am(0)
        do is = 0, ns-1
          if (zn(is).ne.zn_temp .and. am(is).ne.am_temp) then
           inucl=inucl+1 
           znnucl(inucl-1)=zn(is)
           amnucl(inucl-1)=am(is)  
           zn_temp=zn(is)
           am_temp=am(is)
          endif 
*   ..correspondence matrix
          nucl2s(inucl-1,nint(zamax(is)))=is 
          znnucl_int(inucl-1)=nint(znnucl(inucl-1))
        enddo        
!       print *, 'checkpoint_32'
*   ..calculate the N-coefficients and all mass related variables 
        call b2trcn (nnucl, amnucl, N2, N5, N6, N8, N9, N10, 
     & N11, N12, mu, rmu, mum, mdm)
!     endif
       
!       print *, 'checkpoint_33'
*   ..main part  
      do iy = -1, ny
        do ix = -1, nx
*   ..calculate the average Z**2*n
          do inucl = 0, nnucl-1
            z2n(inucl) = 0
            do iz= 1, znnucl_int(inucl)
!             print * 
!       print *,'inucl', inucl
!       print *,'znnucl_int(inucl)', znnucl_int(inucl)  
!       print *,'iz', iz
!       print *,'nucl2s(inucl,iz)', nucl2s(inucl,iz)
              z2n(inucl)=z2n(inucl)+iz**2*na(ix, iy, nucl2s(inucl,iz))
            enddo
          enddo

          amti=amnucl*mp/ti(nx,ny)
          zeffne=sum(z2n) 
*    ..normalized average Z**2*n
          nz2n=z2n/zeffne
*    ..calculate the inversed matrix for heat flux/r-moment
!       print *, 'checkpoint_34'
          call b2trca(nnucl, amti, N5, N6, N9,  
     & N10, N11, N12, nz2n, Ainv)

!       print *, 'checkpoint_35'
*    ..re-normalize part of the inversed matrix  
          do j = 0, 2*nnucl-1
            Ainv_mdf(0:nnucl-1,j)=
     & amti(0:nnucl-1)*Ainv(nnucl:2*nnucl-1,j)
          enddo
*    ..G-objects for the thermal and friction forces	 
          do j = 0, nnucl-1
            G2(:,j)=nz2n(:)*nz2n(j)*N2(:,j)
            G8(:,j)=nz2n(:)*nz2n(j)*N8(:,j)
          enddo

!       print *, 'checkpoint_36'
*    ..calculate the thermal force transport coefficient 
          call b2trtf(nx, ny, nnucl, G2, G8, mum, mdm,
     & Ainv, Ainv_mdf, ix, iy)
!       print *, 'checkpoint_37'
*    ..calculate the friction force transport coefficient 
          call b2trfr(nx, ny, nnucl, G2, G8, mum, rmu,
     & z2n, Ainv, Ainv_mdf, ix, iy)
        enddo
      enddo

!     print *, 'checkpoint_38'
      ncall = ncall+1
      call subend ()
      return


      end subroutine b2trzh
