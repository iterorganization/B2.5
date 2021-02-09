*-----------------------------------------------------------------------
*.specification

      subroutine b2trzh (nx, ny, ns, nnucl, zmaxmax, na, ti)          !som 05.02.2021  
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_constants     
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nnucl, zmaxmax
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
     & nucl2s(0:nnucl-1,0:zmaxmax), zmnucl_int(0:nnucl-1) 
      real (kind=R8) :: zamax_temp, am_temp, amti(0:nnucl-1),
     & amnucl(0:nnucl-1), zmnucl(0:nnucl-1), z2n(0:nnucl-1),
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
      if (ncall.eq.0) then

*   ..calculate the atom number and the atomic mass of each type of ions
        inucl=1
        zmnucl(0)=zamax(0)
        amnucl(0)=am(0)         
        zamax_temp=zamax(0)
        am_temp=am(0)
        do is = 0, ns-1
          if (zamax(is).ne.zamax_temp .and. am(is).ne.am_temp) then
           inucl=inucl+1 
           zmnucl(inucl-1)=zamax(is)
           amnucl(inucl-1)=am(is)  
           zamax_temp=zamax(is)
           am_temp=am(is)
          endif 
*   ..correspondence matrix
          nucl2s(inucl-1,nint(zn(is)))=is 
          zmnucl_int(inucl-1)=nint(zmnucl(inucl-1))
        enddo        
*   ..calculate the N-coefficients and all mass related variables 
        call b2trcn (nnucl, amnucl, N2, N5, N6, N8, N9, N10, 
     & N11, N12, mu, rmu, mum, mdm)
      endif
       
*   ..main part  
      do iy = -1, ny
        do ix = -1, nx
*   ..calculate the average Z**2*n
          do inucl = 0, nnucl-1
            z2n(inucl) = 0
            do iz= 1, zmnucl_int(inucl)
              z2n(inucl)=z2n(inucl)+iz**2*na(ix, iy, nucl2s(inucl,iz))
            enddo
          enddo

          amti=amnucl*mp/ti(nx,ny)
          zeffne=sum(z2n) 
*    ..normalized average Z**2*n
          nz2n=z2n/zeffne
*    ..calculate the inversed matrix for heat flux/r-moment
          call b2trca(nnucl, amti, N5, N6, N9,  
     & N10, N11, N12, nz2n, Ainv)

*    ..re-normalize part of the inversed matrix  
          do j = 0, 2*nnucl-1
            Ainv_mdf(0:nnucl,j)=amti(0:nnucl)*Ainv(nnucl:2*nnucl-1,j)
          enddo
*    ..G-objects for the thermal and friction forces	 
          do j = 0, nnucl-1
            G2(:,j)=nz2n(:)*nz2n(j)*N2(:,j)
            G8(:,j)=nz2n(:)*nz2n(j)*N8(:,j)
          enddo

*    ..calculate the thermal force transport coefficient 
          call b2trtf(nx, ny, nnucl, G2, G8, mum, mdm,
     & Ainv, Ainv_mdf, ix, iy)
*    ..calculate the friction force transport coefficient 
          call b2trfr(nx, ny, nnucl, G2, G8, mum, rmu,
     & z2n, Ainv, Ainv_mdf, ix, iy)
        enddo
      enddo

      ncall = ncall+1
      call subend ()
      return


      end subroutine b2trzh
