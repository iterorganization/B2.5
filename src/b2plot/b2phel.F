      subroutine b2phel(lun)
c
c----------------------------------------------------------------------
c***  Estimates the helium fraction in the pumped gas
c***  and the neutral pressure in PFR
c***  Both estimates are crude and not universal: hard-wired data!
c***
c*** Output files:
c***
c***  fraction.he   helium fraction in the pumped gas
c***  pressure.ntr  neutral pressure in PFR
c----------------------------------------------------------------------
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
      implicit none
      integer lun
c
c*** Local variables
c
      integer i,j,k,l,m,n,is
      real (kind=R8) ::
     . hlp_tot,hlp_he,hlp,hlp_p,hlp_s,hlp_pi,hlp_si,hlp_po,hlp_so
      integer n_wll_hlp
c
c*** The following lines are case-dependent
c<<<
      parameter (n_wll_hlp=14)
      integer ind_wll_hlp(n_wll_hlp)
c
c*** Wall segments representing the duct region:
      data ind_wll_hlp /4,5,6,7,8,9,10,11,12,13,14,64,65,66/
c>>>
c----------------------------------------------------------------------
c*** Estimate the helium fraction.
c
      hlp_tot=0.0_R8
      hlp_he=0.0_R8
      do i=1,n_wll_hlp
        k=ind_wll_hlp(i)
        hlp_he=hlp_he+wldna(k,j_he_at,0)
        do is=1,3
          if (j_h_at(is).eq.0) cycle
          hlp_tot=hlp_tot+wldna(k,j_h_at(is),0)
          do j=1,nmoli
            if(mlcmp(j_h_at(is),j).gt.0) then
              hlp_tot=hlp_tot+wldnm(k,j,0)*mlcmp(j_h_at(is),j)
            end if
          end do
        end do
        do j=1,nmoli
          if(mlcmp(j_he_at,j).gt.0) then
            hlp_he=hlp_he+wldnm(k,j,0)*mlcmp(j_he_at,j)
          end if
        end do
      end do
      hlp_tot=hlp_tot+hlp_he
      hlp=0.0_R8
      if(hlp_tot.gt.0.0_R8) hlp=hlp_he/hlp_tot
      open(lun,file='fraction.he')
      rewind(lun)
      write(lun,'(1p,e13.3)') hlp
      close(lun)
c----------------------------------------------------------------------
c*** Now estimate the neutral pressure in PFR - as the average over the
c*** the plasma cells adjacent to the PFR, weighted with the surface
c*** area
c
      hlp_p=0.0_R8
      hlp_s=0.0_R8
c
c*** Inner divertor
      hlp_pi=0.0_R8
      hlp_si=0.0_R8
      if (lpfrb_i.ge.0) then
        m=lpfrb_i
      else
        m=leftcut(1)+lpfrb_i
      end if
      if (lpfrt_i.ge.0) then
        n=lpfrt_i
      else
        n=leftcut(1)+lpfrt_i
      end if
#ifdef WG_TODO
      do j=m,n
        hlp_s=hlp_s+gs(j,0,1)
        hlp_si=hlp_si+gs(j,0,1)
        hlp=0.0_R8
        do k=1,natmi
          hlp=hlp+dab2(j+1,1,k,1)*tab2(j+1,1,k,1)
        end do
        do k=1,nmoli
          hlp=hlp+dmb2(j+1,1,k,1)*tmb2(j+1,1,k,1)
        end do
        hlp_p=hlp_p+hlp*gs(j,0,1)
        hlp_pi=hlp_pi+hlp*gs(j,0,1)
      end do
#endif
c
c*** Outer divertor
      hlp_po=0.0_R8
      hlp_so=0.0_R8
      if (lpfrt_o.ge.0) then
        m=lpfrt_o
      else
        m=rightcut(1)-lpfrt_o
      endif
      if (lpfrb_o.ge.0) then
        n=lpfrb_o
      else
        n=nx+lpfrb_o
      endif
#ifdef WG_TODO
      do j=m,n
        hlp_s=hlp_s+gs(j,0,1)
        hlp_so=hlp_so+gs(j,0,1)
        hlp=0.0_R8
        do k=1,natmi
          hlp=hlp+dab2(j+1,1,k,1)*tab2(j+1,1,k,1)
        end do
        do k=1,nmoli
          hlp=hlp+dmb2(j+1,1,k,1)*tmb2(j+1,1,k,1)
        end do
        hlp_p=hlp_p+hlp*gs(j,0,1)
        hlp_po=hlp_po+hlp*gs(j,0,1)
      end do
#endif
c
c*** Averaging and output
c
      if(hlp_s.gt.0.0_R8) hlp_p=hlp_p/hlp_s
      if(hlp_si.gt.0.0_R8) hlp_pi=hlp_pi/hlp_si
      if(hlp_so.gt.0.0_R8) hlp_po=hlp_po/hlp_so
      open(lun,file='pressure.ntr')
      write(lun,'(1p,3e14.3)') hlp_p,hlp_pi,hlp_po
      close(lun)
c----------------------------------------------------------------------
      return
      end
