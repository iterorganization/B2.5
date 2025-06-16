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
      use b2mod_b2plot
      use b2mod_b2cmpa
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      use b2us_data
      use b2mod_dimensions
      implicit none
      integer lun
c
c*** Local variables
c
      integer i,j,k,is,iCv
      real (kind=R8) ::
     . hlp,hlp_tot,hlp_he,hlp_s,hlp_si,hlp_so,hlp_a,hlp_m
      real (kind=R8) ::
     . hlp_p(0:DEF_NATM+DEF_NMOL),
     . hlp_pi(0:DEF_NATM+DEF_NMOL), hlp_po(0:DEF_NATM+DEF_NMOL)
      integer n_wll_hlp
      logical lpfr_i, lpfr_o
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
        if (j_he_at.gt.0) hlp_he=hlp_he+wldna(k,j_he_at,0)
        do is=1,3
          if (j_h_at(is).eq.0) cycle
          hlp_tot=hlp_tot+wldna(k,j_h_at(is),0)
          do j=1,nmoli
            if(mlcmp(j_h_at(is),j).gt.0) then
              hlp_tot=hlp_tot+wldnm(k,j,0)*mlcmp(j_h_at(is),j)
            end if
          end do
        end do
        if (j_he_at.gt.0) then
          do j=1,nmoli
            if(mlcmp(j_he_at,j).gt.0) then
              hlp_he=hlp_he+wldnm(k,j,0)*mlcmp(j_he_at,j)
            end if
          end do
        end if
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
      hlp_a=0.0_R8
      hlp_m=0.0_R8
      hlp_p=0.0_R8
      hlp_s=0.0_R8
c
c*** PFR regions
      hlp_pi=0.0_R8
      hlp_si=0.0_R8
      hlp_po=0.0_R8
      hlp_so=0.0_R8
      if (npfr_cvs.gt.0) then
        do j=1,npfr_cvs
          hlp=0.0_R8
          iCv=pfr_cvs(j)
          lpfr_i=mod(mpg%cvReg(iCv),4).eq.3
          lpfr_o=mod(mpg%cvReg(iCv),4).eq.0.and.mpg%cvReg(iCv).ne.0
          hlp_s=hlp_s+geo%cvVol(iCv)
          if (lpfr_i) hlp_si=hlp_si+geo%cvVol(iCv)
          if (lpfr_o) hlp_so=hlp_so+geo%cvVol(iCv)
          if (switch%use_eirene.ne.0) then
            do k=1,natmi
              hlp_a=dab2(iCv,k,1)*tab2(iCv,k,1)
              hlp_p(k)=hlp_p(k)+hlp_a*geo%cvVol(iCv)
              if (lpfr_i) hlp_pi(k)=hlp_pi(k)+hlp_a*geo%cvVol(iCv)
              if (lpfr_o) hlp_po(k)=hlp_po(k)+hlp_a*geo%cvVol(iCv)
              hlp=hlp+hlp_a
            end do
            do k=1,nnmoli
              hlp_m=dmb2(iCv,k,1)*tmb2(iCv,k,1)
              hlp_p(natmi+k)=hlp_p(natmi+k)+hlp_m*geo%cvVol(iCv)
              if (lpfr_i) hlp_pi(natmi+k)=hlp_pi(natmi+k)+
     .         hlp_m*geo%cvVol(iCv)
              if (lpfr_o) hlp_po(natmi+k)=hlp_po(natmi+k)+
     .         hlp_m*geo%cvVol(iCv)
              hlp=hlp+hlp_m
            end do
          else
            do is=0,state%ns-1
              if (is_neutral(is)) then
                if (nint(zn(is)).eq.1) then
                  hlp_a=hlp_a+state%pl%na(iCv,is)*state%pl%tn(iCv)
                else
                  hlp_a=hlp_a+state%pl%na(iCv,is)*state%pl%ti(iCv)
                end if
                hlp_p(b2espcr(is))=
     .           hlp_p(b2espcr(is))+hlp_a*geo%cvVol(iCv)
                if (lpfr_i) hlp_pi(b2espcr(is))=
     .           hlp_pi(b2espcr(is))+hlp_a*geo%cvVol(iCv)
                if (lpfr_o) hlp_po(b2espcr(is))=
     .           hlp_po(b2espcr(is))+hlp_a*geo%cvVol(iCv)
                hlp=hlp+hlp_a
              end if
            end do
          end if
          hlp_p(0)=hlp_p(0)+hlp*geo%cvVol(iCv)
          if (lpfr_i) hlp_pi(0)=hlp_pi(0)+hlp*geo%cvVol(iCv)
          if (lpfr_o) hlp_po(0)=hlp_po(0)+hlp*geo%cvVol(iCv)
        end do
      end if
c
c*** Averaging and output
c
      if(hlp_s.gt.0.0_R8) hlp_p(:)=hlp_p(:)/hlp_s
      if(hlp_si.gt.0.0_R8) hlp_pi(:)=hlp_pi(:)/hlp_si
      if(hlp_so.gt.0.0_R8) hlp_po(:)=hlp_po(:)/hlp_so
      open(lun,file='pressure.ntr')
      write(lun,'(1p,3e14.3,a)') hlp_p(0),hlp_pi(0),hlp_po(0),' Total'
      if (switch%use_eirene.ne.0) then
        do k=1,natmi
          write(lun,'(1p,3e14.3,a)') hlp_p(k),hlp_pi(k),hlp_po(k),
     .                               textan(k-1)
        end do
        do k=1,nmoli
          j=k+natmi
          write(lun,'(1p,3e14.3,a)') hlp_p(j),hlp_pi(j),hlp_po(j),
     .                               textmn(k-1)
        end do
      else
        do k=1,nspecies
          write(lun,'(1p,3e14.3,a)') hlp_p(k),hlp_pi(k),hlp_po(k),
     .                               species_list(k)
        end do
      end if
      close(lun)
c----------------------------------------------------------------------
      return
      end subroutine b2phel
