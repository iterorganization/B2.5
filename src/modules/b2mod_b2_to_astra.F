csw 09jan2014 modifications as suggested by J.Terry/MIT
csw 09jan2014 extra module to define some arrays
      module b2mod_b2_to_astra_com
      use b2mod_types
      integer :: np
      real(kind=R8),allocatable :: p(:),x(:),cb(:),cc(:),cd(:)
      end module b2mod_b2_to_astra_com

      module b2mod_b2_to_astra                          !srv 28.07.11 {
      use b2mod_types
      use b2mod_astra_to_b2
#ifdef B25_EIRENE
      use b2mod_neutrals_namelist
     , ,only: b2espcr, latmscl, mlcmp
#ifndef SOLPS4_3
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use eirdiag ! IGNORE
     , ,only: nnatmi, nnmoli, dab2, dmb2, tab2, tmb2, rfluxa, rfluxm
#endif
#else
      use b2mod_eirdiag
#endif
      implicit none
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif

csw 09jan2014 declare fun to be private
       private fun

       real (kind=R8), save :: neb                        ! boundary electron density      [1/m**3]
       real (kind=R8), save :: teb                        ! boundary electron temperature  [keV]
       real (kind=R8), save :: tib                        ! boundary ion temperature       [keV]
       real (kind=R8), allocatable, save :: se(:)         ! integrated electron source     [1/m**3/s]
       real (kind=R8), allocatable, save :: seh(:)        ! integrated electron heat sink  []
       real (kind=R8), allocatable, save :: sih(:)        ! integrated ion heat sink       []
       real (kind=R8), allocatable, save :: svol(:)       ! integrated volume              [1/m**3]
       real (kind=R8), allocatable, save :: R_midpl(:)    ! radial surface coordinates     [m]
       real (kind=R8), allocatable, save :: Z_midpl(:)    ! Z              coordinates     [m]
       real (kind=R8), allocatable, save :: Rc_midpl(:)   ! radial cell centre coordinates [m]
       real (kind=R8), allocatable, save :: chci_koeff(:,:) !
       real (kind=R8), allocatable, save :: sfna(:,:)     ! integrated partical flux        [1/s]
       real (kind=R8), allocatable, save :: sfna_mdf(:,:) ! integrated modified partical flux    [1/s]
       real (kind=R8), allocatable, save :: sfne(:)       ! integrated electron flux        [1/s]
       real (kind=R8), allocatable, save :: sfhe_mdf(:)   ! integrated electron heat flux   [MW]
       real (kind=R8), allocatable, save :: sfhi_mdf(:)   ! integrated ion heat flux        [MW]
       integer, save :: iy_Ratb, iy_crb_astra
       contains

       subroutine alloc_b2mod_b2_to_astra(nxd, nyd, nsd)
        implicit none
        integer nxd, nyd, nsd
        allocate (se(-1:nyd))
        allocate (seh(-1:nyd))                            !srv 11.06.13
        allocate (sih(-1:nyd))                            !srv 11.06.13
        allocate (svol(-1:nyd))
        allocate (R_midpl(-1:nyd))
        allocate (Z_midpl(-1:nyd))
        allocate (Rc_midpl(-1:nyd))
        allocate (chci_koeff(-1:nxd,-1:nyd))
        allocate (sfna(-1:nyd,0:nsd))
        allocate (sfna_mdf(-1:nyd,0:nsd))
        allocate (sfne(-1:nyd))
        allocate (sfhe_mdf(-1:nyd))
        allocate (sfhi_mdf(-1:nyd))
        se         = 0.0_R8
        seh        = 0.0_R8
        sih        = 0.0_R8
        R_midpl    = 0.0_R8
        Rc_midpl   = 0.0_R8
        chci_koeff = 0.0_R8
        write (*,*) 'Allocated space for b2mod_b2_to_astra'
       end subroutine alloc_b2mod_b2_to_astra

       subroutine dealloc_b2mod_b2_to_astra
        deallocate(se)
        deallocate(seh)
        deallocate(sih)
        deallocate (svol)
        deallocate (R_midpl)
        deallocate (Z_midpl)
        deallocate (Rc_midpl)
        deallocate (chci_koeff)
        deallocate (sfna)
        deallocate (sfna_mdf)
        deallocate (sfne)
        deallocate (sfhe_mdf)
        deallocate (sfhi_mdf)
       end subroutine dealloc_b2mod_b2_to_astra

       subroutine write_b2_to_astra(nxd, nyd, nsd, ismain, ismain0,
     &   nwrite)
        use b2mod_anomalous_transport
        use b2mod_indirect
        use b2mod_geo
        use b2mod_plasma
        use b2mod_constants
        use b2mod_b2cmpa
csw 09jan2014 use extra little module
        use b2mod_b2_to_astra_com

        implicit none
        integer :: nxd, nyd, nsd, ismain, ismain0, nwrite
#ifdef B25_EIRENE
        integer isa, ism
        character chns*3
        real (kind=R8) :: wrk0(-1:nxd,-1:nyd)
#endif

csw 09jan2014 remove np from declaration here
c        integer :: mintopcut, ix, iy, is, iout, np
        integer :: mintopcut, ix, iy, is, iout
        integer :: maxtopcut
        integer, save :: ncall = 0, use_eirene = 0

csw 09jan2014 remove declaration p,x,cb,cc,cd from here
c        real (kind=R8) :: p(0:nxd-1), x(0:nxd-1), cb(0:nxd-1),
c     &  cc(0:nxd-1), cd(0:nxd-1), DSEVAL, DFMIN, arg, y(0:nxd-1)
        real (kind=R8) :: DSEVAL, DFMIN, arg, y(0:nxd-1)

csw 09jan2014 add extra allocates
        allocate(p(0:nxd-1), x(0:nxd-1), cb(0:nxd-1),
     &           cc(0:nxd-1),cd(0:nxd-1))

        if (ncall.eq.0) then
         call ipgeti('write_b2_to_astra_iout', iout)                        !srv 14.03.12
         call ipgeti('b2mndr_eirene', use_eirene)
        endif

        do iy = nyd,-1,-1                                                   !srv 21.05.12 {
         if(nncut.eq.2) then
          np=-1
          do ix = rightcut(nncut), rightcut(nncut-1)-1
           np=np+1
           p(np)=ix
           x(np) = crx(ix,iy,0)
           y(np) = cry(ix,iy,0)
          enddo
         elseif(nncut.eq.1) then
          np=-1
          do ix = leftcut(nncut), rightcut(nncut)-1
           np=np+1
           p(np)=ix
           x(np) = crx(ix,iy,0)
           y(np) = cry(ix,iy,0)
          enddo
         else
          np=-1
          do ix = 0, nxd-1
           np=np+1
           p(np)=ix
           x(np) = crx(ix,iy,0)
           y(np) = cry(ix,iy,0)
          enddo
         endif
!        DO IX=0,NP
!         WRITE(*,*) p(IX),x(IX)
!        ENDDO

         CALL DSPLINE(np+1, p, x, CB, CC, CD)
!        DO IX=0,2*np-2
!         WRITE(*,*) IX/2.0_R8,DSEVAL(np,IX/2.0_R8,p,x,CB,CC,CD)
!        ENDDO
         arg = DFMIN(real(nmdpl-5,R8),real(nmdpl+5,R8),fun,1.0e-7_R8)
         R_midpl(iy) = DSEVAL(np+1,arg,p,x,CB,CC,CD)
!        write (*,"(' arg=',es13.6,' R_midpl=',es13.6)") arg, R_midpl(iy)

         CALL DSPLINE(np+1, p, y, CB, CC, CD)
         Z_midpl(iy)=DSEVAL(np+1,arg,p,y,CB,CC,CD)
         write (*,'(a5,es13.6,a9,es13.6,a9,es13.6)')
     &    ' arg=',arg,' R_midpl=',R_midpl(iy),' Z_midpl=',Z_midpl(iy)
        enddo                                                        !srv 21.05.12 }

        maxtopcut = nyd                                              !srv 25.02.16 {
        if (nncut > 0) maxtopcut = maxval(topcut(1:nncut))

        mintopcut = nyd
        if (nncut > 0) then
          if (nnreg(0).eq.7) then
            mintopcut=topcut(1)
          else
            mintopcut=minval(topcut(1:nncut))
          endif
        endif

        iy_Ratb = -1
        do iy = -1,nyd ! was mintopcut    ! for inner surfaces + one open surface
!         R_midpl(iy) = crx(nmdpl,iy,1)
         if(R_midpl(iy).le.Ratb) iy_Ratb = iy
        enddo

        iy_crb_astra = -1
        do iy = -1,nyd-1 ! was mintopcut-1    ! for inner surfaces
         Rc_midpl(iy) = 0.5_R8*(R_midpl(iy)+R_midpl(iy+1))
         if(Rc_midpl(iy).le.crb_astra) iy_crb_astra = iy
        enddo

        open(70,file='Radial_coordinates.dat')
        write(70,*) ' Number','   Coord. cell face      ',
     &              '   Coord. cell centre    '
!        iy = -1
!        write(70,"(i6,1x,2es25.15)")  iy, R_midpl(iy), R_midpl(iy)
        do iy = -1,mintopcut
          write(70,"(i6,1x,2es25.15)") iy, R_midpl(iy), Rc_midpl(iy)
        enddo
!         write(70,"(i6,1x,2es25.15)") iy, R_midpl(mintopcut+1)
        close(70)

!       neb = ne(nmdpl, iy_Ratb)/1.0e19
!       teb = te(nmdpl, iy_Ratb)/qe/1.0e3_R8
!       tib = ti(nmdpl, iy_Ratb)/qe/1.0e3_R8
!       write(nwrite,"(es25.15, 28x,a,f6.4,a)") neb,
!     &   'boundary electron density in Ratb=',Ratb,'      [1/m**3/1e19]'
!       write(nwrite,"(es25.15, 28x,a,f6.4,a)") teb,
!     &   'boundary electron temperature in Ratb=',Ratb,'  [keV]'
!       write(nwrite,"(es25.15, 28x,a,f6.4,a)") tib,
!     &   'boundary ion temperature in Ratb=',Ratb,'       [keV]'

        do is = 0, nsd-1
         do iy = 0,mintopcut   ! for closed surfaces and through separatrix
          sfna(iy,is)     = 0.0_R8
          sfna_mdf(iy,is) = 0.0_R8
          if(nnreg(0).eq.1.and.periodic_bc.eq.0) cycle
          if(.not.is_neutral(is)) then
           do ix=-1,nxd
            if(mod(region(ix,iy,0),4).eq.1 .or.
     &         mod(region(ix,iy,0),4).eq.2) then
              sfna(iy,is) = sfna(iy,is) + fna(ix,iy,1,is)
              sfna_mdf(iy,is) = sfna_mdf(iy,is) + fna_mdf(ix,iy,1,is)
            endif
           enddo
          endif
         enddo

         do iy = mintopcut+1, maxtopcut  ! for open surfaces of DDN geometry  !srv 20.12.12 {
          sfna(iy,is)     = 0.0_R8
          sfna_mdf(iy,is) = 0.0_R8
          do ix=-1,nxd
            if((topcut(1).lt.topcut(2) .and.                    ! DDN Bottom X-point is active
     &         (region(ix,iy,0).eq.3 .or. region(ix,iy,0).eq.2 .or.
     &          region(ix,iy,0).eq.6 .or. region(ix,iy,0).eq.8)) .or.
     &         (topcut(1).gt.topcut(2) .and.                    ! DDN TOP X-point is active
     &         (region(ix,iy,0).eq.2 .or. region(ix,iy,0).eq.4 .or.
     &          region(ix,iy,0).eq.7 .or. region(ix,iy,0).eq.6))) then
             if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nxd+1) then
               sfna(iy,is) = sfna(iy,is) + fna(ix,iy,1,is)
               sfna_mdf(iy,is) = sfna_mdf(iy,is) + fna_mdf(ix,iy,1,is)
             endif
            endif
          enddo
         enddo

         do iy = maxtopcut+1, nyd        ! for open surfaces
          sfna(iy,is)     = 0.0_R8
          sfna_mdf(iy,is) = 0.0_R8
          do ix=-1,nxd
           if(nncut.eq.1 .or.
     &       (nncut.eq.2 .and. (region(ix,iy,0).eq.7 .or.
     &        region(ix,iy,0).eq.6 .or. region(ix,iy,0).eq.8))) then
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nxd+1) then
             sfna(iy,is) = sfna(iy,is) + fna(ix,iy,1,is)
             sfna_mdf(iy,is) = sfna_mdf(iy,is) + fna_mdf(ix,iy,1,is)
            endif
           endif
          enddo                                                          !srv 20.12.12 }
         enddo

#ifdef B25_EIRENE
         if (is_neutral(is).and.use_eirene.ne.0) then
          do isa = 1,nnatmi
           if (b2espcr(is).ne.latmscl(isa)) cycle
           do iy = 0,mintopcut   ! for closed surfaces and through separatrix
            if (nnreg(0).eq.1 .and. periodic_bc.eq.0) cycle
            do ix=-1,nxd
             if(mod(region(ix,iy,0),4).eq.1 .or.
     &          mod(region(ix,iy,0),4).eq.2) then
              sfna(iy,is) = sfna(iy,is) +
     &                      rfluxa(ix+1,iy+1,isa,1)*sy(ix,iy)
              sfna_mdf(iy,is) = sfna_mdf(iy,is) +
     &                      rfluxa(ix+1,iy+1,isa,1)*sy(ix,iy)
              do ism = 1,nnmoli
               if (mlcmp(isa,ism).gt.0) then
                sfna(iy,is) = sfna(iy,is) + mlcmp(isa,ism)*
     &                        rfluxm(ix+1,iy+1,ism,1)*sy(ix,iy)
                sfna_mdf(iy,is) = sfna_mdf(iy,is) + mlcmp(isa,ism)*
     &                        rfluxm(ix+1,iy+1,ism,1)*sy(ix,iy)
               endif
              enddo
             endif
            enddo
           enddo

           do iy = mintopcut+1, maxtopcut  ! for open surfaces of DDN geometry  !srv 20.12.12 {
            do ix=-1,nxd
             if((topcut(1).lt.topcut(2) .and.                    ! DDN Bottom X-point is active
     &         (region(ix,iy,0).eq.3 .or. region(ix,iy,0).eq.2 .or.
     &          region(ix,iy,0).eq.6 .or. region(ix,iy,0).eq.8)) .or.
     &         (topcut(1).gt.topcut(2) .and.                    ! DDN TOP X-point is active
     &         (region(ix,iy,0).eq.2 .or. region(ix,iy,0).eq.4 .or.
     &          region(ix,iy,0).eq.7 .or. region(ix,iy,0).eq.6))) then
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nxd+1)
     >         then
               sfna(iy,is) = sfna(iy,is) +
     &                       rfluxa(ix+1,iy+1,isa,1)*sy(ix,iy)
               sfna_mdf(iy,is) = sfna_mdf(iy,is) +
     &                       rfluxa(ix+1,iy+1,isa,1)*sy(ix,iy)
               do ism = 1,nnmoli
                if (mlcmp(isa,ism).gt.0) then
                 sfna(iy,is) = sfna(iy,is) + mlcmp(isa,ism)*
     &                         rfluxm(ix+1,iy+1,ism,1)*sy(ix,iy)
                 sfna_mdf(iy,is) = sfna_mdf(iy,is) + mlcmp(isa,ism)*
     &                         rfluxm(ix+1,iy+1,ism,1)*sy(ix,iy)
                endif
               enddo
              endif
             endif
            enddo
           enddo

           do iy = maxtopcut+1, nyd        ! for open surfaces
            do ix=-1,nxd
             if(nncut.eq.1 .or.
     &         (nncut.eq.2 .and. (region(ix,iy,0).eq.7 .or.
     &          region(ix,iy,0).eq.6 .or. region(ix,iy,0).eq.8))) then
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nxd+1) then
               sfna(iy,is) = sfna(iy,is) +
     &                       rfluxa(ix+1,iy+1,isa,1)*sy(ix,iy)
               do ism = 1,nnmoli
                if (mlcmp(isa,ism).gt.0) then
                 sfna(iy,is) = sfna(iy,is) + mlcmp(isa,ism)*
     &                         rfluxm(ix+1,iy+1,ism,1)*sy(ix,iy)
                endif
               enddo
              endif
             endif
            enddo
           enddo
          enddo
         endif                                                           !srv 20.12.12 }
#endif
        enddo

        do iy = 0,mintopcut   ! for closed surfaces and through separatrix
          sfne(iy) = 0.0_R8
          sfhe_mdf(iy) = 0.0_R8
          sfhi_mdf(iy) = 0.0_R8
          if (nnreg(0).eq.1.and.periodic_bc.eq.0) cycle
          do ix=-1,nxd
            if(mod(region(ix,iy,0),4).eq.1 .or.
     &         mod(region(ix,iy,0),4).eq.2) then
             sfhe_mdf(iy) = sfhe_mdf(iy) + fhe_mdf(ix,iy,1)
     &       +fac_ExB(ix,iy)*                                    !srv 11.04.13 {
     &      ((po(ix,iy)-po(leftix(ix,iy),iy))/
     &       (hx(ix,iy)+hx(leftix(ix,iy),iy))+
     &       (po(rightix(ix,iy),iy)-po(ix,iy))/
     &       (hx(rightix(ix,iy),iy)+hx(ix,iy)))*
     &        ne(ix,iy)*te(ix,iy)*sy(ix,iy)*
     &        bb(ix,iy,2)/bb(ix,iy,3)**2                         !srv 11.04.13 }

             sfhi_mdf(iy) = sfhi_mdf(iy) + fhi_mdf(ix,iy,1)
     &       +fac_ExB(ix,iy)*                                    !srv 11.04.13 {
     &      ((po(ix,iy)-po(leftix(ix,iy),iy))/
     &       (hx(ix,iy)+hx(leftix(ix,iy),iy))+
     &       (po(rightix(ix,iy),iy)-po(ix,iy))/
     &       (hx(rightix(ix,iy),iy)+hx(ix,iy)))*
     &        ni(ix,iy,1)*ti(ix,iy)*sy(ix,iy)*
     &        bb(ix,iy,2)/bb(ix,iy,3)**2                         !srv 11.04.13 }

             sfne(iy) = sfne(iy) + fne(ix,iy,1)
            endif
          enddo
        enddo

        do iy = mintopcut+1, maxtopcut  ! for open surfaces of DDN geometry  !srv 20.12.12 {
          sfhe_mdf(iy) = 0.0_R8
          sfhi_mdf(iy) = 0.0_R8
          sfne(iy) = 0.0_R8
          do ix=-1,nxd
            if((topcut(1).lt.topcut(2) .and.                    ! DDN Bottom X-point is active
     &         (region(ix,iy,0).eq.3 .or. region(ix,iy,0).eq.2 .or.
     &          region(ix,iy,0).eq.6 .or. region(ix,iy,0).eq.8)) .or.
     &         (topcut(1).gt.topcut(2) .and.                    ! DDN TOP X-point is active
     &         (region(ix,iy,0).eq.2 .or. region(ix,iy,0).eq.4 .or.
     &          region(ix,iy,0).eq.7 .or. region(ix,iy,0).eq.6))) then
             if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nxd+1) then
              sfhe_mdf(iy) = sfhe_mdf(iy) + fhe_mdf(ix,iy,1)
     &       +fac_ExB(ix,iy)*                                    !srv 11.04.13 {
     &      ((po(ix,iy)-po(leftix(ix,iy),iy))/
     &       (hx(ix,iy)+hx(leftix(ix,iy),iy))+
     &       (po(rightix(ix,iy),iy)-po(ix,iy))/
     &       (hx(rightix(ix,iy),iy)+hx(ix,iy)))*
     &        ne(ix,iy)*te(ix,iy)*sy(ix,iy)*
     &        bb(ix,iy,2)/bb(ix,iy,3)**2                         !srv 11.04.13 }

              sfhi_mdf(iy) = sfhi_mdf(iy) + fhi_mdf(ix,iy,1)
     &       +fac_ExB(ix,iy)*                                    !srv 11.04.13 {
     &      ((po(ix,iy)-po(leftix(ix,iy),iy))/
     &       (hx(ix,iy)+hx(leftix(ix,iy),iy))+
     &       (po(rightix(ix,iy),iy)-po(ix,iy))/
     &       (hx(rightix(ix,iy),iy)+hx(ix,iy)))*
     &        ni(ix,iy,1)*ti(ix,iy)*sy(ix,iy)*
     &        bb(ix,iy,2)/bb(ix,iy,3)**2                         !srv 11.04.13 }

              sfne(iy) = sfne(iy) + fne(ix,iy,1)
             endif
            endif
          enddo
        enddo

        do iy = maxtopcut+1, nyd        ! for open surfaces
          sfhe_mdf(iy) = 0.0_R8
          sfhi_mdf(iy) = 0.0_R8
          sfne(iy) = 0.0_R8
          do ix=-1,nxd
            if(nncut.eq.1 .or.
     &         (nncut.eq.2 .and. (region(ix,iy,0).eq.7 .or.
     &            region(ix,iy,0).eq.6 .or. region(ix,iy,0).eq.8))) then
             if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nxd+1) then
              sfhe_mdf(iy) = sfhe_mdf(iy) + fhe_mdf(ix,iy,1)
     &       +fac_ExB(ix,iy)*                                    !srv 11.04.13 {
     &      ((po(ix,iy)-po(leftix(ix,iy),iy))/
     &       (hx(ix,iy)+hx(leftix(ix,iy),iy))+
     &       (po(rightix(ix,iy),iy)-po(ix,iy))/
     &       (hx(rightix(ix,iy),iy)+hx(ix,iy)))*
     &        ne(ix,iy)*te(ix,iy)*sy(ix,iy)*
     &        bb(ix,iy,2)/bb(ix,iy,3)**2                         !srv 11.04.13 }

              sfhi_mdf(iy) = sfhi_mdf(iy) + fhi_mdf(ix,iy,1)
     &       +fac_ExB(ix,iy)*                                    !srv 11.04.13 {
     &      ((po(ix,iy)-po(leftix(ix,iy),iy))/
     &       (hx(ix,iy)+hx(leftix(ix,iy),iy))+
     &       (po(rightix(ix,iy),iy)-po(ix,iy))/
     &       (hx(rightix(ix,iy),iy)+hx(ix,iy)))*
     &        ni(ix,iy,1)*ti(ix,iy)*sy(ix,iy)*
     &        bb(ix,iy,2)/bb(ix,iy,3)**2                         !srv 11.04.13 }

              sfne(iy) = sfne(iy) + fne(ix,iy,1)
             endif
            endif
          enddo
        enddo                                                            !srv 20.12.12 }
*   ..the simplest linear interpolation
        iy =  iy_crb_astra
        neb = (ne(nmdpl,iy)+(crb_astra-Rc_midpl(iy))*
     &        (ne(nmdpl,iy+1)-ne(nmdpl,iy))/
     &        (Rc_midpl(iy+1)-Rc_midpl(iy)))/1.0e19
        teb = (te(nmdpl,iy)+(crb_astra-Rc_midpl(iy))*
     &        (te(nmdpl,iy+1)-te(nmdpl,iy))/
     &        (Rc_midpl(iy+1)-Rc_midpl(iy)))/qe/1.0e3_R8

        tib = (ti(nmdpl,iy)+(crb_astra-Rc_midpl(iy))*
     &        (ti(nmdpl,iy+1)-ti(nmdpl,iy))/
     &        (Rc_midpl(iy+1)-Rc_midpl(iy)))/qe/1.0e3_R8

        write(nwrite,"(es25.15, 28x,a,f6.4,a)") neb,
     &  'boundary electron density in crb_astra=',crb_astra,
     &                                             '     [1/m**3/1e19]'
        write(nwrite,"(es25.15, 28x,a,f6.4,a)") teb,
     &  'boundary electron temperature in crb_astra=',crb_astra,' [keV]'
        write(nwrite,"(es25.15, 28x,a,f6.4,a)") tib,
     &  'boundary ion temperature in crb_astra=',crb_astra,'      [keV]'

        is = ismain0
        write(nwrite,"(es25.15, 28x,a,i4,a)")
     &     -(sfna(0,is)+sfna(1,is))*0.5_R8/1.0e19_R8,
     &    '- integrated neutral flux',is,
     &    '  in cell centre for inner surface [1e19/s]'
        iy = 0                  ! it is written value in 0-cell instead of -1 cell
        write(nwrite,"(i6,1x,5es25.15)") iy-1, se(iy)/svol(iy)/1.0e19_R8 ! electron source [1e19/s]
     &   ,(sfna_mdf(iy+1,ismain)-sfna_mdf(iy,ismain))/svol(iy)/1.0e19_R8
     &   ,(sfna_mdf(iy+1,is    )-sfna_mdf(iy,is    ))/svol(iy)/1.0e19_R8
     &   , seh(iy)/svol(iy), sih(iy)/svol(iy)                            !srv 11.06.13
        do iy = 0,mintopcut-1                                            ! for Core region
         write(nwrite,"(i6,1x,5es25.15)") iy, se(iy)/svol(iy)/1.0e19_R8  ! electron source [1e19/s]
     &   ,(sfna_mdf(iy+1,ismain)-sfna_mdf(iy,ismain))/svol(iy)/1.0e19_R8
     &   ,(sfna_mdf(iy+1,is    )-sfna_mdf(iy,is    ))/svol(iy)/1.0e19_R8
     &   , seh(iy)/svol(iy), sih(iy)/svol(iy)                            !srv 11.06.13
        enddo

        open(70,file='table_01.dat')                                     !srv 19.09.11 {
        write(70,10)
   10   format(4x,'Rc(m)',5x,2x,'ne(1/m^3)',3x,3x,'Te(keV)',4x,
     &         3x,'Ti(keV)',4x,3x,'D(m^2/s)',3x,3x,'xe(m^2/s)',2x,
     &         3x,'xi(m^2/s)',2x,3x,'Vy_AN(m/s)',1x,3x,'Se(1/m^3/s)',
     &         3x,'Seh(      )',3x,'Sih(      )')                        !srv 11.06.13
        ix=nmdpl
        do iy= -1,nyd-1
         write(70,"(11es14.6)")
     &    Rc_midpl(iy), ne(ix,iy), te(ix,iy)/ev/1.0e3_R8,
     &    ti(ix,iy)/ev/1.0e3_R8, dna0(ix,iy,ismain),
     &    hce0(ix,iy)/ne(ix,iy), hcib(ix,iy,ismain)/na(ix,iy,ismain),
     &    vla0(ix,iy,1,ismain), se(iy)/svol(iy),
     &    seh(iy)/svol(iy), sih(iy)/svol(iy)                             !srv 11.06.13

        enddo

        close(70)                                                        !srv 19.09.11 }

        open(70,file='table_02.dat')                                     !srv 19.09.11 {
        write(70,20)
   20   format(4x,'Rf(m)',5x,4x,'Fle(1/s)',2x,4x,'qe(MW)',5x,
     &         4x,'qi(MW)',5x,3x,'-Fln(1/s)',2x)
        ix=nmdpl
        do iy= 0,nyd
         write(70,"(5es14.6)")
     &    R_midpl(iy), sfna_mdf(iy,ismain), sfhe_mdf(iy)/1.0e6_R8,
     &    sfhi_mdf(iy)/1.0e6_R8,-sfna(iy,0)
        enddo

        close(70)

        open(70,file='check_inv_hi_astra.dat')
        write(70,*) ' Number    Coord    inv_hi_astra'
        do iy=-1,nyd
         write(70,40) iy, Rmaj_astra(iy),inv_hi_astra(iy)
        enddo
   40   format(i6,1x,2e13.5)
        close(70)

        write (*,*) 'Ratb, iy_Ratb', Ratb, iy_Ratb
        write (*,*) 'data for b2mod_b2_to_astra has been written'

#ifdef B25_EIRENE
        if(iout.ne.0) then                                         !srv 14.03.12 {
         do iy = -1,nyd
          do ix = -1,nxd
           wrk0(ix,iy) = sy(ix,iy)
          enddo
         enddo
         call my_out(70,nxd,nyd,gs(-1,-1,1),'write_b2_to_astra_gsy')
         call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_sy')
         do isa = 1,nnatmi
          write (chns,'(i3.3)') isa
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = dab2(ix+1,iy+1,isa,1)
           enddo
          enddo
          call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_dab2'//chns)
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = tab2(ix+1,iy+1,isa,1)
           enddo
          enddo
          call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_tab2'//chns)
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = rfluxa(ix+1,iy+1,isa,1)
           enddo
          enddo
          call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_rfluxa'//chns)
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = rfluxa(ix+1,iy+1,isa,1)*sy(ix,iy)
           enddo
          enddo
          call my_out(70,nxd,nyd,wrk0,
     .     'write_b2_to_astra_rfluxasy'//chns)
!         call my_out(70,nxd,nyd,refluxa(0:nxd+1,0:nyd+1,isa,1),
!     &   'write_b2_to_astra_refluxa'//chns)
!         call my_out(70,nxd,nyd,pfluxa(0:nxd+1,0:nyd+1,isa,1),
!     &   'write_b2_to_astra_pfluxa'//chns)
!         call my_out(70,nxd,nyd,pefluxa(0:nxd+1,0:nyd+1,isa,1),
!     &   'write_b2_to_astra_pefluxa'//chns)
         enddo
         do ism = 1,nnmoli
          write (chns,'(i3.3)') ism
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = dmb2(ix+1,iy+1,ism,1)
           enddo
          enddo
          call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_dmb2'//chns)
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = tmb2(ix+1,iy+1,ism,1)
           enddo
          enddo
          call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_tmb2'//chns)
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = rfluxm(ix+1,iy+1,ism,1)
           enddo
          enddo
          call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_rfluxm'//chns)
          do iy = -1,nyd
           do ix = -1,nxd
            wrk0(ix,iy) = rfluxm(ix+1,iy+1,ism,1)*sy(ix,iy)
           enddo
          enddo
         call my_out(70,nxd,nyd,wrk0,'write_b2_to_astra_rfluxmsy'//chns)
!         call my_out(70,nxd,nyd,refluxm(0:nxd+1,0:nyd+1,ism,1),
!     &   'write_b2_to_astra_refluxm'//chns)
!         call my_out(70,nxd,nyd,pfluxm(0:nxd+1,0:nyd+1,ism,1),
!     &   'write_b2_to_astra_pfluxm'//chns)
!         call my_out(70,nxd,nyd,pefluxm(0:nxd+1,0:nyd+1,ism,1),
!     &   'write_b2_to_astra_pefluxm'//chns)
        enddo
       endif                                                       !srv 14.03.12 }
#endif
       ncall = ncall +1

csw 09jan2014 remove contains: fun
c        contains
c
c         function fun(arg)
c          real (kind=R8) :: fun, arg
c          fun = -DSEVAL(np+1,arg,p,x,CB,CC,CD)
c         end function fun

csw 09jan2014 add dealloc
       deallocate(p,x,cb,cc,cd)
       end subroutine write_b2_to_astra

csw 09jan2014 add 'fun' function
       function fun(arg)
       use b2mod_b2_to_astra_com
       real(kind=R8) :: fun,arg,DSEVAL
       fun=-DSEVAL(np+1,arg,p,x,cb,cc,cd)
       end function fun

       subroutine b2_to_astra_src(nxd, nyd, vol, source, svol, src)   !srv 11.06.13
       use b2mod_indirect
       implicit none
!       integrated source   [1/s]
       integer :: nxd, nyd
       real (kind=R8) :: src(-1:nyd), source(-1:nxd,-1:nyd),
     &                   svol(-1:nyd), vol(-1:nxd,-1:nyd)
       integer :: ix, iy
       integer, save :: ncall = 0, mintopcut


       if(ncall.eq.0) then
         mintopcut = nyd                                    !srv 25.02.16
         if (nncut > 0) then
          if (nnreg(0).eq.7) then
            mintopcut=topcut(1)
          else
            mintopcut=minval(topcut(1:nncut))
          endif
         endif
         do iy = 0, mintopcut-1    ! closed surfaces
          svol(iy) = 0.0_R8
          do ix = -1, nxd
           if(on_closed_surface(ix,iy)) then
            svol(iy) = svol(iy) + vol(ix,iy)
           endif
          enddo
         enddo
         do iy = mintopcut,nyd    ! open surfaces
          svol(iy) = 0.0_R8
          do ix = -1, nxd
            svol(iy) = svol(iy) + vol(ix,iy)
          enddo
         enddo
!         write (*,"(a, i5)") '    svol ncall=', ncall
!         write (*,"(es15.7)") svol
       endif

       do iy = 0, mintopcut-1
         do ix = -1, nxd
          if(on_closed_surface(ix,iy)) then   ! closed surfaces
!          src(iy) = src(iy) + vol(ix,iy)*source(ix,iy)
           src(iy) = src(iy) + source(ix,iy)   ! source has been multiplied by vol
          endif
         enddo
       enddo
       do iy = mintopcut,nyd-1   ! open surfaces
         do ix = 0,nxd-1          !   -1, nxd
!          src(iy) = src(iy) + vol(ix,iy)*source(ix,iy)
           src(iy) = src(iy) + source(ix,iy)   ! source has been multiplied by vol
         enddo
       enddo
!         write (*,"(a, i5)") '       src         svol      ncall=', ncall
!         do iy = 0, nyd-1
!          write (*,"(2es15.7)") src(iy),svol(iy)
!        enddo

       ncall = ncall + 1
       end subroutine b2_to_astra_src

      end module b2mod_b2_to_astra                                  !srv 28.07.11 }
