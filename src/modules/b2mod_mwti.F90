Module b2mod_mwti
  use b2mod_types , only : R8
  Implicit None
  Private
  Public :: b2mwti, rwcdf, output_ds, rwcdf_settime
Contains
  
  subroutine b2mwti (itim, tim, nx, ny, ns, ismain, ismain0, BoRiS)
    use b2mod_geo
    use b2mod_plasma
    use b2mod_rates
    use b2mod_residuals
    use b2mod_sources
    use b2mod_transport
    use b2mod_anomalous_transport
    use b2mod_neutrals_namelist
    use b2mod_work
    use b2mod_indirect
    use b2mod_constants
    use b2mod_tallies
    use b2mod_wall
    use b2mod_b2cmpa
    use b2mod_external
#ifdef B25_EIRENE
    use eirmod_extrab25
#endif
    implicit none
    !   ..input arguments (unchanged on exit)
    integer, Intent(In) :: itim, nx, ny, ns, ismain, ismain0
    real (kind=R8), Intent(In) :: tim, BoRiS
    !   ..output arguments (unspecified on entry)
    !     (none)
    !   ..common blocks
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
    !-----------------------------------------------------------------------
    !.documentation
    !
    !  1. purpose
    !
    !     B2MWTI writes to an un*formatted file a selection of data
    !     describing the state of the plasma.
    !
    !
    !  3. description (see also routine b2cdca)
    !
    !     This routine writes to unit (nput) data that are suitable for
    !     producing movie output. Unit nput must be connected to an
    !     un*formatted file. Presently, only the changes in the plasma
    !     state are written out.
    !
    !
    !-----------------------------------------------------------------------
    !.declarations

    !   ..local variables
    integer ncall, ntstep
    real (kind=R8) :: &
         fnixip(nncutmax), feexip(nncutmax), feixip(nncutmax), &
         fnixap(nncutmax), feexap(nncutmax), feixap(nncutmax), &
         pwmxip(nncutmax), fchxip(nncutmax), fchxap(nncutmax), &
         fetxip(nncutmax), fetxap(nncutmax)
    real (kind=R8) :: &
         fniyip(nncutmax), feeyip(nncutmax), feiyip(nncutmax), &
         fniyap(nncutmax), feeyap(nncutmax), feiyap(nncutmax), &
         pwmxap(nncutmax), fchyip(nncutmax), fchyap(nncutmax), &
         fetyip(nncutmax), fetyap(nncutmax)
    real (kind=R8) :: &
         nemxip(nncutmax), temxip(nncutmax), timxip(nncutmax), &
         nemxap(nncutmax), temxap(nncutmax), timxap(nncutmax), &
         pomxip(nncutmax), pomxap(nncutmax), tpmxip(nncutmax), &
         tpmxap(nncutmax)
    real (kind=R8) :: &
         fnisip(nncutmax), feesip(nncutmax), feisip(nncutmax), &
         fnisap(nncutmax), feesap(nncutmax), feisap(nncutmax), &
         fchsip(nncutmax), fchsap(nncutmax), &
         fetsip(nncutmax), fetsap(nncutmax)
    real (kind=R8) :: &
         fnisipp(nncutmax), feesipp(nncutmax), feisipp(nncutmax), &
         fnisapp(nncutmax), feesapp(nncutmax), feisapp(nncutmax), &
         fchsipp(nncutmax), fchsapp(nncutmax), &
         fetsipp(nncutmax), fetsapp(nncutmax)
    real (kind=R8) :: &
         tmne(1),tmte(1),tmti(1),tmvol

    integer iy, ix, ic, is, ixtl, ixtr, jsep
    integer jxi, jxa, target_offset, ix_off
    integer iyastrt, iyistrt, iylstrt, iyrstrt, iytlstrt, iytrstrt, &
         iyaend,  iyiend,  iylend,  iyrend,  iytlend,  iytrend, &
         nybl, nybr, nytl, nytr, nya, nyi, nc
    !   ..procedures
    external subini, subend, xertst, ipgeti
    Real(kind=R8) :: fnitmp,feetmp,feitmp,fchtmp,fettmp,pwrtmp

#ifndef NO_CDF
    integer imap(maxvdims), iret, ncid
    real (kind=R8) :: &
         nesepi(nncutmax), tesepi(nncutmax), tisepi(nncutmax), &
         nesepm(nncutmax), tesepm(nncutmax), tisepm(nncutmax), &
         nesepa(nncutmax), tesepa(nncutmax), tisepa(nncutmax), &
         posepi(nncutmax), posepm(nncutmax), posepa(nncutmax), &
         dnsepm(nncutmax), dpsepm(nncutmax), kesepm(nncutmax), &
         kisepm(nncutmax), vxsepm(nncutmax), vysepm(nncutmax), &
         vssepm(nncutmax), tpsepi(nncutmax), tpsepa(nncutmax)
    real (kind=R8) :: &
         tmhacore(1), tmhasol(1), tmhadiv(1), slice(-1:ny), tstepn(1)
    real (kind=R8) :: &
         timesa(1)
    character*5 rw
#endif
    !   ..initialisation
    save ncall, ntstep, jxi, jxa, jsep, ixtl, ixtr, target_offset, &
         iyastrt, iyistrt, iylstrt, iyrstrt, iytlstrt, iytrstrt, &
         iyaend,  iyiend,  iylend,  iyrend,  iytlend,  iytrend, &
         nybl, nybr, nytl, nytr, nya, nyi, nc
    data ncall/0/,target_offset/1/

    ! Statement functions
    real (kind=R8) :: &
         hy1
    hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)            

    !-----------------------------------------------------------------------
    !.computation

    ! ..preliminaries
    !   ..subprogram start-up calls
    call subini ('b2mwti')
    !     ..test nx, ny
    call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
    call xertst (1.le.ns, 'faulty argument ns')
    call xertst (0.le.ismain.and.ismain.lt.ns, &
         'invalid main plasma species index ismain')
    !   ..extensive tests on first few calls
    if (ncall.eq.0) then
      !   ..test state
      call get_jsep(nx,ny,jxi,jxa,jsep)
      call ipgeti ('b2mwti_target_offset',target_offset)
      call xertst (0.le.target_offset.and.target_offset.le.1, &
           'faulty internal parameter target_offset')
      write(*,*) 'target_offset ', target_offset
      call output_ds(crx,cry,nx,ny, -1,+target_offset, &
           jsep,iylstrt,iylend,'../dsl')
      call output_ds(crx,cry,nx,ny,jxi,0, &
           jsep,iyistrt,iyiend,'../dsi')
      call output_ds(crx,cry,nx,ny,jxa,0, &
           jsep,iyastrt,iyaend,'../dsa')
      call output_ds(crx,cry,nx,ny, nx,-target_offset, &
           jsep,iyrstrt,iyrend,'../dsr')
      if (nnreg(0).eq.8) then
        ixtl = 0
        do while (rightix(ixtl,max(topcut(1),topcut(2))).ne.nx+1 &
             .and.ixtl.lt.nx)
          ixtl=ixtl+1
        enddo
        ixtr = ixtl
        do while (leftix(ixtr,max(topcut(1),topcut(2))).ne.-2 .and. &
             ixtr.lt.nx)
          ixtr=ixtr+1
        enddo
        call output_ds(crx,cry,nx,ny,ixtl,-target_offset, &
             jsep,iytlstrt,iytlend,'../dstl')
        call output_ds(crx,cry,nx,ny,ixtr,+target_offset, &
             jsep,iytrstrt,iytrend,'../dstr')
        nytl = iytlend - iytlstrt + 1
        nytr = iytrend - iytrstrt + 1
      else
        nytl = 0
        nytr = 0
      endif
      nybl = iylend  - iylstrt  + 1
      nybr = iyrend  - iyrstrt  + 1
      nya  = iyaend  - iyastrt  + 1
      nyi  = iyiend  - iyistrt  + 1
      nc = max(nncut,1)
      ! Target areas
      open(99,file='../dsL')
      do iy=-1,ny
        if(region(-1,iy,0).ne.0) &
             write(99,*) gs(rightix(-1,iy),rightiy(-1,iy),0)
      enddo
      close(99)
      open(99,file='../dsR')
      do iy=-1,ny
        if(region(nx,iy,0).ne.0) &
             write(99,*) gs(nx,iy,0)
      enddo
      close(99)
      if (nnreg(0).eq.8) then
        open(99,file='../dsTL')
        do iy=iytlstrt,iytlend
          write(99,*) gs(ixtl,iy,0)
        enddo
        close(99)
        open(99,file='../dsTR')
        do iy=iytrstrt,iytrend
          write(99,*) gs(rightix(ixtr,iy),rightiy(ixtr,iy),0)
        enddo
        close(99)
      endif
      ! Poloidal contact areas
      open(99,file='../dsLP')
      do iy=-1,ny
        if(region(-1,iy,0).ne.0) &
             write(99,*) gs(rightix(-1,iy),rightiy(-1,iy),0) &
             *qc(rightix(-1,iy),rightiy(-1,iy))
      enddo
      close(99)
      open(99,file='../dsRP')
      do iy=-1,ny
        if(region(nx,iy,0).ne.0) &
             write(99,*) gs(nx,iy,0)*qc(nx,iy)
      enddo
      close(99)
      if (nnreg(0).eq.8) then
        open(99,file='../dsTLP')
        do iy=iytlstrt,iytlend
          write(99,*) gs(ixtl,iy,0)*qc(ixtl,iy)
        enddo
        close(99)
        open(99,file='../dsTRP')
        do iy=iytrstrt,iytrend
          write(99,*) gs(rightix(ixtr,iy),rightiy(ixtr,iy),0) &
               *qc(rightix(ixtr,iy),rightiy(ixtr,iy))
        enddo
        close(99)
      endif
#ifndef NO_CDF
      call b2crtimecdf(nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, iret)
      rw='write'
      iret = nf_open('b2time.nc',NCWRITE,ncid)
      ntstep=0
      imap(1)=1
      tstepn(1) = ntstep
      call rwcdf (rw, ncid, 'ntstep', imap, tstepn, iret)
      iret = nf_close(ncid)
#endif
    endif! ncall == 0
    if(ncall.lt.3) then
      call xertst (0.le.itim, 'faulty parameter itim')
      call xertst (0.0_R8.le.tim, 'faulty parameter tim')
    endif
    !   ..compute change in plasma state
    !
    ntstep = ntstep + 1
    !     write(*,*) 'ntstep = ',ntstep
#ifndef NO_CDF
    rw = 'write'
    iret = nf_open('b2time.nc', NCWRITE, ncid)
    call rwcdf_settime ('time', ntstep)
    timesa(1) = tim
#endif
    !
    !    total flows to the divertor plates
    !
            
    fnixip = 0.0_R8; feexip = 0.0_R8; feixip = 0.0_R8; fchxip = 0.0_R8; fetxip = 0.0_R8
    nemxip = 0.0_R8; temxip = 0.0_R8; timxip = 0.0_R8; pomxip = 0.0_R8; pwmxip = 0.0_R8; tpmxip = 0.0_R8   
    ix = -1 ! 1
    ix_off  = ix + target_offset
    Do iy = iylstrt,iylend
      Call calc_fet(ix,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
      fnixip(1) = fnixip(1) + fnitmp
      feexip(1) = feexip(1) + feetmp
      feixip(1) = feixip(1) + feitmp
      fchxip(1) = fchxip(1) + fchtmp
      fetxip(1) = fetxip(1) + fettmp
      nemxip(1) = max(nemxip(1), ne(ix_off,iy))
      temxip(1) = max(temxip(1), te(ix_off,iy))
      timxip(1) = max(timxip(1), ti(ix_off,iy))
      pomxip(1) = max(pomxip(1), po(ix_off,iy))
      pwmxip(1) = max(pwmxip(1), pwrtmp)
      if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) Then
        tpmxip(1) = max(tpmxip(1), target_temp(xymap(ix,iy),1))
      Endif      
    Enddo
    
    fnixap = 0.0_R8; feexap = 0.0_R8; feixap = 0.0_R8; fchxap = 0.0_R8; fetxap = 0.0_R8
    nemxap = 0.0_R8; temxap = 0.0_R8; timxap = 0.0_R8; pomxap = 0.0_R8; pwmxap = 0.0_R8; tpmxap = 0.0_R8
    ix = nx ! 2
    ix_off  = ix - target_offset
    Do iy = iyrstrt,iyrend
      Call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
      fnixap(1) = fnixap(1) + fnitmp
      feexap(1) = feexap(1) + feetmp
      feixap(1) = feixap(1) + feitmp
      fchxap(1) = fchxap(1) + fchtmp
      fetxap(1) = fetxap(1) + fettmp
      nemxap(1) = max(nemxap(1), ne(ix_off,iy))
      temxap(1) = max(temxap(1), te(ix_off,iy))
      timxap(1) = max(timxap(1), ti(ix_off,iy))
      pomxap(1) = max(pomxap(1), po(ix_off,iy))
      pwmxap(1) = max(pwmxap(1), pwrtmp)
      if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) Then
        tpmxap(1) = max(tpmxap(1), target_temp(xymap(ix,iy),1))
      Endif      
    Enddo
    
    if(nncut.ge.2) then
      ix = ixtr ! 3
      ix_off  = ix + target_offset
      Do iy = iytrstrt,iytrend
        Call calc_fet(ix,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
        fnixap(2) = fnixap(2) + fnitmp
        feexap(2) = feexap(2) + feetmp
        feixap(2) = feixap(2) + feitmp
        fchxap(2) = fchxap(2) + fchtmp
        fetxap(2) = fetxap(2) + fettmp
        nemxap(2) = max(nemxap(2), ne(ix_off,iy))
        temxap(2) = max(temxap(2), te(ix_off,iy))
        timxap(2) = max(timxap(2), ti(ix_off,iy))
        pomxap(2) = max(pomxap(2), po(ix_off,iy))
        pwmxap(2) = max(pwmxap(2), pwrtmp)
        if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) Then
          tpmxap(2) = max(tpmxap(2), target_temp(xymap(ix,iy),1))
        Endif
      Enddo
      
      ix = ixtl ! 4
      ix_off  = ix - target_offset
      Do iy = iytlstrt,iytlend
        Call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp,pwrtmp)
        fnixip(2) = fnixip(2) + fnitmp
        feexip(2) = feexip(2) + feetmp
        feixip(2) = feixip(2) + feitmp
        fchxip(2) = fchxip(2) + fchtmp
        fetxip(2) = fetxip(2) + fettmp
        nemxip(2) = max(nemxip(2), ne(ix_off,iy))
        temxip(2) = max(temxip(2), te(ix_off,iy))
        timxip(2) = max(timxip(2), ti(ix_off,iy))
        pomxip(2) = max(pomxip(2), po(ix_off,iy))
        pwmxip(2) = max(pwmxip(2), pwrtmp)
        if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and. xymap(ix,iy).ne.0) Then
          tpmxip(2) = max(tpmxip(2), target_temp(xymap(ix,iy),1))
        Endif
      Enddo      
    endif

    fnisip = 0.0_R8; feesip = 0.0_R8; feisip = 0.0_R8; fchsip = 0.0_R8; fetsip = 0.0_R8       
    fnisap = 0.0_R8; feesap = 0.0_R8; feisap = 0.0_R8; fchsap = 0.0_R8; fetsap = 0.0_R8    
    fnisipp = 0.0_R8; feesipp = 0.0_R8; feisipp = 0.0_R8; fetsipp = 0.0_R8; fchsipp = 0.0_R8
    fnisapp = 0.0_R8; feesapp = 0.0_R8; feisapp = 0.0_R8; fetsapp = 0.0_R8; fchsapp = 0.0_R8    
    if(nnreg(0).ge.3) then
      Do ic = 1, nncut
        Do iy = -1,jsep
          ix = leftcut(ic) ! 5
          Call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisipp(ic) = fnisipp(ic) + fnitmp
          feesipp(ic) = feesipp(ic) + feetmp
          feisipp(ic) = feisipp(ic) + feitmp
          fchsipp(ic) = fchsipp(ic) + fchtmp
          fetsipp(ic) = fetsipp(ic) + fettmp
          ix = rightcut(ic) ! 6
          Call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisapp(ic) = fnisapp(ic) + fnitmp
          feesapp(ic) = feesapp(ic) + feetmp
          feisapp(ic) = feisapp(ic) + feitmp
          fchsapp(ic) = fchsapp(ic) + fchtmp
          fetsapp(ic) = fetsapp(ic) + fettmp
        Enddo
        Do iy = jsep+1,ny
          ix = leftcut(ic) ! 7
          Call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisip(ic) = fnisip(ic) + fnitmp
          feesip(ic) = feesip(ic) + feetmp
          feisip(ic) = feisip(ic) + feitmp
          fchsip(ic) = fchsip(ic) + fchtmp
          fetsip(ic) = fetsip(ic) + fettmp
          ix = rightcut(ic) ! 8
          Call calc_fet(ix,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fnisap(ic) = fnisap(ic) + fnitmp
          feesap(ic) = feesap(ic) + feetmp
          feisap(ic) = feisap(ic) + feitmp
          fchsap(ic) = fchsap(ic) + fchtmp
          fetsap(ic) = fetsap(ic) + fettmp
        Enddo
      Enddo
    endif

    fniyip = 0.0_R8; feeyip = 0.0_R8; feiyip = 0.0_R8; fetyip = 0.0_R8; fchyip = 0.0_R8
    fniyap = 0.0_R8; feeyap = 0.0_R8; feiyap = 0.0_R8; fetyap = 0.0_R8; fchyap = 0.0_R8

    if(nnreg(0).eq.4) then
      Do ix = -1,nx
        if(region(ix,ny,0).eq.2) then
          iy = ny ! 9
          Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp
        endif
        if(region(ix,ny,0).eq.3.or.region(ix,ny,0).eq.4) then
          iy = ny ! 10
          Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,-1,0).eq.3.or.region(ix,-1,0).eq.4) then
          iy = -1 ! 11
          Call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
      Enddo
    elseif(nnreg(0).eq.5) then
      Do ix = -1,nx
        if(region(ix,ny,0).eq.5) then
          iy = ny ! 12
          Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp          
        endif
        if(region(ix,-1,0).eq.3.or.region(ix,-1,0).eq.4) then
          iy = -1 ! 13
          Call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
      Enddo
    elseif (nnreg(0).eq.8) then
      Do ix = -1,nx
        if(mod(region(ix,ny,0),4).eq.2) then
          iy = ny ! 14
          Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(region(ix,ny,0)/4+1) = fniyip(region(ix,ny,0)/4+1) + fnitmp
          feeyip(region(ix,ny,0)/4+1) = feeyip(region(ix,ny,0)/4+1) + feetmp
          feiyip(region(ix,ny,0)/4+1) = feiyip(region(ix,ny,0)/4+1) + feitmp
          fchyip(region(ix,ny,0)/4+1) = fchyip(region(ix,ny,0)/4+1) + fchtmp
          fetyip(region(ix,ny,0)/4+1) = fetyip(region(ix,ny,0)/4+1) + fettmp          
        endif
        if(region(ix,ny,0).eq.3 .or. region(ix,ny,0).eq.8) then
          iy = ny ! 15
          Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,-1,0).eq.3 .or. region(ix,-1,0).eq.8) then
          iy = -1 ! 16
          Call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(1) = fniyap(1) + fnitmp
          feeyap(1) = feeyap(1) + feetmp
          feiyap(1) = feiyap(1) + feitmp
          fchyap(1) = fchyap(1) + fchtmp
          fetyap(1) = fetyap(1) + fettmp
        endif
        if(region(ix,ny,0).eq.4 .or. region(ix,ny,0).eq.7) then
          iy = ny ! 17
          Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(2) = fniyap(2) + fnitmp
          feeyap(2) = feeyap(2) + feetmp
          feiyap(2) = feiyap(2) + feitmp
          fchyap(2) = fchyap(2) + fchtmp
          fetyap(2) = fetyap(2) + fettmp
        endif
        if(region(ix,-1,0).eq.4 .or. region(ix,-1,0).eq.7) then
          iy = -1 ! 18
          Call calc_fet(ix,iy,'B',-1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyap(2) = fniyap(2) + fnitmp
          feeyap(2) = feeyap(2) + feetmp
          feiyap(2) = feiyap(2) + feitmp
          fchyap(2) = fchyap(2) + fchtmp
          fetyap(2) = fetyap(2) + fettmp
        endif
      Enddo
    else if (nnreg(0).eq.2) then
      Do ix = -1,nx
        if(region(ix,ny,0).eq.2) then
          iy = ny ! 19
          Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
          fniyip(1) = fniyip(1) + fnitmp
          feeyip(1) = feeyip(1) + feetmp
          feiyip(1) = feiyip(1) + feitmp
          fchyip(1) = fchyip(1) + fchtmp
          fetyip(1) = fetyip(1) + fettmp          
        endif
      Enddo
    else
      Do ix = -1,nx
        iy = ny ! 20
        Call calc_fet(ix,iy,'T',1._R8,nx,ny,ns,ismain,BoRiS,fettmp,fnitmp,feetmp,feitmp,fchtmp)
        fniyip(1) = fniyip(1) + fnitmp
        feeyip(1) = feeyip(1) + feetmp
        feiyip(1) = feiyip(1) + feitmp
        fchyip(1) = fchyip(1) + fchtmp
        fetyip(1) = fetyip(1) + fettmp          
      Enddo
    endif ! nnreg check
    
    !
    !    other quantities related to the target plates
    !
#ifndef NO_CDF
    if(nnreg(0).ne.2) then
      nesepi(1) = 0.5_R8 * &
           (ne(-1+target_offset,jsep)+ &
           ne(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
      tesepi(1) = 0.5_R8/ev * &
           (te(-1+target_offset,jsep)+ &
           te(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
      tisepi(1) = 0.5_R8/ev * &
           (ti(-1+target_offset,jsep)+ &
           ti(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
      if(xymap(-1,jsep).gt.0 .and. &
           xymap(topix(-1,jsep),topiy(-1,jsep)).gt.0) then
        tpsepi(1) = 0.5_R8 * &
             (target_temp(xymap(-1,jsep),1)+ &
             target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1))
      else
        tpsepi(1) = 0.0_R8
      endif
      posepi(1) = 0.5_R8 * &
           (po(-1+target_offset,jsep)+ &
           po(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
    else
      nesepi(1) = &
           ne(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
      tesepi(1) = 1.0_R8/ev * &
           te(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
      tisepi(1) = 1.0_R8/ev * &
           ti(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
      if(xymap(topix(-1,jsep),topiy(-1,jsep)).gt.0) then
        tpsepi(1) = &
             target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1)
      else
        tpsepi(1) = 0.0_R8
      endif
      posepi(1) = &
           po(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
    endif
    nesepm(1) = 0.5_R8 * (ne(jxa,jsep)+ &
         ne(topix(jxa,jsep),topiy(jxa,jsep)))
    tesepm(1) = 0.5_R8 * (te(jxa,jsep)+ &
         te(topix(jxa,jsep),topiy(jxa,jsep)))/ev
    tisepm(1) = 0.5_R8 * (ti(jxa,jsep)+ &
         ti(topix(jxa,jsep),topiy(jxa,jsep)))/ev
    posepm(1) = 0.5_R8 * (po(jxa,jsep)+ &
         po(topix(jxa,jsep),topiy(jxa,jsep)))
    dnsepm(1) = 0.5_R8 * (dna0(jxa,jsep,ismain)+ &
         dna0(topix(jxa,jsep),topiy(jxa,jsep),ismain))
    dpsepm(1) = 0.5_R8 * ( &
         dpa0(jxa,jsep,ismain0)*( &
         rza(jxa,jsep,ismain0)*te(jxa,jsep)+ti(jxa,jsep))+ &
         dpa0(topix(jxa,jsep),topiy(jxa,jsep),ismain0)*( &
         rza(topix(jxa,jsep),topiy(jxa,jsep),ismain0)* &
         te(topix(jxa,jsep),topiy(jxa,jsep))+ &
         ti(topix(jxa,jsep),topiy(jxa,jsep))))
    kesepm(1) = 0.5_R8 * (hce0(jxa,jsep)/ne(jxa,jsep)+ &
         hce0(topix(jxa,jsep),topiy(jxa,jsep))/ &
         ne(topix(jxa,jsep),topiy(jxa,jsep)))
    kisepm(1) = 0.5_R8 * (hci0(jxa,jsep)/ni(jxa,jsep,0)+ &
         hci0(topix(jxa,jsep),topiy(jxa,jsep))/ &
         ni(topix(jxa,jsep),topiy(jxa,jsep),0))
    vxsepm(1) = 0.5_R8 * (vla0(jxa,jsep,0,ismain)+ &
         vla0(topix(jxa,jsep),topiy(jxa,jsep),0,ismain))
    vysepm(1) = 0.5_R8 * (vla0(jxa,jsep,1,ismain)+ &
         vla0(topix(jxa,jsep),topiy(jxa,jsep),1,ismain))
    vssepm(1) = 0.5_R8 * ( &
         vsa0(jxa,jsep,ismain)/(mp*am(ismain)*na(jxa,jsep,ismain))+ &
         vsa0(topix(jxa,jsep),topiy(jxa,jsep),ismain)/ &
         (mp*am(ismain)*na(topix(jxa,jsep),topiy(jxa,jsep),ismain)))
    if(nnreg(0).ne.2) then
      nesepa(1) = 0.5_R8 * &
           (ne(nx-target_offset,jsep)+ &
           ne(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
      tesepa(1) = 0.5_R8/ev * &
           (te(nx-target_offset,jsep)+ &
           te(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
      tisepa(1) = 0.5_R8/ev * &
           (ti(nx-target_offset,jsep)+ &
           ti(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
      if(xymap(nx,jsep).gt.0 .and. &
           xymap(topix(nx,jsep),topiy(nx,jsep)).ge.0) then
        tpsepa(1) = 0.5_R8 * &
             (target_temp(xymap(nx,jsep),1)+ &
             target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1))
      else
        tpsepa(1) = 0.0_R8
      endif
      posepa(1) = 0.5_R8 * &
           (po(nx-target_offset,jsep)+ &
           po(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
    else
      nesepa(1) = &
           ne(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
      tesepa(1) = 1.0_R8/ev * &
           te(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
      tisepa(1) = 1.0_R8/ev * &
           ti(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
      if(xymap(topix(nx,jsep),topiy(nx,jsep)).gt.0) then
        tpsepa(1) = &
             target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1)
      else
        tpsepa(1) = 0.0
      endif
      posepa(1) = &
           po(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
    endif
    if(nncut.eq.2) then
      nesepa(2) = 0.5_R8 * (ne(ixtr+target_offset,jsep)+ &
           ne(topix(ixtr+target_offset,jsep), &
           topiy(ixtr+target_offset,jsep)))
      tesepa(2) = 0.5_R8 * (te(ixtr+target_offset,jsep)+ &
           te(topix(ixtr+target_offset,jsep), &
           topiy(ixtr+target_offset,jsep)))/ev
      tisepa(2) = 0.5_R8 * (ti(ixtr+target_offset,jsep)+ &
           ti(topix(ixtr+target_offset,jsep), &
           topiy(ixtr+target_offset,jsep)))/ev
      if(xymap(ixtr,jsep).gt.0 .and. &
           xymap(topix(ixtr,jsep),topiy(ixtr,jsep)).gt.0) then
        tpsepa(2) = 0.5_R8 * &
             (target_temp(xymap(ixtr,jsep),1)+ &
             target_temp(xymap(topix(ixtr,jsep),topiy(ixtr,jsep)),1))
      else
        tpsepa(2) = 0.0_R8
      endif
      posepa(2) = 0.5_R8 * (po(ixtr+target_offset,jsep)+ &
           po(topix(ixtr+target_offset,jsep), &
           topiy(ixtr+target_offset,jsep)))
      nesepm(2) = 0.5_R8 * (ne(jxi,jsep)+ &
           ne(topix(jxi,jsep),topiy(jxi,jsep)))
      tesepm(2) = 0.5_R8 * (te(jxi,jsep)+ &
           te(topix(jxi,jsep),topiy(jxi,jsep)))/ev
      tisepm(2) = 0.5_R8 * (ti(jxi,jsep)+ &
           ti(topix(jxi,jsep),topiy(jxi,jsep)))/ev
      posepm(2) = 0.5_R8 * (po(jxi,jsep)+ &
           po(topix(jxi,jsep),topiy(jxi,jsep)))
      dnsepm(2) = 0.5_R8 * (dna0(jxi,jsep,ismain)+ &
           dna0(topix(jxi,jsep),topiy(jxi,jsep),ismain))
      dpsepm(2) = 0.5_R8 * ( &
           dpa0(jxi,jsep,ismain0)*( &
           rza(jxi,jsep,ismain0)*te(jxi,jsep)+ti(jxi,jsep))+ &
           dpa0(topix(jxi,jsep),topiy(jxi,jsep),ismain0)*( &
           rza(topix(jxi,jsep),topiy(jxi,jsep),ismain0)* &
           te(topix(jxi,jsep),topiy(jxi,jsep))+ &
           ti(topix(jxi,jsep),topiy(jxi,jsep))))
      kesepm(2) = 0.5_R8 * (hce0(jxi,jsep)/ne(jxi,jsep)+ &
           hce0(topix(jxi,jsep),topiy(jxi,jsep))/ &
           ne(topix(jxi,jsep),topiy(jxi,jsep)))
      kisepm(2) = 0.5_R8 * (hci0(jxi,jsep)/ni(jxi,jsep,0)+ &
           hci0(topix(jxi,jsep),topiy(jxi,jsep))/ &
           ni(topix(jxi,jsep),topiy(jxi,jsep),0))
      vxsepm(2) = 0.5_R8 * (vla0(jxi,jsep,0,ismain)+ &
           vla0(topix(jxi,jsep),topiy(jxi,jsep),0,ismain))
      vysepm(2) = 0.5_R8 * (vla0(jxi,jsep,1,ismain)+ &
           vla0(topix(jxi,jsep),topiy(jxi,jsep),1,ismain))
      vssepm(2) = 0.5_R8 * ( &
           vsa0(jxi,jsep,ismain)/(mp*am(ismain)*na(jxi,jsep,ismain))+ &
           vsa0(topix(jxi,jsep),topiy(jxi,jsep),ismain)/ &
           (mp*am(ismain)*na(topix(jxi,jsep),topiy(jxi,jsep),ismain)))
      nesepi(2) = 0.5_R8 * (ne(ixtl-target_offset,jsep)+ &
           ne(topix(ixtl-target_offset,jsep), &
           topiy(ixtl-target_offset,jsep)))
      tesepi(2) = 0.5_R8 * (te(ixtl-target_offset,jsep)+ &
           te(topix(ixtl-target_offset,jsep), &
           topiy(ixtl-target_offset,jsep)))/ev
      tisepi(2) = 0.5_R8 * (ti(ixtl-target_offset,jsep)+ &
           ti(topix(ixtl-target_offset,jsep), &
           topiy(ixtl-target_offset,jsep)))/ev
      if(xymap(ixtl,jsep).gt.0 .and. &
           xymap(topix(ixtl,jsep),topiy(ixtl,jsep)).gt.0) then
        tpsepi(2) = 0.5_R8 * &
             (target_temp(xymap(ixtl,jsep),1)+ &
             target_temp(xymap(topix(ixtl,jsep),topiy(ixtl,jsep)),1))
      else
        tpsepi(2) = 0.0_R8
      endif
      posepi(2) = 0.5_R8 * (po(ixtl-target_offset,jsep)+ &
           po(topix(ixtl-target_offset,jsep), &
           topiy(ixtl-target_offset,jsep)))
    endif
#endif
    !
    temxip(1:nc) = temxip(1:nc)/ev
    timxip(1:nc) = timxip(1:nc)/ev
    temxap(1:nc) = temxap(1:nc)/ev
    timxap(1:nc) = timxap(1:nc)/ev

    tmne(1)=0.0_R8
    tmte(1)=0.0_R8
    tmti(1)=0.0_R8
    tmvol=0.0_R8
    do iy=-1,ny
      do ix=-1,nx
        if(region(ix,iy,0).ne.0) then
          tmne(1)=tmne(1)+ne(ix,iy)*vol(ix,iy)
          tmte(1)=tmte(1)+te(ix,iy)*ne(ix,iy)*vol(ix,iy)
          tmti(1)=tmti(1)+ti(ix,iy)*ni(ix,iy,0)*vol(ix,iy)
          tmvol=tmvol+vol(ix,iy)
        endif
      enddo
    enddo
    tmne(1)=tmne(1)
    tmte(1)=tmte(1)/ev
    tmti(1)=tmti(1)/ev

#ifndef NO_CDF      
    tmhacore(1)=0.0_R8
    tmhasol(1)=0.0_R8
    tmhadiv(1)=0.0_R8
#ifdef B25_EIRENE
    ! note no emission in guard cells and offset of array by 1
    do iy=-1,ny
      do ix=-1,nx
        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and. &
             bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 ) then
          if(mod(region(ix,iy,0),4).eq.1) then
            tmhacore(1)=tmhacore(1)+ &
                 (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
          elseif(mod(region(ix,iy,0),4).eq.3 .or. &
               (mod(region(ix,iy,0),4).eq.0 .and. &
               region(ix,iy,0).ne.0)) then
            tmhadiv(1)=tmhadiv(1)+ &
                 (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
          elseif(mod(region(ix,iy,0),4).eq.2) then
            tmhasol(1)=tmhasol(1)+ &
                 (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
          elseif(region(ix,iy,0).ne.0) then
            write(*,*) 'b2mwti: unknown region @ ', &
                 ix,iy,region(ix,iy,0)
          endif
        endif
      enddo
    enddo
#endif

    imap(1)=1
    tstepn(1) = ntstep
    call rwcdf(rw,ncid,'ntstep',imap,tstepn,iret)
    call rwcdf(rw,ncid,'timesa',imap,timesa,iret)

    imap(1)=1
    imap(2)=1
    call rwcdf(rw,ncid,'fnixip',imap,fnixip,iret)
    call rwcdf(rw,ncid,'feexip',imap,feexip,iret)
    call rwcdf(rw,ncid,'feixip',imap,feixip,iret)
    call rwcdf(rw,ncid,'fetxip',imap,fetxip,iret)
    call rwcdf(rw,ncid,'fchxip',imap,fchxip,iret)
    call rwcdf(rw,ncid,'fnixap',imap,fnixap,iret)
    call rwcdf(rw,ncid,'feexap',imap,feexap,iret)
    call rwcdf(rw,ncid,'feixap',imap,feixap,iret)
    call rwcdf(rw,ncid,'fetxap',imap,fetxap,iret)
    call rwcdf(rw,ncid,'fchxap',imap,fchxap,iret)

    call rwcdf(rw,ncid,'nesepi',imap,nesepi,iret)
    call rwcdf(rw,ncid,'tesepi',imap,tesepi,iret)
    call rwcdf(rw,ncid,'tisepi',imap,tisepi,iret)
    call rwcdf(rw,ncid,'posepi',imap,posepi,iret)
    call rwcdf(rw,ncid,'nesepm',imap,nesepm,iret)
    call rwcdf(rw,ncid,'tesepm',imap,tesepm,iret)
    call rwcdf(rw,ncid,'tisepm',imap,tisepm,iret)
    call rwcdf(rw,ncid,'posepm',imap,posepm,iret)
    call rwcdf(rw,ncid,'dnsepm',imap,dnsepm,iret)
    call rwcdf(rw,ncid,'dpsepm',imap,dpsepm,iret)
    call rwcdf(rw,ncid,'kesepm',imap,kesepm,iret)
    call rwcdf(rw,ncid,'kisepm',imap,kisepm,iret)
    call rwcdf(rw,ncid,'vxsepm',imap,vxsepm,iret)
    call rwcdf(rw,ncid,'vysepm',imap,vysepm,iret)
    call rwcdf(rw,ncid,'vssepm',imap,vssepm,iret)
    call rwcdf(rw,ncid,'nesepa',imap,nesepa,iret)
    call rwcdf(rw,ncid,'tesepa',imap,tesepa,iret)
    call rwcdf(rw,ncid,'tisepa',imap,tisepa,iret)
    call rwcdf(rw,ncid,'posepa',imap,posepa,iret)
    call rwcdf(rw,ncid,'nemxip',imap,nemxip,iret)
    call rwcdf(rw,ncid,'temxip',imap,temxip,iret)
    call rwcdf(rw,ncid,'timxip',imap,timxip,iret)
    call rwcdf(rw,ncid,'pomxip',imap,pomxip,iret)
    call rwcdf(rw,ncid,'nemxap',imap,nemxap,iret)
    call rwcdf(rw,ncid,'temxap',imap,temxap,iret)
    call rwcdf(rw,ncid,'timxap',imap,timxap,iret)
    call rwcdf(rw,ncid,'pomxap',imap,pomxap,iret)
    call rwcdf(rw,ncid,'fniyip',imap,fniyip,iret)
    call rwcdf(rw,ncid,'feeyip',imap,feeyip,iret)
    call rwcdf(rw,ncid,'feiyip',imap,feiyip,iret)
    call rwcdf(rw,ncid,'fetyip',imap,fetyip,iret)
    call rwcdf(rw,ncid,'fchyip',imap,fchyip,iret)
    call rwcdf(rw,ncid,'fniyap',imap,fniyap,iret)
    call rwcdf(rw,ncid,'feeyap',imap,feeyap,iret)
    call rwcdf(rw,ncid,'feiyap',imap,feiyap,iret)
    call rwcdf(rw,ncid,'fetyap',imap,fetyap,iret)
    call rwcdf(rw,ncid,'fchyap',imap,fchyap,iret)
    call rwcdf(rw,ncid,'pwmxip',imap,pwmxip,iret)
    call rwcdf(rw,ncid,'pwmxap',imap,pwmxap,iret)

    imap(1)=1
    call rwcdf(rw,ncid,'tmne',imap,tmne,iret)
    call rwcdf(rw,ncid,'tmte',imap,tmte,iret)
    call rwcdf(rw,ncid,'tmti',imap,tmti,iret)
    call rwcdf(rw,ncid,'tmhacore',imap,tmhacore,iret)
    call rwcdf(rw,ncid,'tmhasol',imap,tmhasol,iret)
    call rwcdf(rw,ncid,'tmhadiv',imap,tmhadiv,iret)

    imap(1)=1
    imap(2)=1
    call rwcdf(rw,ncid,'fnisip',imap,fnisip,iret)
    call rwcdf(rw,ncid,'feesip',imap,feesip,iret)
    call rwcdf(rw,ncid,'feisip',imap,feisip,iret)
    call rwcdf(rw,ncid,'fetsip',imap,fetsip,iret)
    call rwcdf(rw,ncid,'fchsip',imap,fchsip,iret)
    call rwcdf(rw,ncid,'fnisap',imap,fnisap,iret)
    call rwcdf(rw,ncid,'feesap',imap,feesap,iret)
    call rwcdf(rw,ncid,'feisap',imap,feisap,iret)
    call rwcdf(rw,ncid,'fetsap',imap,fetsap,iret)
    call rwcdf(rw,ncid,'fchsap',imap,fchsap,iret)
    call rwcdf(rw,ncid,'fnisipp',imap,fnisipp,iret)
    call rwcdf(rw,ncid,'feesipp',imap,feesipp,iret)
    call rwcdf(rw,ncid,'feisipp',imap,feisipp,iret)
    call rwcdf(rw,ncid,'fetsipp',imap,fetsipp,iret)
    call rwcdf(rw,ncid,'fchsipp',imap,fchsipp,iret)
    call rwcdf(rw,ncid,'fnisapp',imap,fnisapp,iret)
    call rwcdf(rw,ncid,'feesapp',imap,feesapp,iret)
    call rwcdf(rw,ncid,'feisapp',imap,feisapp,iret)
    call rwcdf(rw,ncid,'fetsapp',imap,fetsapp,iret)
    call rwcdf(rw,ncid,'fchsapp',imap,fchsapp,iret)
    !
    imap(1)=nx+2     ! bl
    imap(2)=1
    call rwcdf(rw,ncid,'ne3dl',imap,ne(-1+target_offset,iylstrt),iret)
    call rwcdf(rw,ncid,'te3dl',imap,te(-1+target_offset,iylstrt),iret)
    call rwcdf(rw,ncid,'ti3dl',imap,ti(-1+target_offset,iylstrt),iret)
    call rwcdf(rw,ncid,'po3dl',imap,po(-1+target_offset,iylstrt),iret)
    call rwcdf(rw,ncid,'fn3dl',imap,fna(0,iylstrt,0,ismain),iret)
    call rwcdf(rw,ncid,'fe3dl',imap,fhe(0,iylstrt,0),iret)
    call rwcdf(rw,ncid,'fi3dl',imap,fhi(0,iylstrt,0),iret)
    call rwcdf(rw,ncid,'fc3dl',imap,fch(0,iylstrt,0),iret)
    call rwcdf(rw,ncid,'fl3dl',imap,fne(0,iylstrt,0),iret)
    call rwcdf(rw,ncid,'fo3dl',imap,fni(0,iylstrt,0),iret)
    imap(1)=nx+2     ! i
    imap(2)=1
    call rwcdf(rw,ncid,'ne3di',imap,ne(jxi,iyistrt),iret)
    call rwcdf(rw,ncid,'te3di',imap,te(jxi,iyistrt),iret)
    call rwcdf(rw,ncid,'ti3di',imap,ti(jxi,iyistrt),iret)
    call rwcdf(rw,ncid,'po3di',imap,po(jxi,iyistrt),iret)
    imap(1)=nx+2     ! a
    imap(2)=1
    call rwcdf(rw,ncid,'ne3da',imap,ne(jxa,iyastrt),iret)
    call rwcdf(rw,ncid,'te3da',imap,te(jxa,iyastrt),iret)
    call rwcdf(rw,ncid,'ti3da',imap,ti(jxa,iyastrt),iret)
    call rwcdf(rw,ncid,'po3da',imap,po(jxa,iyastrt),iret)
    imap(1)=nx+2     ! br
    imap(2)=1
    call rwcdf(rw,ncid,'ne3dr',imap,ne(nx-target_offset,iyrstrt),iret)
    call rwcdf(rw,ncid,'te3dr',imap,te(nx-target_offset,iyrstrt),iret)
    call rwcdf(rw,ncid,'ti3dr',imap,ti(nx-target_offset,iyrstrt),iret)
    call rwcdf(rw,ncid,'po3dr',imap,po(nx-target_offset,iyrstrt),iret)
    call rwcdf(rw,ncid,'fn3dr',imap,fna(nx,iyrstrt,0,ismain),iret)
    call rwcdf(rw,ncid,'fe3dr',imap,fhe(nx,iyrstrt,0),iret)
    call rwcdf(rw,ncid,'fi3dr',imap,fhi(nx,iyrstrt,0),iret)
    call rwcdf(rw,ncid,'fc3dr',imap,fch(nx,iyrstrt,0),iret)
    call rwcdf(rw,ncid,'fl3dr',imap,fne(nx,iyrstrt,0),iret)
    call rwcdf(rw,ncid,'fo3dr',imap,fni(nx,iyrstrt,0),iret)
    if (nnreg(0).ge.8) then
      imap(1)=nx+2     ! tr
      imap(2)=1
      call rwcdf(rw,ncid,'ne3dtr',imap,ne(ixtr+target_offset,iytrstrt),iret)
      call rwcdf(rw,ncid,'te3dtr',imap,te(ixtr+target_offset,iytrstrt),iret)
      call rwcdf(rw,ncid,'ti3dtr',imap,ti(ixtr+target_offset,iytrstrt),iret)
      call rwcdf(rw,ncid,'po3dtr',imap,po(ixtr+target_offset,iytrstrt),iret)
      call rwcdf(rw,ncid,'fn3dtr',imap,fna(ixtr+1,iytrstrt,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dtr',imap,fhe(ixtr+1,iytrstrt,0),iret)
      call rwcdf(rw,ncid,'fi3dtr',imap,fhi(ixtr+1,iytrstrt,0),iret)
      call rwcdf(rw,ncid,'fc3dtr',imap,fch(ixtr+1,iytrstrt,0),iret)
      call rwcdf(rw,ncid,'fl3dtr',imap,fne(ixtr+1,iytrstrt,0),iret)
      call rwcdf(rw,ncid,'fo3dtr',imap,fni(ixtr+1,iytrstrt,0),iret)
      imap(1)=nx+2     ! tl
      imap(2)=nx+2
      call rwcdf(rw,ncid,'ne3dtl',imap,ne(ixtl-target_offset,iytlstrt),iret)
      call rwcdf(rw,ncid,'te3dtl',imap,te(ixtl-target_offset,iytlstrt),iret)
      call rwcdf(rw,ncid,'ti3dtl',imap,ti(ixtl-target_offset,iytlstrt),iret)
      call rwcdf(rw,ncid,'po3dtl',imap,po(ixtl-target_offset,iytlstrt),iret)
      call rwcdf(rw,ncid,'fn3dtl',imap,fna(ixtl,iytlstrt,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dtl',imap,fhe(ixtl,iytlstrt,0),iret)
      call rwcdf(rw,ncid,'fi3dtl',imap,fhi(ixtl,iytlstrt,0),iret)
      call rwcdf(rw,ncid,'fc3dtl',imap,fch(ixtl,iytlstrt,0),iret)
      call rwcdf(rw,ncid,'fl3dtl',imap,fne(ixtl,iytlstrt,0),iret)
      call rwcdf(rw,ncid,'fo3dtl',imap,fni(ixtl,iytlstrt,0),iret)
    endif
    !
    imap(1)=1
    imap(2)=1
    slice=0.0_R8
    if (ismain0.ne.ismain) slice(iylstrt:iylend)=na(-1+target_offset,iylstrt:iylend,ismain0)
    slice(0:ny-1)=slice(0:ny-1)+dab2(1,1:ny,1,1)
    call rwcdf(rw,ncid,'an3dl',imap,slice(iylstrt),iret)
    slice=0.0_R8
    if (ismain0.ne.ismain) slice(iyistrt:iyiend)=na(jxi,iyistrt:iyiend,ismain0)
    slice(0:ny-1)=slice(0:ny-1)+dab2(jxi+1,1:ny,1,1)
    call rwcdf(rw,ncid,'an3di',imap,slice(iyistrt),iret)
    slice=0.0_R8
    if (ismain0.ne.ismain) slice(iyastrt:iyaend)=na(jxa,iyastrt:iyaend,ismain0)
    slice(0:ny-1)=slice(0:ny-1)+dab2(jxa+1,1:ny,1,1)
    call rwcdf(rw,ncid,'an3da',imap,slice(iyastrt),iret)
    slice=0.0_R8
    if (ismain0.ne.ismain) slice(iyrstrt:iyrend)=na(nx-target_offset,iyrstrt:iyrend,ismain0)
    slice(0:ny-1)=slice(0:ny-1)+dab2(nx,1:ny,1,1)
    call rwcdf(rw,ncid,'an3dr',imap,slice(iyrstrt),iret)
    if (nnreg(0).ge.8) then
      slice=0.0_R8
      if (ismain0.ne.ismain) slice(iytlstrt:iytlend)= na(ixtl-target_offset,iytlstrt:iytlend,ismain0)
      slice(0:ny-1)=slice(0:ny-1)+dab2(ixtl,1:ny,1,1)
      call rwcdf(rw,ncid,'an3dtl',imap,slice(iytlstrt),iret)
      slice=0.0_R8
      if (ismain0.ne.ismain) slice(iytrstrt:iytrend)= na(ixtr+target_offset,iytrstrt:iytrend,ismain0)
      slice(0:ny-1)=slice(0:ny-1)+dab2(ixtr+1,1:ny,1,1)
      call rwcdf(rw,ncid,'an3dtr',imap,slice(iytrstrt),iret)
    endif
    slice=0.0_R8
    slice(0:ny-1)=dmb2(1,1:ny,1,1)
    call rwcdf(rw,ncid,'mn3dl',imap,slice,iret)
    slice(0:ny-1)=dmb2(jxi+1,1:ny,1,1)
    call rwcdf(rw,ncid,'mn3di',imap,slice,iret)
    slice(0:ny-1)=dmb2(jxa+1,1:ny,1,1)
    call rwcdf(rw,ncid,'mn3da',imap,slice,iret)
    slice(0:ny-1)=dmb2(nx,1:ny,1,1)
    call rwcdf(rw,ncid,'mn3dr',imap,slice,iret)
    if (nnreg(0).ge.8) then
      slice(0:ny-1)=dmb2(ixtl,1:ny,1,1)
      call rwcdf(rw,ncid,'mn3dtl',imap,slice,iret)
      slice(0:ny-1)=dmb2(ixtr+1,1:ny,1,1)
      call rwcdf(rw,ncid,'mn3dtr',imap,slice,iret)
    endif
    !
    slice=0.0_R8
    do iy = iylstrt, iylend
      Call calc_fet(-1,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
    enddo
    call rwcdf(rw,ncid,'ft3dl',imap,slice(iylstrt),iret)
    
    slice=0.0_R8
    do iy = iyrstrt, iyrend
      Call calc_fet(nx,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
    enddo
    call rwcdf(rw,ncid,'ft3dr',imap,slice(iyrstrt),iret)
    if (nnreg(0).ge.8) then
      slice=0.0_R8
      do iy = iytlstrt, iytlend
        Call calc_fet(ixtl,iy,'R',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
      enddo
      call rwcdf(rw,ncid,'ft3dtl',imap,slice(iytlstrt),iret)
      
      slice=0.0_R8
      do iy = iytrstrt, iytrend
        Call calc_fet(ixtr,iy,'L',1._R8,nx,ny,ns,ismain,BoRiS,slice(iy))
      enddo
      call rwcdf(rw,ncid,'ft3dtr',imap,slice(iytrstrt),iret)
    endif
    slice(-1:ny)=dna0(jxi,-1:ny,ismain)
    call rwcdf(rw,ncid,'dn3di',imap,slice,iret)
    slice(-1:ny)=dna0(jxa,-1:ny,ismain)
    call rwcdf(rw,ncid,'dn3da',imap,slice,iret)
    slice(-1:ny)=dpa0(jxi,-1:ny,ismain0)* &
         (rza(jxi,-1:ny,ismain0)*te(jxi,-1:ny)+ti(jxi,-1:ny))
    call rwcdf(rw,ncid,'dp3di',imap,slice,iret)
    slice(-1:ny)=dpa0(jxa,-1:ny,ismain0)* &
         (rza(jxa,-1:ny,ismain0)*te(jxa,-1:ny)+ti(jxa,-1:ny))
    call rwcdf(rw,ncid,'dp3da',imap,slice,iret)
    slice(-1:ny)=fllim0fhi(jxi,-1:ny,1,ismain0)
    call rwcdf(rw,ncid,'lh3di',imap,slice,iret)
    slice(-1:ny)=fllim0fhi(jxa,-1:ny,1,ismain0)
    call rwcdf(rw,ncid,'lh3da',imap,slice,iret)
    slice(-1:ny)=fllim0fna(jxi,-1:ny,1,ismain0)
    call rwcdf(rw,ncid,'ln3di',imap,slice,iret)
    slice(-1:ny)=fllim0fna(jxa,-1:ny,1,ismain0)
    call rwcdf(rw,ncid,'ln3da',imap,slice,iret)
    slice(-1:ny)=hce0(jxi,-1:ny)/ne(jxi,-1:ny)
    call rwcdf(rw,ncid,'ke3di',imap,slice,iret)
    slice(-1:ny)=hce0(jxa,-1:ny)/ne(jxa,-1:ny)
    call rwcdf(rw,ncid,'ke3da',imap,slice,iret)
    slice(-1:ny)=hci0(jxi,-1:ny)/ni(jxi,-1:ny,0)
    call rwcdf(rw,ncid,'ki3di',imap,slice,iret)
    slice(-1:ny)=hci0(jxa,-1:ny)/ni(jxa,-1:ny,0)
    call rwcdf(rw,ncid,'ki3da',imap,slice,iret)
    slice(-1:ny)=vla0(jxi,-1:ny,0,ismain)
    call rwcdf(rw,ncid,'vx3di',imap,slice,iret)
    slice(-1:ny)=vla0(jxa,-1:ny,0,ismain)
    call rwcdf(rw,ncid,'vx3da',imap,slice,iret)
    slice(-1:ny)=vla0(jxi,-1:ny,1,ismain)
    call rwcdf(rw,ncid,'vy3di',imap,slice,iret)
    slice(-1:ny)=vla0(jxa,-1:ny,1,ismain)
    call rwcdf(rw,ncid,'vy3da',imap,slice,iret)
    slice(-1:ny)=vsa0(jxi,-1:ny,ismain)/(mp*am(ismain)*na(jxi,-1:ny,ismain))
    call rwcdf(rw,ncid,'vs3di',imap,slice,iret)
    slice(-1:ny)=vsa0(jxa,-1:ny,ismain)/(mp*am(ismain)*na(jxa,-1:ny,ismain))
    call rwcdf(rw,ncid,'vs3da',imap,slice,iret)
    !
    imap(1)=1
    imap(2)=1
    call rwcdf(rw,ncid,'tpsepi',imap,tpsepi,iret)
    call rwcdf(rw,ncid,'tpsepa',imap,tpsepa,iret)
    call rwcdf(rw,ncid,'tpmxip',imap,tpmxip,iret)
    call rwcdf(rw,ncid,'tpmxap',imap,tpmxap,iret)
    imap(1)=1
    imap(2)=1
    slice=0.0_R8
    if(minval(xymap(-1,0:ny-1)).gt.0) then
      slice(-1)=target_temp(xymap(-1,0),1)
      slice(0:ny-1)=target_temp(xymap(-1,0:ny-1),1)
      slice(ny)=target_temp(xymap(-1,ny-1),1)
    endif
    call rwcdf(rw,ncid,'tp3dl',imap,slice,iret)
    slice=0.0_R8
    if(minval(xymap(nx,0:ny-1)).gt.0) then
      slice(-1)=target_temp(xymap(nx,0),1)
      slice(0:ny-1)=target_temp(xymap(nx,0:ny-1),1)
      slice(ny)=target_temp(xymap(nx,ny-1),1)
    endif
    call rwcdf(rw,ncid,'tp3dr',imap,slice,iret)
    if (nnreg(0).ge.8) then
      slice=0.0_R8
      if(minval(xymap(ixtl,0:ny-1)).gt.0) then
        slice(-1)=target_temp(xymap(ixtl,0),1)
        slice(0:ny-1)=target_temp(xymap(ixtl,0:ny-1),1)
        slice(ny)=target_temp(xymap(ixtl,ny-1),1)
      endif
      call rwcdf(rw,ncid,'tp3dtl',imap,slice,iret)
      slice=0.0_R8
      if(minval(xymap(ixtr,0:ny-1)).gt.0) then
        slice(-1)=target_temp(xymap(ixtr,0),1)
        slice(0:ny-1)=target_temp(xymap(ixtr,0:ny-1),1)
        slice(ny)=target_temp(xymap(ixtr,ny-1),1)
      endif
      call rwcdf(rw,ncid,'tp3dtr',imap,slice,iret)
    endif
    !      
    iret = nf_close(ncid)
#endif

    ! ..return
    ncall = ncall+1
    call subend ()
    return

    !-----------------------------------------------------------------------
    !.end b2mwti

  end subroutine b2mwti

#ifndef NO_CDF
  subroutine b2crtimecdf(nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, iret)
    use b2mod_constants
#     include <netcdf.inc>
    integer nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, iret
    ! netcdf id
    integer  ncid
    ! dimension ids
    integer  nxdim, nydim, nsdim, timedim, &
         nybldim,nytldim,nytrdim,nybrdim,nyadim,nyidim,ncdim
    ! variable ids
    integer  ntstepid, timesaid, fnixipid, feexipid, feixipid, &
         fnixapid, feexapid, feixapid, nesepiid, tesepiid, tisepiid, &
         nesepmid, tesepmid, tisepmid, nesepaid, tesepaid, tisepaid, &
         nemxipid, temxipid, timxipid, nemxapid, temxapid, timxapid, &
         fniyipid, feeyipid, feiyipid, fniyapid, feeyapid, feiyapid, &
         pwmxipid, pwmxapid, tmneid, tmteid, tmtiid, tmhacoreid, &
         tmhasolid, tmhadivid, fnisipid, feesipid, feisipid, fnisapid, &
         feesapid, feisapid, fnisippid, feesippid, feisippid, fnisappid, &
         feesappid, feisappid, ne3dlid, te3dlid, ti3dlid, ne3diid, &
         te3diid, ti3diid, ne3daid, te3daid, ti3daid, ne3drid, te3drid, &
         ti3drid, an3dlid, mn3dlid, an3diid, mn3diid, an3daid, mn3daid, &
         an3drid, mn3drid, fn3dlid, fe3dlid, fi3dlid, fn3drid, fe3drid, &
         fi3drid, ne3dtlid, te3dtlid, ti3dtlid, ne3dtrid, te3dtrid, &
         ti3dtrid, an3dtlid, mn3dtlid, &
         an3dtrid, mn3dtrid, fn3dtlid, fe3dtlid, fi3dtlid, fn3dtrid, &
         fe3dtrid, fi3dtrid, fetxipid, fetxapid, fetyipid, fetyapid, &
         fetsipid, fetsapid, fetsippid, fetsappid
    integer  fchxipid, fchxapid, posepiid, posepmid, posepaid, &
         pomxipid, pomxapid, fchyipid, fchyapid, &
         fchsipid, fchsapid, fchsippid, fchsappid, po3dlid, &
         po3diid, po3daid, po3drid, fc3dlid, fc3drid, &
         fl3dlid, fl3drid, fo3dlid, fo3drid, &
         ft3dlid, ft3drid, po3dtlid, po3dtrid, fc3dtlid, fc3dtrid, &
         fl3dtlid, fl3dtrid, fo3dtlid, fo3dtrid, &
         ft3dtlid, ft3dtrid, &
         dn3diid, dn3daid, dp3diid, dp3daid, ke3diid, ke3daid, &
         ki3diid, ki3daid, vx3diid, vx3daid, vy3diid, vy3daid, &
         vs3diid, vs3daid, lh3diid, lh3daid, ln3diid, ln3daid, &
         dnsepmid, dpsepmid, kesepmid, kisepmid, &
         vxsepmid, vysepmid, vssepmid, &
         tpmxipid, tpmxapid, tp3drid, tp3dlid, tp3dtlid, tp3dtrid, &
         tpsepiid, tpsepaid
    ! variable shapes
    integer :: dims(2)
    ! Create and enter define mode
    iret = nf_create('b2time.nc', ncclob, ncid)
    ! define dimensions
    iret = nf_def_dim(ncid, 'nx', nx+2, nxdim)
    iret = nf_def_dim(ncid, 'ny', ny+2, nydim)
    iret = nf_def_dim(ncid, 'nybl', nybl, nybldim)
    if(nytl.gt.0) iret = nf_def_dim(ncid, 'nytl', nytl, nytldim)
    if(nytr.gt.0) iret = nf_def_dim(ncid, 'nytr', nytr, nytrdim)
    iret = nf_def_dim(ncid, 'nybr', nybr, nybrdim)
    iret = nf_def_dim(ncid, 'nyi', nyi, nyidim)
    iret = nf_def_dim(ncid, 'nya', nya, nyadim)
    iret = nf_def_dim(ncid, 'nc', nc, ncdim)
    iret = nf_def_dim(ncid, 'ns', ns, nsdim)
    iret = nf_def_dim(ncid, 'time', ncunlim, timedim)
    ! define variables
    dims(1) = 0
    iret = nf_def_var(ncid, 'ntstep', NCDOUBLE, 0, dims, ntstepid)
    dims(1) = timedim
    iret = nf_def_var(ncid, 'timesa', NCDOUBLE, 1, dims, timesaid)
    dims(1) = ncdim
    dims(2) = timedim
    iret = nf_def_var(ncid, 'fnixip', NCDOUBLE, 2, dims, fnixipid)
    iret = nf_def_var(ncid, 'feexip', NCDOUBLE, 2, dims, feexipid)
    iret = nf_def_var(ncid, 'feixip', NCDOUBLE, 2, dims, feixipid)
    iret = nf_def_var(ncid, 'fetxip', NCDOUBLE, 2, dims, fetxipid)
    iret = nf_def_var(ncid, 'fchxip', NCDOUBLE, 2, dims, fchxipid)
    iret = nf_def_var(ncid, 'fnixap', NCDOUBLE, 2, dims, fnixapid)
    iret = nf_def_var(ncid, 'feexap', NCDOUBLE, 2, dims, feexapid)
    iret = nf_def_var(ncid, 'feixap', NCDOUBLE, 2, dims, feixapid)
    iret = nf_def_var(ncid, 'fetxap', NCDOUBLE, 2, dims, fetxapid)
    iret = nf_def_var(ncid, 'fchxap', NCDOUBLE, 2, dims, fchxapid)
    iret = nf_def_var(ncid, 'nesepi', NCDOUBLE, 2, dims, nesepiid)
    iret = nf_def_var(ncid, 'tesepi', NCDOUBLE, 2, dims, tesepiid)
    iret = nf_def_var(ncid, 'tisepi', NCDOUBLE, 2, dims, tisepiid)
    iret = nf_def_var(ncid, 'tpsepi', NCDOUBLE, 2, dims, tpsepiid)
    iret = nf_def_var(ncid, 'posepi', NCDOUBLE, 2, dims, posepiid)
    iret = nf_def_var(ncid, 'nesepm', NCDOUBLE, 2, dims, nesepmid)
    iret = nf_def_var(ncid, 'tesepm', NCDOUBLE, 2, dims, tesepmid)
    iret = nf_def_var(ncid, 'tisepm', NCDOUBLE, 2, dims, tisepmid)
    iret = nf_def_var(ncid, 'posepm', NCDOUBLE, 2, dims, posepmid)
    iret = nf_def_var(ncid, 'dnsepm', NCDOUBLE, 2, dims, dnsepmid)
    iret = nf_def_var(ncid, 'dpsepm', NCDOUBLE, 2, dims, dpsepmid)
    iret = nf_def_var(ncid, 'kesepm', NCDOUBLE, 2, dims, kesepmid)
    iret = nf_def_var(ncid, 'kisepm', NCDOUBLE, 2, dims, kisepmid)
    iret = nf_def_var(ncid, 'vxsepm', NCDOUBLE, 2, dims, vxsepmid)
    iret = nf_def_var(ncid, 'vysepm', NCDOUBLE, 2, dims, vysepmid)
    iret = nf_def_var(ncid, 'vssepm', NCDOUBLE, 2, dims, vssepmid)
    iret = nf_def_var(ncid, 'nesepa', NCDOUBLE, 2, dims, nesepaid)
    iret = nf_def_var(ncid, 'tesepa', NCDOUBLE, 2, dims, tesepaid)
    iret = nf_def_var(ncid, 'tisepa', NCDOUBLE, 2, dims, tisepaid)
    iret = nf_def_var(ncid, 'tpsepa', NCDOUBLE, 2, dims, tpsepaid)
    iret = nf_def_var(ncid, 'posepa', NCDOUBLE, 2, dims, posepaid)
    iret = nf_def_var(ncid, 'nemxip', NCDOUBLE, 2, dims, nemxipid)
    iret = nf_def_var(ncid, 'temxip', NCDOUBLE, 2, dims, temxipid)
    iret = nf_def_var(ncid, 'timxip', NCDOUBLE, 2, dims, timxipid)
    iret = nf_def_var(ncid, 'tpmxip', NCDOUBLE, 2, dims, tpmxipid)
    iret = nf_def_var(ncid, 'pomxip', NCDOUBLE, 2, dims, pomxipid)
    iret = nf_def_var(ncid, 'nemxap', NCDOUBLE, 2, dims, nemxapid)
    iret = nf_def_var(ncid, 'temxap', NCDOUBLE, 2, dims, temxapid)
    iret = nf_def_var(ncid, 'timxap', NCDOUBLE, 2, dims, timxapid)
    iret = nf_def_var(ncid, 'tpmxap', NCDOUBLE, 2, dims, tpmxapid)
    iret = nf_def_var(ncid, 'pomxap', NCDOUBLE, 2, dims, pomxapid)
    iret = nf_def_var(ncid, 'fniyip', NCDOUBLE, 2, dims, fniyipid)
    iret = nf_def_var(ncid, 'feeyip', NCDOUBLE, 2, dims, feeyipid)
    iret = nf_def_var(ncid, 'feiyip', NCDOUBLE, 2, dims, feiyipid)
    iret = nf_def_var(ncid, 'fetyip', NCDOUBLE, 2, dims, fetyipid)
    iret = nf_def_var(ncid, 'fchyip', NCDOUBLE, 2, dims, fchyipid)
    iret = nf_def_var(ncid, 'fniyap', NCDOUBLE, 2, dims, fniyapid)
    iret = nf_def_var(ncid, 'feeyap', NCDOUBLE, 2, dims, feeyapid)
    iret = nf_def_var(ncid, 'feiyap', NCDOUBLE, 2, dims, feiyapid)
    iret = nf_def_var(ncid, 'fetyap', NCDOUBLE, 2, dims, fetyapid)
    iret = nf_def_var(ncid, 'fchyap', NCDOUBLE, 2, dims, fchyapid)
    iret = nf_def_var(ncid, 'pwmxip', NCDOUBLE, 2, dims, pwmxipid)
    iret = nf_def_var(ncid, 'pwmxap', NCDOUBLE, 2, dims, pwmxapid)
    dims(1) = timedim
    iret = nf_def_var(ncid, 'tmne', NCDOUBLE, 1, dims, tmneid)
    iret = nf_def_var(ncid, 'tmte', NCDOUBLE, 1, dims, tmteid)
    iret = nf_def_var(ncid, 'tmti', NCDOUBLE, 1, dims, tmtiid)
    iret = nf_def_var(ncid, 'tmhacore', NCDOUBLE, 1, dims, tmhacoreid)
    iret = nf_def_var(ncid, 'tmhasol', NCDOUBLE, 1, dims, tmhasolid)
    iret = nf_def_var(ncid, 'tmhadiv', NCDOUBLE, 1, dims, tmhadivid)
    dims(1) = ncdim
    dims(2) = timedim
    iret = nf_def_var(ncid, 'fnisip', NCDOUBLE, 2, dims, fnisipid)
    iret = nf_def_var(ncid, 'feesip', NCDOUBLE, 2, dims, feesipid)
    iret = nf_def_var(ncid, 'feisip', NCDOUBLE, 2, dims, feisipid)
    iret = nf_def_var(ncid, 'fetsip', NCDOUBLE, 2, dims, fetsipid)
    iret = nf_def_var(ncid, 'fchsip', NCDOUBLE, 2, dims, fchsipid)
    iret = nf_def_var(ncid, 'fnisap', NCDOUBLE, 2, dims, fnisapid)
    iret = nf_def_var(ncid, 'feesap', NCDOUBLE, 2, dims, feesapid)
    iret = nf_def_var(ncid, 'feisap', NCDOUBLE, 2, dims, feisapid)
    iret = nf_def_var(ncid, 'fetsap', NCDOUBLE, 2, dims, fetsapid)
    iret = nf_def_var(ncid, 'fchsap', NCDOUBLE, 2, dims, fchsapid)
    iret = nf_def_var(ncid, 'fnisipp', NCDOUBLE, 2, dims, fnisippid)
    iret = nf_def_var(ncid, 'feesipp', NCDOUBLE, 2, dims, feesippid)
    iret = nf_def_var(ncid, 'feisipp', NCDOUBLE, 2, dims, feisippid)
    iret = nf_def_var(ncid, 'fetsipp', NCDOUBLE, 2, dims, fetsippid)
    iret = nf_def_var(ncid, 'fchsipp', NCDOUBLE, 2, dims, fchsippid)
    iret = nf_def_var(ncid, 'fnisapp', NCDOUBLE, 2, dims, fnisappid)
    iret = nf_def_var(ncid, 'feesapp', NCDOUBLE, 2, dims, feesappid)
    iret = nf_def_var(ncid, 'feisapp', NCDOUBLE, 2, dims, feisappid)
    iret = nf_def_var(ncid, 'fetsapp', NCDOUBLE, 2, dims, fetsappid)
    iret = nf_def_var(ncid, 'fchsapp', NCDOUBLE, 2, dims, fchsappid)
    dims(2) = timedim
    dims(1) = nybldim
    iret = nf_def_var(ncid, 'fn3dl', NCDOUBLE, 2, dims, fn3dlid)
    iret = nf_def_var(ncid, 'fe3dl', NCDOUBLE, 2, dims, fe3dlid)
    iret = nf_def_var(ncid, 'fi3dl', NCDOUBLE, 2, dims, fi3dlid)
    iret = nf_def_var(ncid, 'ft3dl', NCDOUBLE, 2, dims, ft3dlid)
    iret = nf_def_var(ncid, 'fc3dl', NCDOUBLE, 2, dims, fc3dlid)
    iret = nf_def_var(ncid, 'fl3dl', NCDOUBLE, 2, dims, fl3dlid)
    iret = nf_def_var(ncid, 'fo3dl', NCDOUBLE, 2, dims, fo3dlid)
    iret = nf_def_var(ncid, 'ne3dl', NCDOUBLE, 2, dims, ne3dlid)
    iret = nf_def_var(ncid, 'te3dl', NCDOUBLE, 2, dims, te3dlid)
    iret = nf_def_var(ncid, 'ti3dl', NCDOUBLE, 2, dims, ti3dlid)
    iret = nf_def_var(ncid, 'po3dl', NCDOUBLE, 2, dims, po3dlid)
    iret = nf_def_var(ncid, 'an3dl', NCDOUBLE, 2, dims, an3dlid)
    iret = nf_def_var(ncid, 'mn3dl', NCDOUBLE, 2, dims, mn3dlid)
    iret = nf_def_var(ncid, 'tp3dl', NCDOUBLE, 2, dims, tp3dlid)
    if(nytl.gt.0) then
      dims(2) = timedim
      dims(1) = nytldim
      iret = nf_def_var(ncid, 'fn3dtl', NCDOUBLE, 2, dims, fn3dtlid)
      iret = nf_def_var(ncid, 'fe3dtl', NCDOUBLE, 2, dims, fe3dtlid)
      iret = nf_def_var(ncid, 'fi3dtl', NCDOUBLE, 2, dims, fi3dtlid)
      iret = nf_def_var(ncid, 'ft3dtl', NCDOUBLE, 2, dims, ft3dtlid)
      iret = nf_def_var(ncid, 'fc3dtl', NCDOUBLE, 2, dims, fc3dtlid)
      iret = nf_def_var(ncid, 'fl3dtl', NCDOUBLE, 2, dims, fl3dtlid)
      iret = nf_def_var(ncid, 'fo3dtl', NCDOUBLE, 2, dims, fo3dtlid)
      iret = nf_def_var(ncid, 'ne3dtl', NCDOUBLE, 2, dims, ne3dtlid)
      iret = nf_def_var(ncid, 'te3dtl', NCDOUBLE, 2, dims, te3dtlid)
      iret = nf_def_var(ncid, 'ti3dtl', NCDOUBLE, 2, dims, ti3dtlid)
      iret = nf_def_var(ncid, 'po3dtl', NCDOUBLE, 2, dims, po3dtlid)
      iret = nf_def_var(ncid, 'an3dtl', NCDOUBLE, 2, dims, an3dtlid)
      iret = nf_def_var(ncid, 'mn3dtl', NCDOUBLE, 2, dims, mn3dtlid)
      iret = nf_def_var(ncid, 'tp3dtl', NCDOUBLE, 2, dims, tp3dtlid)
    endif
    dims(2) = timedim
    dims(1) = nyidim
    iret = nf_def_var(ncid, 'ne3di', NCDOUBLE, 2, dims, ne3diid)
    iret = nf_def_var(ncid, 'te3di', NCDOUBLE, 2, dims, te3diid)
    iret = nf_def_var(ncid, 'ti3di', NCDOUBLE, 2, dims, ti3diid)
    iret = nf_def_var(ncid, 'po3di', NCDOUBLE, 2, dims, po3diid)
    iret = nf_def_var(ncid, 'an3di', NCDOUBLE, 2, dims, an3diid)
    iret = nf_def_var(ncid, 'mn3di', NCDOUBLE, 2, dims, mn3diid)
    iret = nf_def_var(ncid, 'dn3di', NCDOUBLE, 2, dims, dn3diid)
    iret = nf_def_var(ncid, 'dp3di', NCDOUBLE, 2, dims, dp3diid)
    iret = nf_def_var(ncid, 'lh3di', NCDOUBLE, 2, dims, lh3diid)
    iret = nf_def_var(ncid, 'ln3di', NCDOUBLE, 2, dims, ln3diid)
    iret = nf_def_var(ncid, 'ke3di', NCDOUBLE, 2, dims, ke3diid)
    iret = nf_def_var(ncid, 'ki3di', NCDOUBLE, 2, dims, ki3diid)
    iret = nf_def_var(ncid, 'vx3di', NCDOUBLE, 2, dims, vx3diid)
    iret = nf_def_var(ncid, 'vy3di', NCDOUBLE, 2, dims, vy3diid)
    iret = nf_def_var(ncid, 'vs3di', NCDOUBLE, 2, dims, vs3diid)
    dims(2) = timedim
    dims(1) = nyadim
    iret = nf_def_var(ncid, 'ne3da', NCDOUBLE, 2, dims, ne3daid)
    iret = nf_def_var(ncid, 'te3da', NCDOUBLE, 2, dims, te3daid)
    iret = nf_def_var(ncid, 'ti3da', NCDOUBLE, 2, dims, ti3daid)
    iret = nf_def_var(ncid, 'po3da', NCDOUBLE, 2, dims, po3daid)
    iret = nf_def_var(ncid, 'an3da', NCDOUBLE, 2, dims, an3daid)
    iret = nf_def_var(ncid, 'mn3da', NCDOUBLE, 2, dims, mn3daid)
    iret = nf_def_var(ncid, 'dn3da', NCDOUBLE, 2, dims, dn3daid)
    iret = nf_def_var(ncid, 'dp3da', NCDOUBLE, 2, dims, dp3daid)
    iret = nf_def_var(ncid, 'lh3da', NCDOUBLE, 2, dims, lh3daid)
    iret = nf_def_var(ncid, 'ln3da', NCDOUBLE, 2, dims, ln3daid)
    iret = nf_def_var(ncid, 'ke3da', NCDOUBLE, 2, dims, ke3daid)
    iret = nf_def_var(ncid, 'ki3da', NCDOUBLE, 2, dims, ki3daid)
    iret = nf_def_var(ncid, 'vx3da', NCDOUBLE, 2, dims, vx3daid)
    iret = nf_def_var(ncid, 'vy3da', NCDOUBLE, 2, dims, vy3daid)
    iret = nf_def_var(ncid, 'vs3da', NCDOUBLE, 2, dims, vs3daid)
    dims(2) = timedim
    dims(1) = nybrdim
   iret = nf_def_var(ncid, 'ne3dr', NCDOUBLE, 2, dims, ne3drid)
    iret = nf_def_var(ncid, 'te3dr', NCDOUBLE, 2, dims, te3drid)
    iret = nf_def_var(ncid, 'ti3dr', NCDOUBLE, 2, dims, ti3drid)
    iret = nf_def_var(ncid, 'po3dr', NCDOUBLE, 2, dims, po3drid)
    iret = nf_def_var(ncid, 'an3dr', NCDOUBLE, 2, dims, an3drid)
    iret = nf_def_var(ncid, 'mn3dr', NCDOUBLE, 2, dims, mn3drid)
    iret = nf_def_var(ncid, 'fn3dr', NCDOUBLE, 2, dims, fn3drid)
    iret = nf_def_var(ncid, 'fe3dr', NCDOUBLE, 2, dims, fe3drid)
    iret = nf_def_var(ncid, 'fi3dr', NCDOUBLE, 2, dims, fi3drid)
    iret = nf_def_var(ncid, 'ft3dr', NCDOUBLE, 2, dims, ft3drid)
    iret = nf_def_var(ncid, 'fc3dr', NCDOUBLE, 2, dims, fc3drid)
    iret = nf_def_var(ncid, 'fl3dr', NCDOUBLE, 2, dims, fl3drid)
    iret = nf_def_var(ncid, 'fo3dr', NCDOUBLE, 2, dims, fo3drid)
    iret = nf_def_var(ncid, 'tp3dr', NCDOUBLE, 2, dims, tp3drid)
    if(nytr.gt.0) then
      dims(2) = timedim
      dims(1) = nytrdim
      iret = nf_def_var(ncid, 'ne3dtr', NCDOUBLE, 2, dims, ne3dtrid)
      iret = nf_def_var(ncid, 'te3dtr', NCDOUBLE, 2, dims, te3dtrid)
      iret = nf_def_var(ncid, 'ti3dtr', NCDOUBLE, 2, dims, ti3dtrid)
      iret = nf_def_var(ncid, 'po3dtr', NCDOUBLE, 2, dims, po3dtrid)
      iret = nf_def_var(ncid, 'an3dtr', NCDOUBLE, 2, dims, an3dtrid)
      iret = nf_def_var(ncid, 'mn3dtr', NCDOUBLE, 2, dims, mn3dtrid)
      iret = nf_def_var(ncid, 'fn3dtr', NCDOUBLE, 2, dims, fn3dtrid)
      iret = nf_def_var(ncid, 'fe3dtr', NCDOUBLE, 2, dims, fe3dtrid)
      iret = nf_def_var(ncid, 'fi3dtr', NCDOUBLE, 2, dims, fi3dtrid)
      iret = nf_def_var(ncid, 'ft3dtr', NCDOUBLE, 2, dims, ft3dtrid)
      iret = nf_def_var(ncid, 'fc3dtr', NCDOUBLE, 2, dims, fc3dtrid)
      iret = nf_def_var(ncid, 'fl3dtr', NCDOUBLE, 2, dims, fl3dtrid)
      iret = nf_def_var(ncid, 'fo3dtr', NCDOUBLE, 2, dims, fo3dtrid)
      iret = nf_def_var(ncid, 'tp3dtr', NCDOUBLE, 2, dims, tp3dtrid)
    endif
    ! assign attributes
    iret = nf_put_att_text(ncid, timesaid, 'long_name', 4, 'time')
    iret = nf_put_att_text(ncid, timesaid, 'units', 2, 's ')
    iret = nf_put_att_text(ncid, fnixipid, 'long_name', 51, 'integrated poloidal particle flux, inboard divertor')
    iret = nf_put_att_text(ncid, fnixipid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feexipid, 'long_name', 58, 'integrated poloidal electron energy flux, inboard divertor')
    iret = nf_put_att_text(ncid, feexipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feixipid, 'long_name', 53, 'integrated poloidal ion energy flux, inboard divertor')
    iret = nf_put_att_text(ncid, feixipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetxipid, 'long_name', 53, 'integrated poloidal tot energy flux, inboard divertor')
    iret = nf_put_att_text(ncid, fetxipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchxipid, 'long_name', 45, 'integrated poloidal current, inboard divertor')
    iret = nf_put_att_text(ncid, fchxipid, 'units', 2, 'A ')
    iret = nf_put_att_text(ncid, fnixapid, 'long_name', 52, 'integrated poloidal particle flux, outboard divertor')
    iret = nf_put_att_text(ncid, fnixapid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feexapid, 'long_name', 59, 'integrated poloidal electron energy flux, outboard divertor')
    iret = nf_put_att_text(ncid, feexapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feixapid, 'long_name', 54, 'integrated poloidal ion energy flux, outboard divertor')
    iret = nf_put_att_text(ncid, feixapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetxapid, 'long_name', 54, 'integrated poloidal tot energy flux, outboard divertor')
    iret = nf_put_att_text(ncid, fetxapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchxapid, 'long_name', 46, 'integrated poloidal current, outboard divertor')
    iret = nf_put_att_text(ncid, fchxapid, 'units', 2, 'A ')
    iret = nf_put_att_text(ncid, nesepiid, 'long_name', 45, 'separatrix electron density, inboard divertor')
    iret = nf_put_att_text(ncid, nesepiid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, tesepiid, 'long_name', 49, 'separatrix electron temperature, inboard divertor')
    iret = nf_put_att_text(ncid, tesepiid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tisepiid, 'long_name', 44, 'separatrix ion temperature, inboard divertor')
    iret = nf_put_att_text(ncid, tisepiid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tpsepiid, 'long_name', 46, 'separatrix plate temperature, inboard divertor')
    iret = nf_put_att_text(ncid, tpsepiid, 'units', 2, 'K ')
    iret = nf_put_att_text(ncid, posepiid, 'long_name', 38, 'separatrix potential, inboard divertor')
    iret = nf_put_att_text(ncid, posepiid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, nesepmid, 'long_name', 43, 'separatrix electron density, outer midplane')
    iret = nf_put_att_text(ncid, nesepmid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, tesepmid, 'long_name', 47, 'separatrix electron temperature, outer midplane')
    iret = nf_put_att_text(ncid, tesepmid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tisepmid, 'long_name', 42, 'separatrix ion temperature, outer midplane')
    iret = nf_put_att_text(ncid, tisepmid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, posepmid, 'long_name', 36, 'separatrix potential, outer midplane')
    iret = nf_put_att_text(ncid, posepmid, 'units', 2, 'V ')

    iret = nf_put_att_text(ncid, dnsepmid, 'long_name', 37, 'diffusion coefficient, outer midplane')
    iret = nf_put_att_text(ncid, dnsepmid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, dpsepmid, 'long_name', 46, 'pressure diffusion coefficient, outer midplane')
    iret = nf_put_att_text(ncid, dpsepmid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, kesepmid, 'long_name', 44, 'electron thermal diffusivity, outer midplane')
    iret = nf_put_att_text(ncid, kesepmid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, kisepmid, 'long_name', 39, 'ion thermal diffusivity, outer midplane')
    iret = nf_put_att_text(ncid, kisepmid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, vxsepmid, 'long_name', 39, 'poloidal pinch velocity, outer midplane')
    iret = nf_put_att_text(ncid, vxsepmid, 'units', 6, 'm.s^-1')
    iret = nf_put_att_text(ncid, vysepmid, 'long_name', 37, 'radial pinch velocity, outer midplane')
    iret = nf_put_att_text(ncid, vysepmid, 'units', 6, 'm.s^-1')
    iret = nf_put_att_text(ncid, vssepmid, 'long_name', 37, 'viscosity coefficient, outer midplane')
    iret = nf_put_att_text(ncid, vssepmid, 'units', 12, 'm.kg^-1.s^-1')

    iret = nf_put_att_text(ncid, nesepaid, 'long_name', 46, 'separatrix electron density, outboard divertor')
    iret = nf_put_att_text(ncid, nesepaid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, tesepaid, 'long_name', 50, 'separatrix electron temperature, outboard divertor')
    iret = nf_put_att_text(ncid, tesepaid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tisepaid, 'long_name', 45, 'separatrix ion temperature, outboard divertor')
    iret = nf_put_att_text(ncid, tisepaid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tpsepaid, 'long_name', 47, 'separatrix plate temperature, outboard divertor')
    iret = nf_put_att_text(ncid, tpsepaid, 'units', 2, 'K ')
    iret = nf_put_att_text(ncid, posepaid, 'long_name', 43, 'separatrix potential, outboard divertor')
    iret = nf_put_att_text(ncid, posepaid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, nemxipid, 'long_name', 42, 'maximum electron density, inboard divertor')
    iret = nf_put_att_text(ncid, nemxipid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, temxipid, 'long_name', 46, 'maximum electron temperature, inboard divertor')
    iret = nf_put_att_text(ncid, temxipid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, timxipid, 'long_name', 41, 'maximum ion temperature, inboard divertor')
    iret = nf_put_att_text(ncid, timxipid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tpmxipid, 'long_name', 43, 'maximum plate temperature, inboard divertor')
    iret = nf_put_att_text(ncid, tpmxipid, 'units', 2, 'K ')
    iret = nf_put_att_text(ncid, pomxipid, 'long_name', 35, 'maximum potential, inboard divertor')
    iret = nf_put_att_text(ncid, pomxipid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, nemxapid, 'long_name', 43, 'maximum electron density, outboard divertor')
    iret = nf_put_att_text(ncid, nemxapid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, temxapid, 'long_name', 47, 'maximum electron temperature, outboard divertor')
    iret = nf_put_att_text(ncid, temxapid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, timxapid, 'long_name', 43, 'maximum ion temperature, outboard divertor')
    iret = nf_put_att_text(ncid, timxapid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tpmxapid, 'long_name', 45, 'maximum plate temperature, outboard divertor')
    iret = nf_put_att_text(ncid, tpmxapid, 'units', 2, 'K ')
    iret = nf_put_att_text(ncid, pomxapid, 'long_name', 37, 'maximum potential, outboard divertor')
    iret = nf_put_att_text(ncid, pomxapid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, fniyipid, 'long_name', 45, 'integrated radial particle flux, main chamber')
    iret = nf_put_att_text(ncid, fniyipid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feeyipid, 'long_name', 52, 'integrated radial electron energy flux, main chamber')
    iret = nf_put_att_text(ncid, feeyipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feiyipid, 'long_name', 47, 'integrated radial ion energy flux, main chamber')
    iret = nf_put_att_text(ncid, feiyipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetyipid, 'long_name', 47, 'integrated radial tot energy flux, main chamber')
    iret = nf_put_att_text(ncid, fetyipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchyipid, 'long_name', 39, 'integrated radial current, main chamber')
    iret = nf_put_att_text(ncid, fchyipid, 'units', 2, 'A ')
    iret = nf_put_att_text(ncid, fniyapid, 'long_name', 48, 'integrated radial particle flux, divertor region')
    iret = nf_put_att_text(ncid, fniyapid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feeyapid, 'long_name', 55, 'integrated radial electron energy flux, divertor region')
    iret = nf_put_att_text(ncid, feeyapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feiyapid, 'long_name', 50, 'integrated radial ion energy flux, divertor region')
    iret = nf_put_att_text(ncid, feiyapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetyapid, 'long_name', 50, 'integrated radial tot energy flux, divertor region')
    iret = nf_put_att_text(ncid, fetyapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchyapid, 'long_name', 42, 'integrated radial current, divertor region')
    iret = nf_put_att_text(ncid, fchyapid, 'units', 2, 'A ')
    iret = nf_put_att_text(ncid, pwmxipid, 'long_name', 42, 'maximum total power flux, inboard divertor')
    iret = nf_put_att_text(ncid, pwmxipid, 'units', 6, 'W.m^-2')
    iret = nf_put_att_text(ncid, pwmxapid, 'long_name', 43, 'maximum total power flux, outboard divertor')
    iret = nf_put_att_text(ncid, pwmxapid, 'units', 6, 'W.m^-2')
    iret = nf_put_att_text(ncid, tmneid, 'long_name', 25, 'total number of particles')
    iret = nf_put_att_text(ncid, tmneid, 'units', 2, '  ')
    iret = nf_put_att_text(ncid, tmteid, 'long_name', 21, 'total electron energy')
    iret = nf_put_att_text(ncid, tmteid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tmtiid, 'long_name', 16, 'total ion energy')
    iret = nf_put_att_text(ncid, tmtiid, 'units', 2, 'eV')
    iret = nf_put_att_text(ncid, tmhacoreid, 'long_name', 24, 'H-alpha emissivity, core')
    iret = nf_put_att_text(ncid, tmhacoreid, 'units', 19, 'photons.m^-2.sr^-1?')
    iret = nf_put_att_text(ncid, tmhasolid, 'long_name', 41, 'H-alpha emissivity, SOL above the x-point')
    iret = nf_put_att_text(ncid, tmhasolid, 'units', 19, 'photons.m^-2.sr^-1?')
    iret = nf_put_att_text(ncid, tmhadivid, 'long_name', 35, 'H-alpha emissivity, divertor region')
    iret = nf_put_att_text(ncid, tmhadivid, 'units', 19, 'photons.m^-2.sr^-1?')
    iret = nf_put_att_text(ncid, fnisipid, 'long_name', 54, 'poloidal particle flux, into inboard separatrix throat')
    iret = nf_put_att_text(ncid, fnisipid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feesipid, 'long_name', 61, 'poloidal electron energy flux, into inboard separatrix throat')
    iret = nf_put_att_text(ncid, feesipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feisipid, 'long_name', 56, 'poloidal ion energy flux, into inboard separatrix throat')
    iret = nf_put_att_text(ncid, feisipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetsipid, 'long_name', 58, 'poloidal total energy flux, into inboard separatrix throat')
    iret = nf_put_att_text(ncid, fetsipid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchsipid, 'long_name', 48, 'poloidal current, into inboard separatrix throat')
    iret = nf_put_att_text(ncid, fchsipid, 'units', 2, 'A ')
    iret = nf_put_att_text(ncid, fnisapid, 'long_name', 55, 'poloidal particle flux, into outboard separatrix throat')
    iret = nf_put_att_text(ncid, fnisapid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feesapid, 'long_name', 62, 'poloidal electron energy flux, into outboard separaatrix throat')
    iret = nf_put_att_text(ncid, feesapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feisapid, 'long_name', 57, 'poloidal ion energy flux, into outboard separatrix throat')
    iret = nf_put_att_text(ncid, feisapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetsapid, 'long_name', 59, 'poloidal total energy flux, into outboard separatrix throat')
    iret = nf_put_att_text(ncid, fetsapid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchsapid, 'long_name', 49, 'poloidal current, into outboard separatrix throat')
    iret = nf_put_att_text(ncid, fchsapid, 'units', 2, 'A ')
    iret = nf_put_att_text(ncid, fnisippid, 'long_name', 40, 'poloidal particle flux, core x-pt region')
    iret = nf_put_att_text(ncid, fnisippid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feesippid, 'long_name', 52, 'poloidal electron energy flux, core x-pt flux region')
    iret = nf_put_att_text(ncid, feesippid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feisippid, 'long_name', 47, 'poloidal ion energy flux, core x-pt flux region')
    iret = nf_put_att_text(ncid, feisippid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetsippid, 'long_name', 49, 'poloidal total energy flux, core x-pt flux region')
    iret = nf_put_att_text(ncid, fetsippid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchsippid, 'long_name', 39, 'poloidal current, core x-pt flux region')
    iret = nf_put_att_text(ncid, fchsippid, 'units', 2, 'A ')
    iret = nf_put_att_text(ncid, fnisappid, 'long_name', 48, 'poloidal particle flux, x-pt private flux region')
    iret = nf_put_att_text(ncid, fnisappid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, feesappid, 'long_name', 55, 'poloidal electron energy flux, x-pt private flux region')
    iret = nf_put_att_text(ncid, feesappid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, feisappid, 'long_name', 50, 'poloidal ion energy flux, x-pt private flux region')
    iret = nf_put_att_text(ncid, feisappid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fetsappid, 'long_name', 52, 'poloidal total energy flux, x-pt private flux region')
    iret = nf_put_att_text(ncid, fetsappid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fchsappid, 'long_name', 42, 'poloidal current, x-pt private flux region')
    iret = nf_put_att_text(ncid, fchsappid, 'units', 2, 'A ')

    ! inboard divertor quantities
    iret = nf_put_att_text(ncid, ne3dlid, 'long_name', 34, 'electron density, inboard divertor')
    iret = nf_put_att_text(ncid, ne3dlid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, te3dlid, 'long_name', 38, 'electron temperature, inboard divertor')
    iret = nf_put_att_text(ncid, te3dlid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, te3dlid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, ti3dlid, 'long_name', 33, 'ion temperature, inboard divertor')
    iret = nf_put_att_text(ncid, ti3dlid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, ti3dlid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, tp3dlid, 'long_name', 35, 'plate temperature, inboard divertor')
    iret = nf_put_att_text(ncid, tp3dlid, 'units', 2, 'K ')
    iret = nf_put_att_text(ncid, po3dlid, 'long_name', 26, 'potential, inboard divertor')
    iret = nf_put_att_text(ncid, po3dlid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, an3dlid, 'long_name', 30, 'atom density, inboard divertor')
    iret = nf_put_att_text(ncid, an3dlid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, mn3dlid, 'long_name', 34, 'molecule density, inboard divertor')
    iret = nf_put_att_text(ncid, mn3dlid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, fn3dlid, 'long_name', 44, 'poloidal main species flux, inboard divertor')
    iret = nf_put_att_text(ncid, fn3dlid, 'units', 4, 's^-1')
    iret = nf_put_att_double(ncid, fn3dlid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    iret = nf_put_att_text(ncid, fl3dlid, 'long_name', 40, 'poloidal electron flux, inboard divertor')
    iret = nf_put_att_text(ncid, fl3dlid, 'units', 4, 's^-1')
    iret = nf_put_att_double(ncid, fl3dlid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    iret = nf_put_att_text(ncid, fo3dlid, 'long_name', 35, 'poloidal ion flux, inboard divertor')
    iret = nf_put_att_text(ncid, fo3dlid, 'units', 4, 's^-1')
    iret = nf_put_att_double(ncid, fo3dlid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    iret = nf_put_att_text(ncid, fe3dlid, 'long_name', 47, 'poloidal electron energy flux, inboard divertor')
    iret = nf_put_att_text(ncid, fe3dlid, 'units', 2, 'W ')
    iret = nf_put_att_double(ncid, fe3dlid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    iret = nf_put_att_text(ncid, fi3dlid, 'long_name', 42, 'poloidal ion energy flux, inboard divertor')
    iret = nf_put_att_text(ncid, fi3dlid, 'units', 2, 'W ')
    iret = nf_put_att_double(ncid, fi3dlid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    iret = nf_put_att_text(ncid, ft3dlid, 'long_name', 44, 'poloidal total energy flux, inboard divertor')
    iret = nf_put_att_text(ncid, ft3dlid, 'units', 2, 'W ')
    iret = nf_put_att_double(ncid, ft3dlid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    iret = nf_put_att_text(ncid, fc3dlid, 'long_name', 35, 'poloidal current, inboard divertor')
    iret = nf_put_att_text(ncid, fc3dlid, 'units', 2, 'A ')
    iret = nf_put_att_double(ncid, fc3dlid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    ! inboard midplane quantities
    iret = nf_put_att_text(ncid, ne3diid, 'long_name', 34, 'electron density, inboard midplane')
    iret = nf_put_att_text(ncid, ne3diid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, te3diid, 'long_name', 38, 'electron temperature, inboard midplane')
    iret = nf_put_att_text(ncid, te3diid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, te3diid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, ti3diid, 'long_name', 33, 'ion temperature, inboard midplane')
    iret = nf_put_att_text(ncid, ti3diid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, ti3diid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, po3diid, 'long_name', 27, 'potential, inboard midplane')
    iret = nf_put_att_text(ncid, po3diid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, an3diid, 'long_name', 30, 'atom density, inboard midplane')
    iret = nf_put_att_text(ncid, an3diid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, mn3diid, 'long_name', 34, 'molecule density, inboard midplane')
    iret = nf_put_att_text(ncid, mn3diid, 'units', 4, 'm^-3')

    iret = nf_put_att_text(ncid, dn3diid, 'long_name', 39, 'diffusion coefficient, inboard midplane')
    iret = nf_put_att_text(ncid, dn3diid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, dp3diid, 'long_name', 48, 'pressure diffusion coefficient, inboard midplane')
    iret = nf_put_att_text(ncid, dp3diid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, lh3diid, 'long_name', 50, 'radial neutral heat flux limiter, inboard midplane')
    iret = nf_put_att_text(ncid, lh3diid, 'units', 1, ' ')
    iret = nf_put_att_text(ncid, ln3diid, 'long_name', 50, 'radial neutral part flux limiter, inboard midplane')
    iret = nf_put_att_text(ncid, ln3diid, 'units', 1, ' ')
    iret = nf_put_att_text(ncid, ke3diid, 'long_name', 46, 'electron thermal diffusivity, inboard midplane')
    iret = nf_put_att_text(ncid, ke3diid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, ki3diid, 'long_name', 41, 'ion thermal diffusivity, inboard midplane')
    iret = nf_put_att_text(ncid, ki3diid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, vx3diid, 'long_name', 41, 'poloidal pinch velocity, inboard midplane')
    iret = nf_put_att_text(ncid, vx3diid, 'units', 6, 'm.s^-1')
    iret = nf_put_att_text(ncid, vy3diid, 'long_name', 39, 'radial pinch velocity, inboard midplane')
    iret = nf_put_att_text(ncid, vy3diid, 'units', 6, 'm.s^-1')
    iret = nf_put_att_text(ncid, vs3diid, 'long_name', 39, 'viscosity coefficient, inboard midplane')
    iret = nf_put_att_text(ncid, vs3diid, 'units', 12, 'm.kg^-1.s^-1')

    ! upper inboard divertor quantities
    if(nytl.gt.0) then
      iret = nf_put_att_text(ncid, ne3dtlid, 'long_name', 40, 'electron density, upper inboard divertor')
      iret = nf_put_att_text(ncid, ne3dtlid, 'units', 4, 'm^-3')
      iret = nf_put_att_text(ncid, te3dtlid, 'long_name', 44, 'electron temperature, upper inboard divertor')
      iret = nf_put_att_text(ncid, te3dtlid, 'units', 2, 'eV')
      iret = nf_put_att_double(ncid, te3dtlid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
      iret = nf_put_att_text(ncid, ti3dtlid, 'long_name', 39, 'ion temperature, upper inboard divertor')
      iret = nf_put_att_text(ncid, ti3dtlid, 'units', 2, 'eV')
      iret = nf_put_att_double(ncid, ti3dtlid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
      iret = nf_put_att_text(ncid, tp3dtlid, 'long_name', 41, 'plate temperature, upper inboard divertor')
      iret = nf_put_att_text(ncid, tp3dtlid, 'units', 2, 'K ')
      iret = nf_put_att_text(ncid, po3dtlid, 'long_name', 32, 'potential, upper inboard divertor')
      iret = nf_put_att_text(ncid, po3dtlid, 'units', 2, 'V ')
      iret = nf_put_att_text(ncid, an3dtlid, 'long_name', 36, 'atom density, upper inboard divertor')
      iret = nf_put_att_text(ncid, an3dtlid, 'units', 4, 'm^-3')
      iret = nf_put_att_text(ncid, mn3dtlid, 'long_name', 40, 'molecule density, upper inboard divertor')
      iret = nf_put_att_text(ncid, mn3dtlid, 'units', 4, 'm^-3')
      iret = nf_put_att_text(ncid, fn3dtlid, 'long_name', 50, 'poloidal main species flux, upper inboard divertor')
      iret = nf_put_att_text(ncid, fn3dtlid, 'units', 4, 's^-1')
      iret = nf_put_att_text(ncid, fl3dtlid, 'long_name', 46, 'poloidal electron flux, inboard divertor')
      iret = nf_put_att_text(ncid, fl3dtlid, 'units', 4, 's^-1')
      iret = nf_put_att_text(ncid, fo3dtlid, 'long_name', 41, 'poloidal ion flux, upper inboard divertor')
      iret = nf_put_att_text(ncid, fo3dtlid, 'units', 4, 's^-1')
      iret = nf_put_att_text(ncid, fe3dtlid, 'long_name', 53, 'poloidal electron energy flux, upper inboard divertor')
      iret = nf_put_att_text(ncid, fe3dtlid, 'units', 2, 'W ')
      iret = nf_put_att_text(ncid, fi3dtlid, 'long_name', 48, 'poloidal ion energy flux, upper inboard divertor')
      iret = nf_put_att_text(ncid, fi3dtlid, 'units', 2, 'W ')
      iret = nf_put_att_text(ncid, ft3dtlid, 'long_name', 50, 'poloidal total energy flux, upper inboard divertor')
      iret = nf_put_att_text(ncid, ft3dtlid, 'units', 2, 'W ')
      iret = nf_put_att_text(ncid, fc3dtlid, 'long_name', 41, 'poloidal current, upper inboard divertor')
      iret = nf_put_att_text(ncid, fc3dtlid, 'units', 2, 'A ')
    endif
    ! outboard midplane quantities
    iret = nf_put_att_text(ncid, ne3daid, 'long_name', 35, 'electron density, outboard midplane')
    iret = nf_put_att_text(ncid, ne3daid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, te3daid, 'long_name', 39, 'electron temperature, outboard midplane')
    iret = nf_put_att_text(ncid, te3daid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, te3daid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, ti3daid, 'long_name', 34, 'ion temperature, outboard midplane')
    iret = nf_put_att_text(ncid, ti3daid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, ti3daid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, po3daid, 'long_name', 28, 'potential, outboard midplane')
    iret = nf_put_att_text(ncid, po3daid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, an3daid, 'long_name', 31, 'atom density, outboard midplane')
    iret = nf_put_att_text(ncid, an3daid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, mn3daid, 'long_name', 35, 'molecule density, outboard midplane')
    iret = nf_put_att_text(ncid, mn3daid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, dn3daid, 'long_name', 40, 'diffusion coefficient, outboard midplane')
    iret = nf_put_att_text(ncid, dn3daid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, dp3daid, 'long_name', 49, 'pressure diffusion coefficient, outboard midplane')
    iret = nf_put_att_text(ncid, dp3daid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, lh3daid, 'long_name', 51, 'radial neutral heat flux limiter, outboard midplane')
    iret = nf_put_att_text(ncid, lh3daid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, ln3daid, 'long_name', 51, 'radial neutral part flux limiter, outboard midplane')
    iret = nf_put_att_text(ncid, ln3daid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, ke3daid, 'long_name', 47, 'electron thermal diffusivity, outboard midplane')
    iret = nf_put_att_text(ncid, ke3daid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, ki3daid, 'long_name', 42, 'ion thermal diffusivity, outboard midplane')
    iret = nf_put_att_text(ncid, ki3daid, 'units', 8, 'm^2.s^-1')
    iret = nf_put_att_text(ncid, vx3daid, 'long_name', 42, 'poloidal pinch velocity, outboard midplane')
    iret = nf_put_att_text(ncid, vx3daid, 'units', 6, 'm.s^-1')
    iret = nf_put_att_text(ncid, vy3daid, 'long_name', 40, 'radial pinch velocity, outboard midplane')
    iret = nf_put_att_text(ncid, vy3daid, 'units', 6, 'm.s^-1')
    iret = nf_put_att_text(ncid, vs3daid, 'long_name', 40, 'viscosity coefficient, outboard midplane')
    iret = nf_put_att_text(ncid, vs3daid, 'units', 12, 'm.kg^-1.s^-1')
    ! outboard divertor quantities
    iret = nf_put_att_text(ncid, ne3drid, 'long_name', 35, 'electron density, outboard divertor')
    iret = nf_put_att_text(ncid, ne3drid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, te3drid, 'long_name', 39, 'electron temperature, outboard divertor')
    iret = nf_put_att_text(ncid, te3drid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, te3drid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, ti3drid, 'long_name', 34, 'ion temperature, outboard divertor')
    iret = nf_put_att_text(ncid, ti3drid, 'units', 2, 'eV')
    iret = nf_put_att_double(ncid, ti3drid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
    iret = nf_put_att_text(ncid, tp3drid, 'long_name', 36, 'plate temperature, outboard divertor')
    iret = nf_put_att_text(ncid, tp3drid, 'units', 2, 'K ')
    iret = nf_put_att_text(ncid, po3drid, 'long_name', 28, 'potential, outboard divertor')
    iret = nf_put_att_text(ncid, po3drid, 'units', 2, 'V ')
    iret = nf_put_att_text(ncid, an3drid, 'long_name', 31, 'atom density, outboard divertor')
    iret = nf_put_att_text(ncid, an3drid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, mn3drid, 'long_name', 35, 'molecule density, outboard divertor')
    iret = nf_put_att_text(ncid, mn3drid, 'units', 4, 'm^-3')
    iret = nf_put_att_text(ncid, fn3drid, 'long_name', 45, 'poloidal main species flux, outboard divertor')
    iret = nf_put_att_text(ncid, fn3drid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, fl3drid, 'long_name', 41, 'poloidal electron flux, outboard divertor')
    iret = nf_put_att_text(ncid, fl3drid, 'units', 4, 's^-1')
    iret = nf_put_att_text(ncid, fo3drid, 'long_name', 36, 'poloidal ion flux, outboard divertor')
    iret = nf_put_att_text(ncid, fo3drid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fe3drid, 'long_name', 48, 'poloidal electron energy flux, outboard divertor')
    iret = nf_put_att_text(ncid, fe3drid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fi3drid, 'long_name', 43, 'poloidal ion energy flux, outboard divertor')
    iret = nf_put_att_text(ncid, fi3drid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, ft3drid, 'long_name', 45, 'poloidal total energy flux, outboard divertor')
    iret = nf_put_att_text(ncid, ft3drid, 'units', 2, 'W ')
    iret = nf_put_att_text(ncid, fc3drid, 'long_name', 36, 'poloidal current, outboard divertor')
    iret = nf_put_att_text(ncid, fc3drid, 'units', 2, 'A ')
    ! upper outboard divertor quantities
    if(nytr.gt.0) then
      iret = nf_put_att_text(ncid, ne3dtrid, 'long_name', 41, 'electron density, upper outboard divertor')
      iret = nf_put_att_text(ncid, ne3dtrid, 'units', 4, 'm^-3')
      iret = nf_put_att_text(ncid, te3dtrid, 'long_name', 45, 'electron temperature, upper outboard divertor')
      iret = nf_put_att_text(ncid, te3dtrid, 'units', 2, 'eV')
      iret = nf_put_att_double(ncid, te3dtrid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
      iret = nf_put_att_text(ncid, ti3dtrid, 'long_name', 40, 'ion temperature, upper outboard divertor')
      iret = nf_put_att_text(ncid, ti3dtrid, 'units', 2, 'eV')
      iret = nf_put_att_double(ncid, ti3dtrid, 'scale', NCDOUBLE, 1, (/1.0_R8/ev/))
      iret = nf_put_att_text(ncid, tp3dtrid, 'long_name', 42, 'plate temperature, upper outboard divertor')
      iret = nf_put_att_text(ncid, tp3dtrid, 'units', 2, 'K ')
      iret = nf_put_att_text(ncid, po3dtrid, 'long_name', 34, 'potential, upper outboard divertor')
      iret = nf_put_att_text(ncid, po3dtrid, 'units', 2, 'V ')
      iret = nf_put_att_text(ncid, an3dtrid, 'long_name', 37, 'atom density, upper outboard divertor')
      iret = nf_put_att_text(ncid, an3dtrid, 'units', 4, 'm^-3')
      iret = nf_put_att_text(ncid, mn3dtrid, 'long_name', 41, 'molecule density, upper outboard divertor')
      iret = nf_put_att_text(ncid, mn3dtrid, 'units', 4, 'm^-3')
      iret = nf_put_att_text(ncid, fn3dtrid, 'long_name', 51, 'poloidal main species flux, upper outboard divertor')
      iret = nf_put_att_text(ncid, fn3dtrid, 'units', 4, 's^-1')
      iret = nf_put_att_double(ncid, fn3dtrid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
      iret = nf_put_att_text(ncid, fl3dtrid, 'long_name', 47, 'poloidal electron flux, upper outboard divertor')
      iret = nf_put_att_text(ncid, fl3dtrid, 'units', 4, 's^-1')
      iret = nf_put_att_double(ncid, fl3dtrid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
      iret = nf_put_att_text(ncid, fo3dtrid, 'long_name', 42, 'poloidal ion flux, upper outboard divertor')
      iret = nf_put_att_text(ncid, fo3dtrid, 'units', 2, 'W ')
      iret = nf_put_att_double(ncid, fo3dtrid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
      iret = nf_put_att_text(ncid, fe3dtrid, 'long_name', 54, 'poloidal electron energy flux, upper outboard divertor')
      iret = nf_put_att_text(ncid, fe3dtrid, 'units', 2, 'W ')
      iret = nf_put_att_double(ncid, fe3dtrid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
      iret = nf_put_att_text(ncid, fi3dtrid, 'long_name', 49, 'poloidal ion energy flux, upper outboard divertor')
      iret = nf_put_att_text(ncid, fi3dtrid, 'units', 2, 'W ')
      iret = nf_put_att_double(ncid, fi3dtrid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
      iret = nf_put_att_text(ncid, ft3dtrid, 'long_name', 51, 'poloidal total energy flux, upper outboard divertor')
      iret = nf_put_att_text(ncid, ft3dtrid, 'units', 2, 'W ')
      iret = nf_put_att_double(ncid, ft3dtrid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
      iret = nf_put_att_text(ncid, fc3dtrid, 'long_name', 42, 'poloidal current, upper outboard divertor')
      iret = nf_put_att_text(ncid, fc3dtrid, 'units', 2, 'A ')
      iret = nf_put_att_double(ncid, fc3dtrid, 'scale', NCDOUBLE, 1, (/-1.0_R8/))
    endif

    ! leave define mode
    iret = nf_enddef(ncid)
    iret = nf_close(ncid)
    return
  end subroutine b2crtimecdf

  subroutine rwcdf(rw,ncid,data_name,imap,data_set,iret)
#     include <netcdf.inc>
    !
    character*(*) rw,data_name
    integer ncid,imap(*),iret,i,varid,dimlen
    real(kind=R8), Intent(InOut) :: data_set(*)
    character*(maxncnam) dimnam
    integer vartyp,nvdims,start(maxvdims),count(maxvdims), &
         dimids(maxvdims)
    !
    character*(*) timnam
    character*(maxncnam) timsav
    integer ntsav,ntstep
    integer :: istride, imax
    save timsav,ntsav
    data timsav /'!!!! INVALID NAME !!!!'/
    external subini, subend, xerrab
    !
    call subini ('rwcdf')
    iret = nf_inq_varid(ncid,data_name,varid)
    if(iret.ne.0) call xerrab ('Data name not declared')
    iret = nf_inq_varndims(ncid,varid,nvdims)
    iret = nf_inq_vardimid(ncid,varid,dimids)
    count(1) = 1 ! for scalars
    start(1) = 1
    do i=1,nvdims
      iret = nf_inq_dim(ncid,dimids(i),dimnam,dimlen)
      count(i) = dimlen
      if(dimnam.eq.timsav) then
        start(i)=ntsav
        count(i)=1
      else
        start(i)=1
      endif
    enddo
    istride = 1
    imax = 1
    do i=1,nvdims-1
      istride = istride*imap(i)
      imax = imax + imap(i)*count(i)
    enddo
    iret = nf_inq_vartype(ncid,varid,vartyp)
    if(rw.eq.'read') then
      Select Case (vartyp)
      Case (NCDOUBLE)
        iret = nf_get_vara_double(ncid,varid,start,count, &
             data_set(1:imax:istride))
      Case Default
        call xerrab ('Unknown data type in rwcdf read')
      End Select
    elseif(rw.eq.'write') then
      Select Case (vartyp)
      Case (NCDOUBLE)
        iret = nf_put_vara_double(ncid,varid,start,count, &
             data_set(1:imax:istride))
      Case Default
        call xerrab ('Unknown data type in rwcdf write')
      End Select
    else
      write(*,*) 'Either "read" or "write" must be chosen, not', rw
    endif

    call subend ()
    return
    !
    entry rwcdf_settime(timnam,ntstep)
    write(*,*) 'Saving ',Trim(timnam), &
         ' as the time dimension'
    write(*,*) 'ntstep = ',ntstep
    timsav=timnam
    ntsav=ntstep
    return
  end subroutine rwcdf
#endif
!
  subroutine output_ds(crx,cry,nx,ny,iref,target_offset, &
       jsep,iystart,iyend,filename)
    use b2mod_indirect
    implicit none
    integer nx,ny,iref,jsep,iystart,iyend,target_offset
    real (kind=R8) :: &
         crx(-1:nx,-1:ny,0:3),cry(-1:nx,-1:ny,0:3)
    real (kind=R8) :: &
         ds(-1:ny), ds_offset
    character*(*) filename
    integer ix,iy
    external subini, subend, xertst
    intrinsic sqrt
    real (kind=R8) :: &
         cr,cz

    cr(ix,iy)=0.25_R8* &
         (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
    cz(ix,iy)=0.25_R8* &
         (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))

    call subini ('output_ds')
    iystart=-1
    do while (region(iref,iystart,0).eq.0 .and. iystart.lt.ny)
      iystart=iystart+1
    enddo
    call xertst(iystart.le.ny, 'faulty parameter iystart')
    iyend=ny
    do while (region(iref,iyend,0).eq.0 .and. iyend.gt.-1)
      iyend=iyend-1
    enddo
    call xertst(iyend.ge.-1, 'faulty parameter iyend')
    if(iystart.eq.ny.and.iyend.eq.-1) then
      ! special case [DPC]
      iystart=-1
      iyend=ny
      write(*,*) 'special treatment for iystart, iyend for ', &
           trim(filename)
    endif
    call xertst(iyend.ge.iystart, 'faulty parameter iystart & iyend')
    ds(iystart)= &
         sqrt((cr(iref+target_offset,iystart)- &
         0.5_R8*(crx(iref+target_offset,iystart,0)+ &
         crx(iref+target_offset,iystart,1)))**2+ &
         (cz(iref+target_offset,iystart)- &
         0.5_R8*(cry(iref+target_offset,iystart,0)+ &
         cry(iref+target_offset,iystart,1)))**2)
    do iy=iystart+1,iyend
      ds(iy)=ds(iy-1)+ &
           sqrt((cr(iref+target_offset,iy)- &
           cr(iref+target_offset,iy-1))**2+ &
           (cz(iref+target_offset,iy)- &
           cz(iref+target_offset,iy-1))**2)
    enddo
    if(iystart.le.jsep.and.iyend.gt.jsep) then
      ds_offset=(ds(jsep)+ds(jsep+1))/2.0_R8
      do iy=iystart,iyend
        ds(iy)=ds(iy)-ds_offset
      enddo
    endif
    open(99,file=filename)
    do iy=iystart,iyend
      write(99,*) ds(iy)
    enddo
    close(99)
    call subend ()
    return
  end subroutine output_ds


  Subroutine calc_fet(ix,iy,SIDE,fac_flux,nx,ny,ns,ismain,BoRiS,fet,fni0,fee0,fei0,fch0,pwr)
    use b2mod_plasma   , Only : fna, fhe, fhi, fch, fhm, fhp
    use b2mod_indirect , Only : rightix, rightiy, bottomix, bottomiy, topix, topiy, leftix, leftiy
    use b2mod_external , Only : fhi_ext, pt_ext, ua_ext, am_ext, ns_ext, fa_ext
    use b2mod_constants , Only : ev, mp
    Use b2mod_geo , Only : hx, hy, qz, gs
    Implicit None
    Integer, Intent(In) :: ix, iy
    Integer, Intent(In) :: ismain, nx, ny, ns
    Real(kind=R8), Intent(In) :: BoRiS, fac_flux
    Real(kind=R8), Intent(Out) :: fet
    Real(kind=R8), Intent(Out), Optional :: fni0, fee0, fei0, fch0, pwr
    Character(Len=1) :: SIDE
    ! Local vars
    Integer :: ix_adj, iy_adj, is, ix_flux, iy_flux, idir
    Real(kind=R8) :: kintmp, rpttmp
    Real(kind=R8) :: h(-1:nx,-1:ny)
    ! computation

    Select Case (SIDE)
    Case ('l','L')
      ix_flux = rightix(ix,iy) ! Index to cell with flux entering cell
      iy_flux = rightiy(ix,iy)        
      ix_adj  = rightix(ix,iy) ! Index to cell adjacent
      iy_adj  = rightiy(ix,iy)        
      idir = 0                 ! Index in flux variables (x vs y direction)
      h(-1:nx,-1:ny) = hx(-1:nx,-1:ny)
    Case ('r','R')
      ix_flux = ix
      iy_flux = iy
      ix_adj = leftix(ix,iy)
      iy_adj = leftiy(ix,iy)        
      idir = 0
      h(-1:nx,-1:ny) = hx(-1:nx,-1:ny)
    Case ('t','T')
      ix_flux = ix
      iy_flux = iy
      ix_adj = bottomix(ix,iy)
      iy_adj = bottomiy(ix,iy)
      idir = 1
      h(-1:nx,-1:ny) = hy(-1:nx,-1:ny)*qz(-1:nx,-1:ny,1)
    Case ('b','B')        
      ix_flux = topix(ix,iy)
      iy_flux = topiy(ix,iy)
      ix_adj = topix(ix,iy)
      iy_adj = topiy(ix,iy)
      idir = 1
      h(-1:nx,-1:ny) = hy(-1:nx,-1:ny)*qz(-1:nx,-1:ny,1)
    Case Default
      Call xerrab('Unknown SIDE in calc_fet')
    End Select
    If (Present(fni0)) fni0 = fac_flux*fna(ix_flux,iy_flux,idir,ismain)
    If (Present(fee0)) fee0 = fac_flux*fhe(ix_flux,iy_flux,idir)
    If (Present(fei0)) fei0 = fac_flux*fhi(ix_flux,iy_flux,idir)
    If (Present(fch0)) fch0 = fac_flux*fch(ix_flux,iy_flux,idir)
    fet = fac_flux*(fhe(ix_flux,iy_flux,idir) + fhi(ix_flux,iy_flux,idir) + fhi_ext(ix_flux,iy_flux,idir))
    do is=0,ns-1
      fet = fet + fac_flux*fhm(ix_flux,iy_flux,idir,is)*(1.0_R8-BoRiS) + fac_flux*fhp(ix_flux,iy_flux,idir,is)
    enddo
    do is=0,ns_ext-1
      kintmp = 0.5_R8*am_ext(is)*mp*(ua_ext(ix,iy,is)**2 * h(ix_adj,iy_adj)+ &
           ua_ext(ix_adj,iy_adj,is)**2*h(ix,iy))/(h(ix_adj,iy_adj)+h(ix,iy))
      rpttmp = (pt_ext(ix,iy,is)*h(ix_adj,iy_adj)+pt_ext(ix_adj,iy_adj,is)*h(ix,iy))/(h(ix_adj,iy_adj)+h(ix,iy))
      fet = fet + fac_flux*(rpttmp*ev + kintmp*(1.0_R8-BoRiS))*fa_ext(ix_flux,iy_flux,idir,is)
    enddo
    If (Present(pwr)) pwr = Abs(fet)/gs(ix_flux,iy_flux,idir)

  End Subroutine calc_fet
    
  
End Module b2mod_mwti
