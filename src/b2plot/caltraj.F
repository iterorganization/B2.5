C     =======================================================================
C     SUBROUTINES FOR B2PLOT COMMAND TRAJ
C     =======================================================================
C
C     This routine is used to plot trajectories of particles using the
C     velocity fields uu and vv. To call this routine
C     one must specify at least 8 parameters: ixmin, iymin, ixmax, iymax,
C     # of particles, 1 or 0 for movie or lines, species indices (at least
C     one), and last but not least the total number of species indices one
C     has given in the step before. A valid call could be:
C     0 0 95 35 100 0 0 1 traj
C     This should plot the trajectories of 100 particles which start
C     in any of the cells between cell (0,0) (southwest corner) and
C     cell (95,35) (northeast corner) as lines for hydrogen (0=H).
C     Note that the cell and species indices start at zero!
C     The number of time steps (plots) is given by the "t_#p" command.
C     The total flight time of the particles is given by the "t_pt" command.
C     The velocity fields uu and vv can be overwritten using the
C     "uuse" and "vvse" commands. Their defaults values are ua and fna/na/sy
C     respectively. If the number of particles is positive, then the
C     particles are born randomly, with the probability of being born in a
C     cell being proportional to the cell volume. This probability weighting
C     can be changed by calling "traj" with a negative number of particles.
C     Then the weighting function is whatever field is on top of the matrix
C     stack. For example, this set of commands does the following:
C     fnax na m/ psxperp 0 zsel uuse  ! uutraj set to poloidal neutral flow velocity
C     fnay na m/ psy 0 zsel vvse      ! vvtraj set to radial neutral flow velocity
C     1000 t_#p                       ! 1000 time steps
C     1e-4 t_pt                       ! followed for a total time of 1.e-4 s.
C     na 0 zsel                       ! weighted with neutral density
C     0 0 0 0 -1000 0 0 1 traj        ! over the entire physical domain (0 0 0 0)
C               ! 1000 particles using non-default weighting function (-1000)
C               ! trajectories plotted (0)
C               ! for species 0 (=H neutrals) (0)
C               ! total of species considered is 1 (1)
C     =======================================================================

c     -----------------------------------------------------------------------
c     SETTING THE NUMBER OF PLOTS
c     -----------------------------------------------------------------------

      subroutine t_setnplot(n)

      implicit none
      integer n,nplot
      common /t_com/ nplot

      nplot=n
      write(*,*) 'traj: Setting number of plots: ',nplot

      end subroutine t_setnplot

c     -------------------------------------------------------------------------
c     PLOTTING TRAJECTORIES
c     -------------------------------------------------------------------------

      subroutine caltraj(iarg,specnum)

      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2us_data
#ifndef NO_CDF
      use b2mod_cdf
#endif
      implicit none

c     INPUT
      integer iarg(1:ns+7)              !parameters
      integer specnum                   !location of the number of species

c     the traj common
      integer nplot                     !number of plots
      common /t_com/ nplot              !getting the values
      logical ratio
      common /t_ratio/ ratio

      real (kind=R8), allocatable ::
     1 warr(:,:)                        !the stack

      real (kind=R8) ::
     *   globaltime                     !plottime

c                                       !number of particles per cell
      integer ialt,jalt                 !position in computational domain
      integer ipx(1:npart),ipy(1:npart),ips(1:npart)
      integer ipx_start(1:npart),ipy_start(1:npart)
      real (kind=R8) ::
     *   p_r(1:npart),p_z(1:npart)
      real (kind=R8) ::
     *   p_r_start(1:npart),p_z_start(1:npart)
      real (kind=R8) ::
     *   clr(-1:nx,-1:ny,1:4),          !coordinates of centers of cell borders
     +   clz(-1:nx,-1:ny,1:4)

      real (kind=R8) ::
     *   deltat                         !timestep between plots
      logical ls_active(0:ns-1)         !species active?
      integer isclvl(0:ns-1)
      character*16 texttraj(1:ns)

      integer ip,randzaehler,zaehler
#ifdef WG_TODO
      integer zaehler2
#endif

      integer counter

#ifdef NAG22
      integer IGEN
      integer ISEED(4)
      common /nag22/ igen, iseed
      external G05KBF
#endif
      real (kind=R4) ::
     *   cpu_start, cpu_end

      integer destmode                     !mode for destroying

      integer is,np
      real (kind=R8) ::
     *   zrold,rrold

      integer icol,oocol,lnblnk
      external lnblnk
      character argstring*256,tmpstring*256,titlestring*256
      real (kind=R8) ::
     *   currtime

#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
#     include <netcdf.inc>
      character rw*5
      integer ncid,imap(maxvdims),iret
      integer itstep,ntstep,iCv
      real (kind=R8), allocatable :: lwarr(:,:)
      real (kind=R8) ::
     1   starttime, tstepn(1), dumarr(1)
#endif

      oocol=-1
      allocate(warr(mpg%nCv,0:ns-1))
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
      allocate(lwarr(mpg%nCv,0:ns-1))
      if(moviename .eq. ' ') then
#endif
c     Getting the values for globaltime and deltat
        globaltime=values(16)
        deltat=globaltime/nplot
        currtime=0.0_R8
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
      else
        rw='read'
        iret = nf_open(moviename,NCNOWRIT,ncid)
        imap(1)=1
        imap(2)=1
        call rwcdf (rw, ncid, 'ntstep', imap, tstepn, iret)
        ntstep = nint(tstepn(1))
        nplot=ntstep+1
        itstep=0
        dumarr(1)=starttime
        call rwcdf_settime('time',1)
        call rwcdf(rw,ncid,'time',imap,dumarr,iret)
        call rwcdf(rw,ncid,'na',imap,lwarr,iret)
        currtime=starttime
        dumarr(1)=globaltime
        call rwcdf_settime('time',ntstep)
        call rwcdf(rw,ncid,'time',imap,dumarr,iret)
        globaltime=globaltime-starttime
      endif
#endif

c     writing callup information
      if (iarg(6) .eq. 1) then
         if (b2plot_debug.ge.2) write(*,*) 'MOVIE'
      else
         if (b2plot_debug.ge.2) write(*,*) 'TRAJECTORIES'
      endif
      if(iarg(1) .eq. 0 .or. iarg(2) .eq. 0 .or.
     1   iarg(3) .eq. 0 .or. iarg(4) .eq. 0) then
        destmode=2
      else
        destmode=0
      endif
      if (b2plot_debug.ge.8) write(*,*) 'Set destmode to ',destmode
      if(iarg(1) .le. 0) iarg(1)=0
      if(iarg(2) .le. 0) iarg(2)=0
      if(iarg(3) .eq. 0 .or. iarg(3) .ge. nx) iarg(3)=nx-1
      if(iarg(4) .eq. 0 .or. iarg(4) .ge. ny) iarg(4)=ny-1
      if (b2plot_debug.ge.2)
     + write(*,*) 'plottime in seconds (default=5.0e-4): ', globaltime
      if (b2plot_debug.ge.2)
     + write(*,*) 'number of plots (default=10): ', nplot

c     getting the species numbers
      do is=0,ns-1
        ls_active(is)=.false.
      enddo
      do zaehler=0,iarg(specnum)-1
        ls_active(iarg(7+zaehler))=.true.
        texttraj(zaehler+1)=textpl(iarg(7+zaehler))
      enddo

c     find weighting function for the particle positions
      if (iarg(5) .gt. 0) then      !default: weighting with volume
        if (b2plot_debug.ge.1)
     >   write(*,*) 'weighting randomiser with volume'
        call weightarray(warr,geo%cvVol,
     1         ls_active,iarg(1),iarg(3),iarg(2),iarg(4),1)
      else if (iwrk .gt. 0) then    !user-defined weighting
        if (b2plot_debug.ge.1)
     >   write(*,*) 'weighting randomiser user defined'
        call weightarray(warr,work(1,0,iwrk),
     1   ls_active,iarg(1),iarg(3),iarg(2),iarg(4),ndim(2,iwrk))
      endif

c     initialising the randomiser
#ifndef NAG22
      call g05cbf(1000)
#else
      IGEN=0
      ISEED(1)=1000
      call g05kbf(IGEN,ISEED)
#endif

c     open the Plot-Screen
      tmpstring='Trajectories'
      if (iarg(specnum).eq.1)
     > tmpstring=tmpstring(1:lnblnk(tmpstring))//
     +           texttraj(1)(1:lnblnk(texttraj(1)))
      call customise_label(tmpstring,titlestring)
      call init_page(titlestring,values(5),values(6),
     +               values(7),values(8),1,.true.,.true.,.true.)
      call GOPWK(9,1,3)         !open GFLASH

c     generate beginning positions

      np=abs(iarg(5))
      call cpu_time(cpu_start)
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
      if(moviename .ne. ' ') then
        write(*,*) 'moviename = ',moviename
        if (b2plot_debug.ge.1) write(*,*) 'Changing weighting'
        do is=0,ns-1
          warr(:,is)=warr(:,is)/state%pl%na(:,is)*lwarr(:,is)
        enddo
      endif
#endif
      call allocate_particles(warr,ipx,ipy,ips,p_r,p_z,0,np)
      call cpu_time(cpu_end)
      if (b2plot_debug.ge.8) write(*,'(a,f6.3,a)')
     1 'Initial allocate_particles took ',
     2 cpu_end-cpu_start,
     3 ' cpu-seconds'
      do ip=1,np
        ipx_start(ip)=ipx(ip)
        ipy_start(ip)=ipy(ip)
        p_r_start(ip)=p_r(ip)
        p_z_start(ip)=p_z(ip)
      enddo

C-----Producing the arrays of side middle-points. The arrays are named
C     clr (center line radial) and clz (center line z) respectively.

#ifdef WG_TODO
      do zaehler=-1,nx
        do zaehler2=-1,ny
          clr(zaehler,zaehler2,1)=crx(zaehler,zaehler2,1)+
     +            (crx(zaehler,zaehler2,3)-crx(zaehler,zaehler2,1))/2.
          clz(zaehler,zaehler2,1)=cry(zaehler,zaehler2,1)+
     +            (cry(zaehler,zaehler2,3)-cry(zaehler,zaehler2,1))/2.
          clr(zaehler,zaehler2,2)=crx(zaehler,zaehler2,3)+
     +            (crx(zaehler,zaehler2,2)-crx(zaehler,zaehler2,3))/2.
          clz(zaehler,zaehler2,2)=cry(zaehler,zaehler2,3)+
     +            (cry(zaehler,zaehler2,2)-cry(zaehler,zaehler2,3))/2.
          clr(zaehler,zaehler2,3)=crx(zaehler,zaehler2,2)+
     +            (crx(zaehler,zaehler2,0)-crx(zaehler,zaehler2,2))/2.
          clz(zaehler,zaehler2,3)=cry(zaehler,zaehler2,2)+
     +            (cry(zaehler,zaehler2,0)-cry(zaehler,zaehler2,2))/2.
          clr(zaehler,zaehler2,4)=crx(zaehler,zaehler2,0)+
     +            (crx(zaehler,zaehler2,1)-crx(zaehler,zaehler2,0))/2.
          clz(zaehler,zaehler2,4)=cry(zaehler,zaehler2,0)+
     +            (cry(zaehler,zaehler2,1)-cry(zaehler,zaehler2,0))/2.
        enddo
      enddo
#endif

c     initialise graphics

      call GFLAS1(1)
      call gsplci(1)
      if(iarg(specnum) .eq. 1 .and. iarg(6) .eq. 0) then  ! 1 species and trajectories
        call data_copy(iarg(7))
        call additional_elements(mpg%nCv,.true.) !draws mesh/separatrix/etc... if wanted
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
        if(moviename .eq. ' ') then
#endif
          call lblbartraj(globaltime*1e-3_R8,globaltime,tlinlog)
#if( !defined(NO_CDF) && defined(CDF_MOVIE) )
        else
          call lblbartraj(starttime,globaltime+starttime,1)
        endif
#endif
      else
        call data_copy(0)
        call additional_elements(mpg%nCv,.true.) !draws mesh/separatrix/etc... if wanted
        call specieslblbar(iarg(specnum),isclvl,texttraj)
      endif
      call GFLAS2                              ! sv 27.02.01
      call GFLAS3(1)

C Now we start.:) All particles are calculated and plotted until they
C reach the wall or their total time runs out (particle life). Each particle
C is always advanced one timestep deltaT before the next particle is
C calculated. Some particles die, because there occurred an error. If
C this was the case, the new position of the particle in the computational
C domain is set to be ix=-1,iy=-1 (ii,jj), in which case they are not
C processed any further.

c     MOVIE:
c     in each step a new graphics screen is opened
c     after calculating the new positions of the particles, their positions
c     are represented by dots

c     The density of the plotted particles in the cells should be constant.
c     REFILLING: see subroutine calcweigher
c     DESTROYING: mode 1: If particles cross the core boundary, they are
c                         destroyed if the integral number of particles on
c                         the side they entered is bigger than at beginning.
c                 mode 2: If a particle enters a cell in which there are more
c                         particles than at beginning, it is destroyed.

      if (iarg(6) .eq. 1) then
        do ip=1,np        !draw initial positions on screen
          call lastd
          call gsplci(ips(ip)+1)
          call lined_d(p_r(ip),p_z(ip),p_r(ip),p_z(ip))
        enddo
        do counter=2,nplot
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
          if(moviename .eq. ' ') then
#endif
            globaltime=globaltime-deltat
            currtime=currtime+deltat
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
          else
            itstep=counter-1
            imap(1)=1
            dumarr(1)=currtime
            call rwcdf_settime('time',itstep)
            call rwcdf(rw,ncid,'time',imap,dumarr,iret)
            if (b2plot_debug.ge.2)
     .       write(*,*) 'Processing from ',starttime,' to ',currtime
            deltat=currtime-starttime
            starttime=currtime
            globaltime=globaltime-deltat
            if(itstep .eq. ntstep) then
              if (b2plot_debug.ge.2)
     .         write(*,*) 'Last step: globaltime = ',globaltime
              globaltime=0.0_R8
            endif
            imap(1)=1
            imap(2)=1
            call rwcdf(rw,ncid,'ua',imap,state%pl%ua,iret)
            call rwcdf(rw,ncid,'fna',imap,state%dv%fna,iret)
            call rwcdf(rw,ncid,'na',imap,lwarr,iret)
            do is=0,ns-1
              do iCv=1,mpg%nCv
                uutraj(iCv,is)=state%pl%ua(iCv,is)
#ifdef WG_TODO
                vvtraj(ix,iy,is)=fna(ix,iy,1,1,is)
     .             /na(ix,iy,is)/gs(ix,iy,1)
#endif
              end do
            end do
          endif
#endif
          if (b2plot_debug.ge.2)
     .     write(*,*)'plotnumber: ',counter,'(of ',nplot,')'
          call end_page
          call lastd
          call gsplci(1)
          write(argstring,'(f10.3,a)') currtime*1.0e6_R8,' us'
          call init_page('Trajectories ' //
     2     argstring(1:lnblnk(argstring)),values(5),values(6),
     +     values(7),values(8),1,.true.,.true.,.true.)
          call gflas3(1)
          call cpu_time(cpu_start)
          do ip=1,np
            call lastd
            call gsplci(ips(ip)+1)
            ialt=ipx(ip)
            jalt=ipy(ip)
            if(ialt .lt. 0 .or. ialt .gt. nx-1 .or.
     1         jalt .lt. 0 .or. jalt .gt. ny-1) then
              write(*,*)
     1         'Problem with particle ',ip,': ialt,jalt = ',ialt,jalt
            else
              call trajectories(p_r(ip),p_z(ip),ipx(ip),ipy(ip),
     1         ips(ip),deltat,clr,clz)
              call lined_d(p_r(ip),p_z(ip),p_r(ip),p_z(ip))
            endif
          enddo                 ! ip
          call cpu_time(cpu_end)
          if (b2plot_debug.ge.8) write(*,'(a,f6.3,a)')
     1     'Advancing particles took ',
     2      cpu_end-cpu_start,
     3     ' cpu-seconds'
          if(destmode .eq. 2) then
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
            if(moviename .ne. ' ') then
              do is=0,ns-1
                warr(:,is)=warr(:,is)/state%pl%na(:,is)*lwarr(:,is)
              enddo
            endif
#endif
            call cpu_time(cpu_start)
            call allocate_particles(warr,ipx,ipy,ips,p_r,p_z,np,np)
            call cpu_time(cpu_end)
            if (b2plot_debug.ge.8) write(*,'(a,f6.3,a)')
     1       'Subsequent allocate_particles took ',
     2       cpu_end-cpu_start,
     3       ' cpu-seconds'
          else
            do ip=1,np
              if(ipx(ip) .lt. 0 .or. ipx(ip) .gt. nx-1 .or.
     1           ipy(ip) .lt. 0 .or. ipy(ip) .gt. ny-1) then
                ipx(ip)=ipx_start(ip)
                ipy(ip)=ipy_start(ip)
                p_r(ip)=p_r_start(ip)
                p_z(ip)=p_z_start(ip)
              endif
            enddo
          endif
        enddo                   ! number of plots

c     TRAJECTORIES:

      else                      ! Draw Lines
        randzaehler=0
        do while ((randzaehler .lt. np) .and.
     +   (globaltime .gt. 0.0_R8))
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
          if(moviename .eq. ' ') then
#endif
            globaltime=globaltime-deltat
            currtime=currtime+deltat
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
          else
            itstep=itstep+1
            imap(1)=1
            dumarr(1)=currtime
            call rwcdf_settime('time',itstep)
            call rwcdf(rw,ncid,'time',imap,dumarr,iret)
            if (b2plot_debug.ge.2)
     .       write(*,*) 'Processing from ',starttime,' to ',currtime
            deltat=currtime-starttime
            starttime=currtime
            globaltime=globaltime-deltat
            if(itstep .eq. ntstep) then
              if (b2plot_debug.ge.2)
     .         write(*,*) 'Last step: globaltime = ',globaltime
              globaltime=0.0_R8
            endif
            imap(1)=1
            imap(2)=1
            call rwcdf(rw,ncid,'ua',imap,state%pl%ua,iret)
            call rwcdf(rw,ncid,'na',imap,state%pl%na,iret)
            call rwcdf(rw,ncid,'fna',imap,state%dv%fna,iret)
            do is=0,ns-1
              do iCv=1,mpg%nCv
                uutraj(iCv,is)=state%pl%ua(iCv,is)
#ifdef WG_TODO
                vvtraj(ix,iy,is)=fna(ix,iy,1,1,is)
     .           /na(ix,iy,is)/gs(ix,iy,1)
#endif
              end do
            end do
          endif
#endif
          if(iarg(specnum) .eq. 1) then
            call lblbartrajcol(currtime,icol)
            if(icol .ne. oocol .and. color_set) then
              call lastd
              call gsplci(icol)
              oocol=icol
            endif
          endif
          randzaehler=0
          do ip=1,np
            if (ipy(ip).gt.ny-1 .or. ipy(ip).lt.0 .or.
     +          ipx(ip).gt.nx-1 .or. ipx(ip).lt.0) then
              randzaehler=randzaehler+1
            else
              rrold=p_r(ip)
              zrold=p_z(ip)
              call trajectories(p_r(ip),p_z(ip),ipx(ip),
     +         ipy(ip),ips(ip),deltat,clr,clz)
              if(iarg(specnum).ne.1) then
                icol=ips(ip)+1
                if(icol.ne.oocol .and. color_set) then
                  call lastd
                  call gsplci(icol)
                  oocol=icol
                endif
              endif
              call lined_d(rrold,zrold,p_r(ip),p_z(ip))
            endif
          enddo         !ip
        enddo           !randzaehler
      endif

c     close Plot-Screen

      call end_page
      call GCLWK(9)          !closing GFLASH
      call lastd
      call gsplci(1)

#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
      iret = nf_close(ncid)
      deallocate(lwarr)
#endif
      deallocate(warr)
      return
      end subroutine caltraj

c     -----------------------------------------------------------------------
c     CALCULATING NEW POSITIONS
c     -----------------------------------------------------------------------

      subroutine trajectories(rre,zre,ii,jj,l,deltat,clr,clz)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_debug
      implicit none

C      define particle time tpart, and timestep deltat for plotting :
C      rre/zre real coordinates of the particle
C      rxy and zxy intersection points with the walls
C      d variable needed to check through which wall the particle goes
C      wall: holds through which wall it goes (11=corner 1, 13 wall 13..)
C      Neighbourhood arrays to determine into which cell the
C      particle moves.
C      vr,vz reale transformierte (in physical domain) Geschwindigkeiten

      real (kind=R8) ::
     *   tpart,deltat,timerun
      real (kind=R8) ::
     *   vr,vz
      real (kind=R8) ::
     *   rre,zre, radvance, zadvance
      real (kind=R8) ::
     *   tdist
      real (kind=R8) ::
     *   clr(-1:nx,-1:ny,1:4),clz(-1:nx,-1:ny,1:4)
      integer ii,jj,l,iold,jold
      integer wall,schritte
      logical laufende
      logical inside_rectangle

C     Initial configuration. Particle time is counted downwards,
C     therefore is set to be deltaT here.

      vr=0.0_R8
      vz=0.0_R8
      laufende = .false.
      tpart = deltat
      schritte=0

      do while(.not.laufende)
         schritte=schritte+1

C     The velocities need to be transformed from the computational
C     to the real domain.

         call veloctrans(ii,jj,l,vr,vz,clr,clz)

C      OK, first calculate the time the particle needs to hit a "wall"
C      That is, at first one has to know, which wall it is going to hit!
C      :)

         call intersection(wall,tdist,vz,vr,rre,zre,ii,jj,laufende)
#ifdef WG_TODO
         if (laufende .and. b2plot_debug.ge.12) then
            write(*,*) uutraj(ii,jj,l),vvtraj(ii,jj,l),vr,vz
         endif
#endif

C      We now know through which wall it goes (wall=xy) and how
C      much time the particle is going to need to reach the wall.
C      Now let us see if it makes it to the wall before having to
C      be plotted, that is, before the time tpart runs out, which is
C      counted down from deltaT

         timerun=MIN(tdist,tpart)

         if (timerun .eq. tpart) then
            tpart = deltaT
            laufende=.true.
         else
            tpart = tpart - timerun
            iold = ii
            jold = jj
            if (wall .eq. 13) then
               ii = rightix(iold,jold)
               jj = rightiy(iold,jold)
            else if (wall .eq. 32) then
               ii = topix(iold,jold)
               jj = topiy(iold,jold)
            else if (wall .eq. 20) then
               ii = leftix(iold,jold)
               jj = leftiy(iold,jold)
            else if (wall .eq. 01) then
               ii = bottomix(iold,jold)
               jj = bottomiy(iold,jold)
            else if (wall .eq. 11) then
               jj = bottomiy(rightix(iold,jold),rightiy(iold,jold))
               ii = bottomix(rightix(iold,jold),rightiy(iold,jold))
            else if (wall .eq. 33) then
               jj = topiy(rightix(iold,jold),rightiy(iold,jold))
               ii = topix(rightix(iold,jold),rightiy(iold,jold))
            else if (wall .eq. 22) then
               jj = topiy(leftix(iold,jold),leftiy(iold,jold))
               ii = topix(leftix(iold,jold),leftiy(iold,jold))
            else if (wall .eq. 00) then
               jj = bottomiy(leftix(iold,jold),leftiy(iold,jold))
               ii = bottomix(leftix(iold,jold),leftiy(iold,jold))
            endif
         endif

C        calculate new particle coordinates, advancing them a little
C        bit into the next cell.
CXPB     deltar and deltaz are normally zero, except when crossing periodicity
CXPB     boundaries, then they translate the particle to the spatial
CXPB     location of the next computational cell.

         radvance=vr*MAX(timerun,1.0e-6_R8)*0.01_R8
         zadvance=vz*MAX(timerun,1.0e-6_R8)*0.01_R8

         rre = rre + vr*timerun + radvance
         zre = zre + vz*timerun + zadvance

         do while(
     1    .not. inside_rectangle(rre,zre,crx(ii,jj,1),cry(ii,jj,1),
     2         crx(ii,jj,3),cry(ii,jj,3),crx(ii,jj,2),cry(ii,jj,2),
     3         crx(ii,jj,0),cry(ii,jj,0)))
           radvance=radvance*0.5_R8
           zadvance=zadvance*0.5_R8
           if (Max(Abs(radvance),Abs(zadvance)) .lt. 1.0e-20_R8) then
             if (b2plot_debug.ge.1)
     .        write(*,*)'Particle had to be killed.'
             if (b2plot_debug.ge.1) write(*,*)'No Problem. (?)'
             ii=-1
             jj=-1
             goto 1000
           endif
           rre = rre - radvance
           zre = zre - zadvance
         enddo

 1000    continue

         if(schritte .eq. 10000) then
           if (b2plot_debug.ge.1)
     .      write(*,*) 'Particle had to be killed to save time.'
           if (b2plot_debug.ge.1)
     .      write(*,*) 'A new one will be created to take its place.'
           ii=-1
           jj=-1
         endif

C     Check if particle reached the border.

         if(((jj .gt. ny-1) .or. (jj .lt. 0)) .or.
     1      ((ii .gt. nx-1) .or. (ii .lt. 0))) then
           laufende = .true.
         endif
      enddo
      return
      end subroutine trajectories

c     -------------------------------------------------------------------------
c     CALCULATING THE INTERSECTION WITH THE CELL WALL
c     -------------------------------------------------------------------------

      subroutine intersection(wall,tdist,vz,vr,rre,zre,ii,jj,laufende)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot_debug
      implicit none

      real (kind=R8) ::
     *   r13,r32,r20,r01,z13,z32,z20,z01,rre,zre
      real (kind=R8) ::
     *   d,tdist
      real (kind=R8) ::
     *   vz,vr
      integer wall,ii,jj
      logical laufende,error

      error=.true.

C     Calculate the slope, that is the a in y=ax+b

      call schnittpunktzweierlinien(vr,vz,rre,zre,ii,jj,r13,z13,1,3)
      call schnittpunktzweierlinien(vr,vz,rre,zre,ii,jj,r32,z32,3,2)
      call schnittpunktzweierlinien(vr,vz,rre,zre,ii,jj,r20,z20,2,0)
      call schnittpunktzweierlinien(vr,vz,rre,zre,ii,jj,r01,z01,0,1)

C     Now we know the intersection between the path of the particle
C     and the lines that make up the walls. But we do not know if they
C     are inside the walls of the cells or on the outer part of the
C     line making up the cell.


      if(((r13 .le. crx(ii,jj,3))  .and.  (r13 .ge. crx(ii,jj,1))) .or.
     +   ((r13 .ge. crx(ii,jj,3))  .and.  (r13 .le. crx(ii,jj,1)))) then
        d = r13 + vr
        if (((r13 .le. d)  .and.  (r13 .ge. (rre )))  .or.
     +      ((r13 .ge. d)  .and.  (r13 .le. (rre )))) then
          wall=13
          error=.false.
          tdist = MAX( ABS(r13-rre)/ABS(vr+1.0e-20_R8),
     +                 ABS(z13-zre)/ABS(vz+1.0e-20_R8))
        endif
      endif
      if(((r32 .le. crx(ii,jj,3))  .and.  (r32 .ge. crx(ii,jj,2))) .or.
     +   ((r32 .ge. crx(ii,jj,3))  .and.  (r32 .le. crx(ii,jj,2)))) then
        d = z32 + vz
        if (((z32 .le. d)  .and.  (z32 .ge. (zre)))  .or.
     +      ((z32 .ge. d)  .and.  (z32 .le. (zre)))) then
          wall=32
          error=.false.
          tdist = MAX( ABS(r32-rre)/ABS(vr+1.0e-20_R8),
     +                 ABS(z32-zre)/ABS(vz+1.0e-20_R8))
        endif
      endif
      if(((z20 .le. cry(ii,jj,2))  .and.  (z20 .ge. cry(ii,jj,0)))  .or.
     +   ((z20 .ge. cry(ii,jj,2))  .and.  (z20 .le. cry(ii,jj,0)))) then
        d = r20 + vr
        if (((r20 .le. d)  .and.  (r20 .ge. (rre )))  .or.
     +      ((r20 .ge. d)  .and.  (r20 .le. (rre )))) then
          wall=20
          error=.false.
          tdist = MAX( ABS(r20-rre)/ABS(vr+1.0e-20_R8),
     +                 ABS(z20-zre)/ABS(vz+1.0e-20_R8))
        endif
      endif
      if(((r01 .le. crx(ii,jj,1))  .and.  (r01 .ge. crx(ii,jj,0)))  .or.
     +   ((r01 .ge. crx(ii,jj,1))  .and.  (r01 .le. crx(ii,jj,0)))) then
        d = z01 + vz
        if (((z01 .le. d)  .and.  (z01 .ge. (zre)))  .or.
     +      ((z01 .ge. d)  .and.  (z01 .le. (zre)))) then
          wall=01
          error=.false.
          tdist = MAX( ABS(r01-rre)/ABS(vr+1.0e-20_R8),
     +                 ABS(z01-zre)/ABS(vz+1.0e-20_R8))
        endif
      endif

C Now check the corners for those very rare cases when a particle
C passes through one of them.

       if (r13 .eq. r01) then
           d = z01 + vz
           if (((z01 .le. d)  .and.  (z01 .ge. (zre)) ) .or.
     +         ((z01 .ge. d)  .and.  (z01 .le. (zre)))) then
              error=.false.
               wall=11
               tdist = MAX( ABS(r01-rre)/ABS(vr+1.0e-20_R8),
     +                      ABS(z01-zre)/ABS(vz+1.0e-20_R8))
           endif
       endif
       if (r32 .eq. r13) then
          d = z13 + vz
          if (((z13 .le. d)  .and.  (z13 .ge. (zre)))  .or.
     +        ((z13 .ge. d)  .and.  (z13 .le. (zre)))) then
              wall=33
              error=.false.
              tdist = MAX( ABS(r32-rre)/ABS(vr+1.0e-20_R8),
     +                     ABS(z32-zre)/ABS(vz+1.0e-20_R8))
          endif
       endif
       if (r20 .eq. r32) then
          d = z32 + vz
          if (((z32 .le. d)  .and.  (z32 .ge. (zre)))  .or.
     +        ((z32 .ge. d)  .and.  (z32 .le. (zre)))) then
              wall=22
              error=.false.
              tdist = MAX( ABS(r20-rre)/ABS(vr+1.0e-20_R8),
     +                     ABS(z20-zre)/ABS(vz+1.0e-20_R8))
          endif
       endif
       if (r01 .eq. r20) then
          d = z20 + vz
          if (((z20 .le. d)  .and.  (z20 .ge. (zre)))  .or.
     +        ((z20 .ge. d)  .and.  (z20 .le. (zre)))) then
              wall=00
              error=.false.
              tdist = MAX( ABS(r20-rre)/ABS(vr+1.0e-20_R8),
     +                     ABS(z20-zre)/ABS(vz+1.0e-20_R8))
          endif
       endif

       if (error) then
         if (b2plot_debug.ge.1)
     .    write(*,*)'Error while looking for a wall to cross'
         laufende = .true.
         tdist=0.0_R8
       endif

C      Through one of the walls or corners it has to go! :)
C
      return
      end subroutine intersection

c------------------------------------------------------------------------------
c     CALCULATE THE INTERSECTION OF TWO LINES
c------------------------------------------------------------------------------

      subroutine schnittpunktzweierlinien(vr,vz,rre,zre,ii,jj,x,y,e1,e2)
C
C     x and y are the coordinates of the intersection point of the two
C     lines that are defined a) by a point (rre,zre) and a direction
C     (vr,vz) and b) by two points (e1,e2), which are two corners of
C     the cell (ii,jj).

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot_debug
      implicit none

      integer e1,e2,ii,jj
      real (kind=R8) ::
     *      vz,vr,rre,zre,x,y
      real (kind=R8) ::
     *      a1,a2,b1,b2
      logical einssenk, zweisenk

C     Check for vertical lines. I do not check for very large slope which
C     could lead to numerical problems (?).

      einssenk=.false.
      zweisenk=.false.

      if (vr .ne. 0.0_R8) then
         a1 = vz/vr
         b1 = -a1*rre+zre
         if (b2plot_debug.ge.2400) write(*,*)'a1: ',a1,'   b1: ',b1
      else
         einssenk=.true.
      endif

      if (crx(ii,jj,e2) .ne. crx(ii,jj,e1)) then
        a2 = (cry(ii,jj,e2)-cry(ii,jj,e1))/(crx(ii,jj,e2)-crx(ii,jj,e1))
        b2 = -a2*crx(ii,jj,e1)+cry(ii,jj,e1)
      else
        zweisenk=.true.
      endif

      if ((.not.einssenk) .and. (.not.zweisenk)) then
         if (a1 .ne. a2) then
            x = (b2-b1)/(a1-a2)
            y = a1*x + b1
         else
            x = 10000000000.0_R8
            y = 10000000000.0_R8
C           There is no intersection point. Therefore I set
C           it to be very large. The upper program will recognise this
C           as being outside the cell.
         endif
      else if (einssenk .and. (.not.zweisenk)) then
         x = rre
         y = x*a2 + b2
      else if (zweisenk .and. (.not.einssenk)) then
         x = crx(ii,jj,e1)
         y = x*a1 + b1
      else
         x = 2000000000.0_R8
         y = 2000000000.0_R8
C        as above
      endif

      return
      end subroutine schnittpunktzweierlinien

c     -----------------------------------------------------------------------------------
c     TRANSFORM VELOCITIES
c     -----------------------------------------------------------------------------------

      subroutine veloctrans(ii,jj,l,vr,vz,clr,clz)

      use b2mod_types
      use b2mod_b2plot
      implicit none

      integer ii,jj,l
      real (kind=R8) ::
     *   vr,vz
      real (kind=R8) ::
     *   clr(-1:nx,-1:ny,1:4),clz(-1:nx,-1:ny,1:4)
#ifdef WG_TODO
      real (kind=R8) ::
     *   rechvv(1:2),rechuu(1:2),betrag

C Projection of the vvtraj upon the axis of the cell (side middle points)

      betrag=sqrt(((clr(ii,jj,2)-clr(ii,jj,4))**2+
     +             (clz(ii,jj,2)-clz(ii,jj,4))**2))
      rechvv(1)=(1/betrag)*(clr(ii,jj,2)-clr(ii,jj,4))*
     +               vvtraj(ii,jj,l)
      rechvv(2)=(1/betrag)*(clz(ii,jj,2)-clz(ii,jj,4))*
     +               vvtraj(ii,jj,l)

C Projection of the uutraj upon the axis of the cell (side middle points)

      betrag=sqrt((clr(ii,jj,1)-clr(ii,jj,3))**2 +
     +                (clz(ii,jj,1)-clz(ii,jj,3))**2)
      rechuu(1)=(1/betrag)*(clr(ii,jj,1)-clr(ii,jj,3))*
     +              uutraj(ii,jj,l)
      rechuu(2)=(1/betrag)*(clz(ii,jj,1)-clz(ii,jj,3))*
     +              uutraj(ii,jj,l)

      vr=rechvv(1)+rechuu(1)
      vz=rechvv(2)+rechuu(2)
#endif

      return
      end subroutine veloctrans

c     --------------------------------------------------------------------------------------
c     CALCULATE THE WEIGHT ARRAY FROM THE WEIGHER
c     --------------------------------------------------------------------------------------

      subroutine weightarray(warr,warr_in,ls_active,ix1,ix2,iy1,iy2,nd)

      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      implicit none

      integer ix,iy,is,nd,ix1,ix2,iy1,iy2
      logical ls_active(0:ns-1)
      real (kind=R8) ::
     *      warr_in(-1:nx,-1:ny,0:nd-1),warr(-1:nx,-1:ny,0:ns-1)

      do is=0,ns-1
        do ix=-1,nx
          do iy=-1,ny
            warr(ix,iy,is)=0.0_R8
          enddo
        enddo
      enddo

      if (nd .eq. ns) then
        do is=0,ns-1
          if(ls_active(is)) then
            do ix=ix1,ix2
              do iy=iy1,iy2
                warr(ix,iy,is)=warr_in(ix,iy,is)
              enddo
            enddo
          endif
        enddo
      else if (nd .eq. 1) then
        do is=0,ns-1
          if(ls_active(is)) then
            if (b2plot_debug.ge.2)
     .       write(*,*) 'Allocating for species ',is,' for ',
     1       ix1,' <= ix <= ',ix2,' and ',iy1,' <= iy <= ',iy2
            if (b2plot_debug.ge.2) write(*,*) 'nx,ny,ns=',nx,ny,ns
            do ix=ix1,ix2
              do iy=iy1,iy2
                warr(ix,iy,is)=warr_in(ix,iy,0)
              enddo
            enddo
          endif
        enddo
      else
        write(*,*)
     1   'Source function does not have the right size, try again.'
        return
      endif

      return
      end subroutine weightarray

c     -------------------------------------------------------------------------
c     INITIALISATION
c     -------------------------------------------------------------------------

c     initialisation of the traj-parameter t_#p (number of plots)
c     (for initialisation of t_pt (plottime)
c     see file b2plot.F, initialisation of values(...) )
      block data
      integer nplot
      common /t_com/ nplot
      data nplot /10/
      logical ratio
      common /t_ratio/ ratio
      data ratio /.true./
      end block data
c
      subroutine allocate_particles(warr,ipx,ipy,ips,p_r,p_z,nps,np)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      implicit none

      integer nps,np,ipx(np),ipy(np),ips(np)
      real (kind=R8) ::
     *   warr(-1:nx,-1:ny,0:ns-1),swarr
      real (kind=R8) ::
     *   p_r(np),p_z(np)
      integer iwarr(-1:nx,-1:ny,0:ns-1)
      integer index_arr(-1:nx,-1:ny,0:ns-1),
     1 index_part(npart),killed(npart)
      integer ix,iy,is,ip,ikilled,ipt,ipto,iptoo,nadd
      real (kind=R8) ::
     *   killval,addval,rmin,rmax,zmin,zmax,rre,zre
#ifndef NAG22
      real (kind=R8) ::
     *   g05caf
#else
      integer IGEN
      integer ISEED(4)
      common /nag22/ igen, iseed
      real (kind=R8) ::
     *   g05kaf
#endif
      logical phase1,inside_rectangle,gefunden
      if (b2plot_debug.ge.8)
     . write(*,*) 'Starting allocate_particles: nps, np = ',nps,np
c initialise integer weight array and cross-index array
      swarr=0.0_R8
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            iwarr(ix,iy,is)=0
            index_arr(ix,iy,is)=0
            swarr=swarr+warr(ix,iy,is)
          enddo
        enddo
      enddo
      if(swarr .eq. 0) then
        if (b2plot_debug.ge.8)
     .   write(*,*) 'Weighting array == 0, reset to 1'
        do is=0,ns-1
          do iy=-1,ny
            do ix=-1,nx
              warr(ix,iy,is)=1
              swarr=swarr+warr(ix,iy,is)
            enddo
          enddo
        enddo
      endif
      if (b2plot_debug.ge.8)
     . write(*,*) 'allocate_particles: swarr = ',swarr
c evaluate integer weight array and cross-index arrays
      ikilled=0
      do ip=1,nps
        ix=ipx(ip)
        iy=ipy(ip)
        is=ips(ip)
        if(ix .lt. 0 .or. ix .gt. nx-1 .or.
     1     iy .lt. 0 .or. iy .gt. ny-1) then
          ikilled=ikilled+1
          killed(ikilled)=ip
        else
          iwarr(ix,iy,is)=iwarr(ix,iy,is)+1
          if(index_arr(ix,iy,is) .eq. 0) then
            index_arr(ix,iy,is)=ip
          else
            ipt=index_arr(ix,iy,is)
            do while(ipt .ne. 0)
              ipto=ipt
              ipt=index_part(ipto)
            enddo
            index_part(ipto)=ip
          endif
          index_part(ip)=0
        endif
      enddo
c find particles that need to be killed
      nadd=0
      killval=-1.0_R8
      addval=0.0_R8
      phase1=.true.
 1000 continue
      if (b2plot_debug.ge.8)
     . write(*,5000) 1000,killval,addval,phase1,ikilled,nadd
      do ip=nps+1,np
        ikilled=ikilled+1
        killed(ikilled)=ip
      enddo
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            do while
     .       (warr(ix,iy,is)/swarr*np-iwarr(ix,iy,is) .le. killval)
              ipt=index_arr(ix,iy,is)
              if(ipt .eq. 0) then
                write(*,*) 'ERROR: should not reach here'
              endif
              ipto=0
              do while(ipt .ne. 0)
                iptoo=ipto
                ipto=ipt
                ipt=index_part(ipto)
              enddo
              if(iptoo .eq. 0) then
                index_arr(ix,iy,is)=0
              else
                index_part(iptoo)=0
              endif
              iwarr(ix,iy,is)=iwarr(ix,iy,is)-1
              ikilled=ikilled+1
              killed(ikilled)=ipto
              if(.not.phase1 .and. ikilled .eq. nadd) goto 2000     ! we have killed enough particles
            enddo
          enddo
        enddo
      enddo
      if(.not.phase1 .and. ikilled .lt. nadd) then                  ! we did not kill enough particles, try again
        killval=killval+0.1_R8
        goto 1000
      endif
c find where particles need to be added
 2000 continue
      if (b2plot_debug.ge.8)
     . write(*,5000) 2000,killval,addval,phase1,ikilled,nadd
      nadd=0
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            nadd=int(warr(ix,iy,is)/swarr*np)
     .             -iwarr(ix,iy,is)+nint(addval)+nadd
            if(nadd .gt. 0 .and. ikilled .gt. 0) then
              ipt=index_arr(ix,iy,is)
              ipto=0
              do while(ipt .ne. 0)
                ipto=ipt
                ipt=index_part(ipto)
              enddo
              do while(nadd .gt. 0 .and. ikilled .gt. 0)
                if(ipto .eq. 0) then
                  index_arr(ix,iy,is)=killed(ikilled)
                else
                  index_part(ipto)=killed(ikilled)
                endif
                index_part(killed(ikilled))=0
c add particle of species is in cell (ix,iy)
                ipx(killed(ikilled))=ix
                ipy(killed(ikilled))=iy
                ips(killed(ikilled))=is
                zmin=
     1           MIN(cry(ix,iy,1),cry(ix,iy,3),
     2               cry(ix,iy,2),cry(ix,iy,0))
                rmin=
     1           MIN(crx(ix,iy,1),crx(ix,iy,3),
     2               crx(ix,iy,2),crx(ix,iy,0))
                zmax=
     1           MAX(cry(ix,iy,1),cry(ix,iy,3),
     2               cry(ix,iy,2),cry(ix,iy,0))
                rmax=
     1           MAX(crx(ix,iy,1),crx(ix,iy,3),
     2               crx(ix,iy,2),crx(ix,iy,0))
                gefunden=.false.
                do while (.not.gefunden)
#ifndef NAG22
                  rre=g05caf(3.1_R8)*(rmax-rmin)+rmin
                  zre=g05caf(3.1_R8)*(zmax-zmin)+zmin
#else
                  rre=g05kaf(IGEN,ISEED)*(rmax-rmin)+rmin
                  zre=g05kaf(IGEN,ISEED)*(zmax-zmin)+zmin
#endif
                  gefunden=inside_rectangle(rre,zre,
     1             crx(ix,iy,1),cry(ix,iy,1),
     2             crx(ix,iy,3),cry(ix,iy,3),
     3             crx(ix,iy,2),cry(ix,iy,2),
     4             crx(ix,iy,0),cry(ix,iy,0))
                enddo
                p_r(killed(ikilled))=rre
                p_z(killed(ikilled))=zre
                ikilled=ikilled-1
                nadd=nadd-1
                iwarr(ix,iy,is)=iwarr(ix,iy,is)+1
              enddo
            endif
            if(.not.phase1 .and. ikilled .eq. 0) goto 4000
          enddo
        enddo
      enddo
c at this point we might still need to add either "nadd" or "ikilled" particles
 3000 continue
      phase1=.false.
      if (b2plot_debug.ge.8)
     . write(*,5000) 3000,killval,addval,phase1,ikilled,nadd
      if(nadd .gt. 0) then
        killval=killval+0.1_R8
        goto 1000
      endif
      if(ikilled .gt. 0) then
        addval=addval+0.1_R8
        goto 2000
      endif
 4000 continue
      if (b2plot_debug.ge.8)
     . write(*,5000) 4000, killval,addval,phase1,ikilled,nadd
      return
 5000 format('allocate_particles (',i4,'): ',
     1 'killval, addval, phase1, ikilled, nadd = ',
     2 2(f6.2,1x),l1,2(1x,i6))
      end subroutine allocate_particles
c
#ifdef NO_NAG
#define ranf g05caf
#define ranset g05cbf

      function ranf (dummy)                       !srv 17.12.07
      use b2mod_types
      implicit none
      real(kind=R8) :: ra, ranf, dummy
      real(kind=R8), external :: local_h1rn
      ra=local_h1rn(dummy)
      ranf=ra
      return
      end function ranf

      subroutine ranset (iseed)
      implicit none
      integer :: iseed, iseed1, iseed2
      iseed1=iseed
      iseed2=iseed+1
      call local_h1rnin(iseed1,iseed2)
      return
      end subroutine ranset

      FUNCTION LOCAL_H1RN(DUMMY)
      use b2mod_types
      use b2mod_b2plot_debug
      implicit none
*
*#**********************************************************************
*# RANDOM NUMBER GENERATOR AS ADVOCATED BY F. JAMES FROM PROPOSAL OF   *
*# MARSAGLIA AND ZAMAN FSU-SCRI-87-50 AND MODIFIED BY F. JAMES 1988 TO *
*# PRODUCE VECTOR OF NUMBERS.                                          *
*# ENTRIES ARE:                                                        *
*#     FUNCTION    H1RN(DUMMY)     SINGLE RANDOM NUMBER                *
*#     SUBROUTINE  H1RNV(VEC,LEN)  VECTOR OF RANDOM NUMBERS            *
*#     SUBROUTINE  H1RNIN(IJ,KL)   INITIALISE WITH SEEDS               *
*#     SUBROUTINE  H1RNIV(VEC)     INITIALISE/RESTART WITH SEED ARRAY  *
*#     SUBROUTINE  H1RNSV(VEC)     SAVE SEED ARRAY VEC(100)            *
*#                                                                     *
*# NOTE: -H1RNIN OR H1RNIV MUST BE CALLED BEFORE GENERATING ANY        *
*#        RANDOM NUMBER(S).                                            *
*#       -H1RNSV SAVES SEED ARRAY INTO VEC(100) ONLY. THE USER HAS TO  *
*#        OUTPUT IT.                                                   *
*#                                                                     *
*# CHANGED BY: G. GRINDHAMMER AT: 90/03/14                             *
*# REASON :                                                            *
*#**********************************************************************
*
c      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(kind=R8) :: dummy,C,CD,CM,U,LOCAL_H1RN
      INTEGER :: I,J,ISEED1,ISEED2
      CHARACTER*16    FLAG,CHECK
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG
*
      LOGICAL FIRST
cdpc split logical and data
      data first /.TRUE./
      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG  .NE.  CHECK) THEN
            IF (B2PLOT_DEBUG .GE. 8)
     .       WRITE(6,*) ' H1RN (RANMAR): INITIALISED WITH DEFAULT SEED'
            ISEED1      = 12345
            ISEED2      = 98765
            CALL LOCAL_H1RNIN(ISEED1,ISEED2)
            FIRST = .FALSE.
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
  100 CONTINUE
      LOCAL_H1RN = U(I)-U(J)
      IF(LOCAL_H1RN  .LT.  0.0_R8) LOCAL_H1RN = LOCAL_H1RN + 1.0_R8
      U(I) = LOCAL_H1RN
      I = I - 1
      IF( I  .EQ.  0) I=97
      J = J - 1
      IF( J  .EQ.  0) J=97
      C = C - CD
      IF( C  .LT.  0) C = C + CM
      LOCAL_H1RN = LOCAL_H1RN-C
      IF(LOCAL_H1RN  .LE.  0.0_R8) LOCAL_H1RN = LOCAL_H1RN + 1.0_R8
      IF(LOCAL_H1RN  .GE.  1.0_R8) GOTO 100
      RETURN
      END FUNCTION LOCAL_H1RN

      SUBROUTINE LOCAL_H1RNIN(IJ,KL)
      use b2mod_types
      use b2mod_b2plot_debug
      implicit none
      INTEGER :: IJ,KL,I,II,J,JJ,IP,JP,K,L,M
      REAL(kind=R8) :: C,CD,CM,T,U,S
      CHARACTER*16    FLAG
      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
      COMMON /RASET2/ FLAG
*
      IJ = IABS(IJ)
      KL = IABS(KL)
      IJ = MOD(IJ,31329)
      KL = MOD(KL,30082)
      I  = MOD(IJ/177, 177) + 2
      J  = MOD(IJ, 177)     + 2
      K  = MOD(KL/169, 178) + 1
      L  = MOD(KL, 169)
*
      DO II= 1, 97
         S= 0.0_R8
         T= 0.5_R8
         DO JJ= 1,24
            M = MOD(MOD(I*J,179)*K, 179)
            I = J
            J = K
            K = M
            L = MOD(53*L+1, 169)
            IF ( MOD(L*M,64)  .GE.  32) S = S + T
            T = 0.5_R8*T
         END DO
         U(II) = S
      END DO
      C  =   362436.0_R8/16777216.0_R8
      CD =  7654321.0_R8/16777216.0_R8
      CM = 16777213.0_R8/16777216.0_R8
      IP = 97
      JP = 33
*
      FLAG = 'H1RN INITIALISED'
      IF (B2PLOT_DEBUG .GE. 8) WRITE(6,610) IJ,KL,I,J,K,L
  610 FORMAT(' H1RNIN: H1RN (RANMAR) INITIALISED WITH 2 SEEDS: ',
     >       2I7,4I4)
*
      RETURN
      END SUBROUTINE LOCAL_H1RNIN
#endif
