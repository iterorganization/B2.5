C     =======================================================================
C     SUBROUTINES FOR B2PLOT COMMAND TRAJ
C     =======================================================================
C
C     This routine is used to plot trajectories of particles using the
C     velocity fields uu and vv. To call this routine
C     one must specify at least 4 parameters:
C     # of particles, 1 or 0 for movie or lines, species indices (at least
C     one), and last but not least the total number of species indices one
C     has given in the step before. A valid call could be:
C     100 0 0 1 traj
C     This should plot the trajectories of 100 particles as lines
C     for hydrogen (0=H).
C     The start locations of the trajectories will be chosen from the set
C     of set defined in the ic_active array (default is all real cells).
C     The location (in a structured grid) of the cells from which to launch
C     particles can be chosen using the "t_xy" command, or from the non-zero
C     cells from the matrix stack with the "icse" command.
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
C     -1000 0 0 1 traj                ! over all ic_active cells
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
      use b2mod_dimensions
      implicit none

c     INPUT
      integer iarg(1:state%ns+3)        !parameters
      integer specnum                   !location of the number of species

c     the traj common
      integer nplot                     !number of plots
      common /t_com/ nplot              !getting the values
      logical ratio
      common /t_ratio/ ratio

      real (kind=R8), allocatable ::
     1   warr(:,:)                      !the stack

      real (kind=R8) ::
     *   globaltime                     !plot time
      integer ialt                      !position in computational domain
      integer ipc(1:npart),ips(1:npart) !number of particles per cell
      integer ipc_start(1:npart)
      real (kind=R8) ::
     *   p_r(1:npart),p_z(1:npart)
      real (kind=R8) ::
     *   p_r_start(1:npart),p_z_start(1:npart)

      real (kind=R8) ::
     *   deltat                         !timestep between plots
      logical ls_active(0:state%ns-1)   !species active?
      integer isclvl(0:state%ns-1),iscol(0:state%ns-1)
      character*19 texttraj(1:state%ns)

      integer ip,randzaehler,zaehler

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

      integer icol,oocol
      character argstring*256,tmpstring*256,titlestring*256
      real (kind=R8) ::
     *   currtime

#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
#     include <netcdf.inc>
      character rw*5
      integer ncid,imap(maxvdims),iret
      integer itstep,ntstep,iCv,iFc
      real (kind=R8), allocatable :: lwarr(:,:)
      real (kind=R8) ::
     1   starttime, tstepn(1), dumarr(1)
      real (kind=R8) ::
     *   flux_density(1:mpg%nFc,0:state%ns-1),
     *   cv_flux(1:mpg%nCv,0:state%ns-1)
#endif

      oocol=-1
      iscol=-1
      allocate(warr(mpg%nCv,0:state%ns-1))
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
      allocate(lwarr(mpg%nCv,0:state%ns-1))
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
      if (iarg(2) .eq. 1) then
         if (b2plot_debug.ge.2) write(*,*) 'MOVIE'
      else
         if (b2plot_debug.ge.2) write(*,*) 'TRAJECTORIES'
      endif
      if(all(ic_active(1:mpg%nCi))) then
        destmode=2
      else
        destmode=0
      endif
      if (b2plot_debug.ge.8) write(*,*) 'Set destmode to ',destmode
      if (b2plot_debug.ge.2)
     + write(*,*) 'plottime in seconds (default=5.0e-4): ', globaltime
      if (b2plot_debug.ge.2)
     + write(*,*) 'number of plots (default=10): ', nplot

c     getting the species numbers
      do is=0,state%ns-1
        ls_active(is)=.false.
      enddo
      do zaehler=0,iarg(specnum)-1
        iscol(iarg(3+zaehler))=zaehler+1
        ls_active(iarg(3+zaehler))=.true.
        texttraj(zaehler+1)=textpl(iarg(3+zaehler))
      enddo

c     find weighting function for the particle positions
      if (iarg(1) .gt. 0) then      !default: weighting with volume
        if (b2plot_debug.ge.1)
     >   write(*,*) 'weighting randomiser with volume'
        call weightarray(warr,geo%cvVol,ls_active,ic_active,1)
      else if (iwrk .gt. 0) then    !user-defined weighting
        if (b2plot_debug.ge.1)
     >   write(*,*) 'weighting randomiser user defined'
        call weightarray(warr,work(1,0,iwrk),ls_active,ic_active,
     .   ndim(2,iwrk))
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
     > tmpstring=trim(tmpstring)//trim(texttraj(1))
      call customise_label(tmpstring,titlestring)
      call init_page(titlestring,values(5),values(6),
     +               values(7),values(8),1,.true.,.true.,.true.)
      call GOPWK(9,1,3)         !open GFLASH

c     generate beginning positions

      np=abs(iarg(1))
      call cpu_time(cpu_start)
#if ( !defined(NO_CDF) && defined(CDF_MOVIE) )
      if(moviename .ne. ' ') then
        write(*,*) 'moviename = ',moviename
        if (b2plot_debug.ge.1) write(*,*) 'Changing weighting'
        do is=0,state%ns-1
          warr(:,is)=warr(:,is)/state%pl%na(:,is)*lwarr(:,is)
        enddo
      endif
#endif
      call allocate_particles(warr,ipc,ips,p_r,p_z,0,np)
      call cpu_time(cpu_end)
      if (b2plot_debug.ge.8) write(*,'(a,f6.3,a)')
     1 'Initial allocate_particles took ',
     2 cpu_end-cpu_start,
     3 ' cpu-seconds'
      do ip=1,np
        ipc_start(ip)=ipc(ip)
        p_r_start(ip)=p_r(ip)
        p_z_start(ip)=p_z(ip)
      enddo

c     initialise graphics

      call GFLAS1(1)
      call gsplci(1)
      if(iarg(specnum) .eq. 1 .and. iarg(2) .eq. 0) then  ! 1 species and trajectories
        call data_copy(iarg(3))
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

      if (iarg(2) .eq. 1) then
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
            do is=0,state%ns-1
              do iFc = 1, mpg%nFc
                if (geo%fcQalf(iFc,1).ne.0.0_R8) then
                  flux_density(iFc,is)=state%dv%fna(iFc,1,is)
     .             /(geo%fcS(iFc)*geo%fcQalf(iFc,1))
                else
                  flux_density(iFc,is)=0.0_R8
                end if
              end do
              call intcell(mpg%nFc,mpg%nCv,mpg,mpg%intcellR,
     .         flux_density(1,is),cv_flux(1,is))
              do iCv=1,mpg%nCv
                uutraj(iCv,is)=state%pl%ua(iCv,is)
                vvtraj(iCv,is)=cv_flux(iCv,is)/state%pl%na(iCv,is)
              end do
            end do
          endif
#endif
          if (b2plot_debug.ge.2)
     .     write(*,*)'plotnumber: ',counter,'(of ',nplot,')'
          call end_page
          call lastd
          call gsplci(1)
          write(argstring,'(f10.3,a)') currtime*1.0e6_R8,' \Gm\Rs'
          call init_page('Trajectories '//trim(argstring),
     1     values(5),values(6),values(7),values(8),1,
     2     .true.,.true.,.true.)
          call gflas3(1)
          call cpu_time(cpu_start)
          do ip=1,np
            call lastd
            call gsplci(ips(ip)+1)
            ialt=ipc(ip)
            if(ialt.lt.1.or.ialt.gt.mpg%nCi) then
              write(*,*) 'Problem with particle ',ip,': ialt = ',ialt
            else
              call trajectories(mpg,geo,
     1         p_r(ip),p_z(ip),ipc(ip),ips(ip),deltat)
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
              do is=0,state%ns-1
                warr(:,is)=warr(:,is)/state%pl%na(:,is)*lwarr(:,is)
              enddo
            endif
#endif
            call cpu_time(cpu_start)
            call allocate_particles(warr,ipc,ips,p_r,p_z,np,np)
            call cpu_time(cpu_end)
            if (b2plot_debug.ge.8) write(*,'(a,f6.3,a)')
     1       'Subsequent allocate_particles took ',
     2       cpu_end-cpu_start,
     3       ' cpu-seconds'
          else
            do ip=1,np
              if(ipc(ip).lt.1 .or. ipc(ip).gt.mpg%nCi) then
                ipc(ip)=ipc_start(ip)
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
            do is=0,state%ns-1
              do iFc = 1, mpg%nFc
                if (geo%fcQalf(iFc,1).ne.0.0_R8) then
                  flux_density(iFc,is)=state%dv%fna(iFc,1,is)
     .             /(geo%fcS(iFc)*geo%fcQalf(iFc,1))
                else
                  flux_density(iFc,is)=0.0_R8
                end if
              end do
              call intcell(mpg%nFc,mpg%nCv,mpg,mpg%intcellR,
     .         flux_density(1,is),cv_flux(1,is))
              do iCv=1,mpg%nCv
                uutraj(iCv,is)=state%pl%ua(iCv,is)
                vvtraj(iCv,is)=cv_flux(iCv,is)/state%pl%na(iCv,is)
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
            if (ipc(ip).gt.mpg%nCi .or. ipc(ip).lt.1) then
              randzaehler=randzaehler+1
            else
              rrold=p_r(ip)
              zrold=p_z(ip)
              call trajectories(mpg,geo,
     1         p_r(ip),p_z(ip),ipc(ip),ips(ip),deltat)
              if(iarg(specnum).ne.1) then
                icol=iscol(ips(ip))
                if(icol.ne.oocol .and. icol.ne.-1 .and. color_set) then
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
      subroutine trajectories(mpg,geo,rre,zre,iCv,l,deltat)

      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo

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
     *   rre,zre,radvance,zadvance
      real (kind=R8) ::
     *   tdist,udist,rdist
      integer i,iCv,iCn,iFc,iVx,l,iCold
      integer iwall,jwall,schritte
      logical laufende,laufende2
      logical inside_cell

C     Initial configuration. Particle time is counted downwards,
C     therefore is set to be deltaT here.

      vr = 0.0_R8
      vz = 0.0_R8
      laufende = .false.
      tpart = deltat
      schritte = 0

      do while(.not.laufende)
         schritte = schritte+1

C     The velocities need to be transformed from the computational
C     to the real domain.

         call veloctrans(geo,iCv,l,vr,vz)

C      OK, first calculate the time the particle needs to hit a "wall"
C      That is, at first one has to know, which wall it is going to hit!
C      :)

         call intersection(mpg,geo,
     1    iwall,tdist,vz,vr,rre,zre,iCv,laufende)
         if (laufende .and. b2plot_debug.ge.12) then
            write(*,'(i5,4(1e15.8,2x))')
     1       iCv,uutraj(iCv,l),vvtraj(iCv,l),vr,vz
         endif

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
            iCold = iCv
            if (iwall .gt. 0) then ! we cross a face
               iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+iwall-1)
               if (iCold.eq.mpg%fcCv(iFc,1)) then
                 iCv = mpg%fcCv(iFc,2)
               else if (iCold.eq.mpg%fcCv(iFc,2)) then
                 iCv = mpg%fcCv(iFc,1)
               end if
            else if (iwall .lt. 0) then ! we cross at a vertex with index abs(iwall)
               iVx = mpg%cvVx(mpg%cvVxP(iCv,1)+abs(iwall)-1)
CXPB To find which Cv is on the other side of the vertex
CXPB we look for the nearest second intersection
               iCv = 0
               udist = huge(1.0_R8)
               do i = 1, mpg%vxCvP(iVx,2)
                 iCn = mpg%vxCv(mpg%vxCvP(iVx,1)+i-1)
                 if (iCn.eq.iCold) cycle
                 call intersection(mpg,geo,
     1            jwall,rdist,vz,vr,rre,zre,iCv,laufende2)
                 if (jwall.ne.0) then
                   if (rdist.lt.udist) then
                     udist = rdist
                     iCv = iCn
                   end if
                 end if
               end do
            endif
         endif

C        calculate new particle coordinates, advancing them a little
C        bit into the next cell.
CXPB     deltar and deltaz are normally zero, except when crossing periodicity
CXPB     boundaries, then they translate the particle to the spatial
CXPB     location of the next computational cell.

         if (iCv.gt.mpg%nCi) then
           rre = rre + vr*timerun
           zre = zre + vz*timerun
           goto 1000
         end if

         radvance=vr*MAX(timerun,1.0e-6_R8)*0.01_R8
         zadvance=vz*MAX(timerun,1.0e-6_R8)*0.01_R8

         rre = rre + vr*timerun + radvance
         zre = zre + vz*timerun + zadvance

         do while(.not.inside_cell(mpg,geo,rre,zre,iCv))
           radvance=radvance*0.5_R8
           zadvance=zadvance*0.5_R8
           if (Max(Abs(radvance),Abs(zadvance)) .lt. 1.0e-20_R8) then
             if (b2plot_debug.ge.1)
     .        write(*,*) 'Particle had to be killed.'
             if (b2plot_debug.ge.1) write(*,*) 'No Problem. (?)'
             iCv=0
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
           iCv=0
         endif

C     Check if particle reached the border.

         if(iCv.le.0.or.iCv.gt.mpg%nCi) then
           laufende = .true.
         endif
      enddo
      return
      end subroutine trajectories

c     -------------------------------------------------------------------------
c     CALCULATING THE INTERSECTION WITH THE CELL WALL
c     -------------------------------------------------------------------------

      subroutine intersection(mpg,geo,
     * iwall,tdist,vz,vr,rre,zre,iCv,laufende)

      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_b2plot_debug
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo

      real (kind=R8) ::
     *   r,z,rre,zre
      real (kind=R8) ::
     *   d,tdist
      real (kind=R8) ::
     *   vz,vr
      integer i,j,iwall,iCv,iFc,iVx
      logical laufende,error

      error=.true.
      iwall=0

C     Calculate the slope, that is the a in y=ax+b

      do i=1,mpg%cvFcP(iCv,2)
        iFc=mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
        call schnittpunktzweierlinien(mpg,geo,vr,vz,rre,zre,iFc,r,z)

C     Now we know the intersection between the path of the particle
C     and the lines that make up the walls. But we do not know if they
C     are inside the walls of the cells or on the outer part of the
C     line making up the cell.

        if((r.le.max(geo%vxX(mpg%fcVx(iFc,1)),geo%vxX(mpg%fcVx(iFc,2))))
     .      .and.
     .     (r.ge.min(geo%vxX(mpg%fcVx(iFc,1)),geo%vxX(mpg%fcVx(iFc,2))))
     .      .and.
     .     (z.le.max(geo%vxY(mpg%fcVx(iFc,1)),geo%vxY(mpg%fcVx(iFc,2))))
     .      .and.
     .     (z.ge.min(geo%vxY(mpg%fcVx(iFc,1)),geo%vxY(mpg%fcVx(iFc,2))))
     .      ) then
          d = r + vr
          if (((r .le. d)  .and.  (r .ge. (rre )))  .or.
     +        ((r .ge. d)  .and.  (r .le. (rre )))) then
            if (iwall.eq.0) then
              iwall=i
              error=.false.
              tdist = MAX( ABS(r-rre)/ABS(vr+1.0e-20_R8),
     +                     ABS(z-zre)/ABS(vz+1.0e-20_R8))
            else ! we are hitting two walls, i.e. their common vertex
              do j=1,mpg%cvVxP(iCv,2)
                iVx=mpg%cvVx(mpg%cvVxP(iCv,1)+j-1)
                if (r.eq.geo%vxX(iVx).and.z.eq.geo%vxY(iVx)) iwall=-j
              enddo
            endif
          endif
        endif
      enddo

      if (error) then
        if (b2plot_debug.ge.1)
     .   write(*,*) 'Error while looking for a wall to cross'
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

      subroutine schnittpunktzweierlinien(mpg,geo,vr,vz,rre,zre,iFc,x,y)
C
C     x and y are the coordinates of the intersection point of the two
C     lines that are defined a) by a point (rre,zre) and a direction
C     (vr,vz) and b) by a cell face (iFc).

      use b2mod_types
      use b2mod_b2plot_debug
      use b2us_geo
      use b2us_map
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo

      integer iFc
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

      if (geo%vxX(mpg%fcVx(iFc,1)) .ne. geo%vxX(mpg%fcVx(iFc,2))) then
        a2 = (geo%vxY(mpg%fcVx(iFc,2))-geo%vxY(mpg%fcVx(iFc,1)))/
     .       (geo%vxX(mpg%fcVx(iFc,2))-geo%vxX(mpg%fcVx(iFc,1)))
        b2 = -a2*geo%vxX(mpg%fcVx(iFc,1))+geo%vxY(mpg%fcVx(iFc,1))
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
         x = geo%vxX(mpg%fcVx(iFc,1))
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

      subroutine veloctrans(geo,iCv,l,vr,vz)

      use b2mod_types
      use b2mod_b2plot
      use b2us_geo
      implicit none
      type(geometry), intent(in) :: geo

      integer iCv,l
      real (kind=R8) ::
     *   vr,vz
      real (kind=R8) ::
     *   rechvv(1:2),rechuu(1:2),betrag

      betrag=sqrt(geo%cvEb(iCv,0)**2+geo%cvEb(iCv,1)**2)

C Projection of the vvtraj upon the axis of the cell (side middle points)

      rechvv(1)=(signmf/betrag)*  geo%cvEb(iCv,1) *vvtraj(iCv,l)
      rechvv(2)=(signmf/betrag)*(-geo%cvEb(iCv,0))*vvtraj(iCv,l)

C Projection of the uutraj upon the axis of the cell (side middle points)

      rechuu(1)=(signmf/betrag)*geo%cvEb(iCv,0)*uutraj(iCv,l)
      rechuu(2)=(signmf/betrag)*geo%cvEb(iCv,1)*uutraj(iCv,l)

      vr=rechvv(1)+rechuu(1)
      vz=rechvv(2)+rechuu(2)

      return
      end subroutine veloctrans

c     -------------------------------------------------------------------
c     CALCULATE THE WEIGHT ARRAY FROM THE WEIGHER
c     -------------------------------------------------------------------

      subroutine weightarray(warr,warr_in,ls_active,ic_active,nd)

      use b2mod_types
      use b2mod_b2plot_debug
      use b2us_data
      implicit none

      integer iCv,is,nd
      logical ls_active(0:state%ns-1),ic_active(1:mpg%nCv)
      real (kind=R8) ::
     * warr_in(1:mpg%nCv,0:nd-1),warr(1:mpg%nCv,0:state%ns-1)

      do is=0,state%ns-1
        do iCv=1,mpg%nCv
          warr(iCv,is)=0.0_R8
        enddo
      enddo

      if (nd .eq. state%ns) then
        do is=0,state%ns-1
          if(ls_active(is)) then
            do iCv=1,mpg%nCv
              if (ic_active(iCv)) warr(iCv,is)=warr_in(iCv,is)
            enddo
          endif
        enddo
      else if (nd .eq. 1) then
        do is=0,state%ns-1
          if(ls_active(is)) then
            if (b2plot_debug.ge.2)
     .       write(*,*) 'Allocating for species ',is
            if (b2plot_debug.ge.2) write(*,*) 'nCv,ns=',mpg%nCv,state%ns
            do iCv=1,mpg%nCv
              if (ic_active(iCv)) warr(iCv,is)=warr_in(iCv,0)
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
      subroutine allocate_particles(warr,ipc,ips,p_r,p_z,nps,np)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2us_data
      implicit none

      integer nps,np,ipc(np),ips(np)
      real (kind=R8) ::
     *   warr(1:mpg%nCv,0:state%ns-1),swarr
      real (kind=R8) ::
     *   p_r(np),p_z(np)
      integer iwarr(1:mpg%nCv,0:state%ns-1)
      integer index_arr(1:mpg%nCv,0:state%ns-1),
     1 index_part(npart),killed(npart)
      integer iCv,is,ip,ikilled,ipt,ipto,iptoo,nadd
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
      logical phase1,inside_cell,gefunden

      if (b2plot_debug.ge.8)
     . write(*,*) 'Starting allocate_particles: nps, np = ',nps,np
c initialise integer weight array and cross-index array
      swarr=0.0_R8
      do is=0,state%ns-1
        do iCv=1,mpg%nCv
          iwarr(iCv,is)=0
          index_arr(iCv,is)=0
          swarr=swarr+warr(iCv,is)
        enddo
      enddo
      if(swarr .eq. 0) then
        if (b2plot_debug.ge.8)
     .   write(*,*) 'Weighting array == 0, reset to 1'
        do is=0,state%ns-1
          do iCv=1,mpg%nCi
            warr(iCv,is)=1
            swarr=swarr+warr(iCv,is)
          enddo
        enddo
      endif
      if (b2plot_debug.ge.8)
     . write(*,*) 'allocate_particles: swarr = ',swarr
c evaluate integer weight array and cross-index arrays
      ikilled=0
      do ip=1,nps
        iCv=ipc(ip)
        is=ips(ip)
        if(iCv.lt.1 .or. iCv.gt.mpg%nCi) then
          ikilled=ikilled+1
          killed(ikilled)=ip
        else
          iwarr(iCv,is)=iwarr(iCv,is)+1
          if(index_arr(iCv,is) .eq. 0) then
            index_arr(iCv,is)=ip
          else
            ipt=index_arr(iCv,is)
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
      do is=0,state%ns-1
        do iCv=1,mpg%nCv
          do while
     .     (warr(iCv,is)/swarr*np-iwarr(iCv,is) .le. killval)
            ipt=index_arr(iCv,is)
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
              index_arr(iCv,is)=0
            else
              index_part(iptoo)=0
            endif
            iwarr(iCv,is)=iwarr(iCv,is)-1
            ikilled=ikilled+1
            killed(ikilled)=ipto
            if(.not.phase1 .and. ikilled .eq. nadd) goto 2000     ! we have killed enough particles
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
      do is=0,state%ns-1
        do iCv=1,mpg%nCi
          nadd=int(warr(iCv,is)/swarr*np)
     .           -iwarr(iCv,is)+nint(addval)+nadd
          if(nadd .gt. 0 .and. ikilled .gt. 0) then
            ipt=index_arr(iCv,is)
            ipto=0
            do while(ipt .ne. 0)
              ipto=ipt
              ipt=index_part(ipto)
            enddo
            do while(nadd .gt. 0 .and. ikilled .gt. 0)
              if(ipto .eq. 0) then
                index_arr(iCv,is)=killed(ikilled)
              else
                index_part(ipto)=killed(ikilled)
              endif
              index_part(killed(ikilled))=0
c add particle of species is in cell (ix,iy)
              ipc(killed(ikilled))=iCv
              ips(killed(ikilled))=is
              zmin=minval(
     1         geo%vxY(mpg%cvVx(mpg%cvVxP(iCv,1):
     2                          mpg%cvVxP(iCv,1)+mpg%cvVxP(iCv,2)-1)))
              rmin=minval(
     1         geo%vxX(mpg%cvVx(mpg%cvVxP(iCv,1):
     2                          mpg%cvVxP(iCv,1)+mpg%cvVxP(iCv,2)-1)))
              zmax=maxval(
     1         geo%vxY(mpg%cvVx(mpg%cvVxP(iCv,1):
     2                          mpg%cvVxP(iCv,1)+mpg%cvVxP(iCv,2)-1)))
              rmax=maxval(
     1         geo%vxX(mpg%cvVx(mpg%cvVxP(iCv,1):
     2                          mpg%cvVxP(iCv,1)+mpg%cvVxP(iCv,2)-1)))
              gefunden=.false.
              do while (.not.gefunden)
#ifndef NAG22
                rre=g05caf(3.1_R8)*(rmax-rmin)+rmin
                zre=g05caf(3.1_R8)*(zmax-zmin)+zmin
#else
                rre=g05kaf(IGEN,ISEED)*(rmax-rmin)+rmin
                zre=g05kaf(IGEN,ISEED)*(zmax-zmin)+zmin
#endif
                gefunden=inside_cell(mpg,geo,rre,zre,iCv)
              enddo
              p_r(killed(ikilled))=rre
              p_z(killed(ikilled))=zre
              ikilled=ikilled-1
              nadd=nadd-1
              iwarr(iCv,is)=iwarr(iCv,is)+1
            enddo
          endif
          if(.not.phase1 .and. ikilled .eq. 0) goto 4000
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
