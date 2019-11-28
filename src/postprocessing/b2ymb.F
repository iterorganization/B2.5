! ##########################################################################
! This program first reads the formatted file 'b2ymb.dat' and writes the
! unformatted file 'mnb_in.dat' that is needed by the FreeFem++ heat
! transfer model of the monoblock.
! Then, it reads the files 'wlld_trg[i-o].dat' produced by
! b2plot with the "1000 wlld" command and writes the unformatted files
! 'Wtot_tg[i-o].dat', 'Wpls_tg[i-o].dat', 'X_tg[i-o].dat' and
! 'tg[i-o]_size.dat' that are needed by the FreeFem++ heat transfer model
!  of the monoblock.
! In 'Wtot_tg[i-o].dat' is stored the value of the total heat flux for
! each tile of the monoblock.
! In 'Wpls_tg[i-o].dat' is stored the value of the heat flux due only to
! plasma particles for each tile of the monoblock.
! In 'X_tg[i-o].dat' is stored the value of the poloidal coordinate of the
! center of each tile along the monoblock axis (in [m]) (first column), the
! r-coordinate (m) of the center of the tile (second column) and the
! z-coordinate (m) of the center of the tile (third column), interpolated
! using the data in wlly_trg[i-o].dat.
! In 'tg[i-o]_size.dat' is stored the number of points (i.e. tiles) in which
! the inner/outer monoblocks have been subdivided according to the assigned
! tile poloidal length, which is read from the 'b2ymb.dat' file.
! ##########################################################################

      PROGRAM b2ymb

      IMPLICIT NONE

      INTEGER :: NN, NLINE, NTILES, NPOINTS, II, JJ, KK, INUNIT=8888
      real :: HEIGHT, T_LENGTH, CU_TH, CUCRZR_TH, DIAMETER, W_TH,
     . P_LENGTH, ALPHA, BETA, THETA, T_GAP
      real :: T_COOLANT, M_FLOW, CP_COOLANT, dummy
      INTEGER :: DIMENSIONS(4), REVERSE, nskip, innout, nplots
      real :: target_length, LPOINTS(50)=0.0d0
      real, ALLOCATABLE :: MM(:,:),BB(:,:), rr(:),zz(:),rout(:),zout(:),
     . dd(:)
      CHARACTER*(10) :: INFILE0= 'b2ymb.dat',
     . OUTFILE0 = 'mnb_in.dat'
      CHARACTER*(20) :: OUTFILE(4),INFILE
      CHARACTER*(1) :: trg(2), ss
      data trg /'i','o'/
      LOGICAL :: exists
      logical streql
      external streql

!##########################################################################
!Read input file 'b2ymb.dat' and then write 'mnb_in.dat'
!##########################################################################
      NAMELIST /GEOMETRY/
     . HEIGHT, T_LENGTH, CU_TH, CUCRZR_TH,
     . DIAMETER, W_TH, P_LENGTH, ALPHA, BETA,
     . THETA, T_GAP, REVERSE

      NAMELIST /MESH/
     . NN

      NAMELIST /COOLANT/
     . T_COOLANT, M_FLOW, CP_COOLANT

      NAMELIST /OUTPUT/
     . NPOINTS, LPOINTS, NPLOTS
      print*,"############## README ############"
      print*,"Starting b2ymb. Update the input file b2ymb.dat 
     .with the desired inlet temperature and then choose 
     .between inner [i] our outer [o] target"
      read (*,*) ss
      if (streql(trim(ss),trg(1))) then
         KK=1
         innout=1
         print *,'inner target chosen'
      else
         KK=2
         innout=0
         print *,'outer target chosen'
      endif

      open(unit=inunit,file=infile0,status='old',action='read')
      READ(INUNIT,NML=GEOMETRY)
      READ(INUNIT,NML=MESH)
      READ(INUNIT,NML=COOLANT)
      READ(INUNIT,NML=OUTPUT)
      
      close(INUNIT);
      ! Write output file
      INQUIRE (FILE=outfile0, EXIST=exists)
      if (exists) then
         open (unit=inunit,file=outfile0,action='readwrite')
         close (unit=inunit,status='delete')
      endif
      open(unit=inunit,file=outfile0,status='new',action='write')
      write(inunit,*) HEIGHT, T_LENGTH, CU_TH, CUCRZR_TH, DIAMETER,
     . W_TH, P_LENGTH, T_GAP, ALPHA, BETA, THETA
      write(inunit,*) T_COOLANT, M_FLOW, CP_COOLANT
      write(inunit,*) NN, NPOINTS,innout,nplots, LPOINTS(1:NPOINTS)
      close(INUNIT);

! ##########################################################################
! Read input file 'wlld_trg[i-o].dat' to get heat flux profile on targets and
! then write the relative output files 'Wtot_tg[i-o].dat','Wpls_tg[i-o].dat'
! 'X_tg[i-o].dat' and 'tg[i-o]_size.dat'. Read also related wlly_*.dat to get
! the r-z coordinates.
! ##########################################################################
      
      INFILE='wlld_trg'//trg(KK)//'.dat'
      open(unit=inunit,file=infile,status='old',action='read')
      NLINE=0
      NSKIP=0
 5    READ(INUNIT,*,END=15) ss !File is of unknown length, so first count the number of lines
         if (streql(ss,'#')) then 
           NSKIP=NSKIP+1 ! header line, skip this one
         else
           NLINE = NLINE +1
         end if
      GO TO 5
 15   CLOSE(INUNIT)
      ALLOCATE(MM(NLINE,4))     ! then allocate variable size and read again storing the values. First column will be x coordinate,
                             ! second the total  heat flux, third the plasma heat flux and fourth the number of point
      MM=0.0d0
      MM(1,4)=dble(NLINE)
      open(unit=inunit,file=infile,status='old',action='read')
      do ii=1,nskip
         READ(INUNIT,*) ss
      end do
      DO II=1,NLINE
         READ(INUNIT,*) MM(II,1:2), dummy, dummy, MM(II,3) !The data in the file are x, Wtot, Wrad, Wpart, Wpls...
      ENDDO
      CLOSE(INUNIT)
      
      ! read wlly to get coordinates
      NSKIP=0
      INFILE='wlly_trg'//trg(KK)//'.dat'
      open(unit=inunit,file=infile,status='old',action='read')
 8    READ(INUNIT,*,END=18) ss 
         if (streql(ss,'#')) then 
           NSKIP=NSKIP+1 ! header line, skip this one
         end if
      GO TO 8
 18   CLOSE(INUNIT)     
      allocate(rr(nline+1),zz(nline+1),dd(nline+1))
      open(unit=inunit,file=infile,status='old',action='read')
      do ii=1,nskip
         READ(INUNIT,*) ss
      end do
      DO II=1,NLINE+1
         READ(INUNIT,*) rr(ii), zz(ii) 
      ENDDO
      CLOSE(INUNIT)
      dd(1)=0.0
      do ii=2,nline+1
        dd(ii)=dd(ii-1)+sqrt((rr(ii)-rr(ii-1))**2+(zz(ii)-zz(ii-1))**2)
      end do
      
      !Interpolate the data on a coordinate vector with spatial x coordinates with a distance equal to the tile poloidal length
      target_length=mm(nline,1)-mm(1,1)
      ntiles=nint(target_length/p_length)
      allocate(bb(ntiles,5),rout(ntiles),zout(ntiles))
      bb=0.0d0
      bb(1,1)=mm(1,1)
      bb(1,2)=mm(1,2)
      bb(1,3)=mm(1,3)
      bb(ntiles,1)=mm(nline,1)
      bb(ntiles,2)=mm(nline,2)
      bb(ntiles,3)=mm(nline,3)
      bb(1,4)=ntiles
      do ii=2,ntiles-1  ! create the new x vector
         bb(ii,1)=bb(ii-1,1)+p_length
      enddo
      call interp( dd(1:nline+1), rr(1:nline+1), bb(1:ntiles,1),
     .    rout(1:ntiles), nline+1, ntiles)
      call interp( dd(1:nline+1), zz(1:nline+1), bb(1:ntiles,1),
     .    zout(1:ntiles), nline+1, ntiles)
      call interp( mm(1:nline,1), mm(1:nline,2), bb(2:ntiles-1,1),
     .    bb(2:ntiles-1,2), nline, ntiles-2)
      call interp( mm(1:nline,1), mm(1:nline,3), bb(2:ntiles-1,1),
     .    bb(2:ntiles-1,3), nline, ntiles-2)
      dimensions(1:3)=ntiles
      dimensions(4)=1
      if (REVERSE>0) then
         do ii=1,2
            bb(1:ntiles,5)=bb(1:ntiles,ii+1) !flip the vector in the direction of the fluid flow
            do jj=1,ntiles
               bb(jj,ii+1)=bb(ntiles-jj+1,5)
            enddo
         enddo
         bb(1:ntiles,5)=rout(1:ntiles)
         do jj=1,ntiles
           rout(jj)=bb(ntiles-jj+1,5)
         enddo
         bb(1:ntiles,5)=zout(1:ntiles)
         do jj=1,ntiles
           zout(jj)=bb(ntiles-jj+1,5)
         enddo
      endif

      ! now write the output files, but first check if there are already present. If yes, remove them before writing
      OUTFILE(1)='X_tg'//trg(KK)//'.dat'
      OUTFILE(2)='Wtot_tg'//trg(KK)//'.dat'
      OUTFILE(3)='Wpls_tg'//trg(KK)//'.dat'
      OUTFILE(4)='tg'//trg(KK)//'_size.dat'
      
      ! do x_tg file separately to put the x-z coordinates
      INQUIRE (FILE=outfile(1), EXIST=exists)
      if (exists) then
        open (unit=inunit,file=outfile(1),action='readwrite')
        close (unit=inunit,status='delete')
      endif
      open(unit=inunit,file=outfile(1),status='new',
     .    action='write')
      do ii=1,DIMENSIONS(1)
        write(inunit,*) bb(ii,1),rout(ii),zout(ii)
      end do
      close(INUNIT);
      DO JJ=2,4
         INQUIRE (FILE=outfile(JJ), EXIST=exists)
         if (exists) then
            open (unit=inunit,file=outfile(JJ),action='readwrite')
            close (unit=inunit,status='delete')
         endif
         open(unit=inunit,file=outfile(JJ),status='new',
     .    action='write')
         do ii=1,DIMENSIONS(JJ)
            write(inunit,*) bb(ii,jj)
         end do
         close(INUNIT);
      ENDDO
      DEALLOCATE(MM,BB,rr,zz,dd,rout,zout)
      stop
      end

!##########################################################################
      subroutine interp( xData, yData, xVal, yVal, nData, nVal )
! Interpolation subroutine. Works only for sorted x values.
! Inputs: xData = a vector of the x-values of the data to be interpolated
! yData = a vector of the y-values of the data to be interpolated
! xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values
!##########################################################################
      implicit none

      real, intent(in) :: xData(nData), yData(nData), xVal(nVal)
      real, intent(out) :: yVal(nVal)
      integer :: inputIndex, dataIndex
      integer, intent(in) :: nData, nVal
      real :: minXdata, maxXData, weight

      minXData = xData(1)
      maxXData = xData(nData)
      dataIndex=1


      do inputIndex = 1, nVal

 1       if(xData(dataIndex+1).le.xVal(inputIndex)) then
            dataIndex=dataIndex+1
            goto 1
          endif

         weight = (xVal(inputIndex) - xData(dataIndex))/
     .    (xData(dataIndex+1)-xData(dataIndex))

         yVal(inputIndex) = weight*(yData(dataIndex+1)-yData(dataIndex))
     .    +yData(dataIndex)
      enddo

      end subroutine
