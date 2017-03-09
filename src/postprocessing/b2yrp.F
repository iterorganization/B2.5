! ##########################################################################
! This program first reads the formatted file 'b2yrp.dat' and
! writes the unformatted file 'rps_in.dat' that is needed by the FreeFem++
! heat transfer model of the reflector plates.
! Then, it reads the files 'wlld_[i-o]rp.dat' produced by
! b2plot with the "1000 wlld" command and writes the unformatted files
! 'Wtot_[i-o]rp.dat','Wpls_[i-o]rp.dat', 'X_[i-o]rp.dat' and
! '[i-o]rp_size.dat' that are needed  by the FreeFem++ heat transfer model
! of the reflector plates.
! In 'Wtot_[i-o]rp.dat' is stored the value of the total heat flux for each
! tile of the RP.
! In 'Wpls_tg[i-o].dat' is stored the value of the heat flux due only to
! plasma particles for each tile of the hypervapotron.
! In 'X_[i-o]rp.dat' is stored the value of the poloidal coordinate of the
! center of each tile along the RP axis (in [m]).
! In '[i-o]rp_size.dat' is stored the number of points (i.e. tiles) in which
! the inner/outer RPs have been subdivided according to the assigned
! tile poloidal length, which is read from the 'b2yrp.dat' file.
! ##########################################################################

      PROGRAM b2yrp

      IMPLICIT NONE

      INTEGER :: NN, NLINE, NTILES, NPOINTS, II, JJ, KK, INUNIT=8888
      real :: L1, L2, L3, H1, H2, H3 ,H4, H5, H_HOLE, RH,
     . P_LENGTH, ALPHA, BETA, THETA, T_GAP, dummy
      real :: T_COOLANT, M_FLOW, CP_COOLANT
      INTEGER :: DIMENSIONS(4), REVERSE, nskip=18, innout, nplots
      real :: target_length, LPOINTS(50)=0.0d0
      real, ALLOCATABLE :: MM(:,:),BB(:,:)
      CHARACTER*(10) :: INFILE0= 'b2yrp.dat',
     . OUTFILE0 = 'rps_in.dat'
      CHARACTER*(20) :: skip, OUTFILE(4),INFILE
      CHARACTER*(1) :: trg(2),ss
      data trg /'i','o'/
      LOGICAL :: exists
      logical streql
      external streql

!##########################################################################
! Read input file 'b2yrp.dat' and then write 'rps_in.dat'
!##########################################################################
      NAMELIST /GEOMETRY/
     . L1, L2, L3, H1, H2, H3 ,H4,
     . H5, H_HOLE, RH, P_LENGTH, ALPHA, BETA,
     . THETA, T_GAP, REVERSE

      NAMELIST /MESH/
     . NN

      NAMELIST /COOLANT/
     . T_COOLANT, M_FLOW, CP_COOLANT

      NAMELIST /OUTPUT/
     . NPOINTS, LPOINTS, NPLOTS
      
      print*,"############## README ############"
      print*,"Starting b2yrp. Update the input file b2yrp.dat 
     .with the desired inlet temperature and then choose 
     .between inner [i] our outer [o] reflector plate"
      read (*,*) ss
      if (streql(trim(ss),trg(1))) then
         KK=1
         innout=1
         print *,'inner reflector plate chosen'
      else
         KK=2
         innout=0
         print *,'outer reflector plate chosen'
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
      write(inunit,*) L1, L2, L3, H1, H2, H3 ,H4,
     . H5, H_HOLE, RH, P_LENGTH, T_GAP, ALPHA, BETA, THETA
      write(inunit,*) T_COOLANT, M_FLOW, CP_COOLANT
      write(inunit,*) NN, NPOINTS,innout,nplots, LPOINTS(1:NPOINTS)
      close(INUNIT);

! ##########################################################################
! Read input files 'wlld_[i-o]rp.dat' to get heat flux profile on RPs and
! then write the relative output files 'Wtot_[i-o]rp.dat', 'Wpls_[i-o]rp.dat'
! 'X_[i-o]rp.dat', and '[i-o]rp_size.dat'.
! ##########################################################################

      INFILE='wlld_'//trg(KK)//'rp.dat'
      open(unit=inunit,file=infile,status='old',action='read')
      DO II=1,NSKIP
         READ(INUNIT,*) skip      !First nskip lines are text lines so must be skipped
      end do
      NLINE=0
 5    READ(INUNIT,*,END=15) skip !File is of unknown length, so first count the number of lines
         NLINE = NLINE +1
      GO TO 5
 15   CLOSE(INUNIT)

      ALLOCATE(MM(NLINE,4))     ! then allocate variable size and read again storing the values. First column will be x coordinate,
                             ! second the total  heat flux, the third the radiation heat flux and the fourth  the number of points
      MM=0.0d0
      MM(1,4)=dble(NLINE)
      open(unit=inunit,file=infile,status='old',action='read')
      DO II=1,NSKIP
         READ(INUNIT,*) skip      !First nskip lines are text lines so must be skipped
      end do
      DO II=1,NLINE
          READ(INUNIT,*) MM(II,1:2), dummy, dummy, MM(II,3) !The data in the file are x, Wtot, Wrad, Wpart, Wpls...
      ENDDO
      CLOSE(INUNIT)
      !Interpolate the data on a coordinate vector with spatial x coordinates with a distance equal to the tile poloidal length
      target_length=mm(nline,1)-mm(1,1)
      ntiles=nint(target_length/p_length)
      allocate(bb(ntiles,5))
      bb=0.0d0
      bb(1,1)=mm(1,1)
      bb(1,2)=mm(1,2)
      bb(1,3)=mm(1,3)
      bb(ntiles,2)=mm(nline,2)
      bb(ntiles,3)=mm(nline,3)
      bb(1,4)=ntiles
      do ii=2,ntiles            ! create the new x vector
         bb(ii,1)=bb(ii-1,1)+p_length
      enddo
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
      endif
      ! now write the output files, but first check if there are already present. If yes, remove them before writing
      OUTFILE(1)='X_'//trg(KK)//'rp.dat'
      OUTFILE(2)='Wtot_'//trg(KK)//'rp.dat'
      OUTFILE(3)='Wpls_'//trg(KK)//'rp.dat'
      OUTFILE(4)=trg(KK)//'rp_size.dat'
      DO JJ=1,4
         INQUIRE (FILE=outfile(JJ), EXIST=exists)
         if (exists) then
            open (unit=inunit,file=outfile(JJ),action='readwrite')
            close (unit=inunit,status='delete')
         endif
         open(unit=inunit,file=outfile(JJ),status='new',
     .       action='write')
         do ii=1,DIMENSIONS(JJ)
            write(inunit,*) bb(ii,jj)
         end do
         close(INUNIT);
      ENDDO
      DEALLOCATE(MM,BB)
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
