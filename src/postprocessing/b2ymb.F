! ##########################################################################
! This program first reads the formatted file 'b2ymb.dat' and writes the 
! unformatted file 'mnb_in.dat' that is needed by the FreeFem++ heat 
! transfer model of the monoblock.
! Then, it reads the files 'wlld_trg[i-o].dat' produced by
! b2plot with the "1000 wlld" command and writes the unformatted files 
! 'Wtot_tg[i-o].dat', 'Wrad_tg[i-o].dat', 'X_tg[i-o].dat' and 
! 'tg[i-o]_size.dat' that are needed by the FreeFem++ heat transfer model
!  of the monoblock.
! In 'Wtot_tg[i-o].dat' is stored the value of the total heat flux for 
! each tile of the monoblock. 
! In 'Wpls_tg[i-o].dat' is stored the value of the heat flux due only to
! plasma particles for each tile of the monoblock. 
! In 'X_tg[i-o].dat' is stored the value of the poloidal coordinate of the
! center of each tile along the monoblock axis (in [m]). The strike point 
! is at x=0
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
      INTEGER :: DIMENSIONS(4)
      real :: target_length, LPOINTS(50)=0.0d0
      real, ALLOCATABLE :: MM(:,:),BB(:,:)
      CHARACTER*(10) :: INFILE0= 'b2ymb.dat', 
     . OUTFILE0 = 'mnb_in.dat'
      CHARACTER*(20) :: skip, OUTFILE(4),INFILE
      CHARACTER*(1) :: trg(2)
      LOGICAL :: exists

!##########################################################################
!Read input file 'b2ymb.dat' and then write 'mnb_in.dat'
!##########################################################################
      NAMELIST /GEOMETRY/	                
     . HEIGHT, T_LENGTH, CU_TH, CUCRZR_TH, 
     . DIAMETER, W_TH, P_LENGTH, ALPHA, BETA,
     . THETA, T_GAP

      NAMELIST /MESH/                     
     . NN

      NAMELIST /COOLANT/                  
     . T_COOLANT, M_FLOW, CP_COOLANT

      NAMELIST /OUTPUT/
     . NPOINTS, LPOINTS

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
      write(inunit,*) NN, NPOINTS, LPOINTS(1:NPOINTS)
      close(INUNIT);

! ##########################################################################
! Read input file 'wlld_trg[i-o].dat' to get heat flux profile on targets and
! then write the relative output files 'Wtot_tg[i-o].dat','Wpls_tg[i-o].dat' 
! 'X_tg[i-o].dat' and 'tg[i-o]_size.dat'.
! ##########################################################################
      trg(1)='i'
      trg(2)='o'

      DO KK=1,2 !First reads inner target, then outer target
         INFILE='wlld_trg'//trg(KK)//'.dat'
         open(unit=inunit,file=infile,status='old',action='read')
         READ(INUNIT,*) skip      !First two lines are text lines so must be skipped
         READ(INUNIT,*) skip
         NLINE=0
 5       READ(INUNIT,*,END=15) skip !File is of unknown length, so first count the number of lines
            NLINE = NLINE +1
         GO TO 5
 15      CLOSE(INUNIT)       
         ALLOCATE(MM(NLINE,4))     ! then allocate variable size and read again storing the values. First column will be x coordinate, 
                                ! second the total  heat flux, third the plasma heat flux and fourth the number of point
         MM=0.0d0
         MM(1,4)=dble(NLINE)
         open(unit=inunit,file=infile,status='old',action='read')
         READ(INUNIT,*) skip
         READ(INUNIT,*) skip
         DO II=1,NLINE
            READ(INUNIT,*) MM(II,1:2), dummy, dummy, MM(II,3) !The data in the file are x, Wtot, Wrad, Wpart, Wpls...
         ENDDO
         CLOSE(INUNIT)
         !Interpolate the data on a coordinate vector with spatial x coordinates with a distance equal to the tile poloidal length
         target_length=mm(nline,1)-mm(1,1)
         ntiles=nint(target_length/p_length)
         allocate(bb(ntiles,4))
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
         ! now write the output files, but first check if there are already present. If yes, remove them before writing
         OUTFILE(1)='X_tg'//trg(KK)//'.dat'
         OUTFILE(2)='Wtot_tg'//trg(KK)//'.dat'
         OUTFILE(3)='Wpls_tg'//trg(KK)//'.dat'
         OUTFILE(4)='tg'//trg(KK)//'_size.dat'
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
      ENDDO
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
