!     ##########################################################################
!     This program first reads the formatted file 'b2ydm.dat' and 
!     writes the unformatted file 'rps_in.dat' that is needed by the FreeFem++
!     heat transfer model of the reflector plates.
!     Then, it reads the files 'wlld_[i-o]rp.dat' produced by
!     b2plot with the "1000 wlld" command and writes the unformatted files 
!     'Wtot_[i-o]rp.dat','Wpls_[i-o]rp.dat', 'X_[i-o]rp.dat' and 
!     '[i-o]rp_size.dat' that are needed  by the FreeFem++ heat transfer model 
!     of the reflector plates.
!     ##########################################################################

      PROGRAM b2ydm

      IMPLICIT NONE

      INTEGER :: NN, NLINE, NTILES, II, INUNIT=8888
      real :: P_LENGTH, HEIGHT, T_LENGTH, ALPHA, BETA, THETA, T_GAP
      real :: T_COOLANT, M_FLOW, CP_COOLANT, QTOT=0.0d0, T_OUT
      real :: target_length
      real, ALLOCATABLE :: MM(:,:),BB(:,:)
      CHARACTER*(15) :: INFILE0= 'b2ydm.dat', 
     .     OUTFILE0 = 'dome_out.dat'
      CHARACTER*(20) :: skip, INFILE='wlld_dome.dat'
      LOGICAL :: exists

!##########################################################################
!     Read input file 'b2ydm.dat'
!##########################################################################
      NAMELIST /GEOMETRY/	                
     .     HEIGHT, T_LENGTH, P_LENGTH, ALPHA, BETA,
     .     THETA, T_GAP

      NAMELIST /COOLANT/                 
     .     T_COOLANT, M_FLOW, CP_COOLANT

      open(unit=inunit,file=infile0,status='old',action='read')
      READ(INUNIT,NML=GEOMETRY)	
      READ(INUNIT,NML=COOLANT)
      close(INUNIT);
      
!     ##########################################################################
!     Read input files 'wlld_dome.dat' to get heat flux profile on dome and
!     evaluate outlet coolant temperature.
!     ##########################################################################

      open(unit=inunit,file=infile,status='old',action='read')
      READ(INUNIT,*) skip       !First two lines are text lines so must be skipped
      READ(INUNIT,*) skip
      NLINE=0
 5    READ(INUNIT,*,END=15) skip !File is of unknown length, so first count the number of lines
      NLINE = NLINE +1
      GO TO 5
 15   CLOSE(INUNIT)
      
      ALLOCATE(MM(NLINE,2))     ! then allocate variable size and read again storing the values. First column will be x coordinate, 
                                  ! second the total heat flux. 
      MM=0.0d0
      open(unit=inunit,file=infile,status='old',action='read')
      READ(INUNIT,*) skip
      READ(INUNIT,*) skip
      DO II=1,NLINE
         READ(INUNIT,*) MM(II,1:2) !The data in the file are x, Wtot, Wrad, Wpart, Wpls...
      ENDDO
      CLOSE(INUNIT)
      !Interpolate the data on a coordinate vector with spatial x coordinates with a distance equal to the tile poloidal length
      target_length=mm(nline,1)-mm(1,1)
      ntiles=nint(target_length/p_length)
      allocate(bb(ntiles,2))
      bb=0.0d0
      bb(1,1)=mm(1,1)
      bb(1,2)=mm(1,2)
      bb(ntiles,2)=mm(nline,2)
      do ii=2,ntiles            ! create the new x vector
         bb(ii,1)=bb(ii-1,1)+p_length
      enddo
      
      call interp( mm(1:nline,1), mm(1:nline,2), bb(2:ntiles-1,1),
     .     bb(2:ntiles-1,2), nline, ntiles-2)        
     
      DO II=1,NTILES
         QTOT=QTOT+T_LENGTH*P_LENGTH*bb(II,2)
      ENDDO

      T_OUT=QTOT/(M_FLOW*CP_COOLANT)+T_COOLANT

      ! Write output file
      INQUIRE (FILE=outfile0, EXIST=exists)
      if (exists) then 
         open (unit=inunit,file=outfile0,action='readwrite')
         close (unit=inunit,status='delete')
      endif
      open(unit=inunit,file=outfile0,status='new',action='write')
      write(inunit,*) T_OUT
      close(INUNIT);

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
