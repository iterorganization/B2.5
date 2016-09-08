! ##########################################################################
! This program first reads the formatted file 'monoblock_input.dat' and writes the 
! unformatted file 'mnb_in.dat' that is needed by the FreeFem++ heat 
! transfer model of the monoblock.
! Then, it reads the files 'ld_tg_[i-o].dat' produced by
! b2plot with the "100 wlld" command and writes the unformatted files 
! 'W_tg[i-o].dat', 'X_tg[i-o].dat' and 'tg[i-o]_size.dat' that are needed 
! by the FreeFem++ heat transfer model of the monoblock.
! In 'W_tg[i-o].dat' is stored the value of the heat flux for each tile of
! the monoblock. 
! In 'X_tg[i-o].dat' is stored the value of the poloidal coordinate of the
! center of each tile along the monoblock axis (in [m]). The strike point 
! is at x=0
! In 'tg[i-o]_size.dat' is stored the number of points (i.e. tiles) in which
! the inner/outer monoblocks have been subdivided according to the assigned
! tile poloidal length, which is read from the 'monoblock_input.dat' file.
! ##########################################################################

      PROGRAM b2mb

      IMPLICIT NONE

      INTEGER :: NN, NLINE, NTILES, II, JJ, KK, INUNIT=8888
      real :: HEIGTH, T_LENGTH, CU_TH, CUCRZR_TH, DIAMETER, W_TH, 
     . P_LENGTH, BETA
      real :: T_COOLANT, M_FLOW, CP_COOLANT
      real dum(2), DIMENSIONS(3), target_length
      real, ALLOCATABLE :: MM(:,:),BB(:,:)
      CHARACTER*(30) :: INFILE0= 'monoblock_input.dat', 
     . OUTFILE0 = 'mnb_in.dat'
      CHARACTER*(30) :: testo, OUTFILE(3),INFILE
      CHARACTER*(1) :: trg(2)
      LOGICAL :: exists

!##########################################################################
!Read input file 'monoblock_input.dat' and then write 'mnb_in.dat'
!##########################################################################
      NAMELIST /GEOMETRY/	                
     . HEIGTH, T_LENGTH, CU_TH, CUCRZR_TH, 
     . DIAMETER, W_TH, P_LENGTH, BETA

      NAMELIST /MESH/                     
     . NN

      NAMELIST /COOLANT/                  
     . T_COOLANT, M_FLOW, CP_COOLANT
      open(unit=inunit,file=infile0,status='old',action='read')
      READ(INUNIT,NML=GEOMETRY)	
      READ(INUNIT,NML=MESH)	
      READ(INUNIT,NML=COOLANT)
      close(INUNIT);
      ! Write output file
      INQUIRE (FILE=outfile0, EXIST=exists)
      if (exists) then 
         open (unit=inunit,file=outfile0,action='readwrite')
         close (unit=inunit,status='delete')
      endif
      open(unit=inunit,file=outfile0,status='new',action='write')
      write(inunit,*) HEIGTH, T_LENGTH, CU_TH, CUCRZR_TH, DIAMETER,
     . W_TH, P_LENGTH, BETA
      write(inunit,*) T_COOLANT, M_FLOW, CP_COOLANT
      write(inunit,*) NN
      close(INUNIT);

! ##########################################################################
! Read input file 'ld_tg_[i-o].dat' to get heat flux profile on targets and
! then write the relative output files 'W_tg[i-o].dat', 'X_tg[i-o].dat' and
! 'tg[i-o]_size.dat'.
! ##########################################################################
      trg(1)='i'
      trg(2)='o'

      DO KK=1,2 !First reads inner target, then outer target
         INFILE='ld_tg_'//trg(KK)//'.dat'
         open(unit=inunit,file=infile,status='old',action='read')
         READ(INUNIT,*) testo      !First two lines are text lines so must be skipped
         READ(INUNIT,*) testo
         NLINE=0
 5       READ(INUNIT,*,END=15) dum !File is of unknown length, so first count the number of lines
            NLINE = NLINE +1
         GO TO 5
 15      CLOSE(INUNIT)
         
         ALLOCATE(MM(NLINE,3))     ! then allocate variable size and read again storing the values. First column will be x coordinate, 
                                ! second the heat flux and third the number of point
         MM=0.0d0
         MM(1,3)=dble(NLINE)
         open(unit=inunit,file=infile,status='old',action='read')
         READ(INUNIT,*) testo
         READ(INUNIT,*) testo
         DO II=1,NLINE
            READ(INUNIT,*) MM(II,1),dum,MM(II,2) !The data in the file are x, ne, Te, Wtot, ... so two dummy columns must be read
         ENDDO
         CLOSE(INUNIT)
         !Interpolate the data on a coordinate vector with spatial x coordinates with a distance equal to the tile poloidal length
         target_length=mm(nline,1)-mm(1,1)
         ntiles=nint(target_length/p_length)
         allocate(bb(ntiles,3))
         bb=0.0d0
         bb(1,1)=mm(1,1)
         bb(1,2)=mm(1,2)
         bb(ntiles,2)=mm(nline,2)
         bb(1,3)=ntiles
         do ii=2,ntiles            ! create the new x vector
            bb(ii,1)=bb(ii-1,1)+p_length
         enddo
         
         call interp( mm(1:nline,1), mm(1:nline,2), bb(2:ntiles-1,1),
     .    bb(2:ntiles-1,2), nline, ntiles-2)
         
         dimensions(1:2)=ntiles
         dimensions(3)=1
         
         ! now write the output files, but first check if there are already present. If yes, remove them before writing
         OUTFILE(1)='X_tg'//trg(KK)//'.dat'
         OUTFILE(2)='W_tg'//trg(KK)//'.dat'
         OUTFILE(3)='tg'//trg(KK)//'_size.dat'
         DO JJ=1,3
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
