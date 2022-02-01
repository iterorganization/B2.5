      subroutine b2pwldl(BoRiS)
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot_wall_loading
c
c  version : 12.01.2002 00:37
c
c=======================================================================
c*** Calculate and write down the wall loading over the specified zones.
c*** The "plot zones" are specified in b2pwcom.h. They are assumed to
c*** be continous sets of the wall segments. 
c***
c*** The output tables are written to the files "wlld_NN.dat", one file
c*** for each plot zone. NN here marks the plot zone (e.g., 01, 02, ...) 
c*** The geometry of the shadowing structure and of the plot zones is 
c*** written to the files shadow.dat and wlly_NN.dat, the latter files 
c*** are produced in the subroutine "b2pwlprp". They contain also the 
c*** coordinates of the segment centers and the corresponding values of
c*** the distance along the wall.
c=======================================================================
      implicit none
      real (kind=R8) :: BoRiS

c*** Local variables

      integer lun_w, lfn_hlp
      parameter (lun_w=96, lfn_hlp=mlblprtx+9)
      integer i,l,il
      character fn_hlp*(lfn_hlp)
c=======================================================================

      print *,'Starting b2pwldl...' !###

c*** Calculate the radiation load on the specified wall zones
c*** (data input and geometry initialisation also done there)

      call b2pwrldw

c*** Produce the output files - one for a zone

      print *,'  nwldprt:',nwldprt  !### 
      do 20 il=2,nwldprt !{
        l=len_trim(pltzlbl(il))
c        print *,'il,l,pltzlbl(il)=',il,l,' ',pltzlbl(il)    !###
c*** create the filename for the output
        if(l.eq.0) then !{
c*** use the zone index to label the plot zone
          write(fn_hlp,'(a,i2.2,a)') 'wlld_',il-1,'.dat'
        else !}{
          fn_hlp='wlld_'//pltzlbl(il)(1:l)//'.dat'
        end if !}
c*** Open the file and write the heading line
        print *,'About to open : ',fn_hlp,' :'
        open(lun_w,file=fn_hlp,err=20)
        write(lun_w,'(a,i2,2a)') '#Wall loading zone ',il,' : ',label

        do i=1,nwldlis(il) !{
c*** Write the next line to the output file
          call b2pwldw(lun_w,.false.,il,i,BoRiS)
        end do !}
        close(lun_w)
   20 continue !}

      print *,'Leaving b2pwldl...' !###
      return
c=======================================================================
      end subroutine b2pwldl
