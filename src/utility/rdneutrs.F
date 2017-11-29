
*//RDNEUTRS//
*=======================================================================
*         S U B R O U T I N E  R D N E U T R S
*=======================================================================
*
*
*  1. PURPOSE
*
*     RDNEUTRS      READS THE DATA ON NEUTRAL WALL LOADING
*
*
*=======================================================================
*//DECLARATIONS//

      subroutine rdneutrs(kard,dummy,ldmf)
c
c*** reads the surface data
c
      use b2mod_types
      use b2mod_b2plot
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      integer kard,ldmf,is,ifl
      real  (kind=R8) ::
     * dummy(nlimps,*)
      logical end_of_file
      data end_of_file /.false./
      save end_of_file
c
      if(end_of_file) return
#ifdef B25_EIRENE
      call read_title(kard)
#endif
      do ifl=1,ldmf
        if(nlimi.gt.0) read(kard,910,end=10)
     1   (dummy(is,ifl),is=1,nlimi)
        if(nstsi.gt.0) read(kard,910,end=10)
     1   (dummy(is+nlim,ifl),is=1,nstsi)
      enddo
  910 format(5e16.8)
      return
   10 end_of_file=.true.
      write(*,*) 'End of file seen in rdneutrs'
      return
*=======================================================================
*//END RDNEUTRS//
      end subroutine rdneutrs

      subroutine wrneutrs(kard,dummy,ldmf)
      use b2mod_types
      use b2mod_b2plot
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      integer :: kard,ldmf,is,iif
      real (kind=R8) :: dummy(nlimps,*)
      do iif=1,ldmf
        if(nlimi.gt.0) write(kard,'(5(e16.8))') 
     &                 (dummy(is,iif),is=1,nlim1)
        if(nstsi.gt.0) write(kard,'(5(e16.8))') 
     &                 (dummy(is+nlim,iif),is=1,nstsi)
      enddo
      end subroutine wrneutrs

