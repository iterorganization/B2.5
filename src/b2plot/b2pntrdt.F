      subroutine b2pntrdt(lun,indsrf,nindsrf)

c  version : 25.11.2003 18:33

c=======================================================================
c*** Calculate and write down the profiles of impinging neutral fluxes
c*** and energies on the specified surfaces
c=======================================================================
      use b2mod_types
      use b2mod_b2plot
      use b2mod_neutrals_namelist
      implicit none
      integer lun,nindsrf,indsrf(nindsrf)

      integer nn
      integer i,j,k
      real (kind=R8), allocatable :: hlp_flx(:),hlp_enr(:)
      character hlp_fn*16
      character*8, allocatable :: hlp_hd(:)
c=======================================================================
      print *,'Starting b2pntrdt. nindsrf=', nindsrf      !###
      if(nindsrf.le.0) return
      print *,'   indsrf=',(indsrf(i),i=1,nindsrf)        !###

      nn=natmi+nmoli
      allocate(hlp_flx(nn),hlp_enr(nn))
      allocate(hlp_hd(nn))

      do j=1,nindsrf !{
        k=j
        if(k.lt.0) k=nlim-k
        write(hlp_fn,'(a4,i4.4,a4)') 'ndt_',k,'.dat'
        print *,'filename: ',hlp_fn                     !###
        open(lun,file=hlp_fn)
        close(lun)
      end do !}
      print *,'*** b2pntrdt: neutral fluxes output option ***'
      print *,'*** not yet available ***'
      deallocate(hlp_flx,hlp_enr)
      deallocate(hlp_hd)
      return
!####
c=======================================================================
      end subroutine b2pntrdt
