      module b2mod_b2cmrc_mol
      use b2mod_types
*     (/b2cmrc_mol/ contains molecular rate coefficients)
      integer, save :: rtnt, rtnn, rtnms, rtntr, rtn_mcoll
      integer, save :: rtdumm
! old      parameter (rtnt=40, rtnn=16, rtnsmx=75)
      real (kind=R8), allocatable, save ::
     &  rtt(:), rtn(:), rttr(:), rt_mcoll(:), 
     &  rtlt(:), rtln(:), rtltr(:), rtl_mcoll(:), 
     &  rtl_d(:,:,:),
     &  rtl_nd(:,:,:),
     &  rtl_di1(:,:,:),
     &  rtl_ic(:,:,:),
     &  rtl_dr(:,:,:),
     &  rtl_de(:,:,:),
     &  rtl_di2(:,:,:),
     &  rtl_ec_m(:,:,:,:),
     &  rtl_ec_E(:,:,:)

      contains

      subroutine alloc_b2mod_b2cmrc_mol(rtnt_in,rtnn_in,rtnms_in,
     .           rtntr_in, rtn_mcoll_in)
      implicit none
      integer rtnt_in,rtnn_in,rtnms_in, rtntr_in, 
     .        rtn_mcoll_in

      rtnt=rtnt_in
      rtnn=rtnn_in
      rtnms=rtnms_in
      rtntr = rtntr_in
      rtn_mcoll = rtn_mcoll_in

cc values
      allocate(rtt(0:rtnt))
      allocate(rtn(0:rtnn))
      allocate(rttr(0:rtntr))
      allocate(rt_mcoll(0:rtn_mcoll))
cc log of values
      allocate(rtlt(0:rtnt))
      allocate(rtln(0:rtnn))
      allocate(rtltr(0:rtntr))
      allocate(rtl_mcoll(0:rtn_mcoll))
cc logs of rates
      allocate(rtl_d(0:rtnt,0:rtnn,0:rtnms-1))
      allocate(rtl_nd(0:rtnt,0:rtnn,0:rtnms-1))
      allocate(rtl_di1(0:rtnt,0:rtnn,0:rtnms-1))
      allocate(rtl_ic(0:rtnt,0:rtnn,0:rtnms-1))
      allocate(rtl_dr(0:rtnt,0:rtnn,0:rtnms-1))
      allocate(rtl_de(0:rtnt,0:rtnn,0:rtnms-1))
      allocate(rtl_di2(0:rtnt,0:rtnn,0:rtnms-1))
      allocate(rtl_ec_m(0:rtnt,0:rtntr,0:rtn_mcoll,0:rtnms-1))
      allocate(rtl_ec_E(0:rtnt,0:rtntr, 0:rtnms-1))
      
      return
      end subroutine alloc_b2mod_b2cmrc_mol

      subroutine dealloc_b2mod_b2cmrc_mol
      implicit none


      deallocate(rtt)
      deallocate(rtn)
      deallocate(rtlt)
      deallocate(rtln)
      deallocate(rtl_d)
      deallocate(rtl_nd)
      deallocate(rtl_di1)
      deallocate(rtl_ec_m)
      deallocate(rtl_ec_E)
      deallocate(rtl_ic)
      deallocate(rtl_dr)
      deallocate(rtl_de)
      deallocate(rtl_di2)
      
      return
      end subroutine dealloc_b2mod_b2cmrc_mol


      end module b2mod_b2cmrc_mol
