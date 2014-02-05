      module b2mod_b2cmrc
      use b2mod_types
*     (/b2cmrc/ contains atomic rate coefficients)
      integer, save :: rtnt, rtnn, rtns
      integer, save :: rtdumm
! old      parameter (rtnt=40, rtnn=16, rtnsmx=75)
      real (kind=R8), allocatable, save ::
     *  rtzmin(:), rtzmax(:), rtzn(:),
     &  rtt(:), rtn(:), rtlt(:), rtln(:),
     &  rtlsa(:,:,:),
     &  rtlra(:,:,:),
     &  rtlqa(:,:,:),
     &  rtlrd(:,:,:),
     &  rtlbr(:,:,:),
     &  rtlcx(:,:,:),
     &  rtlza(:,:,:),
     &  rtlz2(:,:,:),
     &  rtlpt(:,:,:),
     &  rtlpi(:,:,:)
      integer, allocatable, save ::
     &  rtyr(:)
      logical, save ::  adpak_used

      contains

      subroutine alloc_b2mod_b2cmrc(rtnt_in,rtnn_in,rtns_in)
      implicit none
      integer rtnt_in,rtnn_in,rtns_in

      rtnt=rtnt_in
      rtnn=rtnn_in
      rtns=rtns_in

      allocate(rtzmin(0:rtns-1))
      allocate(rtzmax(0:rtns-1))
      allocate(rtzn(0:rtns-1))
      allocate(rtt(0:rtnt))
      allocate(rtn(0:rtnn))
      allocate(rtlt(0:rtnt))
      allocate(rtln(0:rtnn))
      allocate(rtlsa(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlra(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlqa(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlrd(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlbr(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlcx(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlza(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlz2(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlpt(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtlpi(0:rtnt,0:rtnn,0:rtns-1))
      allocate(rtyr(0:rtns-1))
      
      return
      end subroutine alloc_b2mod_b2cmrc

      end module b2mod_b2cmrc
