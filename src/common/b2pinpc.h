c  VERSION : 20.03.97 13:08
      integer nwldprtx,nwldprt,lwldstrt,
     ,				 lwldlis,nwldlis,lpumplis,npumplis,lppp
      real (kind=R4) :: ppp1(3,nlim),ppp2(3,nlim)
      parameter (nwldprtx=12)
      common/b2pprm/nwldprt,lwldstrt(nwldprtx),nwldlis(nwldprtx),
     ,	  lwldlis(nlim,nwldprtx),npumplis,lpumplis(nlim),lppp,ppp1,ppp2
