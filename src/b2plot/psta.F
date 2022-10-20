c
      subroutine psta()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
#ifdef B25_EIRENE
      use eirmod_ctrig
      implicit none
      integer it,is,lnblnk
      character*256 stmp
      intrinsic min,max
      external lnblnk

      values(38)=eirene_stack(1,0,iewrk)/triangle_area(1)
      values(39)=eirene_stack(1,0,iewrk)/triangle_area(1)
      do is=0,nedim(2,iewrk)-1
        do it=1,ntrii
          eirene_stack(it,is,iewrk)=
     .     eirene_stack(it,is,iewrk)/triangle_area(it)
          values(38)=min(values(38),eirene_stack(it,is,iewrk))
          values(39)=max(values(39),eirene_stack(it,is,iewrk))
        enddo
      enddo
      stmp=
     . eirene_label(iewrk)(1:lnblnk(eirene_label(iewrk)))//'/(m^{2})'
      call customise_label(stmp,eirene_label(iewrk))
#endif
      return
      end subroutine psta
