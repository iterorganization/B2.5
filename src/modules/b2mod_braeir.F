      module b2mod_braeir
      use b2mod_types
      implicit none
c
c*** Plasma data: ni,te,ti,vv,uu,pr,up,rr,fnix,fniy (Braams ---> Eirene)
c
      real (kind=R8), allocatable ::
     &  dnib(:,:), teb(:), tib(:), pob(:), prb(:),
     &  vvb(:,:), uub(:,:), wwb(:,:), upb(:,:),
     &  rrb(:), volb(:), bfeldb(:),
     &  bpolb(:), bradb(:), btorb(:),
     &  fnixb(:,:), fniyb(:,:),
     &  feixb(:), feiyb(:), feexb(:), feeyb(:),
     &  uudiab(:,:), vvdiab(:,:),
     &  delta_sheathxb(:), delta_sheathyb(:)

      contains

      subroutine alloc_b2mod_braeir(nCv,nFc,nsd)
      implicit none
      integer nCv,nFc,nsd

      allocate(dnib(nCv,nsd))
      allocate(teb(nCv))
      allocate(tib(nCv))
      allocate(pob(nCv))
      allocate(prb(nCv))
      allocate(vvb(nCv,nsd))
      allocate(uub(nCv,nsd))
      allocate(wwb(nCv,nsd))
      allocate(upb(nCv,nsd))
      allocate(rrb(nCv))
      allocate(volb(nCv))
      allocate(bfeldb(nCv))
      allocate(bpolb(nCv))
      allocate(bradb(nCv))
      allocate(btorb(nCv))
      allocate(fnixb(nFc,nsd))
      allocate(fniyb(nFc,nsd))
      allocate(feixb(nFc))
      allocate(feiyb(nFc))
      allocate(feexb(nFc))
      allocate(feeyb(nFc))
      allocate(uudiab(nCv,nsd))
      allocate(vvdiab(nCv,nsd))
      allocate(delta_sheathxb(nFc))
      allocate(delta_sheathyb(nFc))

      dnib=0.0_R8
      teb=0.0_R8
      tib=0.0_R8
      pob=0.0_R8
      prb=0.0_R8
      vvb=0.0_R8
      uub=0.0_R8
      wwb=0.0_R8
      upb=0.0_R8
      rrb=0.0_R8
      volb=0.0_R8
      bfeldb=0.0_R8
      bpolb=0.0_R8
      bradb=0.0_R8
      btorb=0.0_R8
      fnixb=0.0_R8
      fniyb=0.0_R8
      feixb=0.0_R8
      feiyb=0.0_R8
      feexb=0.0_R8
      feeyb=0.0_R8
      uudiab=0.0_R8
      vvdiab=0.0_R8
      delta_sheathxb=0.0_R8
      delta_sheathyb=0.0_R8

      return
      end subroutine alloc_b2mod_braeir

      subroutine dealloc_b2mod_braeir
      implicit none

      deallocate(dnib)
      deallocate(teb)
      deallocate(tib)
      deallocate(pob)
      deallocate(vvb)
      deallocate(uub)
      deallocate(wwb)
      deallocate(prb)
      deallocate(upb)
      deallocate(rrb)
      deallocate(fnixb)
      deallocate(fniyb)
      deallocate(feixb)
      deallocate(feiyb)
      deallocate(feexb)
      deallocate(feeyb)
      deallocate(uudiab)
      deallocate(vvdiab)
      deallocate(volb)
      deallocate(bfeldb)
      deallocate(bpolb)
      deallocate(bradb)
      deallocate(btorb)
      deallocate(delta_sheathxb)
      deallocate(delta_sheathyb)

      return
      end subroutine dealloc_b2mod_braeir

      end module b2mod_braeir
