      module b2mod_braeir
      use b2mod_types
      implicit none
c
c*** Plasma data: ni,te,ti,vv,uu,pr,up,rr,fnix,fniy (Braams ---> Eirene)
c
      real (kind=R8), allocatable ::
     &  dnib(:,:,:), teb(:,:), tib(:,:), pob(:,:), prb(:,:),
     &  vvb(:,:,:), uub(:,:,:), wwb(:,:,:), upb(:,:,:),
     &  rrb(:,:), volb(:,:), bfeldb(:,:),
     &  bpolb(:,:), bradb(:,:), btorb(:,:),
     &  fnixb(:,:,:), fniyb(:,:,:),
     &  feixb(:,:), feiyb(:,:), feexb(:,:), feeyb(:,:),
     &  uudiab(:,:,:), vvdiab(:,:,:),
     &  delta_sheathxb(:,:), delta_sheathyb(:,:), zib(:,:,)

      contains

      subroutine alloc_b2mod_braeir(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(dnib(0:nxd+1,0:nyd+1,nsd))
      allocate(teb(0:nxd+1,0:nyd+1))
      allocate(tib(0:nxd+1,0:nyd+1))
      allocate(pob(0:nxd+1,0:nyd+1))
      allocate(prb(0:nxd+1,0:nyd+1))
      allocate(vvb(0:nxd+1,0:nyd+1,nsd))
      allocate(uub(0:nxd+1,0:nyd+1,nsd))
      allocate(wwb(0:nxd+1,0:nyd+1,nsd))
      allocate(upb(0:nxd+1,0:nyd+1,nsd))
      allocate(rrb(0:nxd+1,0:nyd+1))
      allocate(volb(0:nxd+1,0:nyd+1))
      allocate(bfeldb(0:nxd+1,0:nyd+1))
      allocate(bpolb(0:nxd+1,0:nyd+1))
      allocate(bradb(0:nxd+1,0:nyd+1))
      allocate(btorb(0:nxd+1,0:nyd+1))
      allocate(fnixb(0:nxd+1,0:nyd+1,nsd))
      allocate(fniyb(0:nxd+1,0:nyd+1,nsd))
      allocate(feixb(0:nxd+1,0:nyd+1))
      allocate(feiyb(0:nxd+1,0:nyd+1))
      allocate(feexb(0:nxd+1,0:nyd+1))
      allocate(feeyb(0:nxd+1,0:nyd+1))
      allocate(uudiab(0:nxd+1,0:nyd+1,nsd))
      allocate(vvdiab(0:nxd+1,0:nyd+1,nsd))
      allocate(delta_sheathxb(0:nxd+1,0:nyd+1))
      allocate(delta_sheathyb(0:nxd+1,0:nyd+1))
      allocate(zib(0:nxd+1,0:nyd+1,nsd))

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
      zib=0.0_R8

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
      deallocate(zib)

      return
      end subroutine dealloc_b2mod_braeir

      end module b2mod_braeir
