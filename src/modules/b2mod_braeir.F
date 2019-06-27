      module b2mod_braeir
      use b2mod_types
      implicit none
c
c*** Plasma data: ni,te,ti,vv,uu,pr,up,rr,fnix,fniy (Braams ---> Eirene)
c
      real (kind=R8), allocatable ::
     &  dnib(:,:), teb(:), tib(:), pob(:), prb(:), neb(:),
     &  vvb(:,:), uub(:,:), wwb(:,:), upb(:,:),
     &  rrb(:), volb(:), bfeldb(:),
     &  bbxb(:), bbyb(:), bbzb(:),
     &  fnib(:,:), feib(:), feeb(:),
     &  vxb(:,:), vyb(:,:), vzb(:,:),
     &  uudiab(:,:), vvdiab(:,:),
     &  delta_sheathb(:)

      real (kind=R8), allocatable ::
     &  amb(:), znb(:), zab(:)

      REAL(kind=R8), ALLOCATABLE, SAVE, PUBLIC :: 
     &  PUX(:),  PUY(:),  PVX(:),  PVY(:)

      type :: targetdata
        integer :: ntrgdat
        integer, allocatable :: faces(:)
        real (kind=R8), allocatable :: fcori(:), flength(:)
        real (kind=R8), allocatable :: 
     &    tet(:), tit(:), fe(:), fsh(:),
     &    flux(:,:), dnit(:,:), vxt(:,:), vyt(:,:), vzt(:,:),
     &    fi(:,:), fel(:,:), vpart(:,:), mach(:,:), usr(:,:)
      end type targetdata

      type(targetdata), allocatable :: trgt(:)

      contains

      subroutine alloc_b2mod_braeir(nCv,nFc,nsd,nRc)
      implicit none
      integer, intent(in) :: nCv,nFc,nsd,nRc

      allocate(dnib(nCv,nsd))
      allocate(neb(nCv))
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
      allocate(bbxb(nCv))
      allocate(bbyb(nCv))
      allocate(bbzb(nCv))
      allocate(fnib(nFc,nsd))
      allocate(feib(nFc))
      allocate(feeb(nFc))
      allocate(vxb(nFc,nsd))
      allocate(vyb(nFc,nsd))
      allocate(vzb(nFc,nsd))
      allocate(uudiab(nCv,nsd))
      allocate(vvdiab(nCv,nsd))
      allocate(delta_sheathb(nFc))

      allocate(trgt(nrc))

      allocate(amb(nsd))
      allocate(znb(nsd))
      allocate(zab(nsd))

      ALLOCATE (PUX(NCV))
      ALLOCATE (PUY(NCV))
      ALLOCATE (PVX(NCV))
      ALLOCATE (PVY(NCV))

      dnib=0.0_R8
      neb=0.0_R8
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
      bbxb=0.0_R8
      bbyb=0.0_R8
      bbzb=0.0_R8
      fnib=0.0_R8
      feib=0.0_R8
      feeb=0.0_R8
      vxb=0.0_R8
      vyb=0.0_R8
      vzb=0.0_R8
      uudiab=0.0_R8
      vvdiab=0.0_R8
      delta_sheathb=0.0_R8

      amb=0.0_R8
      znb=0.0_R8
      zab=0.0_R8

      PUX  = 0.0_R8
      PUY  = 0.0_R8
      PVX  = 0.0_R8
      PVY  = 0.0_R8

      return
      end subroutine alloc_b2mod_braeir

      subroutine alloc_target_data(ist,ndat,nfl)
      implicit none
      integer, intent(in) :: ist,ndat,nfl

      if (ist < ubound(trgt,1)) then
        trgt(ist)%ntrgdat = ndat
        if (.not.allocated(trgt(ist)%faces)) then
          allocate (trgt(ist)%faces(ndat))
          allocate (trgt(ist)%fcori(ndat))
          allocate (trgt(ist)%flength(ndat))
          allocate (trgt(ist)%tet(ndat))
          allocate (trgt(ist)%tit(ndat))
          allocate (trgt(ist)%fe(ndat))
          allocate (trgt(ist)%fsh(ndat))
          allocate (trgt(ist)%flux(ndat,nfl))
          allocate (trgt(ist)%dnit(ndat,nfl))
          allocate (trgt(ist)%vxt(ndat,nfl))
          allocate (trgt(ist)%vyt(ndat,nfl))
          allocate (trgt(ist)%vzt(ndat,nfl))
          allocate (trgt(ist)%fi(ndat,nfl))
          allocate (trgt(ist)%fel(ndat,nfl))
          allocate (trgt(ist)%vpart(ndat,nfl))
          allocate (trgt(ist)%mach(ndat,nfl))
          allocate (trgt(ist)%usr(ndat,nfl))
        end if
        
        trgt(ist)%faces = 0
        trgt(ist)%fcori = 0.0_R8
        trgt(ist)%flength = 0.0_R8
        trgt(ist)%tet = 0.0_R8
        trgt(ist)%tit = 0.0_R8
        trgt(ist)%fe = 0.0_R8
        trgt(ist)%fsh = 0.0_R8
        trgt(ist)%flux = 0.0_R8
        trgt(ist)%dnit = 0.0_R8
        trgt(ist)%vxt = 0.0_R8
        trgt(ist)%vyt = 0.0_R8
        trgt(ist)%vzt = 0.0_R8
        trgt(ist)%fi = 0.0_R8
        trgt(ist)%fel = 0.0_R8
        trgt(ist)%vpart = 0.0_R8
        trgt(ist)%mach = 0.0_R8
        trgt(ist)%usr = 0.0_R8
        
      end if

      return
      end subroutine alloc_target_data

      subroutine dealloc_b2mod_braeir
      implicit none
      integer :: i

      deallocate(dnib)
      deallocate(neb)
      deallocate(teb)
      deallocate(tib)
      deallocate(pob)
      deallocate(vvb)
      deallocate(uub)
      deallocate(wwb)
      deallocate(prb)
      deallocate(upb)
      deallocate(rrb)
      deallocate(fnib)
      deallocate(feib)
      deallocate(feeb)
      deallocate(uudiab)
      deallocate(vvdiab)
      deallocate(volb)
      deallocate(bfeldb)
      deallocate(bbxb)
      deallocate(bbyb)
      deallocate(bbzb)
      deallocate(vxb)
      deallocate(vyb)
      deallocate(vzb)
      deallocate(delta_sheathb)

      deallocate (amb)
      deallocate (znb)
      deallocate (zab)

      DEALLOCATE (PUX)
      DEALLOCATE (PUY)
      DEALLOCATE (PVX)
      DEALLOCATE (PVY)


      IF (ALLOCATED(TRGT)) THEN
        DO I=LBOUND(TRGT,1), UBOUND(TRGT,1)
          IF (ALLOCATED(TRGT(I)%FACES)) THEN
            DEALLOCATE (TRGT(I)%FACES)
            DEALLOCATE (TRGT(I)%FCORI)
            DEALLOCATE (TRGT(I)%FLENGTH)
            DEALLOCATE (TRGT(I)%TET)
            DEALLOCATE (TRGT(I)%TIT)
            DEALLOCATE (TRGT(I)%FE)
            DEALLOCATE (TRGT(I)%FSH)
            DEALLOCATE (TRGT(I)%FLUX)
            DEALLOCATE (TRGT(I)%DNIT)
            DEALLOCATE (TRGT(I)%VXT)
            DEALLOCATE (TRGT(I)%VYT)
            DEALLOCATE (TRGT(I)%VZT)
            DEALLOCATE (TRGT(I)%FI)
            DEALLOCATE (TRGT(I)%FEL)
            DEALLOCATE (TRGT(I)%VPART)
            DEALLOCATE (TRGT(I)%MACH)
            DEALLOCATE (TRGT(I)%USR)
          END IF
        ENDDO
        DEALLOCATE (TRGT)
      END IF

      return
      end subroutine dealloc_b2mod_braeir

      end module b2mod_braeir
