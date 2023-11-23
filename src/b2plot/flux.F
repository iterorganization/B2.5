c
      subroutine flux(shotinfo,nCv,xlod,xhid,ylod,yhid)
      use b2mod_types
      use b2mod_b2plot_xy
      use b2mod_constants
      implicit none
#include <graphics_vars.h>
      character*(*) shotinfo
      integer nCv
      real (kind=R8) :: xlod,xhid,ylod,yhid

#ifdef WOS

      integer mdim,lpfxdim,lin
      parameter(mdim=256,lpfxdim=4,lin=1)

      integer m,n,ierr,nshot,nedit,lpfx
      integer lnblnk
      external lnblnk
      real (kind=R4) :: tshot
      character*4 expnam
      character*3 dianam
      real (kind=R4) :: ri(0:mdim),zj(0:mdim),pfm(0:mdim,0:mdim),
     1 pfxx(0:lpfxdim),rpfx(0:lpfxdim),zpfx(0:lpfxdim),
     2 rin(lin),zin(lin),br(lin),bz(lin),bt(lin),fpf(lin),fjp(lin)

      integer i,j
      logical main, streql
      external streql

      data nshot /0/

      save nshot

      if(streql(shotinfo,' ')) then
        nshot=0
        return
      endif

      call init_page(
     1 'Flux surfaces for '//shotinfo(1:lnblnk(shotinfo))//' (Wb/rad)',
     1 xlod,xhid,ylod,yhid,1,.true.,.false.,.true.)

      read(shotinfo,*,err=10,end=10) nshot,tshot,dianam,expnam,nedit
      goto 100
   10 continue
      nedit=0
      read(shotinfo,*,err=20,end=20) nshot,tshot,dianam,expnam
      goto 100
   20 continue
      expnam='AUGD'
      read(shotinfo,*) nshot,tshot,dianam
  100 continue

      m=mdim
      n=mdim
      call kkEQPFM(ierr,expnam,dianam,nshot,nedit,tshot,mdim+1,
     1 m,n,ri,zj,pfm)
      lpfx=lpfxdim
      call kkEQpfx(iERR,expnam,dianam,nSHOT,nEDIT,tSHOT,
     1 LPFx,PFxx,RPFx,zPFx)
      rin(1)=rpfx(0)
      zin(1)=zpfx(0)
      call kkrzBrzt (iERR,expnam,dianam,nSHOT,nEDIT,tSHOT,
     1 rin,zin,lin,Br,Bz,Bt,fPF,fJp)

      pfxx(1)=pfxx(1)/(2.0*pi)
      do i=0,m
        do j=0,n
          pfm(i,j)=pfm(i,j)/(2.0*pi)-pfxx(1)
        enddo
      enddo

      flow=pfm(0,0)
      fhgh=pfm(0,0)
      do i=0,m
        do j=0,n
          flow=min(flow,pfm(i,j))
          fhgh=max(fhgh,pfm(i,j))
        enddo
      enddo
      finc=max(abs(flow),abs(fhgh))

      finc=-30
      flow=1
      fhgh=-1

      main=.true.
      goto 200

      entry flux_add(nCv)
      main=.false.

  200 continue

      if(nshot.eq.0) return

      is=0
      do while(ri(is).lt.xlo.and.is.lt.m)
        is=is+1
      enddo
      ie=m
      do while(ri(ie).gt.xhi.and.ie.gt.0)
        ie=ie-1
      enddo
      js=0
      do while(zj(js).lt.ylo.and.js.lt.n)
        js=js+1
      enddo
      je=n
      do while(zj(je).gt.yhi.and.je.gt.0)
        je=je-1
      enddo

      if(is.ge.ie.or.js.ge.je) return

      CALL CPSETI('MAP - MAPPING FLAG',4)
      CALL CPSETR('ILX - INFORMATION LABEL X POSITION',0.5)
      CALL CPSETR('ILY - INFORMATION LABEL Y POSITION',0.95)
      CALL CPSETI('ILP - INFORMATION LABEL POSITION',-3)
      if(.not.main)
     1 CALL CPSETC('ILT - INFORMATION LABEL TEXT STRING',' ')
      ic=0
      do j=js,je
        do i=is,ie
          ic=ic+1
          x(ic)=ri(i)
          y(ic)=zj(j)
        enddo
      enddo
      nnx=ie-is+1
      nny=je-js+1
      call cpcnrc(pfm(is,js),mdim+1,nnx,nny,flow,fhgh,finc,1+3,0,-682)

      if(wall) call data_copy(0)
      if(main) call additional_elements(nCv,.true.)

#else
      write(*,*) 'Not implemented because WOS not defined'
#endif

      return
      end subroutine flux
