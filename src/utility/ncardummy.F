      subroutine curved_d(X,Y,N)
      use b2mod_types
      implicit none
      integer N
      real (kind=R8) ::
     *   X(N),Y(N)
      real (kind=R4) X4(N),Y4(N)
      X4=X
      Y4=Y
#ifndef PLPLOT
      call curved(X4,Y4,N)
#endif
      return
      end
c
      subroutine gpl_d(N,X,Y)
      use b2mod_types
      implicit none
      integer N
      real (kind=R8) ::
     *   X(N),Y(N)
      real (kind=R4) X4(N),Y4(N)
      X4=X
      Y4=Y
#ifndef PLPLOT
      call gpl(N,X4,Y4)
#endif
      return
      end
c
      subroutine gpm_d(N,X,Y)
      use b2mod_types
      implicit none
      integer N
      real (kind=R8) ::
     *   X(N),Y(N)
      real (kind=R4) X4(N),Y4(N)
      X4=X
      Y4=Y
#ifndef PLPLOT
      call gpm(N,X4,Y4)
#endif
      return
      end
c
      subroutine set_d(VL,VR,VB,VT,WL,WR,WB,WT,LF)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   VL,VR,VB,VT,WL,WR,WB,WT
      integer LF
      real (kind=R4) VL4,VR4,VB4,VT4,WL4,WR4,WB4,WT4
      VL4=VL
      VR4=VR
      VB4=VB
      VT4=VT
      WL4=WL
      WR4=WR
      WB4=WB
      WT4=WT
#ifndef PLPLOT
      call set(VL4,VR4,VB4,VT4,WL4,WR4,WB4,WT4,LF)
#endif
      return
      end
c
      subroutine gtx_d(X,Y,STRING)
      use b2mod_types
      real (kind=R8) ::
     *   X,Y
      real (kind=R4) X4,Y4
      character*(*) STRING
      X4=X
      Y4=Y
#ifndef PLPLOT
      call gtx(X4,Y4,STRING)
#endif
      return
      end
c
      subroutine gsmksc_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=X
#ifndef PLPLOT
      call gsmksc(X4)
#endif
      return
      end
c
      subroutine gschxp_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=X
#ifndef PLPLOT
      call gschxp(X4)
#endif
      return
      end
c
      subroutine gschsp_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=X
#ifndef PLPLOT
      call gschsp(X4)
#endif
      return
      end
c
      subroutine gschh_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=X
#ifndef PLPLOT
      call gschh(X4)
#endif
      return
      end
c
      subroutine gschup_d(X,Y)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X,Y
      real (kind=R4) X4,Y4
      X4=X
      Y4=Y
#ifndef PLPLOT
      call gschup(X4,Y4)
#endif
      return
      end
c
      subroutine gscr_d(WKID,CI,CR,CG,CB)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   CR,CG,CB
      integer WKID,CI
      real (kind=R4) CR4,CG4,CB4
      CR4=CR
      CG4=CG
      CB4=CB
#ifndef PLPLOT
      call gscr(WKID,CI,CR4,CG4,CB4)
#endif
      return
      end
c
      subroutine cpsetr_d(STRING,X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      character*(*) STRING
      real (kind=R4) X4
      X4=X
#ifndef PLPLOT
      call cpsetr(STRING,X4)
#endif
      return
      end
c
