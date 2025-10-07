      subroutine curved_d(X,Y,N)
      use b2mod_types
      implicit none
      integer N
      real (kind=R8) ::
     *   X(N),Y(N)
      real (kind=R4) X4(N),Y4(N)
      X4=real(X,R4)
      Y4=real(Y,R4)
#ifndef PLPLOT
      call curved(X4,Y4,N)
#endif
      return
      end subroutine curved_d
c
      subroutine gpl_d(N,X,Y)
      use b2mod_types
      implicit none
      integer N
      real (kind=R8) ::
     *   X(N),Y(N)
      real (kind=R4) X4(N),Y4(N)
      X4=real(X,R4)
      Y4=real(Y,R4)
#ifndef PLPLOT
      call gpl(N,X4,Y4)
#endif
      return
      end subroutine gpl_d
c
      subroutine gpm_d(N,X,Y)
      use b2mod_types
      implicit none
      integer N
      real (kind=R8) ::
     *   X(N),Y(N)
      real (kind=R4) X4(N),Y4(N)
      X4=real(X,R4)
      Y4=real(Y,R4)
#ifndef PLPLOT
      call gpm(N,X4,Y4)
#endif
      return
      end subroutine gpm_d
c
      subroutine set_d(VL,VR,VB,VT,WL,WR,WB,WT,LF)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   VL,VR,VB,VT,WL,WR,WB,WT
      integer LF
      real (kind=R4) VL4,VR4,VB4,VT4,WL4,WR4,WB4,WT4
      VL4=real(VL,R4)
      VR4=real(VR,R4)
      VB4=real(VB,R4)
      VT4=real(VT,R4)
      WL4=real(WL,R4)
      WR4=real(WR,R4)
      WB4=real(WB,R4)
      WT4=real(WT,R4)
#ifndef PLPLOT
      call set(VL4,VR4,VB4,VT4,WL4,WR4,WB4,WT4,LF)
#endif
      return
      end subroutine set_d
c
      subroutine gtx_d(X,Y,STRING)
      use b2mod_types
      real (kind=R8) ::
     *   X,Y
      real (kind=R4) X4,Y4
      character*(*) STRING
      X4=real(X,R4)
      Y4=real(Y,R4)
#ifndef PLPLOT
      call gtx(X4,Y4,STRING)
#endif
      return
      end subroutine gtx_d
c
      subroutine gsmksc_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=real(X,R4)
#ifndef PLPLOT
      call gsmksc(X4)
#endif
      return
      end subroutine gsmksc_d
c
      subroutine gschxp_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=real(X,R4)
#ifndef PLPLOT
      call gschxp(X4)
#endif
      return
      end subroutine gschxp_d
c
      subroutine gschsp_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=real(X,R4)
#ifndef PLPLOT
      call gschsp(X4)
#endif
      return
      end subroutine gschsp_d
c
      subroutine gschh_d(X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      real (kind=R4) X4
      X4=real(X,R4)
#ifndef PLPLOT
      call gschh(X4)
#endif
      return
      end subroutine gschh_d
c
      subroutine gschup_d(X,Y)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X,Y
      real (kind=R4) X4,Y4
      X4=real(X,R4)
      Y4=real(Y,R4)
#ifndef PLPLOT
      call gschup(X4,Y4)
#endif
      return
      end subroutine gschup_d
c
      subroutine gscr_d(WKID,CI,CR,CG,CB)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   CR,CG,CB
      integer WKID,CI
      real (kind=R4) CR4,CG4,CB4
      CR4=real(CR,R4)
      CG4=real(CG,R4)
      CB4=real(CB,R4)
#ifndef PLPLOT
      call gscr(WKID,CI,CR4,CG4,CB4)
#endif
      return
      end subroutine gscr_d
c
      subroutine cpsetr_d(STRING,X)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   X
      character*(*) STRING
      real (kind=R4) X4
      X4=real(X,R4)
#ifndef PLPLOT
      call cpsetr(STRING,X4)
#endif
      return
      end subroutine cpsetr_d
