      module b2mod_adpak
      use b2mod_types
      implicit none

c ...ADPAK variables from ADSDAT common
      integer nucz, nspc, nvalnc(100)
      real (kind=R8) ::
     *  apn(100,10), aqn(100,10), sigma(10,10), qn(100,10),
     &  en(100,10), ein(100,5), enm(100,5,10),
     &  fnm(100,5,10), fnn(100), enn(100)
c ...ADPAK variables from ADELEC common
      real (kind=R8) ::
     *  rclion(100), rrarec(100), rdirec(100), cizlos(100),
     &  radrrc(100), raddrc(100), radclx(100), radbrm(100)
c ...ADPAK variables from ADNEUT common
      real (kind=R8) ::
     *  rcxrec(100), radbcx(100)
c ...ADPAK variables from ADPARM common
      integer laden, ladtip, ltipok, leci, leciok
c ...ADPAK variables from ADDIEL common
      integer ldrmlt
      real (kind=R8) ::
     *  cdnn(100), cdnm(100),
     &  yhnna, yhnnb, yhnnc, yhnma, yhnmb, yhnmc
c ...ADPAK variables from ADCEQ common
      real (kind=R8) ::
     *  cefrac(100), ceavgz, cerad

      end module b2mod_adpak
