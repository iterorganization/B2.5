      module b2mod_b2cmpt
      use b2mod_types
      use b2mod_b2cmpa
*     (/b2cmpt/ contains physics parameters related to transport)
      real (kind=R8), save ::
     *  cfdf0(0:7,0:nsdecl-1), cfdna(0:7,0:nsdecl-1),
     &  cfdpa(0:7,0:nsdecl-1), cfvla(0:7,0:nsdecl-1),
     &  cfvsa(0:7,0:nsdecl-1), cfhci(0:7,0:nsdecl-1)
      real (kind=R8), save ::
     *  cfhce(0:7), cfsig(0:7), cfalf(0:7), cflim(0:7)
      end module b2mod_b2cmpt
