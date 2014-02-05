c  version : 19.02.97 22:24
c
c***  (Braams ---> Eirene)
c
#     include <BRAEIR.H>
c
c*** Plasma data: ni,te,ti,vv,uu,pr,up,rr,fnix,fniy
c
      real  (kind=R8) ::
     * dnib,teb,tib,pob,vvb,uub,prb,upb,rrb,fnixb,fniyb,
     r feixb,feiyb,feexb,feeyb,volb,bfeldb,uudiab,vvdiab
#ifndef NO_B2_CHEM_SPUT
c
c*** Data for chemical sputtering input from B2 - for compatibility!
c
      real  (kind=R8) ::
     * chemical_sputter_yield,fchar_chemical
      integer igass_chemical,itsput_chemical,issput_chemical
#endif
