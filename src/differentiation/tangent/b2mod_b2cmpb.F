      module b2mod_b2cmpb
      use b2mod_types
      use b2mod_b2cmpa
*     (/b2cmpb/ contains physics parameters for boundary conditions)
      integer cbnrmx
      parameter (cbnrmx=20)
      integer, save ::
     1  cbirso, cbirno, cbirwe, cbirea, cbnrso, cbnrno, cbnrwe,
     &  cbnrea
      real (kind=R8), save ::
     *  cbrbrk(0:cbnrmx-1),
     &  cbsna(0:7,0:nsdecl-1,0:cbnrmx-1),
     &  cbsmo(0:7,0:nsdecl-1,0:cbnrmx-1),
     &  cbshi(0:7,0:nsdecl-1,0:cbnrmx-1),
     &  cbshe(0:7,0:cbnrmx-1), cbsch(0:7,0:cbnrmx-1)
      real (kind=R8), save ::
     *  cbrec(0:7,0:nsdecl-1,0:cbnrmx-1),
     &  cbmsa(0:7,0:cbnrmx-1), cbmsb(0:7,0:nsdecl-1,0:cbnrmx-1)
      end module b2mod_b2cmpb
