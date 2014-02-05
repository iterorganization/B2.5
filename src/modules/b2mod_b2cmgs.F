      module b2mod_b2cmgs
*     (/b2cmgs/ contains a description of singular regions)
      integer nlrmx
      integer, save :: nlreg
      parameter (nlrmx=50)
      integer, save ::  nlxlo(0:nlrmx-1), nlxhi(0:nlrmx-1),
     &  nlylo(0:nlrmx-1), nlyhi(0:nlrmx-1), nlloc(0:nlrmx-1)
      end module b2mod_b2cmgs
