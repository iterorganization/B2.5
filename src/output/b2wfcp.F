*-----------------------------------------------------------------------
*.specification

      subroutine b2wfcp (nprt, istyl, ns)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, istyl, ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WFCP writes to a formatted file the parameters in the physics
*     common blocks.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit nprt the species data and other
*     parameters stored in the physics common blocks: /b2cmpa/,
*     /b2cmpb/ and /b2cmpt/. Unit nprt must be connected to a
*     formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nprt - integer, input.
*     nprt specifies a unit number for formatted output.
*
*     istyl - integer, input.
*     istyl specifies the desired style of the output. If istyl.eq.0
*     then the output will appear in a compact format that is suitable
*     for input into the routine b2rfcp. If istyl.eq.1 then the output
*     will appear in a more elaborate format.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     (parameters in common:)
*
*     /b2cmpa/ - input.
*     /b2cmpa/ specifies basic physics parameters.
*     See the calling routine for description.
*
*     /b2cmpb/ - input.
*     /b2cmpb/ specifies physics parameters for boundary conditions.
*     See the calling routine for description.
*
*     /b2cmpt/ - input.
*     /b2cmpt/ specifies physics parameters for transport.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, i, ireg
      character*(*) fmt0, fmt1
      parameter (fmt0='(0p,2x,'''''''',i3,'''''''',1x,4(f9.3,1x))')
      parameter (fmt1='(1p,2x,'''''''',i3,'''''''',1x,8(g9.2,1x))')
*     (The sequence of eight quotes is not a joke. Because it occurs
*     inside a character string it represents four quotes, which in
*     turn represents one quoted quote in the output format.)
*   ..procedures
      external xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wfcp')
*   ..test nprt, istyl
      call xertst (1.le.nprt, 'faulty argument nprt')
      call xertst (istyl.eq.0.or.istyl.eq.1, 'faulty argument istyl')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')

* ..write data
*  ..use short style?
      if (istyl.eq.0) then
*   ..write heading for /b2cmpa/
       write (nprt,'(a)')
     &  '*b2cmpa         basic parameters'
*   ..write species data
       write (nprt,'(a)')
     &  '*specs      zamin    zamax    zn      am'
       do is = 0, ns-1
        write (nprt,fmt0) is,zamin(is),zamax(is),zn(is),am(is)
       enddo
*   ..write heading for /b2cmpb/
       if (cbnrso.gt.0.and.cbnrno.gt.0.and.
     &     cbnrwe.gt.0.and.cbnrea.gt.0) then
        write (nprt,'(a)')
     &   '*b2cmpb         boundary conditions'
*    ..write region specs
        write (nprt,'(a)')
     &   '*cbregs    nr    breakpoints'
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''south''', cbnrso, (cbrbrk(ireg),ireg=cbirso,cbirso+cbnrso-1)
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''north''', cbnrno, (cbrbrk(ireg),ireg=cbirno,cbirno+cbnrno-1)
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''west''', cbnrwe, (cbrbrk(ireg),ireg=cbirwe,cbirwe+cbnrwe-1)
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''east''', cbnrea, (cbrbrk(ireg),ireg=cbirea,cbirea+cbnrea-1)
*    ..loop over regions
        do ireg = 0, cbnrso+cbnrno+cbnrwe+cbnrea-1
         if (cbirso.le.ireg.and.ireg.lt.cbirso+cbnrso) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''south'''
         else if (cbirno.le.ireg.and.ireg.lt.cbirno+cbnrno) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''north'''
         else if (cbirwe.le.ireg.and.ireg.lt.cbirwe+cbnrwe) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''west'''
         else if (cbirea.le.ireg.and.ireg.lt.cbirea+cbnrea) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''east'''
         endif
*    ..write cbsna
         write (nprt,'(a)')
     &    '*cbsna    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbsna(i,is,ireg),i=0,7)
         enddo
*    ..write cbsmo
         write (nprt,'(a)')
     &    '*cbsmo    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbsmo(i,is,ireg),i=0,7)
         enddo
*    ..write cbshi
         write (nprt,'(a)')
     &    '*cbshi    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbshi(i,is,ireg),i=0,7)
         enddo
*    ..write cbshe
         write (nprt,'(a)')
     &    '*cbshe    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         write (nprt,fmt1) -1, (cbshe(i,ireg),i=0,7)
*    ..write cbsch
         write (nprt,'(a)')
     &    '*cbsch    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         write (nprt,fmt1) -1, (cbsch(i,ireg),i=0,7)
*    ..write cbrec
         write (nprt,'(a)')
     &    '*cbrec    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbrec(i,is,ireg),i=0,7)
         enddo
*    ..write cbmsa
         write (nprt,'(a)')
     &    '*cbmsa    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         write (nprt,fmt1) -1, (cbmsa(i,ireg),i=0,7)
*    ..write cbmsb
         write (nprt,'(a)')
     &    '*cbmsb    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbmsb(i,is,ireg),i=0,7)
         enddo
        enddo
       endif
*   ..write heading for /b2cmpt/
       write (nprt,'(a)')
     &  '*b2cmpt         transport coefficients'
*   ..write cfdf0
       write (nprt,'(a)')
     &  '*cfdf0    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       do is = 0, ns-1
        write (nprt,fmt1) is,(cfdf0(i,is),i=0,7)
       enddo
*   ..write cfdna
       write (nprt,'(a)')
     &  '*cfdna    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       do is = 0, ns-1
        write (nprt,fmt1) is,(cfdna(i,is),i=0,7)
       enddo
*   ..write cfdpa
       write (nprt,'(a)')
     &  '*cfdpa    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       do is = 0, ns-1
        write (nprt,fmt1) is,(cfdpa(i,is),i=0,7)
       enddo
*   ..write cfvla
       write (nprt,'(a)')
     &  '*cfvla    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       do is = 0, ns-1
        write (nprt,fmt1) is,(cfvla(i,is),i=0,7)
       enddo
*   ..write cfvsa
       write (nprt,'(a)')
     &  '*cfvsa    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       do is = 0, ns-1
        write (nprt,fmt1) is,(cfvsa(i,is),i=0,7)
       enddo
*   ..write cfhci
       write (nprt,'(a)')
     &  '*cfhci    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       do is = 0, ns-1
        write (nprt,fmt1) is,(cfhci(i,is),i=0,7)
       enddo
*   ..write cfhce
       write (nprt,'(a)')
     &  '*cfhce    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       write (nprt,fmt1) -1, (cfhce(i),i=0,7)
*   ..write cfsig
       write (nprt,'(a)')
     &  '*cfsig    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       write (nprt,fmt1) -1, (cfsig(i),i=0,7)
*   ..write cfalf
       write (nprt,'(a)')
     &  '*cfalf    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       write (nprt,fmt1) -1, (cfalf(i),i=0,7)
*   ..write cflim
       write (nprt,'(a)')
     &  '*cflim    (0)       (1)       (2)       (3)       (4)'//
     &     '       (5)       (6)       (7)'
       write (nprt,fmt1) -1, (cflim(i),i=0,7)
*  ..write end label
       write (nprt,'(a)')
     &  '*endphy'
*  ..use long style?
      else if (istyl.eq.1) then
*!!   Note: nice 'long style' is yet to be produced.
*   ..write species data
       write (nprt,101) ns, 
     &  (is,zamin(is),zamax(is),zn(is),am(is),is=0,ns-1)
*   ..write blank line
       write (nprt,'()')
*   ..write region specs
       if (cbnrso.gt.0.and.cbnrno.gt.0.and.
     &     cbnrwe.gt.0.and.cbnrea.gt.0) then
        write (nprt,'(a)')
     &   '*cbregs    nr    breakpoints'
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''south''', cbnrso, (cbrbrk(ireg),ireg=cbirso,cbirso+cbnrso-1)
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''north''', cbnrno, (cbrbrk(ireg),ireg=cbirno,cbirno+cbnrno-1)
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''west''', cbnrwe, (cbrbrk(ireg),ireg=cbirwe,cbirwe+cbnrwe-1)
        write (nprt,'(1x,a,t12,i2,(t17,7f9.3))')
     &   '''east''', cbnrea, (cbrbrk(ireg),ireg=cbirea,cbirea+cbnrea-1)
*   ..loop over regions
        do ireg = 0, cbnrso+cbnrno+cbnrwe+cbnrea-1
         if (cbirso.le.ireg.and.ireg.lt.cbirso+cbnrso) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''south'''
         else if (cbirno.le.ireg.and.ireg.lt.cbirno+cbnrno) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''north'''
         else if (cbirwe.le.ireg.and.ireg.lt.cbirwe+cbnrwe) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''west'''
         else if (cbirea.le.ireg.and.ireg.lt.cbirea+cbnrea) then
          write (nprt,'(a/2x,i2,2x,a)') '*region', ireg, '''east'''
         endif
*    ..write cbsna
         write (nprt,'(a)')
     &    '*cbsna    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbsna(i,is,ireg),i=0,7)
         enddo
*    ..write cbsmo
         write (nprt,'(a)')
     &    '*cbsmo    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbsmo(i,is,ireg),i=0,7)
         enddo
*    ..write cbshi
         write (nprt,'(a)')
     &    '*cbshi    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbshi(i,is,ireg),i=0,7)
         enddo
*    ..write cbshe
         write (nprt,'(a)')
     &    '*cbshe    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         write (nprt,fmt1) -1, (cbshe(i,ireg),i=0,7)
*    ..write cbsch
         write (nprt,'(a)')
     &    '*cbsch    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         write (nprt,fmt1) -1, (cbsch(i,ireg),i=0,7)
*    ..write cbrec
         write (nprt,'(a)')
     &    '*cbrec    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbrec(i,is,ireg),i=0,7)
         enddo
*    ..write cbmsa
         write (nprt,'(a)')
     &    '*cbmsa    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         write (nprt,fmt1) -1, (cbmsa(i,ireg),i=0,7)
*    ..write cbmsb
         write (nprt,'(a)')
     &    '*cbmsb    (0)       (1)       (2)       (3)       (4)'//
     &       '       (5)       (6)       (7)'
         do is = 0, ns-1
          write (nprt,fmt1) is,(cbmsb(i,is,ireg),i=0,7)
         enddo
        enddo
       endif
*   ..write transport coefficients
       do is = 0, ns-1
        write (nprt,110) is,(cfdf0(i,is),i=0,7)
       enddo
       do is = 0, ns-1
        write (nprt,111) is,(cfdna(i,is),i=0,7)
       enddo
       do is = 0, ns-1
        write (nprt,112) is,(cfdpa(i,is),i=0,7)
       enddo
       do is = 0, ns-1
        write (nprt,113) is,(cfvla(i,is),i=0,7)
       enddo
       do is = 0, ns-1
        write (nprt,114) is,(cfvsa(i,is),i=0,7)
       enddo
       do is = 0, ns-1
        write (nprt,115) is,(cfhci(i,is),i=0,7)
       enddo
       write (nprt,116) (cfhce(i),i=0,7),(cfsig(i),i=0,7),
     &  (cfalf(i),i=0,7)
       write (nprt,117) (cflim(i),i=0,7)
*   ..write blank line
       write (nprt,'()')
*  ..end branch on style
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  101 format (/'  There are ns species, ns =', i3,
     &  '.  The species properties (zamin,zamax,zn,am) are:'/
     &  ('  species ',i3,',  atomic charge range from ',f5.2,
     &   ' to ',f5.2,',  nuclear charge ',f5.2,
     &   ',  atomic mass ',f6.2))
  110 format (/1p,
     &  '  parameters for reference anomalous diffusivity:',/
     &  ('  species', i3, ', diffusivity coeffs ..  ', 8e11.3))
  111 format (/1p,
     &  '  parameters for anomalous particle transport:',/
     &  ('  species', i3, ', coeff. for grad(na) .. ', 8e11.3))
  112 format (/1p,
     &  '  parameters for anomalous particle transport:',/
     &  ('  species', i3, ', coeff. for grad(pa) .. ', 8e11.3))
  113 format (/1p,
     &  '  parameters for anomalous particle transport:',/
     &  ('  species', i3, ', coeff. for vconv ..    ', 8e11.3))
  114 format (/1p,
     &  '  parameters for anomalous parallel momentum transport:',/
     &  ('  species', i3, ', viscosity coefs ..     ', 8e11.3))
  115 format (/1p,
     &  '  parameters for anomalous atom heat transport:',/
     &  ('  species', i3, ', atom thermal diff ..   ', 8e11.3))
  116 format (/1p,
     &  '  parameters for anomalous electron heat transport',
     &    ' and electric current:',/
     &  '  electron thermal diffusivity ..    ', 8e11.3/
     &  '  electrical conductivity ..         ', 8e11.3/
     &  '  thermo-electric coefficient ..     ', 8e11.3)
  117 format (/1p,
     &  '  parameters for flux limiters:',/
     &  '  electron, atom, viscous ..         ', 8e11.3)

*-----------------------------------------------------------------------
*.end b2wfcp

      end subroutine b2wfcp
