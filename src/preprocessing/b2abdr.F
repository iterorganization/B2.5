*-----------------------------------------------------------------------
*.specification

      subroutine b2abdr (ninpab, noutab)
      use b2mod_types
      use b2mod_math
      use b2mod_time
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2us_prep
      use b2us_io
      use b2mod_ipmain
      use b2mod_subsys
      use b2mod_create_afn_bc_tables
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninpab(1), noutab(1:2)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2ABDR prepares and writes to an external file a database for
*     advanced fluid neutral boundary conditions for a certain combination of
*     incident particle and wall material.
*
*
*  3. description (see also routine b2cdca)
*
*  !! TODO_AFN_BCS
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninpab - (1) integer array, input.
*     ninpab specifies the input unit numbers, as follows.
*     ninpab(1): formatted; provides the incident particle and the wall material, and the desired discretization

*
*     noutab - (1:2) integer array, input.
*     noutab specifies the output unit numbers, as follows.
*     noutab(1): formatted; provides summary printed output.
*     noutab(2): un*formatted; provides the database
*
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      character id*30, incident*8, material*8
      integer ntp_in, nv_in, ntheta_in, nT_in, nM_in,
     &        nSh_in, nphi_in
      real (kind=R8) Tmin_in, Tmax_in, Mmin_in, Mmax_in,
     &               Shmin_in, Shmax_in, Eomin_in
*   ..procedures
      logical streql
      character cnamip*80, cvalip*80
      character line*500, ligne*500
      logical is_comment
      external xertst, xerrab, sfill, streql, cfwuch, cfwuin, is_comment

*-----------------------------------------------------------------------
*.documentation-internal
*
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2abdr')
*   ..test ninpab, noutab
      call xertst (1.le.ninpab(1).and.
     &             1.le.noutab(1).and.1.le.noutab(2),
     &  'faulty argument ninpab, noutab')

* ..set default values for optional input arguments. The default values are those with which the fnn.coarse databases in B2.5/Database/trim.data/refl.data were created.
      nT_in = 134
      nM_in = 11
      nSh_in = 7

      nv_in = 100
      ntp_in = 100000
      ntheta_in = 100
      nphi_in = 180

      Tmin_in = 0.1_R8
      Tmax_in = 1000_R8
      Mmin_in = 0.0_R8
      Mmax_in = 1.1_R8
      Shmin_in = 0.0_R8
      Shmax_in = 6.0_R8
      Eomin_in = 0.000001_R8

* ..extract info from b2ab.dat

*   ..read incident particle
      read (ninpab(1),'(a10)') id
      call xertst (streql(id,'*incident'),
     &  'faulty input id instead of *specs, id='//id)
*   ..read incident species
      read (ninpab(1),*) incident
*   ..read wall material
      read (ninpab(1),'(a10)') id
      call xertst (streql(id,'*material'),
     &  'faulty input id instead of *specs, id='//id)
*   ..read wall material
      read (ninpab(1),*) material

      write (noutab(1),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninpab(1),'(a)',end=2,err=93) line
       if (.not.is_comment(line)) then
         ligne = line
         call strip_spaces(ligne)
         if (ligne(1:1).eq.'''') then
           read (line,*) cnamip, cvalip
           write (noutab(1),'(4x,a,2x,a)') cnamip, trim(cvalip)
           call ipsetc (cnamip, cvalip)
         end if
       endif
       goto 1
    2 continue
      write (noutab(1),'(2x,a)') '(end of list of internal parameters)'

      call ipgeti('b2abdr_ntp', ntp_in)
      call ipgeti('b2abdr_nv', nv_in)
      call ipgeti('b2abdr_ntheta', ntheta_in)
      call ipgeti('b2abdr_nphi', nphi_in)

      call ipgeti('b2abdr_nT', nT_in)
      call ipgeti('b2abdr_nM', nM_in)
      call ipgeti('b2abdr_nSh', nSh_in)

      call ipgetr('b2abdr_Tmin', Tmin_in)
      call ipgetr('b2abdr_Tmax', Tmax_in)
      call ipgetr('b2abdr_Mmin', Mmin_in)
      call ipgetr('b2abdr_Mmax', Mmax_in)
      call ipgetr('b2abdr_Shmin', Shmin_in)
      call ipgetr('b2abdr_Shmax', Shmax_in)
      call ipgetr('b2abdr_Eomin', Eomin_in)

      call SetInput(nT_in, nM_in, nSh_in, ntp_in, nv_in,
     &              ntheta_in, nphi_in, Tmin_in, Tmax_in, Mmin_in,
     &              Mmax_in, Shmin_in, Shmax_in, Eomin_in)

      call read_TRIM(incident,material)
      call SetMass(incident)
      call CalcCoefficientsTRIM()
      call CalcFluidReflCoefficients()
      call write_fluid_TRIM(incident,material)
      call deallocate_all_TRIM()

*   ..close input/output files
      close (ninpab(1))
      close (noutab(1))
      close (noutab(2))

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   93 call xerrab
     1 ('error trying to read non-default internal parameters')

*-----------------------------------------------------------------------
*.end b2abdr

      end subroutine b2abdr
