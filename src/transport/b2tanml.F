*-----------------------------------------------------------------------
*.specification

      subroutine b2tanml (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     & csig_an, po, ne, na, rza, fchanml_a, fchanml)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     & po(nCv), ne(nCv), na(nCv,0:ns-1), rza(nCv,0:ns-1),
     & csig_an(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     & fchanml_a(nFc,0:1,0:ns-1), fchanml(nFc,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TANML computes the poloidal and radial components of the 
*     anomalous current.
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer :: is
      real (kind=R8) :: dpo(nFc,0:1), pov(nVx), wrk(nCv), wrkf(nFc),
     &  wght(nFc,2)
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx, diff_r
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tanml')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')

*   ..init
      wght = 1.0_R8

      if(switch%no_current.ne.0.or.
     &   switch%fch_anomalous.eq.0.0_R8) then                           !srv 27.02.13

        fchanml = 0.0_R8
        fchanml_a = 0.0_R8

      else                                                              !srv 27.02.13

*     ..compute differences of po
        call diff(nCv, nFc, nVx, 0, geo, mpg, po, pov, dpo)

*     ..compute the anomalous current
        fchanml = 0.0_R8
        do is = 0, ns-1
          if (.not.is_neutral(is)) then

*         ..interpolate density ratio to faces
            wrk = na(:,is)*rza(:,is)/ne
            call intface(nCv,nFc,mpg%fcCv,wght,wrk,wrkf)

*         ..compute the current contributions per species
            fchanml_a(:,0,is) = -wrkf*csig_an(:,0)*dpo(:,0)
            fchanml_a(:,1,is) = -wrkf*csig_an(:,1)*dpo(:,1)

*         ..accumulate individual contributions into fchanml
            fchanml = fchanml + fchanml_a(:,:,is)

          endif
        enddo

*     ..fully assign to main ion species if needed
        if (switch%anml_all_ions.eq.0) then
          do is = 0, ns-1
            if (is.eq.ismain) then
              fchanml_a(:,:,is) = fchanml
            else
              fchanml_a(:,:,is) = 0.0_R8
            endif
          enddo
        endif

      endif                                                             !srv 27.02.13

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tanml

      end subroutine b2tanml                                             !srv 20.09.06 }
