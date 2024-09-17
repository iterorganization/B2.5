*-----------------------------------------------------------------------
*.specification
      subroutine b2tfvh(nCv, nFc, ns, geo, mpg, na, ua, floi_vhx)
      use b2mod_types
      use b2mod_subsys
      use b2mod_zhfrtf
      use b2mod_b2zhco
      use b2mod_b2mrel , only: s2nucl
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  ua(nCv,0:ns-1), na(nCv,0:ns-1)
*   ..output arguments
      real (kind=R8) ::
     *  floi_vhx(nFc)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFVH computes the ion velocity-dependent heat flux.
*
*  3. description
*
*     This routine computes:
*     -floi_vhx: the ion velocity-dependent heat flux
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer iC, is, isp
      real (kind=R8) :: cfloi_vh(nCv)
*   ..procedures
      external intface
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2tfvh')
*   ..main
*   ..recalculate velocities and densities
      call b2mod_zhfrtf_prep(nCv,ns,na,ua)

*   ..velocity-dependent heat flux in the cell center
      cfloi_vh = 0.0_R8
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          do isp = 0, nspecies-1
            cfloi_vh = cfloi_vh + na(:,is)*
     &        c_hw_save(:,isp,s2nucl(is)-1)*(ua(:,is)-av_ualpha(:,isp))
          enddo
        endif
      enddo

*   ..velocity-dependent heat flux on the cell faces
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,cfloi_vh,floi_vhx)
      floi_vhx = floi_vhx*geo%fcPbs

      call subend ()
      return

      end subroutine b2tfvh


