
*-----------------------------------------------------------------------
*.specification

      subroutine setwrk0 (phi,ph,v,gg0,np,nv,e,b,vwrk0,gwrk0)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer np, nv
      real (kind=R8) ::
     *  e, b, phi, gg0(0:np,0:nv), v(0:nv), ph(0:np)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: 
     *  vwrk0(0:nv),gwrk0(0:nv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*          
*     Set the values for the initial work arrays vwrk0 and gwrk0.
*     
*     Input parameters:
*     phi    - Incidence angle of the magnetic field with respect to the solid 
*              surface (in degrees).
*     ph     - Incidence angle array (in degrees).
*     v      - Normalised velocity array.
*     gg0    - gg0(ip,iv) is the ratio g/g0 for the (ip,iv) point of the 
*              angle-velocity space.
*     nv     - Number of points in velocity space - 1.
*     np     - Number of points in angle space - 1.
*     e      - Modulus of the electric field. 
*     b      - Modulus of the magnetic field.
*     
*     Output parameters:
*     vwrk0 - Initial work array for v
*     gwrk0 - Initial work array for gg0
*         
*     Array vwrk0 is computed using the formula vwrk0=v*e/b
*     Array gwrk0 is computed by interpolating two successive rows of 
*     matrix gg0, corresponding to phind and phsup. If phi lies outside
*     the ph array, the first or last row of matrix gg0 are just copied 
*     into gwrk0.
*     
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iv, ipinf, ipsup
      real (kind=R8) :: phinf, phsup, gginf, ggsup
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('setwrk0')
*   ..test input parameters
      call xertst (0.le.np, 'faulty parameter np')
      call xertst (0.le.nv, 'faulty parameter nv')
      call xertst (0.lt.b, 'faulty input b')
      call xertst (0.lt.phi, 'faulty input phi')

*   ..computation
*     
*     Assign gwrk0 array. Linear interpolation is used.
*     
      if (phi.le.ph(0)) then
        do iv = 0, nv
          gwrk0(iv) = gg0(0,iv)
        end do
      elseif (phi.gt.ph(np)) then
        do iv = 0, nv
          gwrk0(iv) = gg0(np,iv)
        end do
      else
        ipsup = 0
        phsup = ph(ipsup)
        do while (phsup.lt.phi)
          ipsup = ipsup + 1
          phsup = ph(ipsup)
        end do
        ipinf = ipsup - 1
        phinf = ph(ipinf)
        do iv = 0, nv
          gginf = gg0(ipinf,iv)
          ggsup = gg0(ipsup,iv)
          gwrk0(iv) = gginf + (ggsup-gginf)*
     &            (phi-phinf)/(phsup-phinf)
        end do
      end if

*     
*     Assign vwrk0 array
*     
      do iv = 0, nv
        vwrk0(iv) = v(iv)*e/b
      end do

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end setwrk0

      end subroutine setwrk0
