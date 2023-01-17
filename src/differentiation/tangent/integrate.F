
*-----------------------------------------------------------------------
*.specification

      subroutine integrate (vwrk0,gwrk0,w0,acc,nv,mu,ite)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nv
      real (kind=R8) :: w0, acc, vwrk0(0:nv), gwrk0(0:nv)
*   ..output arguments (unspecified on entry)
      integer ite

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*     
*     Perform the actual integration
*     
*     Input parameters
*     vwrk0 - Initial work array for v
*     gwrk0 - Initial work array for gg0
*     w0    - Work function for the material constituting the given surface.
*     acc   - Required accuracy.
*     nv    - Number of points in velocity space - 1.
*     
*     Output parameters
*     mu     - <g/g0>, averaged value.
*     ite    - Iteration number to reach the desired accuracy.
*     
*     Internal parameters
*     maxite - Maximum iteration number.
*     inv    - Actual number of velocity points considered.
*     vwrk   - Velocity work array.
*     gwrk   - g/g0 work array.
*     dwrk   - Distribution work array.
*     dvwrk  - Velocity increment work array.
*     vth    - Normalisation velocity.
*     eps    - Estimated accuracy.
*     nor    - Normalisation factor.
*     muc    - Control parameter (it should be 1.0 for absolute accuracy).
*     vmx    - Limit in velocity space of the region of fine grid. 
*              Set arbitrarily to 4*vth.
*     
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer inv,iv,ivinf,ivsup,maxite
      parameter (maxite=10)
      real (kind=R8) :: mu,muc,vth,nor,eps,vmx
      real (kind=R8), allocatable :: 
     *  vwrk(:),gwrk(:),dwrk(:),dvwrk(:),vwrkf(:)
*   ..procedures
      intrinsic abs, sqrt, real, exp
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('integrate')
*   ..test input parameters
      call xertst (0.le.nv, 'faulty parameter nv')
      call xertst (0.lt.acc, 'faulty parameter acc')

*   ..computation
      vth = sqrt(2.0_R8*w0/me)
      nor = 2.0_R8/(vth**4)
      vmx = 4.0_R8*vth 

      inv = nv
      eps = 2.0_R8 * acc
      ite = 0

*     
*     Loop to reach the required accuracy.
*     0) Initialise work arrays
*     1) Calculate integral using trapeze rule.
*     2) If the desired accuracy is not attained, double the
*        number of point and repeat.
*     
      allocate (vwrk(0:inv),dwrk(0:inv),dvwrk(0:inv-1))
      
*     
*     Create the initial velocity grid. The last "not infinity" point is
*     arbitrarily chosen to be vwrk(inv-1)=4*vth in order to guarantee an
*     acceptably small number of iterations. It is assumed that velocities
*     above 4*vth does not contribute appreciably to the integral value.
*     
      do iv = 0, inv
        vwrk(iv) = real(iv)/real(inv)*vmx
      end do

*     
*     Iterate to check accuracy
*     
      do while (eps.gt.acc)

*     
*     Increment iterations counter and check value
*       
        ite = ite + 1      
        call xertst (ite.le.maxite,
     .   'The required accuracy has not been reached')
*     
*     Assign the distribution array
*     
        do iv = 0, inv
          dwrk(iv) = (vwrk(iv)**3)*exp(-(vwrk(iv)/vth)**2)
          if (iv.lt.inv) dvwrk(iv) = vwrk(iv+1)-vwrk(iv)
        end do

*     
*     Integrate the distribution only to check accuracy
*     
        muc = 0.0_R8
        do iv = 0, inv-1
          muc = muc + 0.5_R8*(dwrk(iv)+dwrk(iv+1))*dvwrk(iv)
        end do
        muc = muc * nor
        eps = abs(1.0_R8-muc)

*     
*     If accuracy is too low, refine the velocity grid
*     
        if (eps.gt.acc) then
          inv = 2 * inv
          allocate (vwrkf(0:inv))
      
*     
*     Create the refined velocity array
*     
          do iv = 0, inv
            vwrkf(iv) = real(iv)/real(inv)*vmx
          end do

*     
*     Assign the new velocity array
*     
          deallocate (vwrk,dwrk,dvwrk)
          allocate (vwrk(0:inv),dwrk(0:inv),dvwrk(0:inv-1))
          do iv = 0, inv
            vwrk(iv) = vwrkf(iv)
          end do
          deallocate (vwrkf)
        end if
      end do

*     
*     Now interpole g/g0 on the chosen grid.
*     If the chosen velocity values exceed those already stored in VWRK0, 
*     the corresponding angles are set to be constants equals to the last 
*     value stored in GWRK0. This is consistent with data from Tskhakaya 
*     calculations. It is expected this velocity region to be unimportant 
*     except for very low incidence angles phi.
*     
      allocate (gwrk(0:inv))
      gwrk(0) = gwrk0(0)
      do iv = 1, inv
        ivsup = 0
        do while (vwrk0(ivsup).lt.vwrk(iv))
          ivsup = ivsup + 1
          if (ivsup.gt.nv) then
            gwrk(iv) = gwrk0(nv)
            goto 10
          end if
        end do
        ivinf = ivsup-1
        gwrk(iv) = gwrk0(ivinf) + (gwrk0(ivsup)-gwrk0(ivinf))*
     &    (vwrk(iv)-vwrk0(ivinf))/(vwrk0(ivsup)-vwrk0(ivinf))
   10   continue
      end do

*     
*     Finally, perform the integration
*     
      mu = 0.0_R8
      do iv = 0, inv-1
        mu = mu + 0.5_R8*(gwrk(iv)*dwrk(iv)+gwrk(iv+1)*dwrk(iv+1))*
     &                   dvwrk(iv)
      end do
      mu = mu * nor

*     Free memory and return
      deallocate (vwrk,gwrk,dwrk,dvwrk)
                        
*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end integrate

      end subroutine integrate
