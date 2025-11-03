
*-----------------------------------------------------------------------
*.specification

      subroutine ggfill (ph,gg0,np,nv)
      use b2mod_types
      use b2mod_subsys
      implicit none      
*   ..input arguments (unchanged on exit)
      integer np, nv
      real (kind=R8) :: ph(0:np)
*   ..input/output argument 
      real (kind=R8) :: gg0(0:np,0:nv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     Fills the negative values of the matrix gg0 using linear interpolation
*     over columns.
*     
*     Input parameters:
*     nv     - Number of points in velocity space - 1.
*     np     - Number of points in angle space - 1.
*     gg0    - gg0(ip,iv) is the ratio g/g0 for the (ip,iv) point of the 
*              angle-velocity space.
*     ph     - Incidence angle array (in degrees).
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables      
      integer iv, ipinf, ipsup, ipmed
      real (kind=R8) :: gginf, ggsup
*   ..procedures
      external xertst
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('ggfill')
*   ..test input parameters
      call xertst (0.le.np, 'faulty parameter np')
      call xertst (0.le.nv, 'faulty parameter nv')

*   ..check first and last rows
      do iv = 0, nv
        call xertst (0.lt.gg0(0,iv), 
     .   'Error in matrix gg0: Missing value in first row.')
        call xertst (0.lt.gg0(np,iv),
     .   'Error in matrix gg0: Missing value in last row.')
      end do

*   ..fill the matrix, proceed column by column
      do iv = 0, nv
        ipinf = 0
        ipsup = 0
        do while (ipsup.lt.np)
          ipsup = ipsup + 1
          gginf = gg0(ipinf,iv)
          ggsup = gg0(ipsup,iv)          
          do while (ggsup.lt.0.0_R8)
            ipsup = ipsup + 1
            ggsup = gg0(ipsup,iv)
          end do
          if (ipsup.gt.(ipinf+1)) then
            do ipmed = ipinf+1, ipsup-1
              gg0(ipmed,iv) = gginf + (ggsup-gginf)*
     &         (ph(ipmed)-ph(ipinf))/(ph(ipsup)-ph(ipinf)) 
            end do
          end if
          ipinf = ipsup
        end do
      end do

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end ggfill

       end subroutine ggfill
