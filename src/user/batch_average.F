
      subroutine batch_average (n,var,aver,it,batch_size)

      use b2mod_types

      implicit none
 
      integer :: n, it, batch_size, ii
      real (kind=R8) :: var(1:n), aver(1:n), gnew, gold
*-----------------------------------------------------------------------
*
*     Routine updates the batch average "aver" of data "var" in a series
*     indexed by it = 1:nt, i.e. the restarting average over 
*     "batch_size" samples.
*
*     for batch j:
*       aver(:) = mean(var(:),it=(j-1)*batch_size+1:j*batch_size)
*
*     on input:
*       aver(:) = mean(var(:),it=(j-1)*batch_size+1:(j-1)*batch_size+ii)
*
*     on output:
*       aver(:) = mean(var(:),it=(j-1)*batch_size+1:(j-1)*batch_size+ii+1)
*
*     If it = 0, the routine will update the running average
*       aver(:) = mean(var(:),it=1:batch_size)
*
*-----------------------------------------------------------------------

*  .. find index ii of var in the current batch
      ii = mod(it,batch_size)
      if (ii.eq.0) ii = batch_size

*  .. average
      gnew = 1.0_R8 / ii
      gold = (ii - 1) * gnew
      aver = gold * aver + gnew * var

      return
      end subroutine batch_average


      subroutine batch_average_sq (n,var,aver,it,batch_size)

      use b2mod_types

      implicit none
 
      real (kind=R8) :: var(1:n), aver(1:n), gnew, gold
      integer :: n, it, batch_size, fac

      fac = mod(it,batch_size)
      if (fac.eq.0) fac = batch_size
      if (fac.eq.1) then
        aver = var
      else
        gnew = 1.0_R8 / fac
        gold = (fac - 1) * gnew
        aver = gold*aver + gnew*var*var
      endif   

      end subroutine batch_average_sq
