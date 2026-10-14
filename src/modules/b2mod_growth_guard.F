      module b2mod_growth_guard
      
c Growth guard: keeps any single cell from changing too much in one
c nonlinear iteration, to stop ionization runaway and density/energy pile-ups
c There are six types of guard: density, velocity, e/i/n heat, total heat

      use b2mod_types, only : R8
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      implicit none
      private
      public :: growth_begin, growth_report, growth_diagonal,
     &          growth_bound, growth_factor, growth_record
      real(kind=R8), save :: step = 0.0_R8
      real(kind=R8), save :: matrix_peak(6), update_peak(6)
      integer, save :: matrix_count(6), update_count(6)
      integer, save :: matrix_cell(6), update_cell(6)
      contains

      subroutine growth_begin(dtim)
      real(kind=R8), intent(in) :: dtim
      call xertst(ieee_is_finite(dtim).and.dtim.gt.0.0_R8,
     &            'growth guard: invalid timestep')
      step = dtim
      matrix_peak = 0.0_R8
      update_peak = 0.0_R8
      matrix_count = 0
      update_count = 0
      matrix_cell = 0
      update_cell = 0
      end subroutine growth_begin

      recursive function growth_diagonal(res,scale,content,frac,
     &  channel,cell) result(diag)
      real(kind=R8), intent(in) :: res, scale, content, frac
      integer, intent(in) :: channel, cell
      real(kind=R8) :: diag, excess, threshold
c Extra matrix diagonal D that damps the correction to this cell, and is applied
c only if the residual would change the cell more than a fraction 'frac' per step:
c    D       = max(0, |R| - frac*content/dt)/(frac*scale).
c    content = particles / thermal energy / momentum in the cell
c    scale   = the unknown's size (1 for dp/p, T for heat)
c D goes on the matrix only, (A+D)delta=R, so the term is D*delta, not a
c source: a converged steady-state solution, with (R=0) D=0, is untouched.
      diag = 0.0_R8
      call xertst(ieee_is_finite(res).and.
     &  ieee_is_finite(scale).and.ieee_is_finite(content),
     &  'growth guard: nonfinite residual or content')
      call xertst(scale.gt.0.0_R8.and.content.gt.0.0_R8.and.
     &  step.gt.0.0_R8,'growth guard: invalid content or step')
      threshold = frac*(content/step)
      excess = abs(res)-threshold
      if (excess.le.0.0_R8) return
      diag = (excess/scale)/frac
      call xertst(ieee_is_finite(diag),
     &  'growth guard: correction diagonal overflow')
      call growth_record(channel,cell,excess/abs(res),.true.)
      end function growth_diagonal

      recursive function growth_factor(delta,scale,frac) result(alpha)
      real(kind=R8), intent(in) :: delta, scale, frac
      real(kind=R8) :: alpha, bound
c Factor alpha<=1 that shrinks an over-large correction
c 'delta' down to at most frac*scale.
      call xertst(ieee_is_finite(delta).and.ieee_is_finite(scale),
     &  'growth guard: nonfinite correction or scale')
      call xertst(scale.gt.0.0_R8,
     &  'growth guard: nonpositive correction scale')
      alpha = 1.0_R8
      bound = frac*scale
      if (abs(delta).gt.bound) alpha = bound/abs(delta)
      call xertst(alpha.gt.0.0_R8,
     &  'growth guard: correction factor underflow')
      end function growth_factor

      recursive function growth_bound(delta,scale,frac,channel,cell)
     &  result(limited)
      real(kind=R8), intent(in) :: delta, scale, frac
      integer, intent(in) :: channel, cell
      real(kind=R8) :: limited, alpha
      alpha = growth_factor(delta,scale,frac)
      limited = delta
      if (alpha.lt.1.0_R8) then
c Set the increment straight to its bound (same sign) instead of scaling
c a huge delta by a tiny alpha - numerically safer, still nonzero.
        limited = sign(frac*scale,delta)
        call growth_record(channel,cell,1.0_R8-alpha,.false.)
      endif
      end function growth_bound

      recursive subroutine growth_record(channel,cell,activity,matrix)
      integer, intent(in) :: channel, cell
      real(kind=R8), intent(in) :: activity
      logical, intent(in) :: matrix
c Accumulate one channel's count and its peak activity + cell.
      if (activity.le.0.0_R8) return
!$omp critical(b2_growth_guard_diagnostics)
      if (matrix) then
        matrix_count(channel) = matrix_count(channel)+1
        if (activity.gt.matrix_peak(channel)) then
          matrix_peak(channel) = activity
          matrix_cell(channel) = cell
        endif
      else
        update_count(channel) = update_count(channel)+1
        if (activity.gt.update_peak(channel)) then
          update_peak(channel) = activity
          update_cell(channel) = cell
        endif
      endif
!$omp end critical(b2_growth_guard_diagnostics)
      end subroutine growth_record

      subroutine growth_report(itim)
      integer, intent(in) :: itim
      integer :: k
      character(len=2), parameter :: names(6) =
     &  (/'na','ua','Te','Ti','Tn','Tt'/)
c One B2GUARD line per channel: activation count, peak activity, peak
c cell. 'count' totals activations (species x inner solves), not unique
c cells; diag peak = excess/|R|, update peak = clipped fraction 1-alpha.
      do k=1,6
        write(*,'(a,i9,1x,a,2(1x,i9,1x,es12.4,1x,i9))')
     &    'B2GUARD step/channel diag(n,max,cell) update(n,max,cell) ',
     &    itim,names(k),matrix_count(k),matrix_peak(k),matrix_cell(k),
     &    update_count(k),update_peak(k),update_cell(k)
      enddo
      end subroutine growth_report

      end module b2mod_growth_guard
