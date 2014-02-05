      subroutine hess_grad(obf, n, x, hess, lda, grad)
	implicit none
	integer n, lda
	real(8) :: x(n), step, hess(lda,n), grad(n)
	external obf

      integer i, j
      real(8) f(1:4), xmax, xmin, xi, xj
!    compute Hesse matrix. 
!     It costs 2*n**2+1 evaluations of objective functional
!     Gradient costs nothing!
	xmax = maxval(abs(x))
	xmin = minval(abs(x))
	if(xmax.eq.0.0d0) then
	 step = 0.1d0
	else
       step=0.0001d0*sqrt(xmax*xmin)
	endif
!$omp parallel
	do j =2,n
	 do i =1,j-1
	  xi = x(i)
	  xj = x(j)
	  x(i) = xi + step
	  x(j) = xj + step
	  call obf(n, x, f(1))
	  x(i) = xi - step
	  x(j) = xj + step 
	  call obf(n, x, f(2))
	  x(i) = xi + step
	  x(j) = xj - step
	  call obf(n, x, f(3))
	  x(i) = xi - step
	  x(j) = xj - step
	  call obf(n, x, f(4))
	  hess(i,j) = (f(1) - f(2) - f(3) + f(4))/(4.0d0*step**2)
	  hess(j,i) = hess(i,j)
	  x(i) = xi
	  x(j) = xj
	 enddo
	enddo
!     compute the diagonal
 	call obf(n, x, f(2))
	do i =1,n
	  xi = x(i)
	  x(i) = xi - 2.0d0*step
	  call obf(n, x, f(1))
	  x(i) = xi + 2.0d0*step
	  call obf(n, x, f(3))
	  hess(i,i) = (f(1) - 2.0d0*f(2) + f(3))/(4.0d0*step**2)
	  grad(i) = (f(3) - f(1))/(4.0d0*step)
	  x(i) = xi
	enddo
!$omp end parallel

	end subroutine hess_grad
