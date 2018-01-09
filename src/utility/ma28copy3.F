c######date   01 jan 1984     copyright ukaea, harwell.
c######alias ma28adcopy3 ma28bdcopy3 ma28cdcopy3
c###### calls   ma30    mc20    mc22    mc23    mc24
      subroutine ma28adcopy3                                            !srv 26.12.17 {
     & (n, nz, a, licn, irn, lirn,icn,u,ikeep,iw, w, iflag)
      use b2mod_types
      implicit none
      external bd1copy3, bd2copy3, bd3copy3
c this subroutine performs the LU factorisation of a.
c
c the parameters are as follows.....
c n     order of matrix  not altered by subroutine.
c nz    number of nonzeros in input matrix  not altered by subroutine.
c a is a  real*8 array  length licn.  holds nonzeros of matrix on entry
c     and nonzeros of factors on exit.  reordered by mc20a/ad and
c     mc23a/ad and altered by ma30a/ad.
c licn  integer  length of arrays a and icn.  not altered by subroutine.
c irn   integer array of length lirn.  holds row indices on input.
c     used as workspace by ma30a/ad to hold column orientation of
c     matrix.
c lirn  integer  length of array irn. not altered by the subroutine.
c icn   integer array of length licn.  holds column indices on entry
c     and column indices of decomposed matrix on exit. reordered by
c     mc20a/ad and mc23a/ad and altered by ma30a/ad.
c u     real*8 variable  set by user to control bias towards numeric or
c     sparsity pivoting.  u=1.0 gives partial pivoting while u=0. does
c     not check multipliers at all.  values of u greater than one are
c     treated as one while negative values are treated as zero.  not
c     altered by subroutine.
c ikeep  integer array of length 5*n  used as workspace by ma28a/ad
c     (see later comments).  it is not required to be set on entry
c     and, on exit, it contains information about the decomposition.
c     it should be preserved between this call and subsequent calls
c     to ma28b/bd or ma28c/cd.
c     ikeep(i,1),i=1,n  holds the total length of the part of row i
c     in the diagonal block.
c     row ikeep(i,2),i=1,n  of the input matrix is the ith row in
c     pivot order.
c     column ikeep(i,3),i=1,n  of the input matrix is the ith column
c     in pivot order.
c     ikeep(i,4),i=1,n  holds the length of the part of row i in
c     the l part of the l/u decomposition.
c     ikeep(i,5),i=1,n  holds the length of the part of row i in the
c     off-diagonal blocks.  if there is only one diagonal block,
c     ikeep(1,5) will be set to -1.
c iw    integer array of length 8*n.  if the option nsrch.le.n is
c     used, then the length of array iw can be reduced to 7*n.
c w real*8 array  length n.  used by mc24a/ad both as workspace and to
c     return growth estimate in w(1).  the use of this array by ma28a/ad
c     is thus optional depending on common block logical variable grow.
c iflag  integer variable  used as error flag by routine.  a positive
c     or zero value on exit indicates success.  possible negative
c     values are -1 through -14.
c
      integer n, nz, licn, lirn, iflag
      integer irn(lirn), icn(licn), ikeep(n,5), iw(n,8)
      real (kind=R8) ::
     *   a(licn), u, w(n)
c
c common and private variables.
c     common block ma28f/fd is used merely
c     to communicate with common block ma30f/fd  so that the user
c     need not declare this common block in his main program.
c the common block variables are as follows ...
c lp,mp  integer  default value 6 (line printer).  unit number
c     for error messages and duplicate element warning resp.
c nlp,mlp  integer  unit number for messages from ma30a/ad and
c     mc23a/ad resp.  set by ma28a/ad to value of lp.
c lblock  logical  default value true.  if true mc23a/ad is used
c     to first permute the matrix to block lower triangular form.
c grow    logical  default value true.  if true then an estimate
c     of the increase in size of matrix elements during l/u
c     decomposition is given by mc24a/ad.
c eps,rmin,resid  real/real*8 variables not referenced
c     by ma28a/ad.
c irncp,icncp  integer  set to number of compresses on arrays irn and
c     icn/a respectively.
c minirn,minicn  integer  minimum length of arrays irn and icn/a
c     respectively, for success on future runs.
c irank  integer   estimated rank of matrix.
c mirncp,micncp,mirank,mirn,micn integer variables.  used to
c     communicate between ma30f/fd and ma28f/fd values of abovenamed
c     variables with somewhat similar names.
c abort1,abort2  logical variables with default value true.  if false
c     then decomposition will be performed even if the matrix is
c     structurally or numerically singular respectively.
c aborta,abortb  logical variables used to communicate values of
c     abort1 and abort2 to ma30a/ad.
c abort  logical  used to communicate value of abort1 to mc23a/ad.
c abort3  logical variable not referenced by ma28a/ad.
c idisp   integer array  length 2.  used to communicate information
c     on decomposition between this call to ma28a/ad and subsequent
c     calls to ma28b/bd and ma28c/cd.  on exit, idisp(1) and
c     idisp(2) indicate position in arrays a and icn of the
c     first and last elements in the l/u decomposition of the
c     diagonal blocks, respectively.
c numnz  integer  structural rank of matrix.
c num    integer  number of diagonal blocks.
c large  integer  size of largest diagonal block.
c
c see block data for further comments on common block variables.
c see code for comments on private variables.
c
      real (kind=R8) ::
     *   tol, themax, big, dxmax, errmax, dres, cgce,
     * tol1, big1, upriv, rmin, eps, resid, zero
      integer lp, mp
      integer irncp, icncp, minirn, minicn, irank
      integer idisp(2)
      integer ndrop, maxit, noiter, nsrch, istart
      integer ndrop1, nsrch1
      integer nlp
      integer mirncp, micncp, mirank, mirn, micn
      integer mlp, numnz, num, large
      integer lpiv(10),lnpiv(10),mapiv,manpiv,iavpiv,ianpiv,kountl
      logical grow, lblock, abort, abort1, abort2, abort3, aborta,
     * abortb, lbig, lbig1
      common /ma28edc3/ lp, mp, lblock, grow
      common /ma28fdc3/ eps, rmin, resid, irncp, icncp, minirn, minicn,
     * irank, abort1, abort2
      common /ma28gdc3/ idisp
      common /ma28hdc3/ tol, themax, big, dxmax, errmax, dres, cgce,
     * ndrop, maxit, noiter, nsrch, istart, lbig
      common /ma30idc3/ tol1, big1, ndrop1, nsrch1, lbig1
      common /ma30edc3/ nlp, aborta, abortb, abort3
      common /ma30fdc3/ mirncp, micncp, mirank, mirn, micn
      common /mc23bdc3/ mlp, numnz, num, large, abort
      common /lpivotc3/ lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
c
      integer i, i1, ii, iend, j, j1, j2, jj, jay, newj1,
     * newpos, knum, length, move
c
c some  initialisation and transfer of information between
c     common blocks (see earlier comments).
      data zero /0.0e0_R8/
      iflag = 0
      aborta = abort1
      abortb = abort2
      abort = abort1
      mlp = lp
      nlp = lp
      tol1 = tol
      lbig1 = lbig
      nsrch1 = nsrch
c upriv private copy of u is used in case it is outside
c     range  zero to one  and  is thus altered by ma30a/ad.
      upriv = u
c simple data check on input variables and array dimensions.
      if (n.gt.0) go to 10
      iflag = -8
      if (lp.ne.0) write (lp,99999) n
      go to 210
   10 if (nz.gt.0) go to 20
      iflag = -9
      if (lp.ne.0) write (lp,99998) nz
      go to 210
   20 if (licn.ge.nz) go to 30
      iflag = -10
      if (lp.ne.0) write (lp,99997) licn
      go to 210
   30 if (lirn.ge.nz) go to 40
      iflag = -11
      if (lp.ne.0) write (lp,99996) lirn
      go to 210
c
c data check to see if all indices lie between 1 and n.
   40 do 50 i=1,nz
        if (irn(i).gt.0 .and. irn(i).le.n .and. icn(i).gt.0 .and.
     *   icn(i).le.n) go to 50
        if (iflag.eq.0 .and. lp.ne.0) write (lp,99995)
        iflag = -12
        if (lp.ne.0) write (lp,99994) i, a(i), irn(i), icn(i)
   50 continue
      if (iflag.lt.0) go to 220
c
c sort matrix into row order.
      call mc20adcopy3(n, nz, a, icn, iw, irn, 0)
c part of ikeep is used here as a work-array.  ikeep(i,2) is
c     the last row to have a nonzero in column i.  ikeep(i,3)
c     is the off-set of column i from the start of the row.
      do 60 i=1,n
        ikeep(i,2) = 0
        ikeep(i,1) = 0
   60 continue
c
c check for duplicate elements .. summing any such entries and
c     printing a warning message on unit mp.
c move is equal to the number of duplicate elements found.
      move = 0
c the loop also calculates the largest element in the matrix, themax.
      themax = zero
c j1 is position in arrays of first nonzero in row.
      j1 = iw(1,1)
      do 130 i=1,n
        iend = nz + 1
        if (i.ne.n) iend = iw(i+1,1)
        length = iend - j1
        if (length.eq.0) go to 130
        j2 = iend - 1
        newj1 = j1 - move
        do 120 jj=j1,j2
          j = icn(jj)
          themax = max(themax,abs(a(jj)))
          if (ikeep(j,2).eq.i) go to 110
c first time column has ocurred in current row.
          ikeep(j,2) = i
          ikeep(j,3) = jj - move - newj1
          if (move.eq.0) go to 120
c shift necessary because of  previous duplicate element.
          newpos = jj - move
          a(newpos) = a(jj)
          icn(newpos) = icn(jj)
          go to 120
c duplicate element.
  110     move = move + 1
          length = length - 1
          jay = ikeep(j,3) + newj1
          if (mp.ne.0) write (mp,99993) i, j, a(jj)
          a(jay) = a(jay) + a(jj)
          themax = max(themax,abs(a(jay)))
  120   continue
        ikeep(i,1) = length
        j1 = iend
  130 continue
c
c knum is actual number of nonzeros in matrix with any multiple
c     entries counted only once.
      knum = nz - move
      if (.not.lblock) go to 140
c
c perform block triangularisation.
      call mc23adcopy3(n, icn, a, licn, ikeep, idisp, ikeep(1,2),
     *ikeep(1,3), ikeep(1,5), iw(1,3), iw)
      if (idisp(1).gt.0) go to 170
      iflag = -7
      if (idisp(1).eq.-1) iflag = -1
      if (lp.ne.0) write (lp,99992)
      go to 210
c
c block triangularisation not requested.
c move structure to end of data arrays in preparation for
c     ma30a/ad.
c also set lenoff(1) to -1 and set permutation arrays.
  140 do 150 i=1,knum
        ii = knum - i + 1
        newpos = licn - i + 1
        icn(newpos) = icn(ii)
        a(newpos) = a(ii)
  150 continue
      idisp(1) = 1
      idisp(2) = licn - knum + 1
      do 160 i=1,n
        ikeep(i,2) = i
        ikeep(i,3) = i
  160 continue
      ikeep(1,5) = -1
  170 if (lbig) big1 = themax
      if (nsrch.le.n) go to 180
c
c perform l/u decomosition on diagonal blocks.
      call ma30adcopy3(n, icn, a, licn, ikeep, ikeep(1,4), idisp,
     *ikeep(1,2), ikeep(1,3), irn, lirn, iw(1,2), iw(1,3), iw(1,4),
     *iw(1,5), iw(1,6), iw(1,7), iw(1,8), iw, upriv, iflag)
      go to 190
c this call if used if nsrch has been set less than or equal n.
c     in this case, two integer work arrays of length can be saved.
  180 call ma30adcopy3(n, icn, a, licn, ikeep, ikeep(1,4), idisp,
     * ikeep(1,2), ikeep(1,3), irn, lirn, iw(1,2), iw(1,3), iw(1,4),
     * iw(1,5), iw, iw, iw(1,6), iw, upriv, iflag)
c
c transfer common block information.
  190 minirn = max0(mirn,nz)
      minicn = max0(micn,nz)
      irncp = mirncp
      icncp = micncp
      irank = mirank
      ndrop = ndrop1
      if (lbig) big = big1
      if (iflag.ge.0) go to 200
      if (lp.ne.0) write (lp,99991)
      go to 210
c
c reorder off-diagonal blocks according to pivot permutation.
  200 i1 = idisp(1) - 1
      if (i1.ne.0) call mc22adcopy3(n, icn, a, i1,ikeep(1,5),ikeep(1,2),
     * ikeep(1,3), iw, irn)
      i1 = idisp(1)
      iend = licn - i1 + 1
c
c optionally calculate element growth estimate.
      if (grow) call mc24adcopy3(n, icn, a(i1), iend,ikeep,ikeep(1,4),w)
c increment growth estimate by original maximum element.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n.gt.1) w(2) = themax
c set flag if the only error is due to duplicate elements.
      if (iflag.ge.0 .and. move.ne.0) iflag = -14
      go to 220
  210 if (lp.ne.0) write (lp,99990)
  220 return
99999 format (36x, 'n out of range = ', i10)
99998 format (36x, 'nz non positive = ', i10)
99997 format (36x, 'licn too small = ', i10)
99996 format (36x, 'lirn too small = ', i10)
99995 format (' error return from ma28a/adc because indices found out',
     * ' of range')
99994 format (1x, i6, 'th element with value ', 1pd22.14, ' is out o',
     * 'f range with indices ', i8, ' ,', i8)
99993 format (' duplicate element in position ', i8, ' ,', i8,
     * ' with value ', 1pd22.14)
99992 format (36x, 'error return from mc23a/adc')
99991 format (36x, 'error return from ma30a/adc')
99990 format ('+error return from ma28a/adc because ')
      end subroutine ma28adcopy3
      subroutine ma28bdcopy3
     &  (n, nz, a, licn, ivect, jvect, icn, ikeep, iw, w, iflag)
      use b2mod_types
      implicit none
      external bd1copy3, bd2copy3
c this subroutine factorises a matrix of a similar sparsity
c     pattern to that previously factorised by ma28a/ad.
c the parameters are as follows ...
c n      integer  order of matrix  not altered by subroutine.
c nz     integer  number of nonzeros in input matrix  not altered
c     by subroutine.
c a      real/real*8 array  length licn.  holds nonzeros of
c     matrix on entry and nonzeros of factors on exit.  reordered by
c     ma28d/dd and altered by subroutine ma30b/bd.
c licn   integer  length of arrays a and icn.  not altered by
c     subroutine.
c ivect,jvect  integer arrays of length nz.  hold row and column
c     indices of nonzeros respectively.  not altered by subroutine.
c icn    integer array of length licn.  same array as output from
c     ma28a/ad.  unchanged by ma28b/bd.
c ikeep  integer array of length 5*n.  same array as output from
c     ma28a/ad.  unchanged by ma28b/bd.
c iw     integer array  length 5*n.  used as workspace by ma28d/dd and
c     ma30b/bd.
c w      real/real*8 array  length n.  used as workspace
c     by ma28d/dd,ma30b/bd and (optionally) mc24a/ad.
c iflag  integer  used as error flag with positive or zero value
c     indicating success.
c
      integer n, nz, licn, iw(n,5), iflag
      integer ikeep(n,5), ivect(nz), jvect(nz), icn(licn)
      real (kind=R8) ::
     *   a(licn), w(n)
c
c private and common variables.
c unless otherwise stated common block variables are as in ma28a/ad.
c     those variables referenced by ma28b/bd are mentioned below.
c lp,mp  integers  used as in ma28a/ad as unit number for error and
c     warning messages, respectively.
c nlp    integer variable used to give value of lp to ma30e/ed.
c eps    real/real*8  ma30b/bd will output a positive value
c     for iflag if any modulus of the ratio of pivot element to the
c     largest element in its row (u part only) is less than eps (unless
c     eps is greater than 1.0 when no action takes place).
c rmin   real/real*8  variable equal to the value of this
c     minimum ratio in cases where eps is less than or equal to 1.0.
c meps,mrmin  real/real*8 variables used by the subroutine
c     to communicate between common blocks ma28f/fd and ma30g/gd.
c idisp  integer array  length 2  the same as that used by ma28a/ad.
c     it is unchanged by ma28b/bd.
c
c see block data or ma28a/ad for further comments on variables
c     in common.
c see code for comments on private variables.
c
      logical grow, lblock, aborta, abortb, abort1, abort2, abort3,
     * lbig, lbig1
      integer lp, mp
      integer irncp, icncp, minirn, minicn, irank
      integer idisp(2)
      integer ndrop, maxit, noiter, nsrch, istart
      integer nlp
      integer ndrop1, nsrch1
      real (kind=R8) ::
     *   eps, meps, rmin, mrmin, resid, tol,
     * themax, big, dxmax, errmax, dres, cgce, tol1, big1
c
      common /ma28edc3/ lp, mp, lblock, grow
      common /ma28fdc3/ eps, rmin, resid, irncp, icncp, minirn, minicn,
     * irank, abort1, abort2
      common /ma28gdc3/ idisp
      common /ma28hdc3/ tol, themax, big, dxmax, errmax, dres, cgce,
     * ndrop, maxit, noiter, nsrch, istart, lbig
      common /ma30edc3/ nlp, aborta, abortb, abort3
      common /ma30gdc3/ meps, mrmin
      common /ma30idc3/ tol1, big1, ndrop1, nsrch1, lbig1
c
      integer i1, idup, iend
c
c check to see if elements were dropped in previous ma28a/ad call.
      if (ndrop.eq.0) go to 10
      iflag = -15
      write (6,99999) iflag, ndrop
      go to 70
   10 iflag = 0
      meps = eps
      nlp = lp
c simple data check on variables.
      if (n.gt.0) go to 20
      iflag = -11
      if (lp.ne.0) write (lp,99998) n
      go to 60
   20 if (nz.gt.0) go to 30
      iflag = -10
      if (lp.ne.0) write (lp,99997) nz
      go to 60
   30 if (licn.ge.nz) go to 40
      iflag = -9
      if (lp.ne.0) write (lp,99996) licn
      go to 60
c
   40 call ma28ddcopy3
     &  (n, a, licn, ivect, jvect, nz, icn, ikeep, ikeep(1,4),
     * ikeep(1,5), ikeep(1,2), ikeep(1,3), iw(1,3), iw, w(1), iflag)
c themax is largest element in matrix.
      themax = w(1)
      if (lbig) big1 = themax
c idup equals one if there were duplicate elements, zero otherwise.
      idup = 0
      if (iflag.eq.(n+1)) idup = 1
      if (iflag.lt.0) go to 60
c
c perform row-gauss elimination on the structure received from ma28d/dd
      call ma30bdcopy3(n, icn, a, licn, ikeep, ikeep(1,4), idisp,
     * ikeep(1,2), ikeep(1,3), w, iw, iflag)
c
c transfer common block information.
      if (lbig) big1 = big
      rmin = mrmin
      if (iflag.ge.0) go to 50
      iflag = -2
      if (lp.ne.0) write (lp,99995)
      go to 60
c
c optionally calculate the growth parameter.
   50 i1 = idisp(1)
      iend = licn - i1 + 1
      if (grow) call mc24adcopy3(n, icn, a(i1), iend,ikeep,ikeep(1,4),w)
c increment estimate by largest element in input matrix.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n.gt.1) w(2) = themax
c set flag if the only error is due to duplicate elements.
      if (idup.eq.1 .and. iflag.ge.0) iflag = -14
      go to 70
   60 if (lp.ne.0) write (lp,99994)
   70 return
99999 format (' error return from ma28b/bdc with iflag=', i4/i7,
     * ' entries dropped from structure by ma28a/adc')
99998 format (36x, 'n out of range = ', i10)
99997 format (36x, 'nz non positive = ', i10)
99996 format (36x, 'licn too small = ', i10)
99995 format (36x, 'error return from ma30b/bdc')
99994 format ('+error return from ma28b/bdc because ')
      end subroutine ma28bdcopy3
      subroutine ma28cdcopy3(n, a, licn, icn, ikeep, rhs, w, mtype)
      use b2mod_types
      implicit none
      external bd1copy3, bd2copy3
c
c this subroutine uses the factors from ma28a/ad or ma28b/bd to
c     solve a system of equations without iterative refinement.
c the parameters are ...
c n   integer  order of matrix  not altered by subroutine.
c a   real/real*8 array  length licn.  the same array as
c     was used in the most recent call to ma28a/ad or ma28b/bd.
c licn  integer  length of arrays a and icn.  not altered by
c     subroutine.
c icn    integer array of length licn.  same array as output from
c     ma28a/ad.  unchanged by ma28c/cd.
c ikeep  integer array of length 5*n.  same array as output from
c     ma28a/ad.  unchanged by ma28c/cd.
c rhs    real/real*8 array  length n.  on entry, it holds the
c     right hand side.  on exit, the solution vector.
c w      real/real*8 array  length n. used as workspace by
c     ma30c/cd.
c mtype  integer  used to tell ma30c/cd to solve the direct equation
c     (mtype=1) or its transpose (mtype.ne.1).
c
      integer n, licn, mtype
      real (kind=R8) ::
     *   a(licn), rhs(n), w(n), resid, mresid, eps, rmin
      integer irncp, icncp, minirn, minicn, irank
      integer idisp(2)
      integer icn(licn), ikeep(n,5)
      logical abort1, abort2
c common block variables.
c unless otherwise stated common block variables are as in ma28a/ad.
c     those variables referenced by ma28c/cd are mentioned below.
c resid  real/real*8  variable returns maximum residual of
c     equations where pivot was zero.
c mresid  real/real*8 variable used by ma28c/cd to
c     communicate between ma28f/fd and ma30h/hd.
c idisp  integer array  length 2  the same as that used by ma28a/ad.
c     it is unchanged by ma28b/bd.
c
c further information on common block variables can be found in block
c     data or ma28a/ad.
      common /ma28fdc3/ eps, rmin, resid, irncp, icncp, minirn, minicn,
     * irank, abort1, abort2
      common /ma28gdc3/ idisp
      common /ma30hdc3/ mresid
c
c this call performs the solution of the set of equations.
      call ma30cdcopy3(n, icn, a, licn, ikeep, ikeep(1,4), ikeep(1,5),
     * idisp, ikeep(1,2), ikeep(1,3), rhs, w, mtype)
c transfer common block information.
      resid = mresid
      return
      end subroutine ma28cdcopy3
      subroutine ma28ddcopy3(n, a, licn, ivect, jvect,nz,icn,lenr,lenrl,
     * lenoff, ip, iq, iw1, iw, w1, iflag)
      use b2mod_types
      implicit none
c this subroutine need never be called by the user directly.
c     it sorts the matrix as given by the user into the structure of the
c     decomposed form and checks for the presence of duplicate entries or
c     nonzeros lying outside the sparsity pattern of the decomposition
c     it also calculates the largest element in the input matrix.
      integer n, licn, nz, iflag
      real (kind=R8) ::
     *   a(licn), zero, w1, aa
      integer iw(n,2), idisp(2)
      integer icn(licn), ivect(nz), jvect(nz), ip(n), iq(n),
     * lenr(n), iw1(n,3), lenrl(n), lenoff(n)
      integer lp, mp
      logical lblock, grow, blockl
      common /ma28edc3/ lp, mp, lblock, grow
      common /ma28gdc3/ idisp
c
      integer ii, i, jj, j1, j2, iold, inew, jnew, jold,
     * idisp2, iblock, jcomp, idummy, jdummy, midpt
c
      data zero /0.0e0_R8/
      blockl = lenoff(1).ge.0
c iw1(i,3)  is set to the block in which row i lies and the
c     inverse permutations to ip and iq are set in iw1(.,1) and
c     iw1(.,2) resp.
c pointers to beginning of the part of row i in diagonal and
c   off-diagonal blocks are set in iw(i,2) and iw(i,1) resp.
      iblock = 1
      iw(1,1) = 1
      iw(1,2) = idisp(1)
      do 10 i=1,n
        iw1(i,3) = iblock
        if (ip(i).lt.0) iblock = iblock + 1
        ii = iabs(ip(i)+0)
        iw1(ii,1) = i
        jj = iq(i)
        jj = iabs(jj)
        iw1(jj,2) = i
        if (i.eq.1) go to 10
        if (blockl) iw(i,1) = iw(i-1,1) + lenoff(i-1)
        iw(i,2) = iw(i-1,2) + lenr(i-1)
   10 continue
c place each nonzero in turn into its correct location
c    in the a/icn array.
      idisp2 = idisp(2)
      do 170 i=1,nz
c necessary to avoid reference to unassigned element of icn.
        if (i.gt.idisp2) go to 20
        if (icn(i).lt.0) go to 170
   20   iold = ivect(i)
        jold = jvect(i)
        aa = a(i)
c this is a dummy loop for following a chain of interchanges.
c   it will be executed nz times in total.
        do 140 idummy=1,nz
c perform some validity checks on iold and jold.
          if (iold.le.n .and. iold.gt.0 .and. jold.le.n .and.
     *     jold.gt.0) go to 30
          if (lp.ne.0) write (lp,99999) i, a(i), iold, jold
          iflag = -12
          go to 180
   30     inew = iw1(iold,1)
          jnew = iw1(jold,2)
c are we in a valid block and is it diagonal or off-diagonal?
          if (iw1(inew,3)-iw1(jnew,3)) 40, 60, 50
   40     iflag = -13
          if (lp.ne.0) write (lp,99998) iold, jold
          go to 180
   50     j1 = iw(inew,1)
          j2 = j1 + lenoff(inew) - 1
          go to 110
c element is in diagonal block.
   60     j1 = iw(inew,2)
          if (inew.gt.jnew) go to 70
          j2 = j1 + lenr(inew) - 1
          j1 = j1 + lenrl(inew)
          go to 110
   70     j2 = j1 + lenrl(inew)
c binary search of ordered list  .. element in l part of row.
          do 100 jdummy=1,n
            midpt = (j1+j2)/2
            jcomp = iabs(icn(midpt)+0)
            if (jnew-jcomp) 80, 130, 90
   80       j2 = midpt
            go to 100
   90       j1 = midpt
  100     continue
          iflag = -13
          if (lp.ne.0) write (lp,99997) iold, jold
          go to 180
c linear search ... element in l part of row or off-diagonal blocks.
  110     do 120 midpt=j1,j2
            if (iabs(icn(midpt)+0).eq.jnew) go to 130
  120     continue
          iflag = -13
          if (lp.ne.0) write (lp,99997) iold, jold
          go to 180
c equivalent element of icn is in position midpt.
  130     if (icn(midpt).lt.0) go to 160
          if (midpt.gt.nz .or. midpt.le.i) go to 150
          w1 = a(midpt)
          a(midpt) = aa
          aa = w1
          iold = ivect(midpt)
          jold = jvect(midpt)
          icn(midpt) = -icn(midpt)
  140   continue
  150   a(midpt) = aa
        icn(midpt) = -icn(midpt)
        go to 170
  160   a(midpt) = a(midpt) + aa
c set flag for duplicate elements.
        iflag = n + 1
  170 continue
c reset icn array  and zero elements in l/u but not in a.
c also calculate maximum element of a.
  180 w1 = zero
      do 200 i=1,idisp2
        if (icn(i).lt.0) go to 190
        a(i) = zero
        go to 200
  190   icn(i) = -icn(i)
        w1 = max(w1,abs(a(i)))
  200 continue
      return
99999 format (' element ', i6, ' with value ', 1pd22.14, ' has indices',
     * i8, ' ,', i8/36x, 'indices out of range')
99998 format (36x, ' nonzero', i7, ' ,', i6, ' in zero off-diagonal b',
     * 'lock')
99997 format (36x, ' element', i6, ' ,', i6, ' was not in l/u pattern')
      end subroutine ma28ddcopy3
      subroutine ma28idcopy3(n, nz, aorg, irnorg, icnorg, licn, a, icn,
     * ikeep, rhs, x, r, w, mtype, prec, iflag)
      use b2mod_types
      implicit none
      external bd1copy3
c this subroutine uses the factors from an earlier call to ma28a/ad
c     or ma28b/bd to solve the system of equations with iterative
c     refinement.
c
c the parameters are...
c
c n is equal to the order of the matrix. it is not altered by the
c     subroutine.
c nz is equal to the number of entries in the original matrix.  it is
c     not altered by the subroutine.
c for this entry the original matrix must have been saved in
c     aorg,irnorg,icnorg where entry aorg(k) is in row irnorg(k) and
c     column icnorg(k), k=1,...nz.  information about the factors of a
c     is communicated to this subroutine via the parameters licn, a, icn
c     and ikeep where:
c aorg is an array of length nz.  not altered by ma28i/id.
c irnorg is an array of length nz.  not altered by ma28i/id.
c icnorg is an array of length nz.  not altered by ma28i/id.
c licn is equal to the length of arrays a and icn. it is not altered by
c     the subroutine.
c a is an array of length licn. it must be unchanged since the last call
c     to ma28a/ad or ma28b/bd. it is not altered by the subroutine.
c icn, ikeep are the arrays (of lengths licn and 5*n, respectively) of
c     the same names as in the previous all to ma28a/ad. they should be
c     unchanged since this earlier call and    they are not altered by
c     ma28i/id.
c the  other parameters are as follows:
c rhs is an array of length n. the user must set rhs(i) to contain the
c     value of the i-th component of the right hand side. it is not
c     altered by ma28i/id.
c x is an array of length n. if an initial guess of the solution is
c     given (istart equal to 1), then the user must set x(i) to contain
c     the value of the i-th component of the estimated solution.  on
c     exit, x(i) contains the i-th component of the solution vector.
c r is an array of length n. it need not be set on entry.  on exit, r(i)
c     contains the i-th component of an estimate of the error if maxit
c     is greater than 0.
c w is an array of length n. it is used as workspace by ma28i/id.
c mtype must be set to determine whether ma28i/id will solve a*x=rhs
c     (mtype equal to 1) or at*x=rhs (mtype ne 1, zero say). it is not
c     altered by ma28i/id.
c prec should be set by the user to the relative accuracy required. the
c     iterative refinement will terminate if the magnitude of the
c     largest component of the estimated error relative to the largest
c     component in the solution is less than prec.  it is not altered by
c     ma28i/id.
c iflag is a diagnostic flag which will be set to zero on successful
c     exit from ma28i/id, otherwise it will have a nonzero value. the
c     nonzero value iflag can have on exit from ma28i/id are ...
c     -16    indicating that more than maxit iterations are required.
c     -17    indicating that convergence was too slow.
c
      integer n, nz, licn, mtype, iflag
      integer icnorg(nz), irnorg(nz), ikeep(n,5), icn(licn)
      real (kind=R8) ::
     *   a(licn), aorg(nz), rhs(n), r(n), x(n), w(n), prec
c
c private and common variables.
c see block data for comments on variables in common.
c see comments in code for use of private variables.
c
      logical lblock, grow, lbig
      integer lp, mp, ndrop, maxit, noiter, nsrch, istart
      real (kind=R8) ::
     *   tol, themax, big, dxmax, errmax, dres,
     * cgce, d, dd, conver, zero
      common /ma28edc3/ lp, mp, lblock, grow
      common /ma28hdc3/ tol, themax, big, dxmax, errmax, dres, cgce,
     * ndrop, maxit, noiter, nsrch, istart, lbig
c
      integer i, iterat, nrow, ncol
c
      data zero /0.0e0_R8/
c
c   initialisation of noiter, errmax and iflag.
c
      noiter = 0
      errmax = zero
      iflag = 0
c
c  jump if a starting vector has been supplied by the user.
c
      if (istart.eq.1) go to 20
c
c  make a copy of the right-hand side vector.
c
      do 10 i=1,n
        x(i) = rhs(i)
   10 continue
c
c  find the first solution.
c
      call ma28cdcopy3(n, a, licn, icn, ikeep, x, w, mtype)
c
c  stop the computations if   maxit=0.
c
   20 if (maxit.eq.0) go to 160
c
c  calculate the max-norm of the first solution.
c
      dd = 0.0e0_R8                                                   !sv 21.11.00
      do 30 i=1,n
        dd = max(dd,abs(x(i)))
   30 continue
      dxmax = dd
c
c  begin the iterative process.
c
      do 120 iterat=1,maxit
        d = dd
c
c  calculate the residual vector.
c
        do 40 i=1,n
          r(i) = rhs(i)
   40   continue
        if (mtype.eq.1) go to 60
        do 50 i=1,nz
          nrow = irnorg(i)
          ncol = icnorg(i)
          r(ncol) = r(ncol) - aorg(i)*x(nrow)
   50   continue
        go to 80
c mtype=1.
   60   do 70 i=1,nz
          nrow = irnorg(i)
          ncol = icnorg(i)
          r(nrow) = r(nrow) - aorg(i)*x(ncol)
   70   continue
   80   dres = 0.0e0_R8                                               !sv 21.11.00
c
c  find the max-norm of the residual vector.
c
        do 90 i=1,n
          dres = max(dres,abs(r(i)))
   90   continue
c
c  stop the calculations if the max-norm of
c  the residual vector is zero.
c
        if (dres.eq.0.0e0_R8) go to 150                               !sv 21.11.00
c
c  calculate the correction vector.
c
        noiter = noiter + 1
        call ma28cdcopy3(n, a, licn, icn, ikeep, r, w, mtype)
c
c  find the max-norm of the correction vector.
c
        dd = 0.0e0_R8                                                 !sv 21.11.00
        do 100 i=1,n
          dd = max(dd,abs(r(i)))
  100   continue
c
c  check the convergence.
c
        if (dd.gt.d*cgce .and. iterat.ge.2) go to 130
        if (dxmax*10.0e0_R8+dd.eq.dxmax*10.0e0_R8) go to 140          !sv 21.11.00
c
c  attempt to improve the solution.
c
        dxmax = 0.0e0_R8                                              !sv 21.11.00
        do 110 i=1,n
          x(i) = x(i) + r(i)
          dxmax = max(dxmax,abs(x(i)))
  110   continue
c
c  check the stopping criterion.
c
        if (dd.lt.prec*dxmax) go to 140
  120 continue
c more than maxit iterations required.
      iflag = -16
      write (lp,99999) iflag, maxit
      go to 140
c convergence rate unacceptably slow.
  130 iflag = -17
      conver = dd/d
      write (lp,99998) iflag, conver, cgce
c
c  the iterative process is terminated.
c
  140 errmax = dd
  150 continue
  160 return
99999 format (' error return from ma28i/idc with iflag = ', i3/' more ',
     * 'than', i5, ' iterations required')
99998 format (' error return from ma28ic with iflag = ', i3/' converge',
     * 'nce rate of ', 1pe9.2, ' too slow'/' maximum acceptable rate',
     * ' set to ', 1pe9.2)
      end subroutine ma28idcopy3
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias ma30adcopy3
      subroutine ma30adcopy3(nn, icn, a, licn, lenr, lenrl,idisp, ip,iq,
     * irn, lirn, lenc, ifirst, lastr, nextr, lastc, nextc, iptr, ipc,
     * u, iflag)
      use b2mod_types
      implicit none
c if  the user requires a more convenient data interface then the ma28
c     package should be used.  the ma28 subroutines call the ma30
c     subroutines after checking the input data from the user and optionally
c     using mc23a/ad to permute the matrix to block triangular form.
c this package of subroutines (ma30a/ad, ma30b/bd, ma30c/cd and
c     ma30d/dd) performs operations pertinent to the solution of a
c     general sparse n by n system of linear equations (i.e. solve
c     ax=b). structually singular matrices are permitted including
c     those with row or columns consisting entirely of zeros (i.e.
c     including rectangular matrices).  it is assumed that the
c     nonzeros of the matrix a do not differ widely in size.  if
c     necessary a prior call of the scaling subroutine mc19a/ad may be
c     made.
c a discussion of the design of these subroutines is given by duff and
c     reid (acm trans math software 5 pp 18-35,1979 (css 48)) while
c     fuller details of the implementation are given in duff (harwell
c     report aere-r 8730,1977).  the additional pivoting option in
c     ma30a/ad and the use of drop tolerances (see common block
c     ma30i/id) were added to the package after joint work with reid,
c     schaumburg, wasniewski and zlatev (duff, reid, schaumburg,
c     wasniewski and zlatev, harwell report css 135, 1983).
c
c ma30a/ad performs the lu decomposition of the diagonal blocks of the
c     permutation paq of a sparse matrix a, where input permutations
c     p1 and q1 are used to define the diagonal blocks.  there may be
c     nonzeros in the off-diagonal blocks but they are unaffected by
c     ma30a/ad. p and p1 differ only within blocks as do q and q1. the
c     permutations p1 and q1 may be found by calling mc23a/ad or the
c     matrix may be treated as a single block by using p1=q1=i. the
c     matrix nonzeros should be held compactly by rows, although it
c     should be noted that the user can supply the matrix by columns
c     to get the lu decomposition of a transpose.
c
c the parameters are...
c this description should also be consulted for further information on
c     most of the parameters of ma30b/bd and ma30c/cd.
c
c n  is an integer variable which must be set by the user to the order
c     of the matrix.  it is not altered by ma30a/ad.
c icn is an integer array of length licn. positions idisp(2) to
c     licn must be set by the user to contain the column indices of
c     the nonzeros in the diagonal blocks of p1*a*q1. those belonging
c     to a single row must be contiguous but the ordering of column
c     indices with each row is unimportant. the nonzeros of row i
c     precede those of row i+1,i=1,...,n-1 and no wasted space is
c     allowed between the rows.  on output the column indices of the
c     lu decomposition of paq are held in positions idisp(1) to
c     idisp(2), the rows are in pivotal order, and the column indices
c     of the l part of each row are in pivotal order and precede those
c     of u. again there is no wasted space either within a row or
c     between the rows. icn(1) to icn(idisp(1)-1), are neither
c     required nor altered. if mc23a/ad been called, these will hold
c     information about the off-diagonal blocks.
c a is a real/real*8 array of length licn whose entries
c     idisp(2) to licn must be set by the user to the  values of the
c     nonzero entries of the matrix in the order indicated by  icn.
c     on output a will hold the lu factors of the matrix where again
c     the position in the matrix is determined by the corresponding
c     values in icn. a(1) to a(idisp(1)-1) are neither required nor
c     altered.
c licn  is an integer variable which must be set by the user to the
c     length of arrays icn and a. it must be big enough for a and icn
c     to hold all the nonzeros of l and u and leave some "elbow
c     room".  it is possible to calculate a minimum value for licn by
c     a preliminary run of ma30a/ad. the adequacy of the elbow room
c     can be judged by the size of the common block variable icncp. it
c     is not altered by ma30a/ad.
c lenr  is an integer array of length n.  on input, lenr(i) should
c     equal the number of nonzeros in row i, i=1,...,n of the
c     diagonal blocks of p1*a*q1. on output, lenr(i) will equal the
c     total number of nonzeros in row i of l and row i of u.
c lenrl  is an integer array of length n. on output from ma30a/ad,
c     lenrl(i) will hold the number of nonzeros in row i of l.
c idisp  is an integer array of length 2. the user should set idisp(1)
c     to be the first available position in a/icn for the lu
c     decomposition while idisp(2) is set to the position in a/icn of
c     the first nonzero in the diagonal blocks of p1*a*q1. on output,
c     idisp(1) will be unaltered while idisp(2) will be set to the
c     position in a/icn of the last nonzero of the lu decomposition.
c ip  is an integer array of length n which holds a permutation of
c     the integers 1 to n.  on input to ma30a/ad, the absolute value of
c     ip(i) must be set to the row of a which is row i of p1*a*q1. a
c     negative value for ip(i) indicates that row i is at the end of a
c     diagonal block.  on output from ma30a/ad, ip(i) indicates the row
c     of a which is the i-th row in paq. ip(i) will still be negative
c     for the last row of each block (except the last).
c iq is an integer array of length n which again holds a
c     permutation of the integers 1 to n.  on input to ma30a/ad, iq(j)
c     must be set to the column of a which is column j of p1*a*q1. on
c     output from ma30a/ad, the absolute value of iq(j) indicates the
c     column of a which is the j-th in paq.  for rows, i say, in which
c     structural or numerical singularity is detected iq(i) is
c     negated.
c irn  is an integer array of length lirn used as workspace by
c     ma30a/ad.
c lirn  is an integer variable. it should be greater than the
c     largest number of nonzeros in a diagonal block of p1*a*q1 but
c     need not be as large as licn. it is the length of array irn and
c     should be large enough to hold the active part of any block,
c     plus some "elbow room", the  a posteriori  adequacy of which can
c     be estimated by examining the size of common block variable
c     irncp.
c lenc,ifirst,lastr,nextr,lastc,nextc are all integer arrays of
c     length n which are used as workspace by ma30a/ad.  if nsrch is
c     set to a value less than or equal to n, then arrays lastc and
c     nextc are not referenced by ma30a/ad and so can be dummied in
c     the call to ma30a/ad.
c iptr,ipc are integer arrays of length n which are used as workspace
c     by ma30a/ad.
c u  is a real/real*8 variable which should be set by the
c     user to a value between 0. and 1.0. if less than zero it is
c     reset to zero and if its value is 1.0 or greater it is reset to
c     0.9999 (0.999999999 in d version).  it determines the balance
c     between pivoting for sparsity and for stability, values near
c     zero emphasising sparsity and values near one emphasising
c     stability. we recommend u=0.1 as a posible first trial value.
c     the stability can be judged by a later call to mc24a/ad or by
c     setting lbig to .true.
c iflag  is an integer variable. it will have a non-negative value if
c     ma30a/ad is successful. negative values indicate error
c     conditions while positive values indicate that the matrix has
c     been successfully decomposed but is singular. for each nonzero
c     value, an appropriate message is output on unit lp.  possible
c     nonzero values for iflag are ...
c
c -1  the matrix is structually singular with rank given by irank in
c     common block ma30f/fd.
c +1  if, however, the user wants the lu decomposition of a
c     structurally singular matrix and sets the common block variable
c     abort1 to .false., then, in the event of singularity and a
c     successful decomposition, iflag is returned with the value +1
c     and no message is output.
c -2  the matrix is numerically singular (it may also be structually
c     singular) with estimated rank given by irank in common block
c     ma30f/fd.
c +2  the  user can choose to continue the decomposition even when a
c     zero pivot is encountered by setting common block variable
c     abort2 to .false.  if a singularity is encountered, iflag will
c     then return with a value of +2, and no message is output if the
c     decomposition has been completed successfully.
c -3  lirn has not been large enough to continue with the
c     decomposition.  if the stage was zero then common block variable
c     minirn gives the length sufficient to start the decomposition on
c     this block.  for a successful decomposition on this block the
c     user should make lirn slightly (say about n/2) greater than this
c     value.
c -4  licn not large enough to continue with the decomposition.
c -5  the decomposition has been completed but some of the lu factors
c     have been discarded to create enough room in a/icn to continue
c     the decomposition. the variable minicn in common block ma30f/fd
c     then gives the size that licn should be to enable the
c     factorisation to be successful.  if the user sets common block
c     variable abort3 to .true., then the subroutine will exit
c     immediately instead of destroying any factors and continuing.
c -6  both licn and lirn are too small. termination has been caused by
c     lack of space in irn (see error iflag= -3), but already some of
c     the lu factors in a/icn have been lost (see error iflag= -5).
c     minicn gives the minimum amount of space required in a/icn for
c     decomposition up to this point.
c
      integer nn, licn, lirn, iflag
      real (kind=R8) ::
     *   a(licn), u, au, umax, amax, zero, pivrat, pivr,
     * tol, big, anew, aanew, scale
      integer iptr(nn), pivot, pivend, dispc, oldpiv, oldend, pivrow,
     * rowi, ipc(nn), idisp(2), colupd
      integer icn(licn), lenr(nn), lenrl(nn), ip(nn), iq(nn),
     * lenc(nn), irn(lirn), ifirst(nn), lastr(nn), nextr(nn),
     * lastc(nn), nextc(nn)
      integer lp
      integer irncp, icncp, irank, minirn, minicn
      integer ndrop, nsrch
      integer lpiv(10),lnpiv(10),mapiv,manpiv,iavpiv,ianpiv,kountl
      logical abort1, abort2, abort3, lbig
c for comments of common block variables see block data subprogram.
      common /ma30edc3/ lp, abort1, abort2, abort3
      common /ma30fdc3/ irncp, icncp, irank, minirn, minicn
      common /ma30idc3/ tol, big, ndrop, nsrch, lbig
      common /lpivotc3/ lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
c
      integer i, j, k, l, n, i1, i2, ii, j1, j2, jj, kk, ll, iii,
     * lc, lr, nc, nr, nz, nz2, nzpc,
     * ipiv, jpiv, ipos, jpos, ibeg, jbeg, ilast, iend, jend,
     * iop, ifill, irows, istart, idrop, ifir,
     * isw, isw1, isrch, itop, iactiv, indrow, ising, idispc,
     * ijfir, ijpos, ijp1, idummy, jdummy, jnpos, jnew, jold, jval,
     * jcount, jdiff, jroom, jmore, jzer, jzero,
     * kdrop, lenpp, lenpiv, msrch, morei,
     * nblock, num, nnm1, nzcol, nzmin, nzrow
      integer(8) :: jcost, kcost                           !srv 26.02.11 12.03.15
#ifdef DBG
      real (kind=R8) :: product
      real (kind=R8), save :: cutlo
#endif
c
      data umax/.999999999e0_R8/
      data zero /0.0e0_R8/
#ifdef DBG
      cutlo = log(tiny(1.0_R8))
#endif
      msrch = nsrch
      ndrop = 0
      do 1272 kk=1,10
        lnpiv(kk)=0
        lpiv(kk)=0
 1272 continue
      mapiv = 0
      manpiv = 0
      iavpiv = 0
      ianpiv = 0
      kountl = 0
      minirn = 0
      minicn = idisp(1) - 1
      morei = 0
      irank = nn
      irncp = 0
      icncp = 0
      iflag = 0
c reset u if necessary.
      u = min(u,umax)
c ibeg is the position of the next pivot row after elimination step
c     using it.
      u = max(u,zero)
      ibeg = idisp(1)
c iactiv is the position of the first entry in the active part of a/icn.
      iactiv = idisp(2)
c nzrow is current number of nonzeros in active and unprocessed part
c     of row file icn.
      nzrow = licn - iactiv + 1
      minicn = nzrow + minicn
c
c count the number of diagonal blocks and set up pointers to the
c     beginnings of the rows.
c num is the number of diagonal blocks.
      num = 1
      iptr(1) = iactiv
      if (nn.eq.1) go to 20
      nnm1 = nn - 1
      do 10 i=1,nnm1
        if (ip(i).lt.0) num = num + 1
        iptr(i+1) = iptr(i) + lenr(i)
   10 continue
c ilast is the last row in the previous block.
   20 ilast = 0
c
c ***********************************************
c ****    lu decomposition of block nblock   ****
c ***********************************************
c
c each pass through this loop performs lu decomposition on one
c     of the diagonal blocks.
      do 1000 nblock=1,num
        istart = ilast + 1
        do 30 irows=istart,nn
          if (ip(irows).lt.0) go to 40
   30   continue
        irows = nn
   40   ilast = irows
c n is the number of rows in the current block.
c istart is the index of the first row in the current block.
c ilast is the index of the last row in the current block.
c iactiv is the position of the first entry in the block.
c itop is the position of the last entry in the block.
        n = ilast - istart + 1
        if (n.ne.1) go to 90
c
c code for dealing with 1x1 block.
        lenrl(ilast) = 0
        ising = istart
        if (lenr(ilast).ne.0) go to 50
c block is structurally singular.
        irank = irank - 1
        ising = -ising
        if (iflag.ne.2 .and. iflag.ne.-5) iflag = 1
        if (.not.abort1) go to 80
        idisp(2) = iactiv
        iflag = -1
        if (lp.ne.0) write (lp,99999)
c     return
        go to 1120
   50   scale = abs(a(iactiv))
        if (scale.eq.zero) go to 60
        if (lbig) big = max(big,scale)
        go to 70
   60   ising = -ising
        irank = irank - 1
        iptr(ilast) = 0
        if (iflag.ne.-5) iflag = 2
        if (.not.abort2) go to 70
        idisp(2) = iactiv
        iflag = -2
        if (lp.ne.0) write (lp,99998)
        go to 1120
   70   a(ibeg) = a(iactiv)
        icn(ibeg) = icn(iactiv)
        iactiv = iactiv + 1
        iptr(istart) = 0
        ibeg = ibeg + 1
        nzrow = nzrow - 1
   80   lastr(istart) = istart
        ipc(istart) = -ising
        go to 1000
c
c non-trivial block.
   90   itop = licn
        if (ilast.ne.nn) itop = iptr(ilast+1) - 1
c
c set up column oriented storage.
        do 100 i=istart,ilast
          lenrl(i) = 0
          lenc(i) = 0
  100   continue
        if (itop-iactiv.lt.lirn) go to 110
        minirn = itop - iactiv + 1
        pivot = istart - 1
        go to 1100
c
c calculate column counts.
  110   do 120 ii=iactiv,itop
          i = icn(ii)
          lenc(i) = lenc(i) + 1
  120   continue
c set up column pointers so that ipc(j) points to position after end
c     of column j in column file.
        ipc(ilast) = lirn + 1
        j1 = istart + 1
        do 130 jj=j1,ilast
          j = ilast - jj + j1 - 1
          ipc(j) = ipc(j+1) - lenc(j+1)
  130   continue
        do 150 indrow=istart,ilast
          j1 = iptr(indrow)
          j2 = j1 + lenr(indrow) - 1
          if (j1.gt.j2) go to 150
          do 140 jj=j1,j2
            j = icn(jj)
            ipos = ipc(j) - 1
            irn(ipos) = indrow
            ipc(j) = ipos
  140     continue
  150   continue
c dispc is the lowest indexed active location in the column file.
        dispc = ipc(istart)
        nzcol = lirn - dispc + 1
        minirn = max0(nzcol,minirn)
        nzmin = 1
c
c initialise array ifirst.  ifirst(i) = +/- k indicates that row/col
c     k has i nonzeros.  if ifirst(i) = 0, there is no row or column
c     with i nonzeros.
        do 160 i=1,n
          ifirst(i) = 0
  160   continue
c
c compute ordering of row and column counts.
c first run through columns (from column n to column 1).
        do 180 jj=istart,ilast
          j = ilast - jj + istart
          nz = lenc(j)
          if (nz.ne.0) go to 170
          ipc(j) = 0
          go to 180
  170     if (nsrch.le.nn) go to 180
          isw = ifirst(nz)
          ifirst(nz) = -j
          lastc(j) = 0
          nextc(j) = -isw
          isw1 = iabs(isw)
          if (isw.ne.0) lastc(isw1) = j
  180   continue
c now run through rows (again from n to 1).
        do 210 ii=istart,ilast
          i = ilast - ii + istart
          nz = lenr(i)
          if (nz.ne.0) go to 190
          iptr(i) = 0
          lastr(i) = 0
          go to 210
  190     isw = ifirst(nz)
          ifirst(nz) = i
          if (isw.gt.0) go to 200
          nextr(i) = 0
          lastr(i) = isw
          go to 210
  200     nextr(i) = isw
          lastr(i) = lastr(isw)
          lastr(isw) = i
  210   continue
c
c
c **********************************************
c ****    start of main elimination loop    ****
c **********************************************
        do 980 pivot=istart,ilast
c
c first find the pivot using markowitz criterion with stability
c     control.
c jcost is the markowitz cost of the best pivot so far,.. this
c     pivot is in row ipiv and column jpiv.
          nz2 = nzmin
          jcost = int(n,8)*int(n,8)    !srv 26.02.11 12.03.15
c
c examine rows/columns in order of ascending count.
          do 340 l=1,2
            pivrat = zero
            isrch = 1
            ll = l
c a pass with l equal to 2 is only performed in the case of singularity.
            do 330 nz=nz2,n
              if (jcost.le.int((nz-1),8)**2) go to 420               !srv 26.02.11 12.03.15
              ijfir = ifirst(nz)
              if (ijfir) 230, 220, 240
  220         if (ll.eq.1) nzmin = nz + 1
              go to 330
  230         ll = 2
              ijfir = -ijfir
              go to 290
  240         ll = 2
c scan rows with nz nonzeros.
              do 270 idummy=1,n
                if (jcost.le.int((nz-1),8)**2) go to 420             !srv 26.02.11 12.03.15
                if (isrch.gt.msrch) go to 420
                if (ijfir.eq.0) go to 280
c row ijfir is now examined.
                i = ijfir
                ijfir = nextr(i)
c first calculate multiplier threshold level.
                amax = zero
                j1 = iptr(i) + lenrl(i)
                j2 = iptr(i) + lenr(i) - 1
                do 250 jj=j1,j2
                  amax = max(amax,abs(a(jj)))
  250           continue
                au = amax*u
                isrch = isrch + 1
c scan row for possible pivots
                do 260 jj=j1,j2
                  if (abs(a(jj)).le.au .and. l.eq.1) go to 260
                  j = icn(jj)
                  kcost = int((nz-1),8)*int((lenc(j)-1),8)           !srv 26.02.11 12.03.15
                  if (kcost.gt.jcost) go to 260
                  pivr = zero
                  if (amax.ne.zero) pivr = abs(a(jj))/amax
                  if (kcost.eq.jcost .and. (pivr.le.pivrat .or.
     *             nsrch.gt.nn+1)) go to 260
c best pivot so far is found.
                  jcost = kcost
                  ijpos = jj
                  ipiv = i
                  jpiv = j
                  if (msrch.gt.nn+1 .and. jcost.le.int((nz-1),8)**2)
     *               go to 420                                        !srv 26.02.11 12.03.15
                  pivrat = pivr
  260           continue
  270         continue
c
c columns with nz nonzeros now examined.
  280         ijfir = ifirst(nz)
              ijfir = -lastr(ijfir)
  290         if (jcost.le.int(nz,8)*int((nz-1),8)) go to 420         !srv 26.02.11 12.03.15
              if (msrch.le.nn) go to 330
              do 320 idummy=1,n
                if (ijfir.eq.0) go to 330
                j = ijfir
                ijfir = nextc(ijfir)
                i1 = ipc(j)
                i2 = i1 + nz - 1
c scan column j.
                do 310 ii=i1,i2
                  i = irn(ii)
                  kcost = int((nz-1),8)*int((lenr(i)-lenrl(i)-1),8)   !srv 26.02.11 12.03.15
                  if (kcost.ge.jcost) go to 310
c pivot has best markowitz count so far ... now check its
c     suitability on numeric grounds by examining the other nonzeros
c     in its row.
                  j1 = iptr(i) + lenrl(i)
                  j2 = iptr(i) + lenr(i) - 1
c we need a stability check on singleton columns because of possible
c     problems with underdetermined systems.
                  amax = zero
                  do 300 jj=j1,j2
                    amax = max(amax,abs(a(jj)))
                    if (icn(jj).eq.j) jpos = jj
  300             continue
                  if (abs(a(jpos)).le.amax*u .and. l.eq.1) go to 310
                  jcost = kcost
                  ipiv = i
                  jpiv = j
                  ijpos = jpos
                  if (amax.ne.zero) pivrat = abs(a(jpos))/amax
                  if (jcost.le.int(nz,8)*int((nz-1),8)) go to 420     !srv 26.02.11 12.03.15
  310           continue
c
  320         continue
c
  330       continue
c in the event of singularity, we must make sure all rows and columns
c are tested.
            msrch = n
c
c matrix is numerically or structurally singular  ... which it is will
c     be diagnosed later.
            irank = irank - 1
  340     continue
c assign rest of rows and columns to ordering array.
c matrix is structurally singular.
          if (iflag.ne.2 .and. iflag.ne.-5) iflag = 1
          irank = irank - ilast + pivot + 1
          if (.not.abort1) go to 350
          idisp(2) = iactiv
          iflag = -1
          if (lp.ne.0) write (lp,99999)
          go to 1120
  350     k = pivot - 1
          do 390 i=istart,ilast
            if (lastr(i).ne.0) go to 390
            k = k + 1
            lastr(i) = k
            if (lenrl(i).eq.0) go to 380
            minicn = max0(minicn,nzrow+ibeg-1+morei+lenrl(i))
            if (iactiv-ibeg.ge.lenrl(i)) go to 360
            call ma30ddcopy3(a, icn, iptr(istart), n,iactiv,itop,.true.)
c check now to see if ma30d/dd has created enough available space.
            if (iactiv-ibeg.ge.lenrl(i)) go to 360
c create more space by destroying previously created lu factors.
            morei = morei + ibeg - idisp(1)
            ibeg = idisp(1)
            if (lp.ne.0) write (lp,99997)
            iflag = -5
            if (abort3) go to 1090
  360       j1 = iptr(i)
            j2 = j1 + lenrl(i) - 1
            iptr(i) = 0
            do 370 jj=j1,j2
              a(ibeg) = a(jj)
              icn(ibeg) = icn(jj)
              icn(jj) = 0
              ibeg = ibeg + 1
  370       continue
            nzrow = nzrow - lenrl(i)
  380       if (k.eq.ilast) go to 400
  390     continue
  400     k = pivot - 1
          do 410 i=istart,ilast
            if (ipc(i).ne.0) go to 410
            k = k + 1
            ipc(i) = k
            if (k.eq.ilast) go to 990
  410     continue
c
c the pivot has now been found in position (ipiv,jpiv) in location
c     ijpos in row file.
c update column and row ordering arrays to correspond with removal
c     of the active part of the matrix.
  420     ising = pivot
          if (a(ijpos).ne.zero) go to 430
c numerical singularity is recorded here.
          ising = -ising
          if (iflag.ne.-5) iflag = 2
          if (.not.abort2) go to 430
          idisp(2) = iactiv
          iflag = -2
          if (lp.ne.0) write (lp,99998)
          go to 1120
  430     oldpiv = iptr(ipiv) + lenrl(ipiv)
          oldend = iptr(ipiv) + lenr(ipiv) - 1
c changes to column ordering.
          if (nsrch.le.nn) go to 460
          colupd = nn + 1
          lenpp = oldend-oldpiv+1
          if (lenpp.lt.4) lpiv(1) = lpiv(1) + 1
          if (lenpp.ge.4 .and. lenpp.le.6) lpiv(2) = lpiv(2) + 1
          if (lenpp.ge.7 .and. lenpp.le.10) lpiv(3) = lpiv(3) + 1
          if (lenpp.ge.11 .and. lenpp.le.15) lpiv(4) = lpiv(4) + 1
          if (lenpp.ge.16 .and. lenpp.le.20) lpiv(5) = lpiv(5) + 1
          if (lenpp.ge.21 .and. lenpp.le.30) lpiv(6) = lpiv(6) + 1
          if (lenpp.ge.31 .and. lenpp.le.50) lpiv(7) = lpiv(7) + 1
          if (lenpp.ge.51 .and. lenpp.le.70) lpiv(8) = lpiv(8) + 1
          if (lenpp.ge.71 .and. lenpp.le.100) lpiv(9) = lpiv(9) + 1
          if (lenpp.ge.101) lpiv(10) = lpiv(10) + 1
          mapiv = max0(mapiv,lenpp)
          iavpiv = iavpiv + lenpp
          do 450 jj=oldpiv,oldend
            j = icn(jj)
            lc = lastc(j)
            nc = nextc(j)
            nextc(j) = -colupd
            if (jj.ne.ijpos) colupd = j
            if (nc.ne.0) lastc(nc) = lc
            if (lc.eq.0) go to 440
            nextc(lc) = nc
            go to 450
  440       nz = lenc(j)
            isw = ifirst(nz)
            if (isw.gt.0) lastr(isw) = -nc
            if (isw.lt.0) ifirst(nz) = -nc
  450     continue
c changes to row ordering.
  460     i1 = ipc(jpiv)
          i2 = i1 + lenc(jpiv) - 1
          do 480 ii=i1,i2
            i = irn(ii)
            lr = lastr(i)
            nr = nextr(i)
            if (nr.ne.0) lastr(nr) = lr
            if (lr.le.0) go to 470
            nextr(lr) = nr
            go to 480
  470       nz = lenr(i) - lenrl(i)
            if (nr.ne.0) ifirst(nz) = nr
            if (nr.eq.0) ifirst(nz) = lr
  480     continue
c
c move pivot to position lenrl+1 in pivot row and move pivot row
c     to the beginning of the available storage.
c the l part and the pivot in the old copy of the pivot row is
c     nullified while, in the strictly upper triangular part, the
c     column indices, j say, are overwritten by the corresponding
c     entry of iq (iq(j)) and iq(j) is set to the negative of the
c     displacement of the column index from the pivot entry.
          if (oldpiv.eq.ijpos) go to 490
          au = a(oldpiv)
          a(oldpiv) = a(ijpos)
          a(ijpos) = au
          icn(ijpos) = icn(oldpiv)
          icn(oldpiv) = jpiv
c check to see if there is space immediately available in a/icn to
c     hold new copy of pivot row.
  490     minicn = max0(minicn,nzrow+ibeg-1+morei+lenr(ipiv))
          if (iactiv-ibeg.ge.lenr(ipiv)) go to 500
          call ma30ddcopy3(a, icn, iptr(istart), n, iactiv,itop,.true.)
          oldpiv = iptr(ipiv) + lenrl(ipiv)
          oldend = iptr(ipiv) + lenr(ipiv) - 1
c check now to see if ma30d/dd has created enough available space.
          if (iactiv-ibeg.ge.lenr(ipiv)) go to 500
c create more space by destroying previously created lu factors.
          morei = morei + ibeg - idisp(1)
          ibeg = idisp(1)
          if (lp.ne.0) write (lp,99997)
          iflag = -5
          if (abort3) go to 1090
          if (iactiv-ibeg.ge.lenr(ipiv)) go to 500
c there is still not enough room in a/icn.
          iflag = -4
          go to 1090
c copy pivot row and set up iq array.
  500     ijpos = 0
          j1 = iptr(ipiv)
c
          do 530 jj=j1,oldend
            a(ibeg) = a(jj)
            icn(ibeg) = icn(jj)
            if (ijpos.ne.0) go to 510
            if (icn(jj).eq.jpiv) ijpos = ibeg
            icn(jj) = 0
            go to 520
  510       k = ibeg - ijpos
            j = icn(jj)
            icn(jj) = iq(j)
            iq(j) = -k
  520       ibeg = ibeg + 1
  530     continue
c
          ijp1 = ijpos + 1
          pivend = ibeg - 1
          lenpiv = pivend - ijpos
          nzrow = nzrow - lenrl(ipiv) - 1
          iptr(ipiv) = oldpiv + 1
          if (lenpiv.eq.0) iptr(ipiv) = 0
c
c remove pivot row (including pivot) from column oriented file.
          do 560 jj=ijpos,pivend
            j = icn(jj)
            i1 = ipc(j)
            lenc(j) = lenc(j) - 1
c i2 is last position in new column.
            i2 = ipc(j) + lenc(j) - 1
            if (i2.lt.i1) go to 550
            do 540 ii=i1,i2
              if (irn(ii).ne.ipiv) go to 540
              irn(ii) = irn(i2+1)
              go to 550
  540       continue
  550       irn(i2+1) = 0
  560     continue
          nzcol = nzcol - lenpiv - 1
c
c go down the pivot column and for each row with a nonzero add
c     the appropriate multiple of the pivot row to it.
c we loop on the number of nonzeros in the pivot column since
c     ma30d/dd may change its actual position.
c
          nzpc = lenc(jpiv)
          if (nzpc.eq.0) go to 900
          do 840 iii=1,nzpc
            ii = ipc(jpiv) + iii - 1
            i = irn(ii)
c search row i for nonzero to be eliminated, calculate multiplier,
c     and place it in position lenrl+1 in its row.
c  idrop is the number of nonzero entries dropped from row    i
c        because these fall beneath tolerance level.
c
            idrop = 0
            j1 = iptr(i) + lenrl(i)
            iend = iptr(i) + lenr(i) - 1
            do 570 jj=j1,iend
              if (icn(jj).ne.jpiv) go to 570
c if pivot is zero, rest of column is and so multiplier is zero.
              au = zero
              if (a(ijpos).ne.zero) au = -a(jj)/a(ijpos)
              if (lbig) big = max(big,abs(au))
              a(jj) = a(j1)
              a(j1) = au
              icn(jj) = icn(j1)
              icn(j1) = jpiv
              lenrl(i) = lenrl(i) + 1
              go to 580
  570       continue
c jump if pivot row is a singleton.
  580       if (lenpiv.eq.0) go to 840
c now perform necessary operations on rest of non-pivot row i.
            rowi = j1 + 1
            iop = 0
c jump if all the pivot row causes fill-in.
            if (rowi.gt.iend) go to 650
c perform operations on current nonzeros in row i.
c innermost loop.
            lenpp = iend-rowi+1
            if (lenpp.lt.4) lnpiv(1) = lnpiv(1) + 1
            if (lenpp.ge.4 .and. lenpp.le.6) lnpiv(2) = lnpiv(2) + 1
            if (lenpp.ge.7 .and. lenpp.le.10) lnpiv(3) = lnpiv(3) + 1
            if (lenpp.ge.11 .and. lenpp.le.15) lnpiv(4) = lnpiv(4) + 1
            if (lenpp.ge.16 .and. lenpp.le.20) lnpiv(5) = lnpiv(5) + 1
            if (lenpp.ge.21 .and. lenpp.le.30) lnpiv(6) = lnpiv(6) + 1
            if (lenpp.ge.31 .and. lenpp.le.50) lnpiv(7) = lnpiv(7) + 1
            if (lenpp.ge.51 .and. lenpp.le.70) lnpiv(8) = lnpiv(8) + 1
            if (lenpp.ge.71 .and. lenpp.le.100) lnpiv(9) = lnpiv(9) + 1
            if (lenpp.ge.101) lnpiv(10) = lnpiv(10) + 1
            manpiv = max0(manpiv,lenpp)
            ianpiv = ianpiv + lenpp
            kountl = kountl + 1
            do 590 jj=rowi,iend
              j = icn(jj)
              if (iq(j).gt.0) go to 590
              iop = iop + 1
              pivrow = ijpos - iq(j)
#ifdef DBG
              if ((au.ge.tiny(1.0_R8).or.au.le.-tiny(1.0_R8)).and.
     &            (a(pivrow).ge. tiny(1.0_R8).or.
     &             a(pivrow).le.-tiny(1.0_R8))) then
                product = log(abs(au)) + log(abs(a(pivrow)))
                if (product.ge.cutlo) a(jj) = a(jj) + au*a(pivrow)
              end if
#else
              a(jj) = a(jj) + au*a(pivrow)
#endif
              if (lbig) big = max(abs(a(jj)),big)
              icn(pivrow) = -icn(pivrow)
              if (abs(a(jj)).lt.tol) idrop = idrop + 1
  590       continue
c
c  jump if no nonzeros in non-pivot row have been removed
c       because these are beneath the drop-tolerance  tol.
c
            if (idrop.eq.0) go to 650
c
c  run through non-pivot row compressing row so that only
c      nonzeros greater than   tol   are stored.  all nonzeros
c      less than   tol   are also removed from the column structure.
c
            jnew = rowi
            do 630 jj=rowi,iend
              if (abs(a(jj)).lt.tol) go to 600
              a(jnew) = a(jj)
              icn(jnew) = icn(jj)
              jnew = jnew + 1
              go to 630
c
c  remove nonzero entry from column structure.
c
  600         j = icn(jj)
              i1 = ipc(j)
              i2 = i1 + lenc(j) - 1
              do 610 ii=i1,i2
                if (irn(ii).eq.i) go to 620
  610         continue
  620         irn(ii) = irn(i2)
              irn(i2) = 0
              lenc(j) = lenc(j) - 1
              if (nsrch.le.nn) go to 630
c remove column from column chain and place in update chain.
              if (nextc(j).lt.0) go to 630
c jump if column already in update chain.
              lc = lastc(j)
              nc = nextc(j)
              nextc(j) = -colupd
              colupd = j
              if (nc.ne.0) lastc(nc) = lc
              if (lc.eq.0) go to 622
              nextc(lc) = nc
              go to 630
  622         nz = lenc(j) + 1
              isw = ifirst(nz)
              if (isw.gt.0) lastr(isw) = -nc
              if (isw.lt.0) ifirst(nz) = -nc
  630       continue
            do 640 jj=jnew,iend
              icn(jj) = 0
  640       continue
c the value of idrop might be different from that calculated earlier
c     because, we may now have dropped some nonzeros which were not
c     modified by the pivot row.
            idrop = iend + 1 - jnew
            iend = jnew - 1
            lenr(i) = lenr(i) - idrop
            nzrow = nzrow - idrop
            nzcol = nzcol - idrop
            ndrop = ndrop + idrop
  650       ifill = lenpiv - iop
c jump is if there is no fill-in.
            if (ifill.eq.0) go to 750
c now for the fill-in.
            minicn = max0(minicn,morei+ibeg-1+nzrow+ifill+lenr(i))
c see if there is room for fill-in.
c get maximum space for row i in situ.
            do 660 jdiff=1,ifill
              jnpos = iend + jdiff
              if (jnpos.gt.licn) go to 670
              if (icn(jnpos).ne.0) go to 670
  660       continue
c there is room for all the fill-in after the end of the row so it
c     can be left in situ.
c next available space for fill-in.
            iend = iend + 1
            go to 750
c jmore spaces for fill-in are required in front of row.
  670       jmore = ifill - jdiff + 1
            i1 = iptr(i)
c we now look in front of the row to see if there is space for
c     the rest of the fill-in.
            do 680 jdiff=1,jmore
              jnpos = i1 - jdiff
              if (jnpos.lt.iactiv) go to 690
              if (icn(jnpos).ne.0) go to 700
  680       continue
  690       jnpos = i1 - jmore
            go to 710
c whole row must be moved to the beginning of available storage.
  700       jnpos = iactiv - lenr(i) - ifill
c jump if there is space immediately available for the shifted row.
  710       if (jnpos.ge.ibeg) go to 730
            call ma30ddcopy3(a, icn, iptr(istart), n,iactiv,itop,.true.)
            i1 = iptr(i)
            iend = i1 + lenr(i) - 1
            jnpos = iactiv - lenr(i) - ifill
            if (jnpos.ge.ibeg) go to 730
c no space available so try to create some by throwing away previous
c     lu decomposition.
            morei = morei + ibeg - idisp(1) - lenpiv - 1
            if (lp.ne.0) write (lp,99997)
            iflag = -5
            if (abort3) go to 1090
c keep record of current pivot row.
            ibeg = idisp(1)
            icn(ibeg) = jpiv
            a(ibeg) = a(ijpos)
            ijpos = ibeg
            do 720 jj=ijp1,pivend
              ibeg = ibeg + 1
              a(ibeg) = a(jj)
              icn(ibeg) = icn(jj)
  720       continue
            ijp1 = ijpos + 1
            pivend = ibeg
            ibeg = ibeg + 1
            if (jnpos.ge.ibeg) go to 730
c this still does not give enough room.
            iflag = -4
            go to 1090
  730       iactiv = min0(iactiv,jnpos)
c move non-pivot row i.
            iptr(i) = jnpos
            do 740 jj=i1,iend
              a(jnpos) = a(jj)
              icn(jnpos) = icn(jj)
              jnpos = jnpos + 1
              icn(jj) = 0
  740       continue
c first new available space.
            iend = jnpos
  750       nzrow = nzrow + ifill
c innermost fill-in loop which also resets icn.
            idrop = 0
            do 830 jj=ijp1,pivend
              j = icn(jj)
              if (j.lt.0) go to 820
#ifdef DBG
              if((au   .ge.tiny(1.0_R8).or.au   .le.-tiny(1.0_R8)).and.
     &           (a(jj).ge.tiny(1.0_R8).or.a(jj).le.-tiny(1.0_R8))) then
                product = log(abs(au)) + log(abs(a(jj)))
                if (product.ge.cutlo) then
                  anew = au*a(jj)
                  aanew = abs(anew)
                else
                  anew = 0.0_R8
                  aanew = 0.0_R8
                end if
              else
                anew = 0.0_R8
                aanew = 0.0_R8
              end if
#else
              anew = au*a(jj)
              aanew = abs(anew)
#endif
              if (aanew.ge.tol) go to 760
              idrop = idrop + 1
              ndrop = ndrop + 1
              nzrow = nzrow - 1
              minicn = minicn - 1
              ifill = ifill - 1
              go to 830
  760         if (lbig) big = max(aanew,big)
              a(iend) = anew
              icn(iend) = j
              iend = iend + 1
c
c put new entry in column file.
              minirn = max0(minirn,nzcol+lenc(j)+1)
              jend = ipc(j) + lenc(j)
              jroom = nzpc - iii + 1 + lenc(j)
              if (jend.gt.lirn) go to 770
              if (irn(jend).eq.0) go to 810
  770         if (jroom.lt.dispc) go to 780
c compress column file to obtain space for new copy of column.
              call ma30ddcopy3(a, irn, ipc(istart),n,dispc,lirn,.false.)
              if (jroom.lt.dispc) go to 780
              jroom = dispc - 1
              if (jroom.ge.lenc(j)+1) go to 780
c column file is not large enough.
              go to 1100
c copy column to beginning of file.
  780         jbeg = ipc(j)
              jend = ipc(j) + lenc(j) - 1
              jzero = dispc - 1
              dispc = dispc - jroom
              idispc = dispc
              do 790 ii=jbeg,jend
                irn(idispc) = irn(ii)
                irn(ii) = 0
                idispc = idispc + 1
  790         continue
              ipc(j) = dispc
              jend = idispc
              do 800 ii=jend,jzero
                irn(ii) = 0
  800         continue
  810         irn(jend) = i
              nzcol = nzcol + 1
              lenc(j) = lenc(j) + 1
c end of adjustment to column file.
              go to 830
c
  820         icn(jj) = -j
  830       continue
            if (idrop.eq.0) go to 834
            do 832 kdrop=1,idrop
            icn(iend) = 0
            iend = iend + 1
  832       continue
  834       lenr(i) = lenr(i) + ifill
c end of scan of pivot column.
  840     continue
c
c
c remove pivot column from column oriented storage and update row
c     ordering arrays.
          i1 = ipc(jpiv)
          i2 = ipc(jpiv) + lenc(jpiv) - 1
          nzcol = nzcol - lenc(jpiv)
          do 890 ii=i1,i2
            i = irn(ii)
            irn(ii) = 0
            nz = lenr(i) - lenrl(i)
            if (nz.ne.0) go to 850
            lastr(i) = 0
            go to 890
  850       ifir = ifirst(nz)
            ifirst(nz) = i
            if (ifir) 860, 880, 870
  860       lastr(i) = ifir
            nextr(i) = 0
            go to 890
  870       lastr(i) = lastr(ifir)
            nextr(i) = ifir
            lastr(ifir) = i
            go to 890
  880       lastr(i) = 0
            nextr(i) = 0
            nzmin = min0(nzmin,nz)
  890     continue
c restore iq and nullify u part of old pivot row.
c    record the column permutation in lastc(jpiv) and the row
c    permutation in lastr(ipiv).
  900     ipc(jpiv) = -ising
          lastr(ipiv) = pivot
          if (lenpiv.eq.0) go to 980
          nzrow = nzrow - lenpiv
          jval = ijp1
          jzer = iptr(ipiv)
          iptr(ipiv) = 0
          do 910 jcount=1,lenpiv
            j = icn(jval)
            iq(j) = icn(jzer)
            icn(jzer) = 0
            jval = jval + 1
            jzer = jzer + 1
  910     continue
c adjust column ordering arrays.
          if (nsrch.gt.nn) go to 920
          do 916 jj=ijp1,pivend
            j = icn(jj)
            nz = lenc(j)
            if (nz.ne.0) go to 914
            ipc(j) = 0
            go to 916
  914       nzmin = min0(nzmin,nz)
  916     continue
          go to 980
  920     jj = colupd
          do 970 jdummy=1,nn
            j = jj
            if (j.eq.nn+1) go to 980
            jj = -nextc(j)
            nz = lenc(j)
            if (nz.ne.0) go to 924
            ipc(j) = 0
            go to 970
  924       ifir = ifirst(nz)
            lastc(j) = 0
            if (ifir) 930, 940, 950
  930       ifirst(nz) = -j
            ifir = -ifir
            lastc(ifir) = j
            nextc(j) = ifir
            go to 970
  940       ifirst(nz) = -j
            nextc(j) = 0
            go to 960
  950       lc = -lastr(ifir)
            lastr(ifir) = -j
            nextc(j) = lc
            if (lc.ne.0) lastc(lc) = j
  960       nzmin = min0(nzmin,nz)
  970     continue
  980   continue
c ********************************************
c ****    end of main elimination loop    ****
c ********************************************
c
c reset iactiv to point to the beginning of the next block.
  990   if (ilast.ne.nn) iactiv = iptr(ilast+1)
 1000 continue
c
c ********************************************
c ****    end of decomposition of block   ****
c ********************************************
c
c record singularity (if any) in iq array.
      if (irank.eq.nn) go to 1020
      do 1010 i=1,nn
        if (ipc(i).lt.0) go to 1010
        ising = ipc(i)
        iq(ising) = -iq(ising)
        ipc(i) = -ising
 1010 continue
c
c run through lu decomposition changing column indices to that of new
c     order and permuting lenr and lenrl arrays according to pivot
c     permutations.
 1020 istart = idisp(1)
      iend = ibeg - 1
      if (iend.lt.istart) go to 1040
      do 1030 jj=istart,iend
        jold = icn(jj)
        icn(jj) = -ipc(jold)
 1030 continue
 1040 do 1050 ii=1,nn
        i = lastr(ii)
        nextr(i) = lenr(ii)
        iptr(i) = lenrl(ii)
 1050 continue
      do 1060 i=1,nn
        lenrl(i) = iptr(i)
        lenr(i) = nextr(i)
 1060 continue
c
c update permutation arrays ip and iq.
      do 1070 ii=1,nn
        i = lastr(ii)
        j = -ipc(ii)
        nextr(i) = iabs(ip(ii)+0)
        iptr(j) = iabs(iq(ii)+0)
 1070 continue
      do 1080 i=1,nn
        if (ip(i).lt.0) nextr(i) = -nextr(i)
        ip(i) = nextr(i)
        if (iq(i).lt.0) iptr(i) = -iptr(i)
        iq(i) = iptr(i)
 1080 continue
      ip(nn) = iabs(ip(nn)+0)
      idisp(2) = iend
      go to 1120
c
c   ***    error returns    ***
 1090 idisp(2) = iactiv
      if (lp.eq.0) go to 1120
      write (lp,99996)
      go to 1110
 1100 if (iflag.eq.-5) iflag = -6
      if (iflag.ne.-6) iflag = -3
      idisp(2) = iactiv
      if (lp.eq.0) go to 1120
      if (iflag.eq.-3) write (lp,99995)
      if (iflag.eq.-6) write (lp,99994)
 1110 pivot = pivot - istart + 1
      write (lp,99993) pivot, nblock, istart, ilast
      if (pivot.eq.0) write (lp,99992) minirn
c
c
 1120 return
99999 format (' error return from ma30a/adc because matrix is structu',
     * 'rally singular')
99998 format (' error return from ma30a/adc because matrix is numeric',
     * 'ally singular')
99997 format (' lu decomposition destroyed to create more space')
99996 format (' error return from ma30a/adc because licn not big',
     * ' enough')
99995 format (' error return from ma30a/adc because lirn not big',
     * ' enough')
99994 format (' error return from ma30a/adc lirn and licn too small')
99993 format (' at stage ', i5, ' in block ', i5, ' with first row ',
     * i5, ' and last row ', i5)
99992 format (' to continue set lirn to at least ', i8)
      end subroutine ma30adcopy3
      subroutine ma30bdcopy3(n, icn, a, licn, lenr, lenrl,idisp,ip,iq,w,
     * iw, iflag)
      use b2mod_types
      implicit none
c ma30b/bd performs the lu decomposition of the diagonal blocks of a
c     new matrix paq of the same sparsity pattern, using information
c     from a previous call to ma30a/ad. the entries of the input
c     matrix  must already be in their final positions in the lu
c     decomposition structure.  this routine executes about five times
c     faster than ma30a/ad.
c
c we now describe the argument list for ma30b/bd. consult ma30a/ad for
c     further information on these parameters.
c n  is an integer variable set to the order of the matrix.
c icn is an integer array of length licn. it should be unchanged
c     since the last call to ma30a/ad. it is not altered by ma30b/bd.
c a  is a real/real*8 array of length licn the user must set
c     entries idisp(1) to idisp(2) to contain the entries in the
c     diagonal blocks of the matrix paq whose column numbers are held
c     in icn, using corresponding positions. note that some zeros may
c     need to be held explicitly. on output entries idisp(1) to
c     idisp(2) of array a contain the lu decomposition of the diagonal
c     blocks of paq. entries a(1) to a(idisp(1)-1) are neither
c     required nor altered by ma30b/bd.
c licn  is an integer variable which must be set by the user to the
c     length of arrays a and icn. it is not altered by ma30b/bd.
c lenr,lenrl are integer arrays of length n. they should be
c     unchanged since the last call to ma30a/ad. they are not altered
c     by ma30b/bd.
c idisp  is an integer array of length 2. it should be unchanged since
c     the last call to ma30a/ad. it is not altered by ma30b/bd.
c ip,iq  are integer arrays of length n. they should be unchanged
c     since the last call to ma30a/ad. they are not altered by
c     ma30b/bd.
c w  is a real/real*8 array of length n which is used as
c     workspace by ma30b/bd.
c iw  is an integer array of length n which is used as workspace by
c     ma30b/bd.
c iflag  is an integer variable. on output from ma30b/bd, iflag has
c     the value zero if the factorisation was successful, has the
c     value i if pivot i was very small and has the value -i if an
c     unexpected singularity was detected at stage i of the
c     decomposition.
c
      integer n, licn, iflag
      real (kind=R8) ::
     *   a(licn), w(n), au, eps, rowmax, zero, one, rmin,
     * tol, big
      integer lp
      integer ndrop, nsrch
      logical abort1, abort2, abort3, stab, lbig
      integer iw(n), idisp(2), pivpos
      integer icn(licn), lenr(n), lenrl(n), ip(n), iq(n)
c see block data for comments on variables in common.
      common /ma30edc3/ lp, abort1, abort2, abort3
      common /ma30idc3/ tol, big, ndrop, nsrch, lbig
      common /ma30gdc3/ eps, rmin
c
      integer i, j, jj, jay, jayjay, ifin, jfin, istart,
     * ilend, ipivj, ising
#ifdef DBG
      real (kind=R8) :: product
      real (kind=R8), save :: cutlo
#endif
c
      data zero /0.0e0_R8/, one /1.0e0_R8/
#ifdef DBG
      cutlo = log(tiny(1.0_R8))
#endif
      stab = eps.le.one
      rmin = eps
      ising = 0
      iflag = 0
      do 10 i=1,n
        w(i) = zero
   10 continue
c set up pointers to the beginning of the rows.
      iw(1) = idisp(1)
      if (n.eq.1) go to 25
      do 20 i=2,n
        iw(i) = iw(i-1) + lenr(i-1)
   20 continue
c
c   ****   start  of main loop    ****
c at step i, row i of a is transformed to row i of l/u by adding
c     appropriate multiples of rows 1 to i-1.
c     .... using row-gauss elimination.
   25 do 160 i=1,n
c istart is beginning of row i of a and row i of l.
        istart = iw(i)
c ifin is end of row i of a and row i of u.
        ifin = istart + lenr(i) - 1
c ilend is end of row i of l.
        ilend = istart + lenrl(i) - 1
        if (istart.gt.ilend) go to 90
c load row i of a into vector w.
        do 30 jj=istart,ifin
          j = icn(jj)
          w(j) = a(jj)
   30   continue
c
c add multiples of appropriate rows of  i to i-1  to row i.
        do 70 jj=istart,ilend
          j = icn(jj)
c ipivj is position of pivot in row j.
          ipivj = iw(j) + lenrl(j)
c form multiplier au.
#ifdef DBG
          if((w(j).ge.tiny(1.0_R8).or.w(j).le.-tiny(1.0_R8)).and.
     &        (a(ipivj).ge. tiny(1.0_R8).or.
     &         a(ipivj).le.-tiny(1.0_R8))) then
            product = log(abs(w(j))) - log(abs(a(ipivj)))
            if (product.ge.cutlo) then
              au = -w(j)/a(ipivj)
            else
              au = 0.0_R8
            end if
          else
            au = 0.0_R8
          end if
#else
          au = -w(j)/a(ipivj)
#endif
          if (lbig) big = max(abs(au),big)
          w(j) = au
c au * row j (u part) is added to row i.
          ipivj = ipivj + 1
          jfin = iw(j) + lenr(j) - 1
          if (ipivj.gt.jfin) go to 70
c innermost loop.
          if (lbig) go to 50
          do 40 jayjay=ipivj,jfin
            jay = icn(jayjay)
#ifdef DBG
            if((au.ge.tiny(1.0_R8).or.au.le.-tiny(1.0_R8)).and.
     &         (a(jayjay).ge. tiny(1.0_R8).or.
     &          a(jayjay).le.-tiny(1.0_R8))) then
              product = log(abs(au)) + log(abs(a(jayjay)))
              if (product.ge.cutlo) w(jay) = w(jay) + au*a(jayjay)
            end if
#else
            w(jay) = w(jay) + au*a(jayjay)
#endif
   40     continue
          go to 70
   50     do 60 jayjay=ipivj,jfin
            jay = icn(jayjay)
            w(jay) = w(jay) + au*a(jayjay)
            big = max(abs(w(jay)),big)
   60     continue
   70   continue
c
c reload w back into a (now l/u)
        do 80 jj=istart,ifin
          j = icn(jj)
          a(jj) = w(j)
          w(j) = zero
   80   continue
c we now perform the stability checks.
   90   pivpos = ilend + 1
        if (iq(i).gt.0) go to 140
c matrix had singularity at this point in ma30a/ad.
c is it the first such pivot in current block ?
        if (ising.eq.0) ising = i
c does current matrix have a singularity in the same place ?
        if (pivpos.gt.ifin) go to 100
        if (a(pivpos).ne.zero) go to 170
c it does .. so set ising if it is not the end of the current block
c check to see that appropriate part of l/u is zero or null.
  100   if (istart.gt.ifin) go to 120
        do 110 jj=istart,ifin
          if (icn(jj).lt.ising) go to 110
          if (a(jj).ne.zero) go to 170
  110   continue
  120   if (pivpos.le.ifin) a(pivpos) = one
        if (ip(i).gt.0 .and. i.ne.n) go to 160
c end of current block ... reset zero pivots and ising.
        do 130 j=ising,i
          if ((lenr(j)-lenrl(j)).eq.0) go to 130
          jj = iw(j) + lenrl(j)
          a(jj) = zero
  130   continue
        ising = 0
        go to 160
c matrix had nonzero pivot in ma30a/ad at this stage.
  140   if (pivpos.gt.ifin) go to 170
        if (a(pivpos).eq.zero) go to 170
        if (.not.stab) go to 160
        rowmax = zero
        do 150 jj=pivpos,ifin
          rowmax = max(rowmax,abs(a(jj)))
  150   continue
        if (abs(a(pivpos))/rowmax.ge.rmin) go to 160
        iflag = i
        rmin = abs(a(pivpos))/rowmax
c   ****    end of main loop    ****
  160 continue
c
      go to 180
c   ***   error return   ***
  170 if (lp.ne.0) write (lp,99999) i
      iflag = -i
c
  180 return
99999 format (' error return from ma30b/bdc singularity detected in ',
     * 'row', i8)
      end subroutine ma30bdcopy3
      subroutine ma30cdcopy3(n, icn, a, licn,lenr,lenrl,lenoff,idisp,ip,
     * iq, x, w, mtype)
      use b2mod_types
      implicit none
c ma30c/cd uses the factors produced by ma30a/ad or ma30b/bd to solve
c     ax=b or a transpose x=b when the matrix p1*a*q1 (paq) is block
c     lower triangular (including the case of only one diagonal
c     block).
c
c we now describe the argument list for ma30c/cd.
c n  is an integer variable set to the order of the matrix. it is not
c     altered by the subroutine.
c icn is an integer array of length licn. entries idisp(1) to
c     idisp(2) should be unchanged since the last call to ma30a/ad. if
c     the matrix has more than one diagonal block, then column indices
c     corresponding to nonzeros in subdiagonal blocks of paq must
c     appear in positions 1 to idisp(1)-1. for the same row those
c     entries must be contiguous, with those in row i preceding those
c     in row i+1 (i=1,...,n-1) and no wasted space between rows.
c     entries may be in any order within each row. it is not altered
c     by ma30c/cd.
c a  is a real/real*8 array of length licn.  entries
c     idisp(1) to idisp(2) should be unchanged since the last call to
c     ma30a/ad or ma30b/bd.  if the matrix has more than one diagonal
c     block, then the values of the nonzeros in subdiagonal blocks
c     must be in positions 1 to idisp(1)-1 in the order given by icn.
c     it is not altered by ma30c/cd.
c licn  is an integer variable set to the size of arrays icn and a.
c     it is not altered by ma30c/cd.
c lenr,lenrl are integer arrays of length n which should be
c     unchanged since the last call to ma30a/ad. they are not altered
c     by ma30c/cd.
c lenoff  is an integer array of length n. if the matrix paq (or
c     p1*a*q1) has more than one diagonal block, then lenoff(i),
c     i=1,...,n should be set to the number of nonzeros in row i of
c     the matrix paq which are in subdiagonal blocks.  if there is
c     only one diagonal block then lenoff(1) may be set to -1, in
c     which case the other entries of lenoff are never accessed. it is
c     not altered by ma30c/cd.
c idisp  is an integer array of length 2 which should be unchanged
c     since the last call to ma30a/ad. it is not altered by ma30c/cd.
c ip,iq are integer arrays of length n which should be unchanged
c     since the last call to ma30a/ad. they are not altered by
c     ma30c/cd.
c x is a real/real*8 array of length n. it must be set by
c     the user to the values of the right hand side vector b for the
c     equations being solved.  on exit from ma30c/cd it will be equal
c     to the solution x required.
c w  is a real/real*8 array of length n which is used as
c     workspace by ma30c/cd.
c mtype is an integer variable which must be set by the user. if
c     mtype=1, then the solution to the system ax=b is returned; any
c     other value for mtype will return the solution to the system a
c     transpose x=b. it is not altered by ma30c/cd.
c
      integer n, licn, mtype
      real (kind=R8) ::
     *   a(licn), x(n), w(n), wii, wi, resid, zero
      logical neg, nobloc
      integer idisp(2)
      integer icn(licn), lenr(n), lenrl(n), lenoff(n), ip(n), iq(n)
c see block data for comments on variables in common.
      common /ma30hdc3/ resid
c
      integer i, ib, ii, iii, j, j1, j2, j3, jj, k, lj1, lj2,
     * iend, iback, ifirst, ilast, iblock, iblend,
     * jpiv, jpivp1, lt, ltend, numblk
#ifdef DBG
      real (kind=R8) :: product
      real (kind=R8), save :: cutlo
#endif
c
      data zero /0.0e0_R8/
c
#ifdef DBG
      cutlo = log(tiny(1.0_R8))
#endif
c the final value of resid is the maximum residual for an inconsistent
c     set of equations.
      resid = zero
c nobloc is .true. if subroutine block has been used previously and
c     is .false. otherwise.  the value .false. means that lenoff
c     will not be subsequently accessed.
      nobloc = lenoff(1).lt.0
      if (mtype.ne.1) go to 140
c
c we now solve   a * x = b.
c neg is used to indicate when the last row in a block has been
c     reached.  it is then set to true whereafter backsubstitution is
c     performed on the block.
      neg = .false.
c ip(n) is negated so that the last row of the last block can be
c     recognised.  it is reset to its positive value on exit.
      ip(n) = -ip(n)
c preorder vector ... w(i) = x(ip(i))
      do 10 ii=1,n
        i = ip(ii)
        i = iabs(i)
        w(ii) = x(i)
   10 continue
c lt holds the position of the first nonzero in the current row of the
c     off-diagonal blocks.
      lt = 1
c ifirst holds the index of the first row in the current block.
      ifirst = 1
c iblock holds the position of the first nonzero in the current row
c     of the lu decomposition of the diagonal blocks.
      iblock = idisp(1)
c if i is not the last row of a block, then a pass through this loop
c     adds the inner product of row i of the off-diagonal blocks and w
c     to w and performs forward elimination using row i of the lu
c     decomposition.   if i is the last row of a block then, after
c     performing these aforementioned operations, backsubstitution is
c     performed using the rows of the block.
      do 120 i=1,n
        wi = w(i)
        if (nobloc) go to 30
        if (lenoff(i).eq.0) go to 30
c operations using lower triangular blocks.
c ltend is the end of row i in the off-diagonal blocks.
        ltend = lt + lenoff(i) - 1
        do 20 jj=lt,ltend
          j = icn(jj)
          wi = wi - a(jj)*w(j)
   20   continue
c lt is set the beginning of the next off-diagonal row.
        lt = ltend + 1
c set neg to .true. if we are on the last row of the block.
   30   if (ip(i).lt.0) neg = .true.
        if (lenrl(i).eq.0) go to 50
c forward elimination phase.
c iend is the end of the l part of row i in the lu decomposition.
        iend = iblock + lenrl(i) - 1
        do 40 jj=iblock,iend
          j = icn(jj)
#ifdef DBG
          if((a(jj).ge.tiny(1.0_R8).or.a(jj).le.-tiny(1.0_R8)).and.
     &       (w(j) .ge.tiny(1.0_R8).or.w(j) .le.-tiny(1.0_R8))) then
            product = log(abs(a(jj))) + log(abs(w(j)))
            if (product.ge.cutlo) wi = wi + a(jj)*w(j)
          end if
#else
          wi = wi + a(jj)*w(j)
#endif
   40   continue
c iblock is adjusted to point to the start of the next row.
   50   iblock = iblock + lenr(i)
        w(i) = wi
        if (.not.neg) go to 120
c back substitution phase.
c j1 is position in a/icn after end of block beginning in row ifirst
c     and ending in row i.
        j1 = iblock
c are there any singularities in this block?  if not, continue with
c     the backsubstitution.
        ib = i
        if (iq(i).gt.0) go to 70
        do 60 iii=ifirst,i
          ib = i - iii + ifirst
          if (iq(ib).gt.0) go to 70
          j1 = j1 - lenr(ib)
          resid = max(resid,abs(w(ib)))
          w(ib) = zero
   60   continue
c entire block is singular.
        go to 110
c each pass through this loop performs the back-substitution
c     operations for a single row, starting at the end of the block and
c     working through it in reverse order.
   70   do 100 iii=ifirst,ib
          ii = ib - iii + ifirst
c j2 is end of row ii.
          j2 = j1 - 1
c j1 is beginning of row ii.
          j1 = j1 - lenr(ii)
c jpiv is the position of the pivot in row ii.
          jpiv = j1 + lenrl(ii)
          jpivp1 = jpiv + 1
c jump if row ii of u has no nonzeros.
          if (j2.lt.jpivp1) go to 90
          wii = w(ii)
          do 80 jj=jpivp1,j2
            j = icn(jj)
#ifdef DBG
            if((a(jj).ge.tiny(1.0_R8).or.a(jj).le.-tiny(1.0_R8)).and.
     &         (w(j) .ge.tiny(1.0_R8).or.w(j) .le.-tiny(1.0_R8))) then
              product = log(abs(a(jj))) + log(abs(w(j)))
              if (product.ge.cutlo) wii = wii - a(jj)*w(j)
            end if
#else
            wii = wii - a(jj)*w(j)
#endif
   80     continue
          w(ii) = wii
   90     w(ii) = w(ii)/a(jpiv)
  100   continue
  110   ifirst = i + 1
        neg = .false.
  120 continue
c
c reorder solution vector ... x(i) = w(iqinverse(i))
      do 130 ii=1,n
        i = iq(ii)
        i = iabs(i)
        x(i) = w(ii)
  130 continue
      ip(n) = -ip(n)
      go to 320
c
c
c we now solve   atranspose * x = b.
c preorder vector ... w(i)=x(iq(i))
  140 do 150 ii=1,n
        i = iq(ii)
        i = iabs(i)
        w(ii) = x(i)
  150 continue
c lj1 points to the beginning the current row in the off-diagonal
c     blocks.
      lj1 = idisp(1)
c iblock is initialised to point to the beginning of the block after
c     the last one ]
      iblock = idisp(2) + 1
c ilast is the last row in the current block.
      ilast = n
c iblend points to the position after the last nonzero in the
c     current block.
      iblend = iblock
c each pass through this loop operates with one diagonal block and
c     the off-diagonal part of the matrix corresponding to the rows
c     of this block.  the blocks are taken in reverse order and the
c     number of times the loop is entered is min(n,no. blocks+1).
      do 290 numblk=1,n
        if (ilast.eq.0) go to 300
        iblock = iblock - lenr(ilast)
c this loop finds the index of the first row in the current block..
c     it is first and iblock is set to the position of the beginning
c     of this first row.
        do 160 k=1,n
          ii = ilast - k
          if (ii.eq.0) go to 170
          if (ip(ii).lt.0) go to 170
          iblock = iblock - lenr(ii)
  160   continue
  170   ifirst = ii + 1
c j1 points to the position of the beginning of row i (lt part) or pivot
        j1 = iblock
c forward elimination.
c each pass through this loop performs the operations for one row of the
c     block.  if the corresponding entry of w is zero then the
c     operations can be avoided.
        do 210 i=ifirst,ilast
          if (w(i).eq.zero) go to 200
c jump if row i singular.
          if (iq(i).lt.0) go to 220
c j2 first points to the pivot in row i and then is made to point to the
c     first nonzero in the u transpose part of the row.
          j2 = j1 + lenrl(i)
          wi = w(i)/a(j2)
          if (lenr(i)-lenrl(i).eq.1) go to 190
          j2 = j2 + 1
c j3 points to the end of row i.
          j3 = j1 + lenr(i) - 1
          do 180 jj=j2,j3
            j = icn(jj)
            w(j) = w(j) - a(jj)*wi
  180     continue
  190     w(i) = wi
  200     j1 = j1 + lenr(i)
  210   continue
        go to 240
c deals with rest of block which is singular.
  220   do 230 ii=i,ilast
          resid = max(resid,abs(w(ii)))
          w(ii) = zero
  230   continue
c back substitution.
c this loop does the back substitution on the rows of the block in
c     the reverse order doing it simultaneously on the l transpose part
c     of the diagonal blocks and the off-diagonal blocks.
  240   j1 = iblend
        do 280 iback=ifirst,ilast
          i = ilast - iback + ifirst
c j1 points to the beginning of row i.
          j1 = j1 - lenr(i)
          if (lenrl(i).eq.0) go to 260
c j2 points to the end of the l transpose part of row i.
          j2 = j1 + lenrl(i) - 1
          do 250 jj=j1,j2
            j = icn(jj)
            w(j) = w(j) + a(jj)*w(i)
  250     continue
  260     if (nobloc) go to 280
c operations using lower triangular blocks.
          if (lenoff(i).eq.0) go to 280
c lj2 points to the end of row i of the off-diagonal blocks.
          lj2 = lj1 - 1
c lj1 points to the beginning of row i of the off-diagonal blocks.
          lj1 = lj1 - lenoff(i)
          do 270 jj=lj1,lj2
            j = icn(jj)
            w(j) = w(j) - a(jj)*w(i)
  270     continue
  280   continue
        iblend = j1
        ilast = ifirst - 1
  290 continue
c reorder solution vector ... x(i)=w(ipinverse(i))
  300 do 310 ii=1,n
        i = ip(ii)
        i = iabs(i)
        x(i) = w(ii)
  310 continue
c
  320 return
      end subroutine ma30cdcopy3
      subroutine ma30ddcopy3(a, icn, iptr, n, iactiv, itop, reals)
      use b2mod_types
      implicit none
c this subroutine performs garbage collection operations on the
c     arrays a, icn and irn.
c iactiv is the first position in arrays a/icn from which the compress
c     starts.  on exit, iactiv equals the position of the first entry
c     in the compressed part of a/icn
c
      integer n, itop, iactiv
      real (kind=R8) ::
     *   a(itop)
      logical reals
      integer iptr(n)
      integer icn(itop)
c see block data for comments on variables in common.
      integer irncp, icncp, irank, minirn, minicn
      common /ma30fdc3/ irncp, icncp, irank, minirn, minicn
c
      integer j, jpos, k, kn, kl
c
      if (reals) icncp = icncp + 1
      if (.not.reals) irncp = irncp + 1
c set the first nonzero entry in each row to the negative of the
c     row/col number and hold this row/col index in the row/col
c     pointer.  this is so that the beginning of each row/col can
c     be recognised in the subsequent scan.
      do 10 j=1,n
        k = iptr(j)
        if (k.lt.iactiv) go to 10
        iptr(j) = icn(k)
        icn(k) = -j
   10 continue
      kn = itop + 1
      kl = itop - iactiv + 1
c go through arrays in reverse order compressing to the back so
c     that there are no zeros held in positions iactiv to itop in icn.
c     reset first entry of each row/col and pointer array iptr.
      do 30 k=1,kl
        jpos = itop - k + 1
        if (icn(jpos).eq.0) go to 30
        kn = kn - 1
        if (reals) a(kn) = a(jpos)
        if (icn(jpos).ge.0) go to 20
c first nonzero of row/col has been located
        j = -icn(jpos)
        icn(jpos) = iptr(j)
        iptr(j) = kn
   20   icn(kn) = icn(jpos)
   30 continue
      iactiv = kn
      return
      end subroutine ma30ddcopy3
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc13dcopy3
      subroutine mc13dcopy3(n,icn,licn,ip,lenr,ior,ib,num,iw)
      use b2mod_types
      implicit none
      integer n, licn, num
      integer ip(n)
      integer icn(licn),lenr(n),ior(n),ib(n),iw(n,3)
      call
     & mc13ecopy3(n,icn,licn,ip,lenr,ior,ib,num,iw(1,1),iw(1,2),iw(1,3))
      return
      end subroutine mc13dcopy3
      subroutine mc13ecopy3
     &  (n,icn,licn,ip,lenr,arp,ib,num,lowl,numb,prev)
      use b2mod_types
      implicit none
      integer n, licn, num
      integer stp,dummy
      integer ip(n)
c
c arp(i) is one less than the number of unsearched edges leaving
c     node i.  at the end of the algorithm it is set to a
c     permutation which puts the matrix in block lower
c     triangular form.
c ib(i) is the position in the ordering of the start of the ith
c     block.  ib(n+1-i) holds the node number of the ith node
c     on the stack.
c lowl(i) is the smallest stack position of any node to which a path
c     from node i has been found.  it is set to n+1 when node i
c     is removed from the stack.
c numb(i) is the position of node i in the stack if it is on
c     it, is the permuted order of node i for those nodes
c     whose final position has been found and is otherwise zero.
c prev(i) is the node at the end of the path when node i was
c     placed on the stack.
      integer icn(licn),lenr(n),arp(n),ib(n),lowl(n),numb(n),
     1prev(n)
c
      integer i, i1, i2, ii, j, k, iv, iw,
     * icnt, lcnt, isn, ist, ist1, nnm1
c
c   icnt is the number of nodes whose positions in final ordering have
c     been found.
      icnt=0
c num is the number of blocks that have been found.
      num=0
      nnm1=n+n-1
c
c initialisation of arrays.
      do 20 j=1,n
      numb(j)=0
      arp(j)=lenr(j)-1
   20 continue
c
c
      do 120 isn=1,n
c look for a starting node
      if (numb(isn).ne.0) go to 120
      iv=isn
c ist is the number of nodes on the stack ... it is the stack pointer.
      ist=1
c put node iv at beginning of stack.
      lowl(iv)=1
      numb(iv)=1
      ib(n)=iv
c
c the body of this loop puts a new node on the stack or backtracks.
      do 110 dummy=1,nnm1
      i1=arp(iv)
c have all edges leaving node iv been searched.
      if (i1.lt.0) go to 60
      i2=ip(iv)+lenr(iv)-1
      i1=i2-i1
c
c look at edges leaving node iv until one enters a new node or
c     all edges are exhausted.
      do 50 ii=i1,i2
      iw=icn(ii)
c has node iw been on stack already.
      if (numb(iw).eq.0) go to 100
c update value of lowl(iv) if necessary.
      lowl(iv)=min0(lowl(iv),lowl(iw))
   50 continue
c
c there are no more edges leaving node iv.
      arp(iv)=-1
c is node iv the root of a block.
   60 if (lowl(iv).lt.numb(iv)) go to 90
c
c order nodes in a block.
      num=num+1
      ist1=n+1-ist
      lcnt=icnt+1
c peel block off the top of the stack starting at the top and
c     working down to the root of the block.
      do 70 stp=ist1,n
      iw=ib(stp)
      lowl(iw)=n+1
      icnt=icnt+1
      numb(iw)=icnt
      if (iw.eq.iv) go to 80
   70 continue
   80 ist=n-stp
      ib(num)=lcnt
c are there any nodes left on the stack.
      if (ist.ne.0) go to 90
c have all the nodes been ordered.
      if (icnt.lt.n) go to 120
      go to 130
c
c backtrack to previous node on path.
   90 iw=iv
      iv=prev(iv)
c update value of lowl(iv) if necessary.
      lowl(iv)=min0(lowl(iv),lowl(iw))
      go to 110
c
c put new node on the stack.
 100  arp(iv)=i2-ii-1
      prev(iw)=iv
      iv=iw
      ist=ist+1
      lowl(iv)=ist
      numb(iv)=ist
      k=n+1-ist
      ib(k)=iv
  110 continue
c
  120 continue
c
c
c put permutation in the required form.
  130 do 140 i=1,n
      ii=numb(i)
      arp(ii)=i
  140 continue
      return
      end subroutine mc13ecopy3
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc20adcopy3 mc20bdcopy3
      subroutine mc20adcopy3(nc,maxa,a,inum,jptr,jnum,jdisp)
      use b2mod_types
      implicit none
c
      integer nc, maxa, jdisp
      integer inum(maxa),jnum(maxa),jptr(nc)
      real (kind=R8) ::
     *   a(maxa),ace,acep
c
      integer i, j, ja, jb, k, kr, ice, icep, jce, jcep, loc, null
c
c     ******************************************************************
c
      null=-jdisp
c**      clear jptr
      do 60 j=1,nc
      jptr(j)=0
   60 continue
c**      count the number of elements in each column.
      do 120 k=1,maxa
      j=jnum(k)+jdisp
      jptr(j)=jptr(j)+1
  120 continue
c**      set the jptr array
      k=1
      do 150 j=1,nc
      kr=k+jptr(j)
      jptr(j)=k
      k=kr
  150 continue
c
c**      reorder the elements into column order.  the algorithm is an
c        in-place sort and is of order maxa.
      do 230 i=1,maxa
c        establish the current entry.
      jce=jnum(i)+jdisp
      if(jce.eq.0) go to 230
      ace=a(i)
      ice=inum(i)
c        clear the location vacated.
      jnum(i)=null
c        chain from current entry to store items.
      do 200 j=1,maxa
c        current entry not in correct position.  determine correct
c        position to store entry.
      loc=jptr(jce)
      jptr(jce)=jptr(jce)+1
c        save contents of that location.
      acep=a(loc)
      icep=inum(loc)
      jcep=jnum(loc)
c        store current entry.
      a(loc)=ace
      inum(loc)=ice
      jnum(loc)=null
c        check if next current entry needs to be processed.
      if(jcep.eq.null) go to 230
c        it does.  copy into current entry.
      ace=acep
      ice=icep
      jce=jcep+jdisp
  200 continue
c
  230 continue
c
c**      reset jptr vector.
      ja=1
      do 250 j=1,nc
      jb=jptr(j)
      jptr(j)=ja
      ja=jb
  250 continue
      return
      end subroutine mc20adcopy3
      subroutine mc20bdcopy3(nc,maxa,a,inum,jptr)
      use b2mod_types
      implicit none
      integer nc, maxa
      real (kind=R8) ::
     *   a(maxa),ace
      integer inum(maxa), jptr(nc)
c
      integer j, k, ik, jj, ice, kor, klo, kmax, kdummy
c
c     ******************************************************************
c
      kmax=maxa
      do 30 jj=1,nc
      j=nc+1-jj
      klo=jptr(j)+1
      if(klo.gt.kmax)go to 30
      kor=kmax
      do 25 kdummy=klo,kmax
c items kor, kor+1, .... ,kmax are in order
      ace=a(kor-1)
      ice=inum(kor-1)
      do 10 k=kor,kmax
      ik=inum(k)
      if(iabs(ice).le.iabs(ik))go to 20
      inum(k-1)=ik
      a(k-1)=a(k)
10    continue
      k=kmax+1
20    inum(k-1)=ice
      a(k-1)=ace
      kor=kor-1
25    continue
c        next column
      kmax=klo-2
30    continue
      return
      end subroutine mc20bdcopy3
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc21acopy3
      subroutine mc21acopy3(n,icn,licn,ip,lenr,iperm,numnz,iw)
      use b2mod_types
      implicit none
      integer n, licn, numnz
      integer ip(n)
      integer icn(licn),lenr(n),iperm(n),iw(n,4)
      call mc21bcopy3(n,icn,licn,ip,lenr,iperm,numnz,iw(1,1),iw(1,2),
     1iw(1,3),iw(1,4))
      return
      end subroutine mc21acopy3
      subroutine mc21bcopy3
     &  (n,icn,licn,ip,lenr,iperm,numnz,pr,arp,cv,out)
      use b2mod_types
      implicit none
      integer n, licn, numnz
      integer ip(n)
c   pr(i) is the previous row to i in the depth first search.
c it is used as a work array in the sorting algorithm.
c   elements (iperm(i),i) i=1, ... n  are nonzero at the end of the
c algorithm unless n assignments have not been made.  in which case
c (iperm(i),i) will be zero for n-numnz entries.
c   cv(i) is the most recent row extension at which column i
c was visited.
c   arp(i) is one less than the number of nonzeros in row i
c which have not been scanned when looking for a cheap assignment.
c   out(i) is one less than the number of nonzeros in row i
c which have not been scanned during one pass through the main loop.
      integer icn(licn),lenr(n),iperm(n),pr(n),cv(n),
     1arp(n),out(n)
c
      integer i, j, k, ii, in1, in2, j1, jord, kk, ioutk
c
c   initialisation of arrays.
      do 10 i=1,n
      arp(i)=lenr(i)-1
      cv(i)=0
      iperm(i)=0
   10 continue
      numnz=0
c
c
c   main loop.
c   each pass round this loop either results in a new assignment
c or gives a row with no assignment.
      do 130 jord=1,n
      j=jord
      pr(j)=-1
      do 100 k=1,jord
c look for a cheap assignment
      in1=arp(j)
      if (in1.lt.0) go to 60
      in2=ip(j)+lenr(j)-1
      in1=in2-in1
      do 50 ii=in1,in2
      i=icn(ii)
      if (iperm(i).eq.0) go to 110
   50 continue
c   no cheap assignment in row.
      arp(j)=-1
c   begin looking for assignment chain starting with row j.
   60 out(j)=lenr(j)-1
c inner loop.  extends chain by one or backtracks.
      do 90 kk=1,jord
      in1=out(j)
      if (in1.lt.0) go to 80
      in2=ip(j)+lenr(j)-1
      in1=in2-in1
c forward scan.
      do 70 ii=in1,in2
      i=icn(ii)
      if (cv(i).eq.jord) go to 70
c   column i has not yet been accessed during this pass.
      j1=j
      j=iperm(i)
      cv(i)=jord
      pr(j)=j1
      out(j1)=in2-ii-1
      go to 100
   70 continue
c
c   backtracking step.
   80 j=pr(j)
      if (j.eq.-1) go to 130
   90 continue
c
  100 continue
c
c   new assignment is made.
  110 iperm(i)=j
      arp(j)=in2-ii-1
      numnz=numnz+1
      do 120 k=1,jord
      j=pr(j)
      if (j.eq.-1) go to 130
      ii=ip(j)+lenr(j)-out(j)-2
      i=icn(ii)
      iperm(i)=j
  120 continue
c
  130 continue
c
c   if matrix is structurally singular, we now complete the
c permutation iperm.
      if (numnz.eq.n) return
      do 140 i=1,n
      arp(i)=0
  140 continue
      k=0
      do 160 i=1,n
      if (iperm(i).ne.0) go to 150
      k=k+1
      out(k)=i
      go to 160
  150 j=iperm(i)
      arp(j)=i
  160 continue
      k=0
      do 170 i=1,n
      if (arp(i).ne.0) go to 170
      k=k+1
      ioutk=out(k)
      iperm(ioutk)=i
  170 continue
      return
      end subroutine mc21bcopy3
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc22adcopy3
      subroutine mc22adcopy3(n,icn,a,nz,lenrow,ip,iq,iw,iw1)
      use b2mod_types
      implicit none
      integer n, nz
      real (kind=R8) ::
     *   a(nz),aval
      integer iw(n,2)
      integer icn(nz),lenrow(n),ip(n),iq(n),iw1(nz)
      integer i, j, j2, jj, iold, ipos, jnum, jval,
     * ichain, length, newpos
      if (nz.le.0) go to 1000
      if (n.le.0) go to 1000
c set start of row i in iw(i,1) and lenrow(i) in iw(i,2)
      iw(1,1)=1
      iw(1,2)=lenrow(1)
      do 10 i=2,n
      iw(i,1)=iw(i-1,1)+lenrow(i-1)
      iw(i,2)=lenrow(i)
   10 continue
c permute lenrow according to ip.  set off-sets for new position
c     of row iold in iw(iold,1) and put old row indices in iw1 in
c     positions corresponding to the new position of this row in a/icn.
      jj=1
      do 20 i=1,n
      iold=ip(i)
      iold=iabs(iold)
      length=iw(iold,2)
      lenrow(i)=length
      if (length.eq.0) go to 20
      iw(iold,1)=iw(iold,1)-jj
      j2=jj+length-1
      do 15 j=jj,j2
      iw1(j)=iold
 15   continue
      jj=j2+1
 20   continue
c set inverse permutation to iq in iw(.,2).
      do 30 i=1,n
      iold=iq(i)
      iold=iabs(iold)
      iw(iold,2)=i
 30   continue
c permute a and icn in place, changing to new column numbers.
c
c ***   main loop   ***
c each pass through this loop places a closed chain of column indices
c     in their new (and final) positions ... this is recorded by
c     setting the iw1 entry to zero so that any which are subsequently
c     encountered during this major scan can be bypassed.
      do 200 i=1,nz
      iold=iw1(i)
      if (iold.eq.0) go to 200
      ipos=i
      jval=icn(i)
c if row iold is in same positions after permutation go to 150.
      if (iw(iold,1).eq.0) go to 150
      aval=a(i)
c **  chain loop  **
c each pass through this loop places one (permuted) column index
c     in its final position  .. viz. ipos.
      do 100 ichain=1,nz
c newpos is the original position in a/icn of the element to be placed
c in position ipos.  it is also the position of the next element in
c     the chain.
      newpos=ipos+iw(iold,1)
c is chain complete ?
      if (newpos.eq.i) go to 130
      a(ipos)=a(newpos)
      jnum=icn(newpos)
      icn(ipos)=iw(jnum,2)
      ipos=newpos
      iold=iw1(ipos)
      iw1(ipos)=0
c **  end of chain loop  **
 100  continue
 130  a(ipos)=aval
 150  icn(ipos)=iw(jval,2)
c ***   end of main loop   ***
 200  continue
c
 1000 return
      end subroutine mc22adcopy3
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc23adcopy3
c###### calls   mc13    mc21
      subroutine mc23adcopy3
     *  (n,icn,a,licn,lenr,idisp,ip,iq,lenoff,iw,iw1)
      use b2mod_types
      implicit none
      integer n, licn
      real (kind=R8) ::
     *   a(licn)
      integer idisp(2),iw1(n,2)
      logical abort
      integer icn(licn),lenr(n),ip(n),iq(n),lenoff(n),iw(n,5)
      integer lp,numnz,num,large
      common /mc23bdc3/ lp,numnz,num,large,abort
c input ... n,icn .. a,icn,lenr ....
      integer i, j, k, i1, i2, ii, jj, ibeg, iend, inew, iold,
     * jold, jnew, ilend, iblock, irowb, irowe, jnpos, leni, nz
c
c set up pointers iw(.,1) to the beginning of the rows and set lenoff
c     equal to lenr.
      iw1(1,1)=1
      lenoff(1)=lenr(1)
      if (n.eq.1) go to 20
      do 10 i=2,n
      lenoff(i)=lenr(i)
      iw1(i,1)=iw1(i-1,1)+lenr(i-1)
   10 continue
c idisp(1) points to the first position in a/icn after the
c     off-diagonal blocks and untreated rows.
   20 idisp(1)=iw1(n,1)+lenr(n)
c
c find row permutation ip to make diagonal zero-free.
      call mc21acopy3(n,icn,licn,iw1,lenr,ip,numnz,iw)
c
c possible error return for structurally singular matrices.
      if (numnz.ne.n.and.abort) go to 170
c
c iw1(.,2) and lenr are permutations of iw1(.,1) and lenr/lenoff
c     suitable for entry
c     to mc13dcopy3 since matrix with these row pointer and length arrays
c     has maximum number of nonzeros on the diagonal.
      do 30 ii=1,n
      i=ip(ii)
      iw1(ii,2)=iw1(i,1)
      lenr(ii)=lenoff(i)
   30 continue
c
c find symmetric permutation iq to block lower triangular form.
      call mc13dcopy3(n,icn,licn,iw1(1,2),lenr,iq,iw(1,4),num,iw)
c
      if (num.ne.1) go to 60
c
c action taken if matrix is irreducible.
c whole matrix is just moved to the end of the storage.
      do 40 i=1,n
      lenr(i)=lenoff(i)
      ip(i)=i
      iq(i)=i
   40 continue
      lenoff(1)=-1
c idisp(1) is the first position after the last element in the
c     off-diagonal blocks and untreated rows.
      nz=idisp(1)-1
      idisp(1)=1
c idisp(2) is the position in a/icn of the first element in the
c     diagonal blocks.
      idisp(2)=licn-nz+1
      large=n
      if (nz.eq.licn) go to 230
      do 50 k=1,nz
      j=nz-k+1
      jj=licn-k+1
      a(jj)=a(j)
      icn(jj)=icn(j)
   50 continue
c 230 = return
      go to 230
c
c data structure reordered.
c
c form composite row permutation ... ip(i) = ip(iq(i)).
   60 do 70 ii=1,n
      i=iq(ii)
      iw(ii,1)=ip(i)
   70 continue
      do 80 i=1,n
      ip(i)=iw(i,1)
   80 continue
c
c run through blocks in reverse order separating diagonal blocks
c     which are moved to the end of the storage.  elements in
c     off-diagonal blocks are left in place unless a compress is
c     necessary.
c
c ibeg indicates the lowest value of j for which icn(j) has been
c     set to zero when element in position j was moved to the
c     diagonal block part of storage.
      ibeg=licn+1
c iend is the position of the first element of those treated rows
c     which are in diagonal blocks.
      iend=licn+1
c large is the dimension of the largest block encountered so far.
      large=0
c
c num is the number of diagonal blocks.
      do 150 k=1,num
      iblock=num-k+1
c i1 is first row (in permuted form) of block iblock.
c i2 is last row (in permuted form) of block iblock.
      i1=iw(iblock,4)
      i2=n
      if (k.ne.1) i2=iw(iblock+1,4)-1
      large=max0(large,i2-i1+1)
c go through the rows of block iblock in the reverse order.
      do 140 ii=i1,i2
      inew=i2-ii+i1
c we now deal with row inew in permuted form (row iold in original
c     matrix).
      iold=ip(inew)
c if there is space to move up diagonal block portion of row go to 110
      if (iend-idisp(1).ge.lenoff(iold)) go to 110
c
c in-line compress.
c moves separated off-diagonal elements and untreated rows to
c     front of storage.
      jnpos=ibeg
      ilend=idisp(1)-1
      if (ilend.lt.ibeg) go to 190
      do 90 j=ibeg,ilend
      if (icn(j).eq.0) go to 90
      icn(jnpos)=icn(j)
      a(jnpos)=a(j)
      jnpos=jnpos+1
   90 continue
      idisp(1)=jnpos
      if (iend-jnpos.lt.lenoff(iold)) go to 190
      ibeg=licn+1
c reset pointers to the beginning of the rows.
      do 100 i=2,n
      iw1(i,1)=iw1(i-1,1)+lenoff(i-1)
  100 continue
c
c row iold is now split into diag. and off-diag. parts.
  110 irowb=iw1(iold,1)
      leni=0
      irowe=irowb+lenoff(iold)-1
c backward scan of whole of row iold (in original matrix).
      if (irowe.lt.irowb) go to 130
      do 120 jj=irowb,irowe
      j=irowe-jj+irowb
      jold=icn(j)
c iw(.,2) holds the inverse permutation to iq.
c     ..... it was set to this in mc13dcopy3.
      jnew=iw(jold,2)
c if (jnew.lt.i1) then ....
c element is in off-diagonal block and so is left in situ.
      if (jnew.lt.i1) go to 120
c element is in diagonal block and is moved to the end of the storage.
      iend=iend-1
      a(iend)=a(j)
      icn(iend)=jnew
      ibeg=min0(ibeg,j)
      icn(j)=0
      leni=leni+1
  120 continue
c
      lenoff(iold)=lenoff(iold)-leni
  130 lenr(inew)=leni
  140 continue
c
      ip(i2)=-ip(i2)
  150 continue
c resets ip(n) to positive value.
      ip(n)=-ip(n)
c idisp(2) is position of first element in diagonal blocks.
      idisp(2)=iend
c
c this compress is used to move all off-diagonal elements to the
c     front of the storage.
      if (ibeg.gt.licn) go to 230
      jnpos=ibeg
      ilend=idisp(1)-1
      do 160 j=ibeg,ilend
      if (icn(j).eq.0) go to 160
      icn(jnpos)=icn(j)
      a(jnpos)=a(j)
      jnpos=jnpos+1
  160 continue
c idisp(1) is first position after last element of off-diagonal blocks.
      idisp(1)=jnpos
      go to 230
c
c
c error return
  170 if (lp.ne.0) write(lp,180) numnz
  180 format(33x,' matrix is structurally singular, rank = ',i6)
      idisp(1)=-1
      go to 210
  190 if (lp.ne.0) write(lp,200) n
  200 format(33x,' licn not big enough increase by ',i6)
      idisp(1)=-2
  210 if (lp.ne.0) write(lp,220)
  220 format('+error return from mc23adc because')
c
  230 return
      end subroutine mc23adcopy3
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc24adcopy3
      subroutine mc24adcopy3(n,icn,a,licn,lenr,lenrl,w)
      use b2mod_types
      implicit none
      integer n, licn
      real (kind=R8) ::
     *   a(licn),w(n),amaxl,wrowl,amaxu,zero
      integer icn(licn),lenr(n),lenrl(n)
      integer i, j, j0, j1, j2, jj
      data zero/0.0e0_R8/
      amaxl=zero
      do 10 i=1,n
      w(i)=zero
 10   continue
      j0=1
      do 100 i=1,n
      if (lenr(i).eq.0) go to 100
      j2=j0+lenr(i)-1
      if (lenrl(i).eq.0) go to 50
c calculation of 1-norm of l.
      j1=j0+lenrl(i)-1
      wrowl=zero
      do 30 jj=j0,j1
      wrowl=wrowl+abs(a(jj))
 30   continue
c amaxl is the maximum norm of columns of l so far found.
      amaxl=max(amaxl,wrowl)
      j0=j1+1
c calculation of norms of columns of u (max-norms).
 50   j0=j0+1
      if (j0.gt.j2) go to 90
      do 80 jj=j0,j2
      j=icn(jj)
      w(j)=max(abs(a(jj)),w(j))
 80   continue
 90   j0=j2+1
 100  continue
c amaxu is set to maximum max-norm of columns of u.
      amaxu=zero
      do 200 i=1,n
      amaxu=max(amaxu,w(i))
 200  continue
c grofac is max u max-norm times max l 1-norm.
      w(1)=amaxl*amaxu
      return
      end subroutine mc24adcopy3
      block data bd1copy3
      use b2mod_types
      implicit none
c
c comments on all the common block variables are given here even
c     though some are not initialised by block data.
c lp,mp are used by the subroutine as the unit numbers for its warning
c     and diagnostic messages. default value for both is 6 (for line
c     printer output). the user can either reset them to a different
c     stream number or suppress the output by setting them to zero.
c     while lp directs the output of error diagnostics from the
c     principal subroutines and internally called subroutines, mp
c     controls only the output of a message which warns the user that he
c     has input two or more nonzeros a(i), . . ,a(k) with the same row
c     and column indices.  the action taken in this case is to proceed
c     using a numerical value of a(i)+...+a(k). in the absence of other
c     errors, iflag will equal -14 on exit.
c lblock is a logical variable which controls an option of first
c     preordering the matrix to block lower triangular form (using
c     harwell subroutine mc23a). the preordering is performed if lblock
c     is equal to its default value of .true. if lblock is set to
c     .false. , the option is not invoked and the space allocated to
c     ikeep can be reduced to 4*n+1.
c grow is a logical variable. if it is left at its default value of
c     .true. , then on return from ma28a/ad or ma28b/bd, w(1) will give
c     an estimate (an upper bound) of the increase in size of elements
c     encountered during the decomposition. if the matrix is well
c     scaled, then a high value for w(1), relative to the largest entry
c     in the input matrix, indicates that the lu decomposition may be
c     inaccurate and the user should be wary of his results and perhaps
c     increase u for subsequent runs.  we would like to emphasise that
c     this value only relates to the accuracy of our lu decomposition
c     and gives no indication as to the singularity of the matrix or the
c     accuracy of the solution.  this upper bound can be a significant
c     overestimate particularly if the matrix is badly scaled. if an
c     accurate value for the growth is required, lbig (q.v.) should be
c     set to .true.
c eps,rmin are real*8 variables. if, on entry to ma28b/bd, eps is less
c     than one, then rmin will give the smallest ratio of the pivot to
c     the largest element in the corresponding row of the upper
c     triangular factor thus monitoring the stability of successive
c     factorisations. if rmin becomes very large and w(1) from
c     ma28b/bd is also very large, it may be advisable to perform a
c     new decomposition using ma28a/ad.
c resid is a real*8 variable which on exit from ma28c/cd gives the value
c     of the maximum residual over all the equations unsatisfied because
c     of dependency (zero pivots).
c irncp,icncp are integer variables which monitor the adequacy of "elbow
c     room" in irn and a/icn respectively. if either is quite large (say
c     greater than n/10), it will probably pay to increase the size of
c     the corresponding array for subsequent runs. if either is very low
c     or zero then one can perhaps save storage by reducing the size of
c     the corresponding array.
c minirn,minicn are integer variables which, in the event of a
c     successful return (iflag ge 0 or iflag=-14) give the minimum size
c     of irn and a/icn respectively which would enable a successful run
c     on an identical matrix. on an exit with iflag equal to -5, minicn
c     gives the minimum value of icn for success on subsequent runs on
c     an identical matrix. in the event of failure with iflag= -6, -4,
c     -3, -2, or -1, then minicn and minirn give the minimum value of
c     licn and lirn respectively which would be required for a
c     successful decomposition up to the point at which the failure
c     occurred.
c irank is an integer variable which gives an upper bound on the rank of
c     the matrix.
c abort1 is a logical variable with default value .true.  if abort1 is
c     set to .false.  then ma28a/ad will decompose structurally singular
c     matrices (including rectangular ones).
c abort2 is a logical variable with default value .true.  if abort2 is
c     set to .false. then ma28a/ad will decompose numerically singular
c     matrices.
c idisp is an integer array of length 2. on output from ma28a/ad, the
c     indices of the diagonal blocks of the factors lie in positions
c     idisp(1) to idisp(2) of a/icn. this array must be preserved
c     between a call to ma28a/ad and subsequent calls to ma28b/bd,
c     ma28c/cd or ma28i/id.
c tol is a real*8 variable.  if it is set to a positive value, then any
c     nonzero whose modulus is less than tol will be dropped from the
c     factorisation.  the factorisation will then require less storage
c     but will be inaccurate.  after a run of ma28a/ad with tol positive
c     it is not possible to use ma28b/bd and the user is recommended to
c     use ma28i/id to obtain the solution.  the default value for tol is
c     0.0.
c themax is a real*8 variable.  on exit from ma28a/ad, it will hold the
c     largest entry of the original matrix.
c big is a real*8 variable. if lbig has been set to .true., big will hold
c     the largest entry encountered during the factorisation by ma28a/ad
c     or ma28b/bd.
c dxmax is a real*8 variable. on exit from ma28i/id, dxmax will be set to
c     the largest component of the solution.
c errmax is a real*8 variable.  on exit from ma28i/id, if maxit is
c     positive, errmax will be set to the largest component in the
c     estimate of the error.
c dres is a real*8 variable.  on exit from ma28i/id, if maxit is positive,
c     dres will be set to the largest component of the residual.
c cgce is a real*8 variable. it is used by ma28i/id to check the
c     convergence rate.  if the ratio of successive corrections is
c     not less than cgce then we terminate since the convergence
c     rate is adjudged too slow.
c ndrop is an integer variable. if tol has been set positive, on exit
c     from ma28a/ad, ndrop will hold the number of entries dropped from
c     the data structure.
c maxit is an integer variable. it is the maximum number of iterations
c     performed by ma28i/id. it has a default value of 16.
c noiter is an integer variable. it is set by ma28i/id to the number of
c     iterative refinement iterations actually used.
c nsrch is an integer variable. if nsrch is set to a value less than n,
c     then a different pivot option will be employed by ma28a/ad.  this
c     may result in different fill-in and execution time for ma28a/ad.
c     if nsrch is less than or equal to n, the workspace array iw can be
c     reduced in length.  the default value for nsrch is 32768.
c istart is an integer variable. if istart is set to a value other than
c     zero, then the user must supply an estimate of the solution to
c     ma28i/id.  the default value for istart is zero.
c lbig is a logical variable. if lbig is set to .true., the value of the
c     largest element encountered in the factorisation by ma28a/ad or
c     ma28b/bd is returned in big.  setting lbig to .true.  will
c     increase the time for ma28a/ad marginally and that for ma28b/bd
c     by about 20%.  the default value for lbig is .false.
c
      real (kind=R8) ::
     *   eps, rmin, resid
      real (kind=R8) :: tol, themax, big, dxmax, errmax, dres, cgce
      integer lp, mp
      integer irncp, icncp, minirn, minicn, irank
      integer ndrop, maxit, noiter, nsrch, istart
      logical lblock, grow
      logical abort1, abort2
      logical lbig
      common /ma28edc3/ lp, mp, lblock, grow
      common /ma28fdc3/ eps, rmin, resid, irncp, icncp, minirn, minicn,
     * irank, abort1, abort2
      common /ma28hdc3/ tol, themax, big, dxmax, errmax, dres, cgce,
     * ndrop, maxit, noiter, nsrch, istart, lbig
c     common /ma28gdc3/ idisp(2)
      data eps /1.0e-4_R8/
      data tol /0.0e0_R8/, cgce /0.5e0_R8/
      data lp /6/, mp /6/
      data maxit /16/, nsrch /32768/, istart /0/, lbig /.false./
      data abort1 /.true./, abort2 /.true./
      data lblock /.true./, grow /.true./
      end block data bd1copy3
      block data bd2copy3
      use b2mod_types
      implicit none
c although all common block variables do not have default values,
c     we comment on all the common block variables here.
c
c common block ma30e/ed holds control parameters ....
c     common /ma30edc3/ lp, abort1, abort2, abort3
c the integer lp is the unit number to which the error messages are
c     sent. lp has a default value of 6.  this default value can be
c     reset by the user, if desired.  a value of 0 suppresses all
c     messages.
c the logical variables abort1,abort2,abort3 are used to control the
c     conditions under which the subroutine will terminate.
c if abort1 is .true. then the subroutine will exit  immediately on
c     detecting structural singularity.
c if abort2 is .true. then the subroutine will exit immediately on
c     detecting numerical singularity.
c if abort3 is .true. then the subroutine will exit immediately when
c     the available space in a/icn is filled up by the previously
c     decomposed, active, and undecomposed parts of the matrix.
c the default values for abort1,abort2,abort3 are set to .true.,.true.
c     and .false. respectively.
c
c the variables in the common block ma30f/fd are used to provide the
c     user with information on the decomposition.
c     common /ma30fdc3/ irncp, icncp, irank, minirn, minicn
c irncp and icncp are integer variables used to monitor the adequacy
c     of the allocated space in arrays irn and a/icn respectively, by
c     taking account of the number of data management compresses
c     required on these arrays. if irncp or icncp is fairly large (say
c     greater than n/10), it may be advantageous to increase the size
c     of the corresponding array(s).  irncp and icncp are initialised
c     to zero on entry to ma30a/ad and are incremented each time the
c     compressing routine ma30d/dd is entered.
c icncp is the number of compresses on a/icn.
c irncp is the number of compresses on irn.
c irank is an integer variable which gives an estimate (actually an
c     upper bound) of the rank of the matrix. on an exit with iflag
c     equal to 0, this will be equal to n.
c minirn is an integer variable which, after a successful call to
c     ma30a/ad, indicates the minimum length to which irn can be
c     reduced while still permitting a successful decomposition of the
c     same matrix. if, however, the user were to decrease the length
c     of irn to that size, the number of compresses (irncp) may be
c     very high and quite costly. if lirn is not large enough to begin
c     the decomposition on a diagonal block, minirn will be equal to
c     the value required to continue the decomposition and iflag will
c     be set to -3 or -6. a value of lirn slightly greater than this
c     (say about n/2) will usually provide enough space to complete
c     the decomposition on that block. in the event of any other
c     failure minirn gives the minimum size of irn required for a
c     successful decomposition up to that point.
c minicn is an integer variable which after a successful call to
c     ma30a/ad, indicates the minimum size of licn required to enable
c     a successful decomposition. in the event of failure with iflag=
c     -5, minicn will, if abort3 is left set to .false., indicate the
c     minimum length that would be sufficient to prevent this error in
c     a subsequent run on an identical matrix. again the user may
c     prefer to use a value of icn slightly greater than minicn for
c     subsequent runs to avoid too many conpresses (icncp). in the
c     event of failure with iflag equal to any negative value except
c     -4, minicn will give the minimum length to which licn could be
c     reduced to enable a successful decomposition to the point at
c     which failure occurred.  notice that, on a successful entry
c     idisp(2) gives the amount of space in a/icn required for the
c     decomposition while minicn will usually be slightly greater
c     because of the need for "elbow room".  if the user is very
c     unsure how large to make licn, the variable minicn can be used
c     to provide that information. a preliminary run should be
c     performed with abort3 left set to .false. and licn about 3/2
c     times as big as the number of nonzeros in the original matrix.
c     unless the initial problem is very sparse (when the run will be
c     successful) or fills in extremely badly (giving an error return
c     with iflag equal to -4), an error return with iflag equal to -5
c     should result and minicn will give the amount of space required
c     for a successful decomposition.
c
c common block ma30g/gd is used by the ma30b/bd entry only.
c     common /ma30gdc3/ eps, rmin
c eps is a real/real*8 variable. it is used to test for
c     small pivots. its default value is 1.0e-4 (1.0d-4 in d version).
c     if the user sets eps to any value greater than 1.0, then no
c     check is made on the size of the pivots. although the absence of
c     such a check would fail to warn the user of bad instability, its
c     absence will enable ma30b/bd to run slightly faster. an  a
c     posteriori  check on the stability of the factorisation can be
c     obtained from mc24a/ad.
c rmin is a real/real*8 variable which gives the user some
c     information about the stability of the decomposition.  at each
c     stage of the lu decomposition the magnitude of the pivot apiv
c     is compared with the largest off-diagonal entry currently in its
c     row (row of u), rowmax say. if the ratio
c                       min (apiv/rowmax)
c     where the minimum is taken over all the rows, is less than eps
c     then rmin is set to this minimum value and iflag is returned
c     with the value +i where i is the row in which this minimum
c     occurs.  if the user sets eps greater than one, then this test
c     is not performed. in this case, and when there are no small
c     pivots rmin will be set equal to eps.
c
c common block ma30h/hd is used by ma30c/cd only.
c     common /ma30hdc3/ resid
c resid is a real/real*8 variable. in the case of singular
c     or rectangular matrices its final value will be equal to the
c     maximum residual for the unsatisfied equations; otherwise its
c     value will be set to zero.
c
c common  block ma30i/id controls the use of drop tolerances, the
c     modified pivot option and the the calculation of the largest
c     entry in the factorisation process. this common block was added
c     to the ma30 package in february, 1983.
c     common /ma30idc3/ tol, big, ndrop, nsrch, lbig
c tol is a real/real*8 variable.  if it is set to a positive
c     value, then ma30a/ad will drop from the factors any nonzero
c     whose modulus is less than tol.  the factorisation will then
c     require less storage but will be inaccurate.  after a run of
c     ma30a/ad where entries have been dropped, ma30b/bd  should not
c     be called.  the default value for tol is 0.0.
c big is a real/real*8 variable.  if lbig has been set to
c     .true., big will be set to the largest entry encountered during
c     the factorisation.
c ndrop is an integer variable. if tol has been set positive, on exit
c     from ma30a/ad, ndrop will hold the number of entries dropped
c     from the data structure.
c nsrch is an integer variable. if nsrch is set to a value less than
c     or equal to n, then a different pivot option will be employed by
c     ma30a/ad.  this may result in different fill-in and execution
c     time for ma30a/ad. if nsrch is less than or equal to n, the
c     workspace arrays lastc and nextc are not referenced by ma30a/ad.
c     the default value for nsrch is 32768.
c lbig is a logical variable. if lbig is set to .true., the value of
c     the largest entry encountered in the factorisation by ma30a/ad
c     is returned in big.  setting lbig to .true.  will marginally
c     increase the factorisation time for ma30a/ad and will increase
c     that for ma30b/bd by about 20%.  the default value for lbig is
c     .false.
c
      real (kind=R8) ::
     *   eps, rmin, tol, big
      integer lp
      integer ndrop, nsrch
      logical abort1, abort2, abort3, lbig
      common /ma30edc3/ lp, abort1, abort2, abort3
      common /ma30gdc3/ eps, rmin
      common /ma30idc3/ tol, big, ndrop, nsrch, lbig
      data eps /1.0e-4_R8/, tol /0.0e0_R8/, big /0.0e0_R8/
      data lp /6/, nsrch /32768/
      data lbig /.false./
      data abort1 /.true./, abort2 /.true./, abort3 /.false./
      end block data bd2copy3
      block data bd3copy3
      use b2mod_types
      implicit none
      integer lp,numnz,num,large
      logical abort
      common /mc23bdc3/ lp,numnz,num,large,abort
      data lp/6/,abort/.false./
      end block data bd3copy3                                           !srv 26.12.17 }
