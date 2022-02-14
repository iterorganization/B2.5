      subroutine nice(xmin,xmax,nomaj,nomin,fmtc,nfmt,nice_axes)
      character*10 fmtc
      real xmin, xmax
      integer nomaj, nomin, nfmt
      logical nice_axes
      real dx, scale
      integer iscale, ir, il
      real ceil, floor
      external ceil, floor
      dx=(xmax-xmin)/10
      iscale=int(-floor(alog10(dx)))
      ir=iscale
      scale=10.0**(iscale)
      dx=ceil(dx*scale)
      if(dx.ge.3.and.dx.le.5.and.nice_axes) dx=5
      if(dx.ge.6.and.dx.le.9.and.nice_axes) dx=10
      nomin=int(dx)
      if(nomin.gt.5.and.nice_axes) nomin=1
      dx=dx/scale
      xmin=floor(xmin/dx)*dx
      xmax=ceil(xmax/dx)*dx
      il=int(ceil(alog10(max(abs(xmin),abs(xmax)))))
      if(10.0**il.le.max(abs(xmin),abs(xmax))) il=il+1
      if(abs(il).le.3) then
        nfmt=max(1,il)+2+max(0,ir)
      else
        nfmt=il+8+ir-1
      endif
      nomaj=nint((xmax-xmin)/dx)
      if((nomaj.gt.5.and.nfmt.gt.5).or.
     1 (nomaj.gt.10.and.nfmt.le.5)) then
        if(mod(nomaj,2).eq.1) then
          xmax=xmin+(xmax-xmin)*float(nomaj+1)/float(nomaj)
          il=int(ceil(alog10(max(abs(xmin),abs(xmax)))))
          if(10.0**il.le.max(abs(xmin),abs(xmax))) il=il+1
          if(abs(il).le.3) then
            nfmt=max(1,il)+2+max(0,ir)
          else
            nfmt=il+8+ir-1
          endif
          nomaj=nomaj+1
        endif
        nomaj=nomaj/2
        nomin=nomin*2
      endif
      if(abs(il).le.3) then
        write(fmtc,'(a,i2,a,i2,a)')
     1   '(f',
     2   nfmt,
     3   '.',
     4   max(0,ir),
     5   ')'
      else
        write(fmtc,'(a,i2,a,i2,a)')
     1   '(1pe',
     2   nfmt,
     3   '.',
     4   il+ir-1,
     5   ')'
      endif
      return
      end
