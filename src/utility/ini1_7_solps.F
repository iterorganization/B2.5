
      subroutine iniilu_7(n1,n2)                                   !sv 11.01.00
      use b2mod_types
      use b2mod_ma28_for_7diag
      use b2mod_indirect
      implicit none
      integer n1,n2

      write(*,*) 
     . 'INIILU_7: Initialising arrays for 7-point-stencil matrix'  !sv 11.01.00
      call iniilu_7_w(n_nonzeros,n1,n2,irn,ish,jsh,ijsh,icn,nc0)   !sv 11.01.00
      write(*,*) 
     . 'INIILU_7: Initialised arrays for 7-point-stencil matrix'   !sv 11.01.00
      return
      end subroutine iniilu_7

      subroutine iniilu_7_w(n_nonzeros,n1,n2,irn,ish,jsh,ijsh,icn,nc0)
      use b2mod_types
      use b2mod_indirect
      implicit none
c here: counting of indices from 1 to nx+2 and from 1 to ny+2
c main code : indices from -1 to nx and from -1 to ny
c therefore: shift of neighbourindices necessary: i,j gets i-2,j-2

      integer n1,n2,n_nonzeros,nc0
      integer irn(n_nonzeros), ish(n_nonzeros),
     1 jsh(n_nonzeros), ijsh(n_nonzeros), icn(n_nonzeros)

*  -- WORKSPACE ARGUMENTS --
      INTEGER N,NZ,I,J,K,IJ,JH,iw,jw                               !sv 26.01.00
      INTEGER I1(N1,N2,-1:5),J1(N1,N2,-1:5)                        !sv 26.01.00
      INTEGER I1H(12),IS(12),JS(12), IJS(12)                       !sv 26.01.00
      INTEGER HOLD,HOLDI,HOLDJ, HOLDIJ

      external subini, subend, xertst

      call subini ('iniilu_7_w')
      n= n1*n2
      i1=0
      j1=0
c
c   ASO1                                                           !sv 26.01.00
c
      do i=1,n1                                                    !sv 26.01.00
         do j=1,n2                                                 !sv 26.01.00
          if (bottomiy(i-2,j-2).ne.-2) then
            iw = bottomix(i-2,j-2)                                 !sv 26.01.00
            jw = bottomiy(i-2,j-2)                                 !sv 26.01.00
            if (bottomiy(iw,jw).ne.-2) then
              i1(i,j,-1)=i+(j-1)*n1                                !sv 26.01.00
              j1(i,j,-1)=2+bottomix(iw,jw)+(2+bottomiy(iw,jw)-1)*n1 
            endif
          endif
        enddo
      enddo                                                        !sv 26.01.00
c
c   ASO
c
      do i=1,n1
        do j=1,n2
          if (bottomiy(i-2,j-2).ne.-2) then
            i1(i,j,0)=i+(j-1)*n1
            j1(i,j,0)=2+bottomix(i-2,j-2)+(2+bottomiy(i-2,j-2)-1)*n1
          endif
        enddo
      enddo
c
c   AWE
c
      do i=1,n1
        do j=1,n2
          if (leftix(i-2,j-2).ne.-2) then
            i1(i,j,1)=i+(j-1)*n1
            j1(i,j,1)=2+leftix(i-2,j-2)+(2+leftiy(i-2,j-2)-1)*n1
          endif
        enddo
      enddo
c
c   ACP
c
      do i=1,n1
         do j=1,n2
          i1(i,j,2)=i+(j-1)*n1
          j1(i,j,2)=i+(j-1)*n1
       enddo
      enddo
c
c   AEA
c
      do i=1,n1
        do j=1,n2
          if (rightix(i-2,j-2).ne.n1-1) then
            i1(i,j,3)=i+(j-1)*n1
            j1(i,j,3)=2+rightix(i-2,j-2)+(2+rightiy(i-2,j-2)-1)*n1
          endif
        enddo
      enddo
c
c   ANO
c
      do i=1,n1
        do j=1,n2
          if (topiy(i-2,j-2).ne.n2-1) then
            i1(i,j,4)=i+(j-1)*n1
            j1(i,j,4)=2+topix(i-2,j-2)+(2+topiy(i-2,j-2)-1)*n1
          endif
        enddo
      enddo
c
c   ANO1                                                           !sv 26.01.00
c
      do i=1,n1                                                    !sv 26.01.00
         do j=1,n2                                                 !sv 26.01.00
          if (topiy(i-2,j-2).ne.n2-1) then
            iw = topix(i-2,j-2)                                    !sv 26.01.00
            jw = topiy(i-2,j-2)                                    !sv 26.01.00
            if (topiy(iw,jw).ne.n2-1) then
              i1(i,j,5)=i+(j-1)*n1                                 !sv 26.01.00
              j1(i,j,5)=2+topix(iw,jw)+(2+topiy(iw,jw)-1)*n1       !sv 26.01.00
            endif
          endif
        enddo
      enddo

      nz=0
      nc0=1
      icn(nc0)=nz+1
      do jh=1,n
          k=0
          nc0=nc0+1
          do i=1,n1
            do j=1,N2
              do ij=-1,5                                           !sv 26.01.00
                if (j1(i,j,ij).eq.jh) then
                     k=k+1
                     i1h(k)=i1(i,j,ij)
                     is(k)=i
                     js(k)=j
                     ijs(k)=ij
                endif
            enddo
          enddo
        enddo

C SORT  
          do i=1,k-1
            do j=i+1,k
              if (i1h(i).gt.i1h(j)) then
               hold=i1h(i)
               holdI=is(i)
               holdJ=js(i)
               holdIJ=ijs(i)
               i1h(i)=i1h(j)
               is(i)=is(j)
               js(i)=js(j)
               ijs(i)=ijs(j)
               i1H(j)=hold
               is(j)=holdi
               js(j)=holdj
               ijs(j)=holdij
              endif
          enddo
        enddo
C   END SORT
          do i=1,k
            nz=nz+1
            irn(nz)=i1h(i)
            ish(nz)=is(i)
            jsh(nz)=js(i)
            ijsh(nz)=ijs(i)
            icn(nz)=jh
        enddo
      enddo
      nc0=nz
      call xertst (nz.le.N_NONZEROS, 'Increase N_NONZEROS in b2ux7p')   !sv 26.01.00

      call subend ()
      return
      end subroutine iniilu_7_w


