      subroutine class(nsdecl,nzdecl,nsl,ncl,t,zsp,xi,la,lab,ds,
     +                 coeffc)
c-----------------------------------------------------------
c     this subroutine calculates the classical contribution
c     to the transport
c 
c     input  nsdecl  : maximum number of species
c            nzdecl  : maximum number of charge states
c            nsl   : actual number of species
c            ncl   : array(nsdecl+1) actual number of charge 
c                   states per species
c            t    : array(nsdecl+1) temperature of every species
c            zsp  : array(nsdecl+1,nzdecl) the charge in units of e
c                   of every charge state
c            xi   : array(nsdecl+1,nzdecl) weigting factor for the 
c                   reduced charge state method
c            la   : array(3,3,nsdecl+1) test particle friction 
c                   coefficients
c            lab  : array(3,3,nsdecl+1,nsdecl+1) friction coefficients
c            ds   : array(nsdecl+1,nzdecl,2) the thermodynamic
c                   forces. for definition see routine 
c                   neoart
c 
c     output coeffc : the tranport coefficients (see neoart
c                     for definition
c-----------------------------------------------------------

      use b2mod_types

      implicit none

      integer ::
     *     nsdecl,nzdecl,nsl,ncl(0:nsdecl),i,j,k,l
      real (kind=R8) ::
     *     xi(0:nsdecl,0:nzdecl-1),la(3,3,0:nsdecl),
     &     lab(3,3,0:nsdecl,0:nsdecl),ds(0:nsdecl,0:nzdecl-1,2),
     &     coeffc(0:nsdecl,0:nzdecl-1,4),v,zsp(0:nsdecl,0:nzdecl-1),
     &     t(0:nsdecl)

      do i = 0, nsl-1
        do j = 0, ncl(i)-1
          
          coeffc(i,j,1) = 0.
          coeffc(i,j,2) = 0.
 
          do k = 0, nsl-1
            do l = 0, ncl(k)-1

              v = t(k)/zsp(k,l)
              coeffc(i,j,1) = coeffc(i,j,1) + v*xi(i,j)*
     +         xi(k,l)*(lab(1,1,i,k)*ds(k,l,1)+lab(1,2,i,k)*
     +           ds(k,l,2)) 
              coeffc(i,j,2) = coeffc(i,j,2) + v*xi(i,j)*
     +         xi(k,l)*(lab(2,1,i,k)*ds(k,l,1)+lab(2,2,i,k)*
     +           ds(k,l,2)) 
              if ((k.eq.i).and.(l.eq.j)) then 
                coeffc(i,j,1) = coeffc(i,j,1) + v*xi(i,j)*(
     +            la(1,1,i)*ds(i,j,1)+la(1,2,i)*ds(i,j,2))
                coeffc(i,j,2) = coeffc(i,j,2) + v*xi(i,j)*(
     +            la(2,1,i)*ds(i,j,1)+la(2,2,i)*ds(i,j,2))
              endif
            enddo
          enddo
        enddo
      enddo

      return
      end
