      module b2us_cellhelper
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none


      public :: isVertInFace, isMember, isQuad,
     .          isTriangle, isAligned



      contains


      logical function isVertInFace(iVx,iFc,m)
        integer :: iVx, iFc, i
        type(mapping) :: m

        if ((m%fcVx(iFc,1) .eq. iVx) .or. (m%fcVx(iFc,2).eq.iVx)) then
           isVertInFace = 1
        else 
            isVertInFace = 0
        endif
   
         
      end function isVertInFace

***********************************************************************


      logical function isMember(a,b,nb)
        integer :: a, nb, i
        integer :: b(1:nb)

        !Check whether a integer 'a' is a member of 
        ! array b with length nb
        isMember = 0

        do i = 1,nb
            if (b(i) .eq. a) then
                 isMember = 1
            endif
        enddo
      end function isMember

***********************************************************************

      logical function isQuad(iCv,m)
        integer :: iCv
        type(mapping) :: m
       
        isQuad = 0

        if (m%cvVxP(iCv,2).eq.4) then
           isQuad = 1
        endif

      end function isQuad

***********************************************************************

      logical function isTriangle(iCv,m)
         integer :: iCv
         type(mapping) :: m

         isTriangle = 0

         if (m%cvVxP(iCv,2).eq.3) then
            isTriangle = 1
         endif
      end function isTriangle

***********************************************************************


      logical function isAligned(iFc,m,g)
        use b2us_geo

        integer :: iFc, verts(1:2)
        real(kind=R8) :: verts_psi(1:2)
        type(mapping) :: m
        type(geometry) :: g
        intrinsic abs
        
        isAligned = 0

        verts = m%fcVx(iFc,:)
        verts_psi = g%vxFpsi(verts)

        if (abs(verts_psi(1) - verts_psi(2)) .lt. 1.0e-4_R8) then
           isAligned = 1
        endif
      end function isAligned

      end module b2us_cellhelper