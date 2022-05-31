      module b2us_cellhelper
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none


      public :: isVertInFace, isMember, isQuad,
     .          isTriangle, isAligned, isPoloidal



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


      logical function isAligned(iFc,m,g,tol)
        !use b2us_geo

        integer :: iFc, verts(1:2)
        real(kind=R8) :: verts_psi(1:2), tol
        type(mapping) :: m
        type(geometry) :: g
        intrinsic abs
        
        isAligned = 0

        verts = m%fcVx(iFc,:)
        verts_psi = g%vxFpsi(verts)

        !tol should be 1.0e-4_R8 as default

        if (abs(verts_psi(1) - verts_psi(2)) .lt. tol) then
           isAligned = 1
        endif
      end function isAligned

***********************************************************************

      logical function isPoloidal(iFc,m,g,cvFpsi,tol)
      
        type(mapping) :: m
        type(geometry) :: g
        integer :: iFc, cells(1:2), faces3(0:2), faces4_1(0:3), 
     &             faces4_2(0:3), i1, i2, i  
        real(kind=R8) :: tol, cvFpsi(m%nCv)


        isPoloidal = 0 

        i1 = 0
        i2 = 0

        if (.not.isAligned(iFc,m,g,tol) ) then ! not aligned - skips aligned radial faces
            cells = m%fcCv(iFc,:)
            !all combinations of cells

            !two quads
            if ((m%cvFcP(cells(1),2).eq.4) .and.
     &           (m%cvFcP(cells(2),2).eq.4)) then
                
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+4-1)
                faces4_2 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+4-1)

                do i = 0,3
                   if (isAligned(faces4_1(i),m,g,tol)) then
                       i1 = i1 + 1
                   endif 
                enddo
                do i = 0,3
                   if (isAligned(faces4_2(i),m,g,tol)) then
                       i2 = i2 + 1
                   endif 
                enddo 
                if  (( i1 .eq. 2 .and. i2 .eq. 2 ) .or. 
     &              ( i1 .eq. 1 .and. i2.eq.2 .and. 
     &              abs(cvFpsi(cells(1)) - cvFpsi(cells(2)) 
     &               .lt.tol )) .or.
     &              ( i1 .eq. 2 .and. i2.eq.1 .and. 
     &              abs(cvFpsi(cells(1)) - cvFpsi(cells(2))
     &               .lt. tol ))) then
                   isPoloidal = 1 
                endif               
        
            !tria + quad
            elseif ((m%cvFcP(cells(1),2).eq.3) .and.
     &           (m%cvFcP(cells(2),2).eq.4)) then

                faces3 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+3-1)
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+4-1)
               do i = 0,2
                 if (isAligned(faces3(i),m,g,tol)) then
                    i1 = i1 + 1
                 endif
               enddo 
               do i = 0,3
                 if (isAligned(faces4_1(i),m,g,tol)) then
                    i2 = i2 + 1
                 endif
               enddo
               if ((i1 .eq. 1).and. (i2.eq.2)) then
                  isPoloidal = 1
               endif                                  
            !quad + tria
            elseif ((m%cvFcP(cells(1),2).eq.4) .and.
     &           (m%cvFcP(cells(2),2).eq.3)) then
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+4-1)
                faces3 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+3-1)
               do i = 0,3
                 if (isAligned(faces4_1(i),m,g,tol)) then
                    i1 = i1 + 1
                 endif
               enddo 
               do i = 0,2
                 if (isAligned(faces3(i),m,g,tol)) then
                    i2 = i2 + 1
                 endif
               enddo
               if ((i1 .eq. 2).and. (i2.eq.1)) then
                  isPoloidal = 1
               endif  


            !quad + guard
            elseif ((m%cvFcP(cells(1),2).eq.4) .and.
     &           (m%cvFcP(cells(2),2).eq.1)) then
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+4-1)
               do i = 0,3
                 if (isAligned(faces4_1(i),m,g,tol)) then
                    i1 = i1 + 1
                 endif
               enddo
               if (i1.eq.2) then
                 isPoloidal = 1
               endif            
        
            !guard + quad
            elseif ((m%cvFcP(cells(1),2).eq.1) .and.
     &           (m%cvFcP(cells(2),2).eq.4)) then
                 faces4_2 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+4-1)     
               do i = 0,3
                 if (isAligned(faces4_2(i),m,g,tol)) then
                    i2 = i2 + 1
                 endif
               enddo
               if (i2.eq.2) then
                 isPoloidal = 1
               endif 

            !tria + guard
            elseif ((m%cvFcP(cells(1),2).eq.3) .and.
     &           (m%cvFcP(cells(2),2).eq.1)) then 
                 faces3 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+3-1)     
               do i = 0,2
                 if (isAligned(faces3(i),m,g,tol)) then
                    i1 = i1 + 1
                 endif
               enddo
               if (i1.eq.1) then
                 isPoloidal = 1
               endif         
            !guard + tria
            elseif ((m%cvFcP(cells(1),2).eq.1) .and.
     &           (m%cvFcP(cells(2),2).eq.3)) then                 
                 faces3 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+3-1)     
               do i = 0,2
                 if (isAligned(faces3(i),m,g,tol)) then
                    i2 = i2 + 1
                 endif
               enddo
               if (i2.eq.1) then
                 isPoloidal = 1
               endif 

            endif          
            


        endif   
        

      end function


      end module b2us_cellhelper