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
           isVertInFace = .true.
        else 
            isVertInFace = .false.
        endif
   
         
      end function isVertInFace

***********************************************************************


      logical function isMember(a,b,nb)
        integer :: a, nb, i
        integer :: b(1:nb)

        !Check whether a integer 'a' is a member of 
        ! array b with length nb
        isMember = .false.

        do i = 1,nb
            if (b(i) .eq. a) then
                 isMember = .true.
            endif
        enddo
      end function isMember

***********************************************************************

      logical function isQuad(iCv,m)
        integer :: iCv
        type(mapping) :: m
       
        isQuad = .false.

        if (m%cvVxP(iCv,2).eq.4) then
           isQuad = .true.
        endif

      end function isQuad

***********************************************************************

      logical function isTriangle(iCv,m)
         integer :: iCv
         type(mapping) :: m

         isTriangle = .false.

         if (m%cvVxP(iCv,2).eq.3) then
            isTriangle = .true.
         endif
      end function isTriangle

***********************************************************************


      logical function isAligned(iFc,m)

        integer :: iFc
        type(mapping) :: m
        
        isAligned = m%fcAligned(iFc)

      end function isAligned

***********************************************************************

      logical function isPoloidal(iFc,m,g,cvFpsi)
      
        type(mapping) :: m
        type(geometry) :: g
        integer :: iFc, cells(1:2), faces3(0:2), faces4_1(0:3), 
     &             faces4_2(0:3), i1, i2, i, faces5(0:4),
     &             faces3_1(0:2), faces3_2(0:2)
        real(kind=R8) :: tol, cvFpsi(m%nCv)


        isPoloidal = .false. 

        i1 = 0
        i2 = 0

        if (.not.isAligned(iFc,m) ) then ! not aligned - skips aligned radial faces
            cells = m%fcCv(iFc,:)
            tol = 1e-1_R8*abs( g%vxFpsi(m%fcVx(iFc,1)) - 
     &             g%vxFpsi(m%fcVx(iFc,2)))
            !all combinations of cells

            !two quads
            if ((m%cvFcP(cells(1),2).eq.4) .and.
     &           (m%cvFcP(cells(2),2).eq.4)) then
                
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+4-1)
                faces4_2 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+4-1)

                do i = 0,3
                   if (isAligned(faces4_1(i),m)) then
                       i1 = i1 + 1
                   endif 
                enddo
                do i = 0,3
                   if (isAligned(faces4_2(i),m)) then
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
                   isPoloidal = .true. 
                endif               
        
            !tria + quad
            elseif ((m%cvFcP(cells(1),2).eq.3) .and.
     &           (m%cvFcP(cells(2),2).eq.4)) then

                faces3 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+3-1)
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+4-1)
               do i = 0,2
                 if (isAligned(faces3(i),m)) then
                    i1 = i1 + 1
                 endif
               enddo 
               do i = 0,3
                 if (isAligned(faces4_1(i),m)) then
                    i2 = i2 + 1
                 endif
               enddo
               if ((i1 .eq. 1).and. (i2.eq.2)) then
                  isPoloidal = .true.
               endif                                  
            !quad + tria
            elseif ((m%cvFcP(cells(1),2).eq.4) .and.
     &           (m%cvFcP(cells(2),2).eq.3)) then
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+4-1)
                faces3 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+3-1)
               do i = 0,3
                 if (isAligned(faces4_1(i),m)) then
                    i1 = i1 + 1
                 endif
               enddo 
               do i = 0,2
                 if (isAligned(faces3(i),m)) then
                    i2 = i2 + 1
                 endif
               enddo
               if ((i1 .eq. 2).and. (i2.eq.1)) then
                  ispoloidal = .true.
               endif  


            !quad + guard
            elseif ((m%cvFcP(cells(1),2).eq.4) .and.
     &           (m%cvFcP(cells(2),2).eq.1)) then
                faces4_1 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+4-1)
               do i = 0,3
                 if (isAligned(faces4_1(i),m)) then
                    i1 = i1 + 1
                 endif
               enddo
               if (i1.eq.2) then
                 ispoloidal = .true.
               endif            
        
            !guard + quad
            elseif ((m%cvFcP(cells(1),2).eq.1) .and.
     &           (m%cvFcP(cells(2),2).eq.4)) then
                 faces4_2 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+4-1)     
               do i = 0,3
                 if (isAligned(faces4_2(i),m)) then
                    i2 = i2 + 1
                 endif
               enddo
               if (i2.eq.2) then
                 ispoloidal = .true.
               endif 

            !tria + guard
            elseif ((m%cvFcP(cells(1),2).eq.3) .and.
     &           (m%cvFcP(cells(2),2).eq.1)) then 
                 faces3 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+3-1)     
               do i = 0,2
                 if (isAligned(faces3(i),m)) then
                    i1 = i1 + 1
                 endif
               enddo
               if (i1.eq.1) then
                 ispoloidal = .true.
               endif         
            !guard + tria
            elseif ((m%cvFcP(cells(1),2).eq.1) .and.
     &           (m%cvFcP(cells(2),2).eq.3)) then                 
                 faces3 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+3-1)     
               do i = 0,2
                 if (isAligned(faces3(i),m)) then
                    i2 = i2 + 1
                 endif
               enddo
               if (i2.eq.1) then
                 ispoloidal = .true.
               endif
            ! two trias   
            elseif  ((m%cvFcP(cells(1),2).eq.3) .and.
     &           (m%cvFcP(cells(2),2).eq.3)) then
                  faces3_1 = 
     &         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+3-1)
                  faces3_2 = 
     &         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+3-1)
               do i = 0,2
                 if (isAligned(faces3_1(i),m)) then
                    i1 = i1 + 1
                 endif   
                 if (isAligned(faces3_2(i),m)) then
                    i2 = i2 + 1
                 endif               
               enddo
               if ((i1.eq.1).and.(i2.eq.1)) then
                 ispoloidal = .true.
               endif
             
            !pent + quad   
            !elseif ((m%cvFcP(cells(1),2).eq.4) .and.
            !&           (m%cvFcP(cells(2),2).eq.5)) then 
            !    faces4_1 =
            !&         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+4-1) 
            !   do i = 0,3
            !     if (isAligned(faces4_1(i),m)) then
            !        i1 = i1 + 1
            !     endif
            !   enddo  
            !   faces5 = 
            !&         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+5-1) 
            !   do i = 0,4
            !     if (isAligned(faces5(i),m)) then
            !        i2 = i2 + 1
            !     endif
            !   enddo 
            !   if ((i1.eq.2) .and. (i2.eq.2)) then 
            !     ispoloidal = .true. 
            !   endif    
            !elseif ((m%cvFcP(cells(1),2).eq.5) .and.
            !&           (m%cvFcP(cells(2),2).eq.4)) then
            !    faces4_1 =
            !&         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+4-1) 
            !   do i = 0,3
            !     if (isAligned(faces4_1(i),m)) then
            !        i1 = i1 + 1
            !     endif
            !   enddo  
            !   faces5 = 
            !&         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+5-1) 
            !   do i = 0,4
            !     if (isAligned(faces5(i),m)) then
            !        i2 = i2 + 1
            !     endif
            !   enddo         
            !   if ((i1.eq.2) .and. (i2.eq.2)) then 
            !     ispoloidal = .true. 
            !   endif   
            !elseif ((m%cvFcP(cells(1),2).eq.3) .and.
            !&           (m%cvFcP(cells(2),2).eq.5)) then
            !     faces3 = 
            !&         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+3-1)
            !   do i = 0,2
            !     if (isAligned(faces3(i),m)) then
            !        i1 = i1 + 1
            !     endif
            !   enddo 
            !   faces5 = 
            !&         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+5-1) 
            !   do i = 0,4
            !     if (isAligned(faces5(i),m)) then
            !        i2 = i2 + 1
            !     endif
            !   enddo         
            !   if ((i1.eq.1) .and. (i2.eq.2)) then 
            !     ispoloidal = .true. 
            !   endif               
            !elseif ((m%cvFcP(cells(1),2).eq.5) .and.
            !&           (m%cvFcP(cells(2),2).eq.3)) then
            !     faces3 = 
            !&         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+3-1)
            !   do i = 0,2
            !     if (isAligned(faces3(i),m)) then
            !        i1 = i1 + 1
            !     endif
            !   enddo 
            !   faces5 = 
            !&         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+5-1) 
            !   do i = 0,4
            !     if (isAligned(faces5(i),m)) then
            !        i2 = i2 + 1
            !     endif
            !   enddo         
            !   if ((i1.eq.1) .and. (i2.eq.2)) then 
            !     ispoloidal = .true. 
            !   endif
            !elseif ((m%cvFcP(cells(1),2).eq.1) .and.
            !&           (m%cvFcP(cells(2),2).eq.5)) then 
            !   faces5 = 
            !&         m%cvFc(m%cvFcP(cells(2),1):m%cvFcP(cells(2),1)+5-1) 
            !   do i = 0,4
            !     if (isAligned(faces5(i),m)) then
            !        i1 = i1 + 1
            !     endif
            !   enddo         
            !   if (i1.eq.2)  then 
            !     ispoloidal = .true. 
            !   endif
            !elseif ((m%cvFcP(cells(1),2).eq.5) .and.
            !&           (m%cvFcP(cells(2),2).eq.1)) then 
            !   faces5 = 
            !&         m%cvFc(m%cvFcP(cells(1),1):m%cvFcP(cells(1),1)+5-1) 
            !   do i = 0,4
            !     if (isAligned(faces5(i),m)) then
            !        i1 = i1 + 1
            !     endif
            !   enddo         
            !   if (i1.eq.2)  then 
            !     ispoloidal = .true. 
            !   endif                                    
            !else
            !  call xerrab('isPoloidal: case not implemented')    
            endif          
            
        endif

        

      end function


      end module b2us_cellhelper
