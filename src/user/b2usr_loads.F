C
CVK CORRECT CALCULATION OF THE HEAT LOADS IN B2
C
      subroutine b2usr_loads(nx,ny,nxtl,nxtr,BoRiS,wtarg_max)
      use b2mod_types
      use b2mod_geo
      use b2mod_rates
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_external
      use b2mod_constants
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist

      IMPLICIT NONE
#     include <DIMENSIONS.F>

      integer nx, ny, nxtl, nxtr
      real (kind=R8) :: BoRiS

      INTEGER,PARAMETER :: NTRGM=DEF_NCUT

      REAL (KIND=R8),INTENT(OUT) :: wtarg_max(NTRGM)

      REAL (KIND=R8),SAVE :: STRG(0:DEF_NYD-1,NTRGM) !AREAS OF SURFACE ELEMENTS
      INTEGER,SAVE :: LTRGPOS(NTRGM) !POLOIDAL INDEX FOR EACH TARGET
      INTEGER,SAVE :: NTRGS,         !NUMBER OF TARGETS
     i                SHIFT(NTRGM)
      REAL (KIND=R8) :: WTARG(0:ny-1,NTRGM),wnutrg(0:ny-1,NTRGM),
     r                  wnutpr(0:ny-1,NTRGM)
      INTEGER :: K,J,L,IY,IS
      INTEGER :: geometryType
      REAL (KIND=R8) :: HLP
      CHARACTER*32 :: FRM
      LOGICAL :: target_east,target_west

      LOGICAL :: IFIRST
      DATA IFIRST /.TRUE./

      if(ifirst) then
C PREPARING GEOMETRY DATA
        geometryType = geometryId(nnreg, isymm, periodic_bc, topcut)
        if(nncut.eq.0) then
          ntrgs=0
          if (geometryType.eq.GEOMETRY_LINEAR .or.
     &        geometryType.eq.GEOMETRY_CYLINDER) then
            target_east = .false.
            target_west = .false.
            do k = 1, nbc
              if(bcchar(k).eq.'E'.and.bcpos(k).eq.-1) then
                target_east = bcene(k).eq.3 .or. bcene(k).eq.12 .or.
     &                        bcene(k).eq.15
              end if
              if(bcchar(k).eq.'W'.and.bcpos(k).eq.nx) then
                target_west = bcene(k).eq.3 .or. bcene(k).eq.12 .or.
     &                        bcene(k).eq.15
              end if
            end do
            if (target_west) then
              ntrgs=1
              ltrgpos(1)=0
            end if
            if (target_east) then
              ntrgs=ntrgs+1
              ltrgpos(ntrgs)=nx-1
            end if
          end if
        else if(nncut.eq.1) then
          ntrgs=2
          ltrgpos(1)=0
          ltrgpos(2)=nx-1
        else
          ntrgs=4
          ltrgpos(1)=0
          ltrgpos(2)=nxtl-1
          ltrgpos(3)=nxtr+1
          ltrgpos(4)=nx-1
        end if
C AREAS OF THE SURFACE ELEMENTS
        l=0
        do j=1,ntrgs !{
          k=ltrgpos(j)
          wtarg_max(j)=0.0_R8
          do iy=0,ny-1
            strg(iy,j)=gs(k+l,iy,0)
          end do
          call b2neut_ind(2,k+l,1,ny,shift(j)) !index range in the Eirene arrays
          if(shift(j).eq.-1) then
            write(*,*) "Error in b2usr_loads: the EIRENE data "//
     w                 "is apparently missing. ",
     w                 "No neutral fluxes will be calculated."
          else if(shift(j).lt.0) then
            write(*,*) "Error in b2usr_loads: the EIRENE data on "//
     w                 "neutral fluxes is incorrect. ",
     w                 "No neutral fluxes will be calculated."
          else
            call b2neut_surfcheck(strg(0,j),j,ny,shift(j))
          end if
          l=1-l
        end do !}
        ifirst=.false.
      end if

      write(*,*) "Calculating incident power flux in b2usr_loads"

      l=0
      do j=1,ntrgs !{
        k=ltrgpos(j)
        call b2neut_eng(wnutrg(0,j),wnutpr(0,j),
     .                  ny,shift(j))
        wtarg_max(j)=0.0_R8
        do iy=0,ny-1
          hlp=abs(fht(k+l,iy,0,0)+fhi_ext(k+l,iy,0,0)+
     .            fht(k+l,iy,0,1)+fhi_ext(k+l,iy,0,1))
          do is=0,ns_ext-1
            hlp=hlp+abs((fa_ext(k+l,iy,0,0,is)
     .                  +fa_ext(k+l,iy,0,1,is))*(pt_ext(k,iy,is)*ev+
     .                 (1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .                  ua_ext(k,iy,is)**2+ta_ext(k,iy,is))))
          end do
          hlp=hlp/strg(iy,j)
          wtarg(iy,j)=hlp+wnutrg(iy,j)-wnutpr(iy,j)
          if(wtarg(iy,j).gt.wtarg_max(j))
     .     wtarg_max(j)=wtarg(iy,j)
csw 02jan2017 corrected flip of l variable (i.e.put it AFTER iy loop)
        end do
        l=1-l
      end do !}

      if(debug_flags(81).gt.0) then
        write(frm,700) ntrgs
  700   format('(',i3,'(2x,e15.7))')
C total power
        open(3,file='targb2.datv')
        do iy=0,ny-1
          write(3,frm) (wtarg(iy,j),j=1,ntrgs)
        end do
        close(3)
C neutrals
        open(3,file='targb2n.datv')
        do iy=0,ny-1
          write(3,frm) (wnutrg(iy,j),j=1,ntrgs)
        end do
        close(3)
C plasma
        do j=1,ntrgs
          do iy=0,ny-1
            wtarg(iy,j)=wtarg(iy,j)-wnutrg(iy,j)
          end do
        end do
        open(3,file='targb2pl.datv')
        do iy=0,ny-1
          write(3,frm) (wtarg(iy,j),j=1,ntrgs)
        end do
        close(3)
      end if
      return
      end subroutine b2usr_loads
