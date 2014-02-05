      subroutine b2prdl
c
c  version : 14.05.96 23:27
c
c=====================================================
c*** b2prdl: post-processor for B2-Eirene
c***
c***  Produces files with
c***   (a) distribution of the radiation source for calculation of the
c***       wall loading (a sum of radiation from neutrals, impurities,
c***       and Bremsstrahlung)
c***   (b) distribution of the power coming to the targets with
c***       the plasma particles (parallel kinetic + potential energy)
c***
c*** Is to be called from b2plot!
c***
c*** Output files:
c***
c***  radsrc.b2e   with source data
c***  targload.dat with loading data
c=====================================================
      use b2mod_types
      use b2mod_b2plot
      use b2mod_sources
      implicit none
      integer lun
      parameter (lun=96)
c*** Local variables
c
      integer ix,iy,is
      real (kind=R8) ::
     *   u
c=====================================================
      print *,'Starting b2prdl...'
c*** Write down the radiation data
c
      open(lun,file='radsrc.b2e')
      write(lun,'(3x,a6,9x,a2,15x,a2,14x,a4)') 'ix  iy','cr','cz','Wrad'
      do ix=1,nx
          do iy=1,ny
              u = 0.0_R8
              do is=0,ns-1
                  u=u+rqrad(ix,iy,is)+rqbrm(ix,iy,is)
              end do
              u=1.e-6_R8*u
              write (lun,700) ix,iy,cr(ix,iy),cz(ix,iy),u
  700         format(1x,2i4,3e17.7)
          end do
      end do
c
      end
