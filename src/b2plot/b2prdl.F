      subroutine b2prdl
c
c  version : 14.05.96 23:27
c
c=====================================================
c*** b2prdl: postprocessor for B2-Eirene
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
c=====================================================
      use b2mod_types
      use b2us_data
#ifdef B25_EIRENE
      use eirmod_comusr
      use eirmod_wneutrals
#endif
      implicit none
      integer lun
      parameter (lun=96)
c*** Local variables
c
      integer iCv,is
      real (kind=R8) ::
     *   u
c=====================================================
      print *,'Starting b2prdl...'
c*** Write down the radiation data
c
      open(lun,file='radsrc.b2e')
      write(lun,'(3x,a3,9x,a2,15x,a2,14x,a4)') 'iCv','cr','cz','Wrad'
      do iCv=1,mpg%nCi
        u = 0.0_R8
        do is=0,state%ns-1
          u=u+state%srw%rqrad(iCv,is)+state%srw%rqbrm(iCv,is)
        end do
#ifdef B25_EIRENE
        do is=1,natmi
          u=u-eneutrad(iCv,is,0)
        end do
        do is=1,nmoli
          u=u-emolrad(iCv,is,0)
        end do
        do is=1,nioni
          u=u-eionrad(iCv,is,0)
        end do
#endif
        u=1.e-6_R8*u/geo%cvVol(iCv)
        write (lun,700) iCv,geo%cvX(iCv),geo%cvY(iCv),u
  700   format(1x,i6,3e17.7)
      end do
      print*,'Wrote file radsrc.b2e'
c
      return
      end subroutine b2prdl
