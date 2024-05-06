      module b2mod_astra_to_b2                                !srv 28.07.11 {
       use b2mod_types
       implicit none
       
       integer, save :: ny_astra = 0                          ! number data points from ASTRA
       real (kind=R8), save :: Ratb   = 0.0_R8                ! ASTRA transport barrier        [m]
       real (kind=R8), save :: RTOR                           ! ASTRA                          [m]
       real (kind=R8), save :: ZRD45                          ! ASTRA                          [m]
       real (kind=R8), save :: crb_astra                      ! coordinate of boundary for ASTRA    [m]
       real (kind=R8), save :: neb_astra                      ! boundary electron density           [1/m**3]
       real (kind=R8), save :: teb_astra                      ! boundary electron temperature       [keV]
       real (kind=R8), save :: tib_astra                      ! boundary electron temperature       [keV]
       real (kind=R8), save :: fneb_astra                     ! boundary electron total flux        [1/s]
       real (kind=R8), save :: fheb_astra                     ! boundary electron total heat flux   [W]
       real (kind=R8), save :: fhib_astra                     ! boundary ion total heat flux        [W]
       real (kind=R8), save :: fSeNBI_astra                   ! boundary electron source due to NBI [1/s]
       real (kind=R8), allocatable, save :: Rmaj_astra(:)     ! major radius                        [m]
       real (kind=R8), allocatable, save :: Rmajc_astra(:)    ! major radius                        [m]
       real (kind=R8), allocatable, save :: a_astra(:)        ! major radius                        [m]
       real (kind=R8), allocatable, save :: rho_astra(:)      ! major radius                        [m]
       real (kind=R8), allocatable, save :: PsiN_astra(:)     ! major radius                         
       real (kind=R8), allocatable, save :: sqrtPsiN_astra(:) ! major radius                         
       real (kind=R8), allocatable, save :: vol_astra(:)      ! volume                              [m**3]
       real (kind=R8), allocatable, save :: ne_astra(:)       ! electron density                    [1/m**3]
       real (kind=R8), allocatable, save :: te_astra(:)       ! electron temperature                [keV]
       real (kind=R8), allocatable, save :: ti_astra(:)       ! ion temperature                     [keV]
       real (kind=R8), allocatable, save :: Ge_astra(:)       !                                     [1/s]
       real (kind=R8), allocatable, save :: qe_astra(:)       ! electron total heat flux            [MW]
       real (kind=R8), allocatable, save :: qi_astra(:)       ! ion total heat flux                 [MW]
       real (kind=R8), allocatable, save :: Se_astra(:)       ! ion total heat flux                 [m^-3s-1]
       real (kind=R8), allocatable, save :: SeNBI_astra(:)    ! ion total heat flux                 [s^-1]
       real (kind=R8), allocatable, save :: Pe_astra(:)       ! ion total heat flux                 [MW/m^3]
       real (kind=R8), allocatable, save :: Pi_astra(:)       ! ion total heat flux                 [MW/m^3]
       real (kind=R8), allocatable, save :: D_astra(:)        ! ion total heat flux                 [m^2/s]
       real (kind=R8), allocatable, save :: chi_e_astra(:)    ! electron heat conduct.              [m**2/s]
       real (kind=R8), allocatable, save :: chi_i_astra(:)    ! ions heat conduct.                  [m**2/s]
       real (kind=R8), allocatable, save :: chi_N_astra(:)    ! neoclass. ions heat conduct.        [m**2/s]
       real (kind=R8), allocatable, save :: qi_neo_astra(:)   ! neoclass. ions heat flux.           [MW]
       real (kind=R8), allocatable, save :: Vcovv_astra(:)    ! convective anomalous velocity       [m/s]
       real (kind=R8), allocatable, save :: PeNB_astra(:)     ! electron NB source                  [MW/m^3]
       real (kind=R8), allocatable, save :: PiNB_astra(:)     ! ion NB source                       [MW/m^3]
       real (kind=R8), allocatable, save :: inv_hi_astra(:)   ! inverse Kaveeva ions heat conduct.  [s/m**2]

       real (kind=R8), allocatable, save :: bbaverage2(:)     ! 
       real (kind=R8), allocatable, save :: sheNB_astra(:,:)
       real (kind=R8), allocatable, save :: shiNB_astra(:,:)

       contains

       subroutine alloc_b2mod_astra_to_b2(mpg, nxd, nyd)
       use b2us_map
       implicit none
       type(mapping), intent(in) :: mpg
       integer, intent(in) :: nxd, nyd

       allocate (Rmaj_astra(-1:nyd))
       allocate (Rmajc_astra(-1:nyd))
       allocate (a_astra(-1:nyd))
       allocate (rho_astra(-1:nyd))
       allocate (PsiN_astra(-1:nyd))
       allocate (sqrtPsiN_astra(-1:nyd))
       allocate (vol_astra(-1:nyd))
       allocate (ne_astra(-1:nyd))
       allocate (te_astra(-1:nyd))
       allocate (ti_astra(-1:nyd))
       allocate (Ge_astra(-1:nyd))
       allocate (qe_astra(-1:nyd))
       allocate (qi_astra(-1:nyd))
       allocate (Se_astra(-1:nyd))
       allocate (SeNBI_astra(-1:nyd))
       allocate (Pe_astra(-1:nyd))
       allocate (Pi_astra(-1:nyd))
       allocate (D_astra(-1:nyd))
       allocate (chi_e_astra(-1:nyd))
       allocate (chi_i_astra(-1:nyd))
       allocate (chi_N_astra(-1:nyd))
       allocate (qi_neo_astra(-1:nyd))
       allocate (Vcovv_astra(-1:nyd))
       allocate (PeNB_astra(-1:nyd))
       allocate (PiNB_astra(-1:nyd))
       allocate (bbaverage2(-1:nyd))
       allocate (inv_hi_astra(-1:nyd))
       allocate (sheNB_astra(mpg%nCv,0:3))
       allocate (shiNB_astra(mpg%nCv,0:3))

       Rmaj_astra     = 0.0_R8
       Rmajc_astra    = 0.0_R8
       a_astra        = 0.0_R8
       rho_astra      = 0.0_R8
       PsiN_astra     = 0.0_R8
       sqrtPsiN_astra = 0.0_R8
       vol_astra    = 0.0_R8
       ne_astra     = 0.0_R8
       te_astra     = 0.0_R8
       ti_astra     = 0.0_R8
       Ge_astra     = 0.0_R8
       qe_astra     = 0.0_R8
       qi_astra     = 0.0_R8
       Se_astra     = 0.0_R8
       SeNBI_astra  = 0.0_R8
       Pe_astra     = 0.0_R8
       Pi_astra     = 0.0_R8
       D_astra      = 0.0_R8
       chi_e_astra  = 0.0_R8
       chi_i_astra  = 0.0_R8
       chi_N_astra  = 0.0_R8
       qi_neo_astra = 0.0_R8
       Vcovv_astra  = 0.0_R8
       PeNB_astra   = 0.0_R8
       PiNB_astra   = 0.0_R8
       inv_hi_astra = 0.0_R8
       write (*,*) 'Allocated space for b2mod_astra_to_b2'
       end subroutine alloc_b2mod_astra_to_b2
       
       subroutine dealloc_b2mod_astra_to_b2

       deallocate(Rmaj_astra)
       deallocate(Rmajc_astra)
       deallocate(a_astra)
       deallocate(rho_astra)
       deallocate(PsiN_astra)
       deallocate(sqrtPsiN_astra)
       deallocate(vol_astra)
       deallocate(ne_astra)
       deallocate(te_astra)
       deallocate(ti_astra)
       deallocate(Ge_astra)
       deallocate(qe_astra)
       deallocate(qi_astra)
       deallocate(Se_astra)
       deallocate(SeNBI_astra)
       deallocate(Pe_astra)
       deallocate(Pi_astra)
       deallocate(D_astra)
       deallocate(chi_e_astra)
       deallocate(chi_i_astra)
       deallocate(chi_N_astra)
       deallocate(qi_neo_astra)
       deallocate (Vcovv_astra)
       deallocate (PeNB_astra)
       deallocate (PiNB_astra)
       deallocate(bbaverage2)
       deallocate (inv_hi_astra)
       deallocate (sheNB_astra)
       deallocate (shiNB_astra)

       end subroutine dealloc_b2mod_astra_to_b2

       subroutine read_astra_to_b2(nyd, nread)
       use b2mod_subsys
       implicit none
       integer :: nyd, nread
#ifdef ASTRA
       integer :: iy
       real (kind=R8) :: value
       logical :: ex, LRTOR, LZRD45
       character(64) :: fname
       character(6)  :: chname
       real (kind=R8) :: b(1:nyd+2),c(1:nyd+2),d(1:nyd+2), dseval
#endif
       external xerrab

       call subini ('read_astra_to_b2')
       call xertst (0.le.nyd, 'faulty argument nyd')
       call xertst (1.le.nread, 'faulty argument nread')
#ifndef ASTRA
       call xerrab ('Compile with -DASTRA option to couple to ASTRA!')
#else
*   .. get data from ASTRA
       read(nread,"(es22.15)")
       do iy=-1,nyd
         read(nread,10,end=100) Rmaj_astra(iy), Rmajc_astra(iy),
     &   a_astra(iy), rho_astra(iy), PsiN_astra(iy), sqrtPsiN_astra(iy),
     &   vol_astra(iy), ne_astra(iy), te_astra(iy), ti_astra(iy),
     &   Ge_astra(iy), qe_astra(iy), qi_astra(iy), Se_astra(iy),
     &   SeNBI_astra(iy), Pe_astra(iy), Pi_astra(iy), D_astra(iy),
     &   chi_e_astra(iy), chi_i_astra(iy), chi_N_astra(iy),
     &   qi_neo_astra(iy), Vcovv_astra(iy), PeNB_astra(iy),
     &   PiNB_astra(iy) 
         ny_astra = ny_astra +1
       enddo

  100  write(*,*) ny_astra,' lines from astra_to_b2 have been read'
*   ..the simplest linear interpolation   
       do iy = 1, ny_astra-1
        b(iy)=(Rmaj_astra(iy)-Rmaj_astra(iy-1))/
     &           (a_astra(iy)-a_astra(iy-1))
       enddo
       b(ny_astra) = b(ny_astra-1)
       c = 0.0_R8
       d = 0.0_R8
!       do i = -1, n 
!        fun(i) = dseval(ny_astra, coord(i), a_astra, Rmaj_astra, b, c, d)
!        write(*,*) coord(i), fun(i)
!       enddo

       neb_astra    = ne_astra(-1)
       teb_astra    = te_astra(-1)
       tib_astra    = ti_astra(-1)
       fneb_astra   = Ge_astra(-1)                                    !srv 04.10.11
       fheb_astra   = qe_astra(-1)*1.0e6_R8
       fhib_astra   = qi_astra(-1)*1.0e6_R8
       fSeNBI_astra = SeNBI_astra(-1)
*   .. get Ratb
       fname = 'datahead.tmpl'
       call find_file(fname,ex)
       if(ex) then
         LRTOR =.false.
         LZRD45=.false.
         open(70,file=fname)
         read(70,30)
         read(70,30)
         do while (.not.LRTOR .or. .not.LZRD45)
          read(70,30,end=200) chname,value
!          write(*,*) chname,value
          if(trim(chname).eq.'RTOR') then
           LRTOR = .true.
           RTOR = value
           Ratb = Ratb + RTOR
!           write(*,30) 'RTOR =',RTOR
          endif
          if(trim(chname).eq.'ZRD45') then
           LZRD45 = .true.
           ZRD45 = value
           Ratb = Ratb + ZRD45
!           write(*,30) 'ZRD45=',ZRD45
          endif
         enddo

  200    if(LRTOR) then
          write(*,*) 'RTOR has been read from file ',
     &                trim(fname),'  ', RTOR
         else
          write(*,*) 'RTOR has not been read from file ',
     &                trim(fname)
          call xerrab ('RTOR not read from ASTRA file')
         endif
         if(LZRD45) then
          write(*,*) 'ZRD45 has been read from file ',
     &                trim(fname),'  ', ZRD45
         else
          write(*,*) 'ZRD45 has not been read from file ',
     &                trim(fname)
          call xerrab ('ZRD45 not read from ASTRA file')
         endif

         close(70)
       else
         write(*,*) ' There is no file ', trim(fname)
         call xerrab ('ASTRA file not found!')
       endif

       crb_astra=dseval(ny_astra, ZRD45, a_astra, Rmaj_astra, b, c, d)

   10  format(25e13.5)
c   10   format(91x,e13.5,13x,e13.5,13x,e13.5,e13.5,13x,e13.5,65x,2e13.5)
c   20   format(117x,e13.5,130x,2e13.5)
   30  format(a6,9x,f7.4)
#endif

       call subend ()
       return
       end subroutine read_astra_to_b2

       subroutine write_to_check_astra_to_b2(nwrite)
       implicit none
       integer :: nwrite
       integer :: iy

       write(nwrite,"(es22.15)")
       write(nwrite,10) (Rmaj_astra(iy), Rmajc_astra(iy),
     &   a_astra(iy), rho_astra(iy), PsiN_astra(iy), sqrtPsiN_astra(iy),
     &   vol_astra(iy), ne_astra(iy), te_astra(iy), ti_astra(iy),
     &   Ge_astra(iy), qe_astra(iy), qi_astra(iy), Se_astra(iy),
     &   SeNBI_astra(iy), Pe_astra(iy), Pi_astra(iy), D_astra(iy),
     &   chi_e_astra(iy), chi_i_astra(iy), chi_N_astra(iy),
     &   qi_neo_astra(iy), Vcovv_astra(iy), PeNB_astra(iy),
     &   PeNB_astra(iy), iy=-1,ny_astra-1) 

       write(nwrite,30) 'RTOR =',RTOR
       write(nwrite,30) 'ZRD45=',ZRD45
       write(nwrite,30) 'Ratb = RTOR + ZRD45 =',Ratb
       write(nwrite,30) 'crb_astra =',crb_astra
   10  format(25e13.5)
c   10   format(91x,e13.5,13x,e13.5,13x,e13.5,e13.5,13x,e13.5,65x,2e13.5)
c   20   format(117x,e13.5,130x,2e13.5)
   30  format(a,9x,f7.4)

      end subroutine write_to_check_astra_to_b2

      end module b2mod_astra_to_b2                           !srv 28.07.11 }
