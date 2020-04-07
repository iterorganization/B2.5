      module b2us_plasma
 
      use b2mod_types
      implicit none
      
 !     private
      public :: B2state, B2Plasma, B2Derivatives, B2coeff, B2rates,
     .          B2ratesWork, B2Source, B2SourceWork, B2PlasmaSnapShot,
     .          B2Diagnostic, B2Update
 
  ! Plasma state derived types

  ! Main plasma state parameters
      type B2Plasma
        integer :: ns ! total number of species (1:ns or 0:ns-1)
        real(R8), dimension(:,:), allocatable :: na, ua  ! second index is species
        real(R8), dimension(:), allocatable :: po, te, ti
      end type B2Plasma

  ! Quantities derived from plasma state
      type B2Derivatives
        ! fluxes/flux densities
        ! dimensions: fc, direction (, species)
        real(R8), allocatable :: fch(:,:), fch_32(:,:), 
     .                           fch_52(:,:), fch_p(:,:),
     .                           fchdia(:,:), fchin(:,:),
     .                           fchvispar(:,:), fchvisper(:,:),
     .                           fchvisq(:,:), fchinert(:,:),
     .                           fchanml(:,:),
     .                           fch_pi_c(:,:), fch_pi_f(:,:)
        real(R8), allocatable :: fni_32(:,:), fni_52(:,:), 
     .                           fni(:,:), fni_he(:,:)
        real(R8), allocatable :: fna(:,:,:), fna_mdf(:,:,:), 
     .                           fna_52(:,:,:), fna_32(:,:,:), 
     .                           fna_52nd(:,:,:), fna_32nd(:,:,:),
     .                           fna_nodrift(:,:,:), fna_he(:,:,:),
     .                           fnaPSch(:,:,:), fna_fcor(:,:,:),
     .                           fna_eir(:,:,:),
     .                           fmo(:,:,:)
        real(R8), allocatable :: fne(:,:), fne_he(:,:), fne_32(:,:), 
     .                           fne_52(:,:), fne_eir(:,:),
     .                           fhe(:,:), fhe_mdf(:,:), fhet(:,:),
     .                           fhePSch(:,:), fhe_eir(:,:)
        real(R8), allocatable :: fhi(:,:), fhi_mdf(:,:), fhit(:,:),          
     .                           fhiPSch(:,:), fhi_eir(:,:)
        real(R8), allocatable :: fhm(:,:,:), fhp(:,:,:), fhj(:,:), 
     .                           fht(:,:)


        ! kinetic energy
        ! dimensions: cv, species
        real(R8), allocatable :: kinrgy(:,:)

        ! numerical coefficients (maybe move to separate object)
        real(R8), dimension(:,:), allocatable :: conc
        real(R8), dimension(:,:), allocatable :: flob, floe, floi, 
     .                                           floe_noc, floi_noc
        real(R8), dimension(:,:,:), allocatable :: conb, cone, coni
        real(R8), dimension(:), allocatable :: fllime, fllimi
        
        ! residuals
        real(R8), dimension(:,:), allocatable :: resmo
        real(R8), dimension(:,:), allocatable :: resco          
        real(R8), dimension(:), allocatable :: respo, reshe, reshi, 
     .                                         resht, resmt

        ! corrections
        real(R8), dimension(:,:), allocatable :: corua, corpa
        real(R8), dimension(:), allocatable :: corut, corpo, cortt, 
     .                                           corte, corti

        ! coupling coefficients
        real(R8), dimension(:,:,:), allocatable :: pcca
        real(R8), dimension(:,:), allocatable :: pccm

        ! derived state quantities
        ! ni(:,0) is ions, ni(:,1) is neutrals, ue is parallel electron velocity
        ! ne2(:) is electron density, ne2 is elect. den. computed w. squared charge rz2
        ! pa(:,0:ns-1) is partial pressure
        real(R8), allocatable :: ne(:), ni(:,:), ue(:), ne2(:), pa(:,:),
     .                           pz(:), lnlam(:)

        ! drift velocities, at cell faces
        ! dimensions: fc, direction, species
        real(R8), allocatable :: uadia(:,:,:), vadia(:,:,:), 
     .                           wadia(:,:,:), vaecrb(:,:,:), 
     .                           vedia(:,:), wedia(:,:), veecrb(:,:),
     .                           facdrift(:), fac_ExB(:), fac_vis(:)


      end type B2Derivatives

  ! Transport coefficients (velocities, diffusivities...)
      type B2Coeff
        ! at faces
        real(R8), dimension(:,:), allocatable :: csig, calf, csig_an,
     .                                           calf_an, csig_cl, 
     .                                           calf_cl
        real(R8), dimension(:,:,:,:), allocatable :: csigin
        real(R8), dimension(:,:), allocatable :: chve, chce
        real(R8), dimension(:,:), allocatable :: chvi, chci
        real(R8), dimension(:), allocatable :: chvemx, chvimx
        real(R8), allocatable :: cvla(:,:,:), cdna(:,:,:), 
     .                           cdpa(:,:,:), cvsa(:,:,:),
     .                           cvlahz(:,:,:), cdnahz(:,:,:), 
     .                           cdpahz(:,:,:), cvsahz(:,:,:),
     .                           cddi(:,:,:), cvsahz_cl(:,:,:),
     .                           chcb(:,:,:), cvsa_an(:,:,:),
     .                           cvmahz(:,:,:), cvsahz_eff(:,:,:)
        real(R8), allocatable :: cthe(:,:), cthi(:,:), cvsa_cl(:,:,:)
        real(R8), allocatable :: ceqp(:)
          
        real(R8), allocatable :: fllim0fhi(:,:,:), fllimvisc(:,:),
     .                           fllim0fna(:,:,:)

        real(R8), allocatable :: vsaf_cl(:,:,:)

        real(R8), allocatable, dimension(:) :: sig0, hce0, hci0, alf0
        real(R8), allocatable, dimension(:,:) :: dpa0, dna0, vsa0, 
     .                                           hcib
        real(R8), allocatable, dimension(:,:,:) :: vla0, vma0

        real(R8), allocatable, dimension(:) :: 
     .                           fllim_ki, fllim_ke, fllim_al, 
     .                           fllim_al_c, f_luc_ke, f_luc_ki,
     .                           f_luc_et, f_luc_sg, f_luc_al,
     .                           alfx_c, sigx_c


cwdk move to dv
!        ! velocities at cell centers
!        ! dimensions: cv, direction, species
!        real(R8), allocatable :: uadia(:,:,:,:), wadia(:,:,:,:), 
!     .                           vaecrb(:,:,:,:)

      end type B2Coeff
  
      type B2Source 
        real(R8), dimension(:,:), allocatable :: sch, she, shi, sne
        real(R8), dimension(:,:,:), allocatable :: smo, smq, sna
      end type B2Source

      type B2SourceEir
        real(R8), dimension(:,:), allocatable :: sch, she, shi, sne
        real(R8), dimension(:,:,:), allocatable :: smo, smq, sna
      end type B2SourceEir

      type B2SourceWork
        real(R8), dimension(:,:), allocatable :: sch0, she0, shi0, sne0
        real(R8), dimension(:,:,:), allocatable :: smo0, smq0, sna0         
        real(R8), dimension(:), allocatable :: 
     .            b2stbm_sch, b2stbm_she, b2stbm_shi, b2stbm_sne
        real(R8), dimension(:,:), allocatable :: 
     .            b2stbm_smo, b2stbm_smq, b2stbm_sna       
      end type B2SourceWork

  ! b2mod_rates
      type B2Rates

        integer :: nscx, nscxmax = 3
        integer,allocatable :: iscx(:)

        ! nCv, dir, ns, 0:nscxmax-1
        real(R8), allocatable :: rlcx(:,:,:,:)
        ! nCv, dir, is
        real(R8), allocatable, dimension(:,:,:) :: rlqa, rlrd, rlbr, 
     .                                             rlra, rlsa, rlza, 
     .                                             rlz2, rlpt, rlpi
        ! nCv, is
        real(R8), allocatable, dimension(:,:) :: rza, rz2, rpt, rpi
      end type B2Rates

  ! Type for data transfer in b2sral
      type B2RatesWork          
        real(R8), allocatable :: rsa(:,:), rra(:,:), rqa(:,:), rrd(:,:),
     .                           rbr(:,:), rcx(:,:) 
      end type B2RatesWork

  ! A plasma state snapshot for timestepping. Includes a combination of quantities from above types
      type B2PlasmaSnapshot
        real(R8), dimension(:,:), allocatable :: na, ua  ! second index is species
        real(R8), dimension(:), allocatable :: po, te, ti
        real(R8), allocatable :: ne(:), ni(:,:)
        real(R8), allocatable :: fch(:,:), fna(:,:,:), 
     .                           fhi(:,:), fhe(:,:)
        real(R8), allocatable :: kinrgy(:,:)
      end type B2PlasmaSnapshot
 
      type B2Diagnostic
        ! evaluated sources for the individual equations
        real(R8), dimension(:,:), allocatable :: srcna, srcmo
        real(R8), dimension(:), allocatable :: srcpo, srche, srchi, 
     .                                         srcmt

        ! b2mod_residuals
        real(R8), allocatable ::  aresco(:), aresmo(:)
        real(R8) :: aresmt, areshe, areshi, arespo         

        ! b2nmndt: corrections
        real(R8), allocatable :: acorpa(:), acorua(:)
        real(R8) :: acorut, acorte, acorti, acorpo

        ! b2mod_diag
        real (R8), allocatable ::  rescoreg(:,:) , resmoreg(:,:), 
     .                             reshereg(:), reshireg(:)
      end type B2Diagnostic

      type B2Update
        real(R8), dimension(:,:), allocatable :: ua, na, pa
        real(R8), dimension(:), allocatable :: po, te, ti
      end type B2Update

  ! Full set of state variables  
      type B2State
        integer :: ns           ! number of species
        integer :: ismain       ! index of main species
        integer :: ismain0      ! index of main neutral species

        type(B2Plasma) :: pl 
        type(B2Coeff) :: co
        type(B2Derivatives) :: dv
        type(B2Source) :: sr
        type(B2SourceEir) :: sr_eir
        type(B2SourceWork) :: srw
        type(B2Rates) :: rt
        type(B2RatesWork) :: rtw
        type(B2PlasmaSnapshot) :: psnl, psnc
        type(B2Diagnostic) :: diag
        type(B2Update) :: update
      end type B2State


 ! External plasma
 ! As simplified type containg all data on possible external species
 ! To be considered: move to a separate module b2us_plasma_ext?
      type B2StateExt

        character*13, dimension(:), allocatable :: text(:)
        integer :: ns
        logical, dimension(:), allocatable :: is_neutral(:)
        logical :: initialised_externals = .false.
        logical :: allocated_externals = .false.

        real(R8), dimension(:), allocatable :: zn, am,
     .                                         ne, ne2, ue
        real(R8), dimension(:,:), allocatable :: za, za2, pt, 
     .                                           na, ni, ua, ta,
     .                                           fhi
        real(R8), dimension(:,:,:), allocatable :: fa

        real(R8), dimension(:), allocatable ::  sne, she, shi, sch
        real(R8), dimension(:,:), allocatable :: sna, smo

      end type B2StateExt

      contains

      subroutine createB2State (nCv, nFc, ns, st)
      implicit none
      integer, intent (in) :: nCv, nFc, ns
      type(B2State), intent(inout) :: st

      call createB2Plasma( nCv, ns, st % pl )
      call createB2Coeff( nCv, nFC, ns, st % co )
      call createB2Derivatives( nCv, nFc, ns, st % dv )
      call createB2Source( nCv, ns, st % sr )   
      call createB2SourceWork( nCv, ns, st % srw )   
      call createB2Rates( nCv, ns, st % rt )
      call createB2RatesWork( nCv, ns, st % rtw )
      call createB2PlasmaSnapshot( nCv, nFc, ns, st % psnl )
      call createB2PlasmaSnapshot( nCv, nFc, ns, st % psnc )
!WG_TODO      call createB2Diagnostic( nCv, ns, nnreg, st % diag )
      call createB2Update( nCv, ns, st % update )

      return
      end subroutine createB2State

      subroutine destroyB2State ( st )
      implicit none
      type (B2State), intent (inout) :: st

      call destroyB2Plasma( st % pl )
      call destroyB2Coeff( st % co )
      call destroyB2Derivatives( st % dv )
      call destroyB2Source( st % sr )   
      call destroyB2SourceWork( st % srw )   
      call destroyB2Rates( st % rt )
      call destroyB2RatesWork( st % rtw )
      call destroyB2PlasmaSnapshot(  st % psnl )
      call destroyB2PlasmaSnapshot( st % psnc )
      call destroyB2Diagnostic( st % diag )
      call destroyB2Update( st % update )

      return
      end subroutine destroyB2State


      subroutine createB2Plasma( nCv, ns, plasma )
      integer, intent(in) :: nCv, ns
      type(B2Plasma), intent(inout) :: plasma

      if (allocated(plasma % na)) return

      plasma % ns = ns

      allocate( plasma % na( nCv, 0 : ns - 1 ) )
      allocate( plasma % ua( nCv, 0 : ns - 1 ) )
      allocate( plasma % po( nCv ) )
      allocate( plasma % te( nCv ) )
      allocate( plasma % ti( nCv ) )

      plasma % na = 0._R8
      plasma % ua = 0._R8
      plasma % po = 0._R8
      plasma % te = 0._R8
      plasma % ti = 0._R8

      return      
      end subroutine createB2Plasma


      subroutine destroyB2Plasma( plasma )
      type(B2Plasma), intent(inout) :: plasma
      
      plasma % ns = 0
      deallocate( plasma % na )   
      deallocate( plasma % ua )   
      deallocate( plasma % po )   
      deallocate( plasma % te )   
      deallocate( plasma % ti )   

      return      
      end subroutine destroyB2Plasma
      
      
      subroutine createB2Coeff( nCv, nFC, ns, coeff )
      integer, intent(in) :: nCv, nFc, ns
      type(B2Coeff), intent(inout) :: coeff

      if (allocated(coeff % calf)) return

      allocate( coeff % calf( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % calf_cl( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % calf_an( nFc, 0:1 ), SOURCE=0._R8 )
      
      allocate( coeff % csig( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % csig_cl( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % csig_an( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % csigin( nFc, 0:1, 0:ns-1, 0:ns-1 ),SOURCE=0._R8)
      
      allocate( coeff % chve( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % chce( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % chvemx( nFc ), SOURCE=0._R8 )
      
      allocate( coeff % chvi( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % chci( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( coeff % chcb( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % chvimx( nFc ), SOURCE=0._R8 )
      
      allocate( coeff % cvla( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cdna( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cdpa( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cvsa( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cvsa_cl( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cvsa_an( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )

      allocate( coeff % cvlahz( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cdnahz( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cdpahz( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cvsahz( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cvsahz_cl( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cvsahz_eff( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cvmahz( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )

      allocate( coeff % vsaf_cl( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )

      allocate( coeff % cddi( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )

      allocate( coeff % cthe( nFc, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % cthi( nFc, 0:ns-1 ), SOURCE=0._R8 )
      
      allocate( coeff % ceqp( nCv ), SOURCE=0._R8 )

      allocate( coeff % fllim0fna( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )      
      allocate( coeff % fllim0fhi( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff % fllimvisc( nFc, 0:ns-1 ), SOURCE=0._R8 )


! moved here from b2coeffwork
      
      allocate( coeff%dna0( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff%dpa0( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff%vla0( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff%vma0( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff%vsa0( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff%hci0( nCv ), SOURCE=0._R8 )
      allocate( coeff%hcib( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( coeff%hce0( nCv ), SOURCE=0._R8 )
      allocate( coeff%sig0( nCv ), SOURCE=0._R8 )
      allocate( coeff%alf0( nCv ), SOURCE=0._R8 )

      allocate( coeff%fllim_ki( nFc ), SOURCE=0._R8 )
      allocate( coeff%fllim_ke( nFc ), SOURCE=0._R8 )
      allocate( coeff%fllim_al( nFc ), SOURCE=0._R8 )
      allocate( coeff%fllim_al_c( nCv ), SOURCE=0._R8 )
      allocate( coeff%f_luc_ke( nFc ), SOURCE=0._R8 )
      allocate( coeff%f_luc_ki( nFc ), SOURCE=0._R8 )
      allocate( coeff%f_luc_et( nFc ), SOURCE=0._R8 )
      allocate( coeff%f_luc_sg( nFc ), SOURCE=0._R8 )
      allocate( coeff%f_luc_al( nFc ), SOURCE=0._R8 )
      allocate( coeff%alfx_c( nCv ), SOURCE=0._R8 )
      allocate( coeff%sigx_c( nCv ), SOURCE=0._R8 )

      return
      end subroutine createB2Coeff

      
      subroutine destroyB2Coeff( coeff )
      type(B2Coeff), intent(inout) :: coeff
      
      deallocate( coeff % csig )
      deallocate( coeff % csigin )
      deallocate( coeff % csig_an )
      deallocate( coeff % csig_cl )
      deallocate( coeff % calf )
      deallocate( coeff % calf_an )
      
      deallocate( coeff % chve )
      deallocate( coeff % chce )
      deallocate( coeff % chvemx )
      
      deallocate( coeff % chvi )
      deallocate( coeff % chci )
      deallocate( coeff % chcb )
      deallocate( coeff % chvimx )
      
      deallocate( coeff % cvla )
      deallocate( coeff % cdna )
      deallocate( coeff % cdpa )
      deallocate( coeff % cvsa )
      deallocate( coeff % cvsa_cl )
      deallocate( coeff % cvsa_an )

      deallocate( coeff % cvlahz )
      deallocate( coeff % cdnahz )
      deallocate( coeff % cdpahz )
      deallocate( coeff % cvsahz )
      deallocate( coeff % cvsahz_cl )
      deallocate( coeff % cvsahz_eff )
      deallocate( coeff % cvmahz )

      deallocate( coeff % vsaf_cl )

      deallocate( coeff % cddi )

      deallocate( coeff % cthe )
      deallocate( coeff % cthi )

      
      deallocate( coeff % ceqp )

      deallocate( coeff % fllim0fna )      
      deallocate( coeff % fllim0fhi )
      deallocate( coeff % fllimvisc )
      
      deallocate( coeff%dna0 )
      deallocate( coeff%dpa0 )
      deallocate( coeff%vla0 )
      deallocate( coeff%vma0 )
      deallocate( coeff%vsa0 )
      deallocate( coeff%hci0 )
      deallocate( coeff%hcib )
      deallocate( coeff%hce0 )
      deallocate( coeff%sig0 )
      deallocate( coeff%alf0 )

      deallocate( coeff%fllim_ki )
      deallocate( coeff%fllim_ke )
      deallocate( coeff%fllim_al )
      deallocate( coeff%fllim_al_c )
      deallocate( coeff%f_luc_ke )
      deallocate( coeff%f_luc_ki )
      deallocate( coeff%f_luc_et )
      deallocate( coeff%f_luc_sg )
      deallocate( coeff%f_luc_al )
      deallocate( coeff%alfx_c )
      deallocate( coeff%sigx_c )
      
      return
      end subroutine destroyB2Coeff
      
      
      
      subroutine createB2Derivatives( nCv, nFc, ns, dv )
      integer, intent(in) :: nCv, nFc, ns
      type(B2Derivatives), intent(inout) :: dv

      if (allocated(dv % fch)) return
      
! fluxes
      allocate( dv % fch( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fch_32( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fch_52( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fch_p( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fchdia( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fchin( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fchvispar( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fchvisper( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fchvisq( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fchinert( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fchanml( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fch_pi_c( nCv, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fch_pi_f( nFc, 0:1 ), SOURCE=0._R8 )                       
      allocate( dv % fne( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fne_he( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fne_32( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fne_52( nFc, 0:1 ), SOURCE=0._R8 )
      
      allocate( dv % fni( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fni_he( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fni_32( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fni_52( nFc, 0:1 ), SOURCE=0._R8 )
      
      allocate( dv % fna( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fna_mdf( nFc, 0:1 , 0:ns-1 ), SOURCE=0._R8 )    
      allocate( dv % fna_he( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fna_32( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fna_52( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fna_32nd( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fna_52nd( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fnaPSch( nFc, 0:1 , 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fna_fcor( nFc, 0:1 , 0:ns-1 ), SOURCE=0._R8 )        
      allocate( dv % fna_eir( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % fne_eir( nFc, 0:1 ), SOURCE=0._R8 )

      allocate( dv % fna_nodrift( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      
      allocate( dv % fmo( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )

      allocate( dv % fhe( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhe_mdf( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhePSch( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhet( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fhe_eir( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhi( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhi_mdf( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhiPSch( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhit( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fhi_eir( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % fhm( nFc, 0:1, 0:ns-1), SOURCE=0._R8 ) 
      allocate( dv % fhp( nFc, 0:1, 0:ns-1), SOURCE=0._R8 ) 
      allocate( dv % fhj( nFc, 0:1 ), SOURCE=0._R8 ) 
      allocate( dv % fht( nFc, 0:1 ), SOURCE=0._R8 ) 
      
! kinetic energy
      allocate( dv % kinrgy( nCv, 0:ns-1 ), SOURCE=0._R8 )    
      
! derived state quantities
      allocate( dv % ne( nCv ), SOURCE=0._R8 )    
      allocate( dv % ne2( nCv ), SOURCE=0._R8 )    
      allocate( dv % ni( nCv, 0:1 ), SOURCE=0._R8 )    
      allocate( dv % ue( nCv ), SOURCE=0._R8 )
      allocate( dv % pa( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % pz( nCv ), SOURCE=0._R8 )
      allocate( dv % lnlam( nCv ), SOURCE=0._R8 ) 

! drift flows
      allocate( dv % uadia ( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % vadia ( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % wadia ( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % vaecrb( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % vedia ( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % wedia ( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % veecrb( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % facdrift (nFc), SOURCE=0._R8)
      allocate( dv % fac_ExB  (nFc), SOURCE=0._R8)
      allocate( dv % fac_vis  (nFc), SOURCE=0._R8)

      
! numerical coefficients
      allocate( dv % flob (nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % conb (nFc, 0:1, 0:4), SOURCE=0._R8 )
      allocate( dv % conc( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % fllime( nFc ), SOURCE=0._R8 )    
      allocate( dv % floe( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % floe_noc( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % cone( nFc, 0:1, 0:2), SOURCE=0._R8 )
      allocate( dv % fllimi( nFc ), SOURCE=0._R8 )    
      allocate( dv % floi( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % floi_noc( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( dv % coni( nFc, 0:1, 0:2 ), SOURCE=0._R8 )
      
! residuals
      allocate( dv % resmo( nCv, 0:ns-1 ) )
      allocate( dv % resco( nCv, 0:ns-1 ) )
      allocate( dv % respo( nCv ) )
      allocate( dv % reshe( nCv ) )
      allocate( dv % reshi( nCv ) )
      allocate( dv % resht( nCv ) )
      allocate( dv % resmt( nCv ) )
      dv % resmo = 0.0_R8                             !srv 22.05.18 {     
      dv % resco = 0.0_R8
      dv % respo = 0.0_R8
      dv % reshe = 0.0_R8
      dv % reshi = 0.0_R8
      dv % resht = 0.0_R8
      dv % resmt = 0.0_R8                             !srv 22.05.18 }    
      
! corrections
      allocate( dv % corua( nCv, 0:ns-1 ) )
      allocate( dv % corpa( nCv, 0:ns-1 ) )
      allocate( dv % corut( nCv ) )
      allocate( dv % corpo( nCv ) )
      allocate( dv % cortt( nCv ) )
      allocate( dv % corte( nCv ) )
      allocate( dv % corti( nCv ) )
      dv % corua = 0.0_R8                             !srv 22.05.18 {     
      dv % corpa = 0.0_R8
      dv % corut = 0.0_R8
      dv % corpo = 0.0_R8
      dv % cortt = 0.0_R8
      dv % corte = 0.0_R8
      dv % corti = 0.0_R8                             !srv 22.05.18 }    
! coupling coefficients
      allocate( dv % pcca( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( dv % pccm( nCv, 0:1 ), SOURCE=0._R8 )

      return      
      end subroutine createB2Derivatives
      
      
      subroutine destroyB2Derivatives( dv )
      type(B2Derivatives), intent(inout) :: dv
      
      deallocate( dv % fne )
      deallocate( dv % fne_he )
      deallocate( dv % fne_32 )
      deallocate( dv % fne_52 )
      deallocate( dv % fni )
      deallocate( dv % fni_he )
      deallocate( dv % fni_32 )
      deallocate( dv % fni_52 )
      deallocate( dv % fch )
      deallocate( dv % fch_32 )
      deallocate( dv % fch_52 )
      deallocate( dv % fch_p )
      deallocate( dv % fchdia )
      deallocate( dv % fchin )
      deallocate( dv % fchvispar )
      deallocate( dv % fchvisper )
      deallocate( dv % fchvisq )
      deallocate( dv % fchinert )
      deallocate( dv % fchanml )
      deallocate( dv % fch_pi_c )
      deallocate( dv % fch_pi_f )
      deallocate( dv % fna )
      deallocate( dv % fna_mdf )
      deallocate( dv % fna_he )
      deallocate( dv % fna_32 )
      deallocate( dv % fna_52 )
      deallocate( dv % fna_32nd )
      deallocate( dv % fna_52nd )
      deallocate( dv % fnaPSch )
      deallocate( dv % fna_fcor )
      deallocate( dv % fna_nodrift )
      deallocate( dv % fna_eir )
      deallocate( dv % fne_eir )
      
      deallocate( dv % fmo )

      deallocate( dv % ne )
      deallocate( dv % ne2 )
      deallocate( dv % ni )
      deallocate( dv % ue )
      deallocate( dv % pa )
      deallocate( dv % pz )
      deallocate( dv % lnlam)

      deallocate( dv % uadia )
      deallocate( dv % vadia )
      deallocate( dv % wadia )
      deallocate( dv % vaecrb )
      deallocate( dv % vedia )
      deallocate( dv % wedia )
      deallocate( dv % veecrb )
      deallocate( dv % facdrift )
      deallocate( dv % fac_ExB )
      deallocate( dv % fac_vis )

      deallocate (dv % flob )
      deallocate (dv % conb )


! kinetic energy
      deallocate( dv % kinrgy )
      
      deallocate( dv % fhi )
      deallocate( dv % fhi_mdf )
      deallocate( dv % fhiPSch )
      deallocate( dv % fhit )
      deallocate( dv % fhi_eir )
      deallocate( dv % fhe )
      deallocate( dv % fhe_mdf )
      deallocate( dv % fhePSch )
      deallocate( dv % fhet )
      deallocate( dv % fhe_eir )

      deallocate( dv % fhm )
      deallocate( dv % fhp )
      deallocate( dv % fhj )
      deallocate( dv % fht )

      deallocate( dv % conc )
      deallocate( dv % fllime )
      deallocate( dv % floe )
      deallocate( dv % floe_noc )
      deallocate( dv % cone )
      deallocate( dv % fllimi )
      deallocate( dv % floi )
      deallocate( dv % floi_noc )
      deallocate( dv % coni )
      
      deallocate( dv % resmo )
      deallocate( dv % resco )
      deallocate( dv % respo )
      deallocate( dv % reshe )
      deallocate( dv % reshi )
      deallocate( dv % resht )
      deallocate( dv % resmt )
      
      deallocate( dv % corua )
      deallocate( dv % corpa )
      deallocate( dv % corut )
      deallocate( dv % corpo )
      deallocate( dv % cortt )
      deallocate( dv % corte )
      deallocate( dv % corti )
      deallocate( dv % pcca )
      deallocate( dv % pccm )

      return           
      end subroutine destroyB2Derivatives
      
      
      subroutine createB2Source( nCv, ns, sr )
      integer, intent(in) :: nCv, ns
      type(B2Source), intent(inout) :: sr

      if (allocated(sr % sch)) return
      
      allocate( sr % sch( nCv, 0:3 ), SOURCE=0._R8 )
      allocate( sr % she( nCv, 0:3 ), SOURCE=0._R8 )
      allocate( sr % shi( nCv, 0:3 ), SOURCE=0._R8 )
      allocate( sr % smo( nCv, 0:3, 0:ns-1 ), SOURCE=0._R8 )
      allocate( sr % smq( nCv, 0:3, 0:ns-1 ), SOURCE=0._R8 )
      allocate( sr % sna( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( sr % sne( nCv, 0:1 ), SOURCE=0._R8 )

      return      
      end subroutine createB2Source

      subroutine destroyB2Source( sr )
      type(B2Source), intent(inout) :: sr
      
      deallocate( sr % sch )
      deallocate( sr % she )
      deallocate( sr % shi )
      deallocate( sr % smo )
      deallocate( sr % smq )
      deallocate( sr % sna )
      deallocate( sr % sne )

      return
      end subroutine destroyB2Source

      subroutine createB2SourceEir( nCv, ns, nstra, sr_eir )
      integer, intent(in) :: nCv, ns, nstra
      type(B2SourceEir), intent(inout) :: sr_eir

      if (allocated(sr_eir % sch)) return
      
      allocate( sr_eir % sch( nCv, nstra), SOURCE=0._R8 )
      allocate( sr_eir % she( nCv, nstra), SOURCE=0._R8 )
      allocate( sr_eir % shi( nCv, nstra), SOURCE=0._R8 )
      allocate( sr_eir % smo( nCv, 0:ns-1 , nstra), SOURCE=0._R8 )
      allocate( sr_eir % sna( nCv, 0:ns-1 , nstra), SOURCE=0._R8 )
      allocate( sr_eir % sne( nCv, nstra), SOURCE=0._R8 )

      return      
      end subroutine createB2SourceEir
      
      subroutine destroyB2SourceEir( sr_eir )
      type(B2SourceEir), intent(inout) :: sr_eir
      
      deallocate( sr_eir % sch )
      deallocate( sr_eir % she )
      deallocate( sr_eir % shi )
      deallocate( sr_eir % smo )
      deallocate( sr_eir % sna )
      deallocate( sr_eir % sne )

      return
      end subroutine destroyB2SourceEir


      subroutine createB2SourceWork( nCv, ns, srw )
      integer, intent(in) :: nCv, ns
      type(B2SourceWork), intent(inout) :: srw

      if (allocated(srw % sch0)) return
      
      allocate( srw % sch0( nCv, 0:3 ), SOURCE=0._R8 )
      allocate( srw % she0( nCv, 0:3 ), SOURCE=0._R8 )
      allocate( srw % shi0( nCv, 0:3 ), SOURCE=0._R8 )
      allocate( srw % smo0( nCv, 0:3, 0:ns-1 ), SOURCE=0._R8 )
      allocate( srw % smq0( nCv, 0:3, 0:ns-1 ), SOURCE=0._R8 )
      allocate( srw % sna0( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( srw % sne0( nCv, 0:1 ), SOURCE=0._R8 )

      allocate( srw % b2stbm_sch( nCv ), SOURCE=0._R8 )
      allocate( srw % b2stbm_she( nCv ), SOURCE=0._R8 )
      allocate( srw % b2stbm_shi( nCv ), SOURCE=0._R8 )
      allocate( srw % b2stbm_smo( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( srw % b2stbm_smq( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( srw % b2stbm_sna( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( srw % b2stbm_sne( nCv ), SOURCE=0._R8 )

      return      
      end subroutine createB2SourceWork

      subroutine destroyB2SourceWork( srw )
      type(B2SourceWork), intent(inout) :: srw
      
      deallocate( srw % sch0 )
      deallocate( srw % she0 )
      deallocate( srw % shi0 )
      deallocate( srw % smo0 )
      deallocate( srw % smq0 )
      deallocate( srw % sna0 )
      deallocate( srw % sne0 )
      
      deallocate( srw % b2stbm_sch )
      deallocate( srw % b2stbm_she )
      deallocate( srw % b2stbm_shi )
      deallocate( srw % b2stbm_smo )
      deallocate( srw % b2stbm_smq )
      deallocate( srw % b2stbm_sna )
      deallocate( srw % b2stbm_sne )

      return
      end subroutine destroyB2SourceWork
      
      
      subroutine createB2Rates( nCv, ns, rt )
      use b2mod_rates
      
      integer, intent(in) :: nCv, ns
      type(B2Rates), intent(inout) :: rt

      if (allocated(rt % iscx)) return
      
      allocate (rt%iscx( 0:rt%nscxmax-1 ), SOURCE=0 )

      allocate( rt%rlcx( nCv, 0:1, 0:ns-1, 0:rt%nscxmax-1),SOURCE=0._R8)
      
      allocate( rt%rlqa( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlrd( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlbr( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlra( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlsa( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlza( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlz2( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlpt( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rlpi( nCv, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      
      allocate( rt%rza( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rz2( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rpt( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rt%rpi( nCv, 0:ns-1 ), SOURCE=0._R8 )

      return
      end subroutine createB2Rates
      
      
      subroutine destroyB2Rates( rt )
      type(B2Rates), intent(inout) :: rt
      
      deallocate (rt%iscx )
      deallocate( rt%rlcx )
      deallocate( rt%rlqa )
      deallocate( rt%rlrd )
      deallocate( rt%rlbr )
      deallocate( rt%rlra )
      deallocate( rt%rlsa )
      deallocate( rt%rlza )
      deallocate( rt%rlz2 )
      deallocate( rt%rlpt )
      deallocate( rt%rlpi )
      deallocate( rt%rza )
      deallocate( rt%rz2 )
      deallocate( rt%rpt )
      deallocate( rt%rpi )

      return      
      end subroutine destroyB2Rates
      
      
      subroutine createB2RatesWork( nCv, ns, rtw )
      integer, intent(in) :: nCv, ns
      type(B2RatesWork), intent(inout) :: rtw

      if (allocated(rtw % rsa)) return
      
      allocate( rtw%rsa( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rtw%rra( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rtw%rqa( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rtw%rrd( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rtw%rbr( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( rtw%rcx( nCv, 0:ns-1 ), SOURCE=0._R8 )

      return
      end subroutine createB2RatesWork
      
      
      subroutine destroyB2RatesWork( rtw )
      type(B2RatesWork), intent(inout) :: rtw
      
      if ( allocated( rtw%rsa ) ) deallocate( rtw%rsa )
      if ( allocated( rtw%rra ) ) deallocate( rtw%rra )
      if ( allocated( rtw%rqa ) ) deallocate( rtw%rqa )
      if ( allocated( rtw%rrd ) ) deallocate( rtw%rrd )
      if ( allocated( rtw%rbr ) ) deallocate( rtw%rbr )
      if ( allocated( rtw%rcx ) ) deallocate( rtw%rcx )
      return
      end subroutine destroyB2RatesWork
      
      subroutine createB2PlasmaSnapshot( nCv, nFc, ns, snap )
      integer, intent(in) :: nCv, nFc, ns
      type(B2PlasmaSnapshot), intent(inout) :: snap

      if (allocated(snap % na)) return
      
      allocate( snap % na( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( snap % ua( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( snap % po( nCv ), SOURCE=0._R8 )
      allocate( snap % te( nCv ), SOURCE=0._R8 )
      allocate( snap % ti( nCv ), SOURCE=0._R8 )
      
      allocate( snap % ne( nCv ), SOURCE=0._R8 ) 
      allocate( snap % ni( nCv, 0:1 ), SOURCE=0._R8 )    
      
      allocate( snap % fch( nFc, 0:1 ), SOURCE=0._R8 )
      allocate( snap % fna( nFc, 0:1, 0:ns-1 ), SOURCE=0._R8 )
      allocate( snap % fhe( nFc, 0:1 ), SOURCE=0._R8 )    
      allocate( snap % fhi( nFc, 0:1 ), SOURCE=0._R8 )    
      
      allocate( snap % kinrgy( nCv, 0:ns-1 ), SOURCE=0._R8 )    

      return
      end subroutine createB2PlasmaSnapshot
      
      
      subroutine destroyB2PlasmaSnapshot( snap )
      type(B2PlasmaSnapshot), intent(inout) :: snap
      
      deallocate( snap % na )   
      deallocate( snap % ua )   
      deallocate( snap % po )   
      deallocate( snap % te )   
      deallocate( snap % ti )   
      deallocate( snap % ne )
      deallocate( snap % ni )
      deallocate( snap % fch )
      deallocate( snap % fna )
      deallocate( snap % fhe )
      deallocate( snap % fhi )
      deallocate( snap % kinrgy )

      return
      end subroutine destroyB2PlasmaSnapshot
      
      subroutine createB2Diagnostic( nCv, ns, nnreg, diag )
      
      integer, intent(in) :: nCv, ns, nnreg(0:1)
      type(B2Diagnostic), intent(inout) :: diag

      if (allocated(diag % srcna)) return
      
      allocate( diag % srcna( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( diag % srcmo( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( diag % srcpo( nCv ), SOURCE=0._R8 )
      allocate( diag % srche( nCv ), SOURCE=0._R8 )
      allocate( diag % srchi( nCv ), SOURCE=0._R8 )
      allocate( diag % srcmt( nCv ), SOURCE=0._R8 )
      
      allocate( diag % aresco( 0:ns-1 ), SOURCE=0._R8 )
      allocate( diag % aresmo( 0:ns-1 ), SOURCE=0._R8 )
      
      allocate( diag % acorpa( 0:ns-1 ), SOURCE=0._R8 )
      allocate( diag % acorua( 0:ns-1 ), SOURCE=0._R8 )
      
      allocate( diag % rescoreg(0:nnreg(0), 0:ns-1 ), SOURCE=0._R8 )
      allocate( diag % resmoreg(0:nnreg(0), 0:ns-1 ), SOURCE=0._R8 )
      allocate( diag % reshereg(0:nnreg(0)), SOURCE=0._R8 )
      allocate( diag % reshireg(0:nnreg(0)), SOURCE=0._R8 )

      return
      end subroutine createB2Diagnostic
      
      
      subroutine destroyB2Diagnostic( diag )
      type(B2Diagnostic), intent(inout) :: diag
      
      deallocate( diag % srcna)   
      deallocate( diag % srcmo )   
      deallocate( diag % srcpo )   
      deallocate( diag % srche )   
      deallocate( diag % srchi )   
      deallocate( diag % srcmt )   
      deallocate( diag % aresco )
      deallocate( diag % aresmo )
      deallocate( diag % acorpa )
      deallocate( diag % acorua )
      deallocate( diag % rescoreg )
      deallocate( diag % resmoreg )
      deallocate( diag % reshereg )
      deallocate( diag % reshireg )

      return
      end subroutine destroyB2Diagnostic
      
      
      subroutine createB2Update( nCv, ns, update )
      integer, intent(in) :: nCv, ns
      type(B2Update), intent(inout) :: update

      if (allocated(update % pa)) return
      
! corrections
      allocate( update % pa( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( update % ua( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( update % na( nCv, 0:ns-1 ), SOURCE=0._R8 )
      allocate( update % po( nCv ), SOURCE=0._R8 )
      allocate( update % te( nCv ), SOURCE=0._R8 )
      allocate( update % ti( nCv ), SOURCE=0._R8 )

      return      
      end subroutine createB2Update
      
      
      subroutine destroyB2Update( update )
      type(B2Update), intent(inout) :: update
      
      deallocate( update % pa )
      deallocate( update % ua )
      deallocate( update % na )
      deallocate( update % po )
      deallocate( update % te )
      deallocate( update % ti )

      return
      end subroutine destroyB2Update
         
      

 ! External plasma
      subroutine createB2StateExt( nCv, nFc, ns_ext, state_ext )
      integer, intent(in) :: nCv, nFc, ns_ext
      type(B2StateExt), intent(inout) :: state_ext

      if (allocated(state_ext % text)) return

      allocate( state_ext % text (0:ns_ext-1) )
      allocate( state_ext % is_neutral (0:ns_ext-1) )
      allocate( state_ext % zn (0:ns_ext-1) )
      allocate( state_ext % am (0:ns_ext-1) )
      allocate( state_ext % za (1:nCv, 0:ns_ext-1) )
      allocate( state_ext % za2 (1:nCv, 0:ns_ext-1) )
      allocate( state_ext % pt (1:nCv, 0:ns_ext-1) )
      allocate( state_ext % ne (1:nCv) )
      allocate( state_ext % ne2 (1:nCv) )
      allocate( state_ext % ue (1:nCv))
      allocate( state_ext % na (1:nCv, 0:ns_ext-1) )
      allocate( state_ext % ni (1:nCv, 0:1) )
      allocate( state_ext % ua (1:nCv, 0:ns_ext-1) )
      allocate( state_ext % ta (1:nCv, 0:ns_ext-1) )
      allocate( state_ext % fhi (1:nFc, 0:1) )
      allocate( state_ext % fa (1:nFc, 0:1, 0:ns_ext-1) )
      allocate( state_ext % sne (1:nCv) )
      allocate( state_ext % she (1:nCv)  )
      allocate( state_ext % shi (1:nCv) )
      allocate( state_ext % sch (1:nCv) )
      allocate( state_ext % sna (1:nCv, 0:ns_ext-1) )
      allocate( state_ext % smo (1:nCv, 0:ns_ext-1) )
      state_ext % text       = ''
      state_ext % is_neutral = .false.
      state_ext % zn         = 0.0_R8
      state_ext % am         = 0.0_R8
      state_ext % za         = 0.0_R8
      state_ext % za2        = 0.0_R8
      state_ext % pt         = 0.0_R8
      state_ext % ne         = 0.0_R8
      state_ext % ne2        = 0.0_R8
      state_ext % ue         = 0.0_R8
      state_ext % na         = 0.0_R8
      state_ext % ni         = 0.0_R8
      state_ext % ua         = 0.0_R8
      state_ext % ta         = 0.0_R8
      state_ext % fhi        = 0.0_R8
      state_ext % fa         = 0.0_R8
      state_ext % sne        = 0.0_R8
      state_ext % she        = 0.0_R8
      state_ext % shi        = 0.0_R8
      state_ext % sch        = 0.0_R8
      state_ext % sna        = 0.0_R8
      state_ext % smo        = 0.0_R8

      return
      end subroutine createB2StateExt


      subroutine destroyB2StateExt( state_ext )
        type(B2StateExt), intent(inout) :: state_ext

        deallocate( state_ext % text )
        deallocate( state_ext % is_neutral )
        deallocate( state_ext % zn )
        deallocate( state_ext % am )
        deallocate( state_ext % za )
        deallocate( state_ext % za2 )
        deallocate( state_ext % pt )
        deallocate( state_ext % ne )
        deallocate( state_ext % ne2 )
        deallocate( state_ext % ue )
        deallocate( state_ext % na )
        deallocate( state_ext % ni )
        deallocate( state_ext % ua )
        deallocate( state_ext % ta )
        deallocate( state_ext % fhi )
        deallocate( state_ext % fa )
        deallocate( state_ext % sne )
        deallocate( state_ext % she )
        deallocate( state_ext % shi )
        deallocate( state_ext % sch )
        deallocate( state_ext % sna )
        deallocate( state_ext % smo )

      return
      end subroutine destroyB2StateExt

     
! some routines borrowed from b2nStateHelper
      subroutine getB2PlasmaSnapshot( pl, dv, snap )
        type(B2Plasma), intent(in) :: pl
        type(B2Derivatives), intent(in) :: dv
        type(B2PlasmaSnapshot), intent(inout) :: snap

        snap % na = pl % na
        snap % ua = pl % ua
        snap % po = pl % po
        snap % te = pl % te
        snap % ti = pl % ti
        snap % ne = dv % ne
        snap % ni = dv % ni
        snap % fna = dv % fna
        snap % fch = dv % fch
        snap % fhe = dv % fhe
        snap % fhi = dv % fhi
        snap % kinrgy = dv % kinrgy

      return
      end subroutine getB2PlasmaSnapshot

      subroutine putB2PlasmaSnapshot( pl, dv, snap )
        type(B2Plasma), intent(inout) :: pl
        type(B2Derivatives), intent(inout) :: dv
        type(B2PlasmaSnapshot), intent(in) :: snap

        pl % na = snap % na
        pl % ua = snap % ua
        pl % po = snap % po
        pl % te = snap % te
        pl % ti = snap % ti
        dv % ne = snap % ne
        dv % ni = snap % ni
        dv % fna = snap % fna
        dv % fch = snap % fch
        dv % fhe = snap % fhe
        dv % fhi = snap % fhi
        dv % kinrgy = snap % kinrgy
    
      return
      end subroutine putB2PlasmaSnapshot

***********************************************************************


      subroutine read_plasma_state (nput, ns, nCv, nFc, st)
      use b2mod_time
      implicit none
      integer, intent(in) :: nput, ns, nCv, nFc
      type(B2State), intent(inout) :: st
      integer :: ix, iy, ic, ifc, nfx, idum(3)
      real (kind=R8) :: time(1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('read_plasma_state')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
 
      call createB2State (nCv, nFc, ns, st)
     
* ..read state and fluxes

*   ..read plasma state
      call cfrure (nput, nCv*ns,   st % pl % na, 'na')
      call cfrure (nput, nCv,      st % dv % ne, 'ne')
      call cfrure (nput, nCv*ns,   st % pl % ua, 'ua')
      call cfrure (nput, nFc*2*ns, st % dv % uadia, 'uadia') 
      call cfrure (nput, nCv,      st % pl % te, 'te')
      call cfrure (nput, nCv,      st % pl % ti, 'ti')
      call cfrure (nput, nCv,      st % pl % po, 'po')

*   ..read fluxes
      call cfrure (nput, nFc*2*ns, st % dv % fna, 'fna')
      call cfrure (nput, nFc*2,    st % dv % fhe, 'fhe')
      call cfrure (nput, nFc*2,    st % dv % fhi, 'fhi')
      call cfrure (nput, nFc*2,    st % dv % fch, 'fch')
      call cfrure (nput, nFc*2,    st % dv % fch_32, 'fch_32')
      call cfrure (nput, nFc*2,    st % dv % fch_52, 'fch_52')
      call cfrure (nput, nCv*ns,   st % dv % kinrgy, 'kinrgy')

*   ..read time
      call cfrure (nput, 1, time, 'time')
      tim = time(1)

*   ..read parallel current
      call cfrure (nput, nFc*2,    st % dv % fch_p, 'fch_p')

*   ..read various fluxes and velocities needed by advanced SOLPS-ITER features
      call cfrure (nput, nFc*2*ns, st % dv % fna_mdf, 'fna_mdf')
      call cfrure (nput, nFc*2,    st % dv % fhe_mdf, 'fhe_mdf')
      call cfrure (nput, nFc*2,    st % dv % fhi_mdf, 'fhi_mdf')
      call cfrure (nput, nFc*2*ns, st % dv % fna_fcor, 'fna_fcor')
      call cfrure (nput, nFc*2*ns, st % dv % fna_nodrift, 'fna_nodrift')
      call cfrure (nput, nFc*2*ns, st % dv % fna_he, 'fna_he')
      call cfrure (nput, nFc*2*ns, st % dv % fnaPSch, 'fnaPSch')
      call cfrure (nput, nFc*2,    st % dv % fhePSch, 'fhePSch')
      call cfrure (nput, nFc*2,    st % dv % fhiPSch, 'fhiPSch')
      call cfrure (nput, nFc*2*ns, st % dv % fna_eir, 'fna_eir')
      call cfrure (nput, nFc*2,    st % dv % fne_eir, 'fne_eir')
      call cfrure (nput, nFc*2,    st % dv % fhe_eir, 'fhe_eir')
      call cfrure (nput, nFc*2,    st % dv % fhi_eir, 'fhi_eir')
      call cfrure (nput, nFc*2*ns, st % dv % fna_32, 'fna_32')
      call cfrure (nput, nFc*2*ns, st % dv % fna_52, 'fna_52')
      call cfrure (nput, nFc*2,    st % dv % fni_32, 'fni_32')
      call cfrure (nput, nFc*2,    st % dv % fni_52, 'fni_52')
      call cfrure (nput, nFc*2,    st % dv % fne_32, 'fne_32')
      call cfrure (nput, nFc*2,    st % dv % fne_52, 'fne_52')
      call cfrure (nput, nFc*2,    st % dv % fchdia, 'fchdia')
      call cfrure (nput, nFc*2,    st % dv % fchin, 'fchin')
      call cfrure (nput, nFc*2,    st % dv % fchvispar, 'fchvispar')
      call cfrure (nput, nFc*2,    st % dv % fchvisper, 'fchvisper')
      call cfrure (nput, nFc*2,    st % dv % fchvisq, 'fchvisq')
      call cfrure (nput, nFc*2,    st % dv % fchinert, 'fchinert')
      call cfrure (nput, nFc*2*ns, st % dv % vaecrb, 'vaecrb')
      call cfrure (nput, nFc*2*ns, st % dv % vadia, 'vadia')
      call cfrure (nput, nFc*2*ns, st % dv % wadia, 'wadia')
      call cfrure (nput, nFc*2,    st % dv % veecrb, 'veecrb')
      call cfrure (nput, nFc*2,    st % dv % vedia, 'vedia')
      call cfrure (nput, nFc*2,    st % dv % floe_noc, 'floe_noc')
      call cfrure (nput, nFc*2,    st % dv % floi_noc, 'floi_noc')

* ..return
      call subend ()
      return

      end subroutine read_plasma_state


***********************************************************************


      subroutine write_plasma_state (nput, ns, nCv, nFc, st)
      use b2mod_time
      implicit none
      integer, intent(in) :: nput, ns, nCv, nFc
      type(B2State), intent(inout) :: st
      integer ::  idum(0:2) 
      real (kind=R8) :: time(1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('write_plasma_state')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
      
* ..save state and fluxes

*   ..save plasma state
      call cfwure (nput, nCv*ns,   st % pl % na, 'na')
      call cfwure (nput, nCv,      st % dv % ne, 'ne')
      call cfwure (nput, nCv*ns,   st % pl % ua, 'ua')
      call cfwure (nput, nFc*2*ns, st % dv % uadia, 'uadia') 
      call cfwure (nput, nCv,      st % pl % te, 'te')
      call cfwure (nput, nCv,      st % pl % ti, 'ti')
      call cfwure (nput, nCv,      st % pl % po, 'po')

*   ..save fluxes
      call cfwure (nput, nFc*2*ns, st % dv % fna, 'fna')
      call cfwure (nput, nFc*2,    st % dv % fhe, 'fhe')
      call cfwure (nput, nFc*2,    st % dv % fhi, 'fhi')
      call cfwure (nput, nFc*2,    st % dv % fch, 'fch')
      call cfwure (nput, nFc*2,    st % dv % fch_32, 'fch_32')
      call cfwure (nput, nFc*2,    st % dv % fch_52, 'fch_52')
      call cfwure (nput, nCv*ns,   st % dv % kinrgy, 'kinrgy')

*   ..save time
      time(1)=tim
      call cfwure (nput, 1, time, 'time')

*   ..save parallel current
      call cfwure (nput, nFc*2,    st % dv % fch_p, 'fch_p')

*   ..saving various fluxes and velocities needed by advanced SOLPS-ITER features
      call cfwure (nput, nFc*2*ns, st % dv % fna_mdf, 'fna_mdf')
      call cfwure (nput, nFc*2,    st % dv % fhe_mdf, 'fhe_mdf')
      call cfwure (nput, nFc*2,    st % dv % fhi_mdf, 'fhi_mdf')
      call cfwure (nput, nFc*2*ns, st % dv % fna_fcor, 'fna_fcor')
      call cfwure (nput, nFc*2*ns, st % dv % fna_nodrift, 'fna_nodrift')
      call cfwure (nput, nFc*2*ns, st % dv % fna_he, 'fna_he')
      call cfwure (nput, nFc*2*ns, st % dv % fnaPSch, 'fnaPSch')
      call cfwure (nput, nFc*2,    st % dv % fhePSch, 'fhePSch')
      call cfwure (nput, nFc*2,    st % dv % fhiPSch, 'fhiPSch')
      call cfwure (nput, nFc*2*ns, st % dv % fna_eir, 'fna_eir')
      call cfwure (nput, nFc*2,    st % dv % fne_eir, 'fne_eir')
      call cfwure (nput, nFc*2,    st % dv % fhe_eir, 'fhe_eir')
      call cfwure (nput, nFc*2,    st % dv % fhi_eir, 'fhi_eir')
      call cfwure (nput, nFc*2*ns, st % dv % fna_32, 'fna_32')
      call cfwure (nput, nFc*2*ns, st % dv % fna_52, 'fna_52')
      call cfwure (nput, nFc*2,    st % dv % fni_32, 'fni_32')
      call cfwure (nput, nFc*2,    st % dv % fni_52, 'fni_52')
      call cfwure (nput, nFc*2,    st % dv % fne_32, 'fne_32')
      call cfwure (nput, nFc*2,    st % dv % fne_52, 'fne_52')
      call cfwure (nput, nFc*2,    st % dv % fchdia, 'fchdia')
      call cfwure (nput, nFc*2,    st % dv % fchin, 'fchin')
      call cfwure (nput, nFc*2,    st % dv % fchvispar, 'fchvispar')
      call cfwure (nput, nFc*2,    st % dv % fchvisper, 'fchvisper')
      call cfwure (nput, nFc*2,    st % dv % fchvisq, 'fchvisq')
      call cfwure (nput, nFc*2,    st % dv % fchinert, 'fchinert')
      call cfwure (nput, nFc*2*ns, st % dv % vaecrb, 'vaecrb')
      call cfwure (nput, nFc*2*ns, st % dv % vadia, 'vadia')
      call cfwure (nput, nFc*2*ns, st % dv % wadia, 'wadia')
      call cfwure (nput, nFc*2,    st % dv % veecrb, 'veecrb')
      call cfwure (nput, nFc*2,    st % dv % vedia, 'vedia')
      call cfwure (nput, nFc*2,    st % dv % floe_noc, 'floe_noc')
      call cfwure (nput, nFc*2,    st % dv % floi_noc, 'floi_noc')

* ..return
      call subend ()
      return

      end subroutine write_plasma_state

      end module b2us_plasma
      
