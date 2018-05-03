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
        real(R8), allocatable :: fch(:,:,:), fch_32(:,:,:), 
     .                           fch_52(:,:,:), fch_p(:,:,:)
        real(R8), allocatable :: fni_32(:,:,:), fni_52(:,:,:), 
     .                           fni(:,:,:)
        real(R8), allocatable :: fna(:,:,:,:), fna_52(:,:,:,:), 
     .                           fna_32(:,:,:,:), fna_nodrift(:,:,:,:) 
!    .                         , fna_52nd(:,:,:), fna_32nd(:,:,:)
        real(R8), allocatable :: fne(:,:,:), fne_32(:,:,:), 
     .                           fne_52(:,:,:), 
     .                           fhe(:,:,:), fhet(:,:,:)
        real(R8), allocatable :: fhi(:,:,:), fhit(:,:,:)          


        ! kinetic energy
        real(R8), allocatable :: kinrgy(:,:)

        ! numerical coefficients (maybe move to separate object)
        real(R8), dimension(:,:,:), allocatable :: conc
        real(R8), dimension(:,:,:), allocatable :: floe, floi
        real(R8), dimension(:,:,:,:), allocatable :: cone, coni
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
        real(R8), dimension(:,:,:,:), allocatable :: pcca
        real(R8), dimension(:,:,:), allocatable :: pccm

        ! derived state quantities
        ! ni(:,0) is ions, ni(:,1) is neutrals, ue is parallel electron velocity
        ! ne2(:) is electron density, ne2 is elect. den. computed w. squared charge rz2
        real(R8), allocatable :: ne(:), ni(:,:), ue(:), ne2(:)

      end type B2Derivatives

  ! Transport coefficients (velocities, diffusivities...)
      type B2Coeff
        ! at faces
        real(R8), dimension(:,:,:), allocatable :: csig, calf, csig_an
        real(R8), dimension(:,:,:,:,:), allocatable :: csigin
        real(R8), dimension(:,:,:), allocatable :: chve, chce
        real(R8), dimension(:,:,:), allocatable :: chvi, chci
        real(R8), dimension(:,:), allocatable :: chvemx, chvimx
        real(R8), allocatable :: cvla(:,:,:,:), cdna(:,:,:,:), 
     .                           cdpa(:,:,:,:), cvsa(:,:,:,:)
        real(R8), allocatable :: cthe(:,:), cthi(:,:), cvsa_cl(:,:,:,:)
        real(R8), allocatable :: ceqp(:)
          
        real(R8), allocatable :: fllim0fhi(:,:,:,:), fllimvisc(:,:) 

        real(R8), allocatable, dimension(:) :: sig0, hce0, hci0, alf0
        real(R8), allocatable, dimension(:,:) :: dpa0, dna0, vsa0, 
     .                                             hcib
        real(R8), allocatable, dimension(:,:,:) :: vla0

        ! velocities at cell centers
        ! dimensions: cv, direction, species
        real(R8), allocatable :: uadia(:,:,:,:), wadia(:,:,:,:), 
     .                           vaecrb(:,:,:,:)

      end type B2Coeff
  
      type B2Source 
        real(R8), dimension(:,:), allocatable :: sch, she, shi
        real(R8), dimension(:,:,:), allocatable :: smo, smq, sna
      end type B2Source

      type B2SourceWork
        real(R8), dimension(:,:), allocatable :: sch0, she0, shi0
        real(R8), dimension(:,:,:), allocatable :: smo0, smq0, sna0         
      end type B2SourceWork

  ! b2mod_rates
      type B2Rates
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
        real(R8), allocatable :: fch(:,:,:), fna(:,:,:,:), 
     .                           fhi(:,:,:), fhe(:,:,:)
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
        integer :: ns ! number of species
        integer :: ismain       ! index of main species

        type(B2Plasma) :: pl 
        type(B2Coeff) :: co
        type(B2Derivatives) :: dv
        type(B2Source) :: sr
        type(B2SourceWork) :: srw
        type(B2Rates) :: rt
        type(B2PlasmaSnapshot) :: psnl, psnc
        type(B2Diagnostic) :: diag
        type(B2Update) :: update
      end type B2State

      contains

      subroutine createB2Plasma( nCv, ns, plasma )
      integer, intent(in) :: nCv, ns
      type(B2Plasma), intent(inout) :: plasma

      plasma % ns = ns

      allocate( plasma % na( nCv, 0 : ns - 1 ) )
      allocate( plasma % ua( nCv, 0 : ns - 1 ) )
      allocate( plasma % po( nCv ) )
      allocate( plasma % te( nCv ) )
      allocate( plasma % ti( nCv ) )
      
      end subroutine createB2Plasma


      subroutine destroyB2Plasma( plasma )
      type(B2Plasma), intent(inout) :: plasma
      
      plasma % ns = 0
      deallocate( plasma % na )   
      deallocate( plasma % ua )   
      deallocate( plasma % po )   
      deallocate( plasma % te )   
      deallocate( plasma % ti )   
      
      end subroutine destroyB2Plasma
      
      
      subroutine createB2Coeff( nCv, nFC, ns, coeff )
      integer, intent(in) :: nCv, nFc, ns
      type(B2Coeff), intent(inout) :: coeff
      
      allocate( coeff % csig_an( nFc, 0:1, 0:1 ) )
      allocate( coeff % calf( nFc, 0:1, 0:1 ) )
      
      allocate( coeff % csig( nFc, 0:1, 0:1 ) )
      allocate( coeff % csigin( nFc, 0:1, 0:1, 0:ns-1, 0:ns-1 ) )
      
      allocate( coeff % chve( nFc, 0:1, 0:1 ) )
      allocate( coeff % chce( nFc, 0:1, 0:1 ) )
      allocate( coeff % chvemx( nFc, 0:1 ) )
      
      allocate( coeff % chvi( nFc, 0:1, 0:1 ) )
      allocate( coeff % chci( nFc, 0:1, 0:1 ) )
      allocate( coeff % chvimx( nFc, 0:1 ) )
      
      allocate( coeff % cvla( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( coeff % cdna( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( coeff % cdpa( nFc, 0:1, 0:1, 0:ns-1 ) )
      
      allocate( coeff % cvsa( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( coeff % cthe( nFc, 0:ns-1 ) )
      allocate( coeff % cthi( nFc, 0:ns-1 ) )
      allocate( coeff % cvsa_cl( nFc, 0:1, 0:1, 0:ns-1 ) )
      
      allocate( coeff % ceqp( nCv ) )
      
      allocate( coeff % fllim0fhi( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( coeff % fllimvisc( nFc, 0:ns-1 ) )
      
      allocate( coeff % uadia( nCv, 0:1, 0:1, 0:ns-1 ) )
      allocate( coeff % wadia( nCv, 0:1, 0:1, 0:ns-1 ) )
      allocate( coeff % vaecrb( nCv, 0:1, 0:1, 0:ns-1 ) )
      
! moved here from b2coeffwork
      
      allocate( coeff%dna0( nCv, 0:ns-1 ) )
      allocate( coeff%dpa0( nCv, 0:ns-1 ) )
      allocate( coeff%vla0( nCv, 0:1, 0:ns-1 ) )
      allocate( coeff%vsa0( nCv, 0:ns-1 ) )
      allocate( coeff%hci0( nCv ) )
      allocate( coeff%hcib( nCv, 0:ns-1 ) )
      allocate( coeff%hce0( nCv ) )
      allocate( coeff%sig0( nCv ) )
      allocate( coeff%alf0( nCv ) )

      end subroutine createB2Coeff

      
      subroutine destroyB2Coeff( coeff )
      type(B2Coeff), intent(inout) :: coeff
      
      deallocate( coeff % csig )
      deallocate( coeff % csigin )
      deallocate( coeff % csig_an )
      deallocate( coeff % calf )
      
      deallocate( coeff % chve )
      deallocate( coeff % chce )
      deallocate( coeff % chvemx )
      
      deallocate( coeff % chvi )
      deallocate( coeff % chci )
      deallocate( coeff % chvimx )
      
      deallocate( coeff % cvla )
      deallocate( coeff % cdna )
      deallocate( coeff % cdpa )
      
      deallocate( coeff % cvsa )
      deallocate( coeff % cthe )
      deallocate( coeff % cthi )
      deallocate( coeff % cvsa_cl )
      
      deallocate( coeff % ceqp )
      
      deallocate( coeff % fllim0fhi )
      deallocate( coeff % fllimvisc )
      
      deallocate( coeff % uadia )
      deallocate( coeff % wadia )
      deallocate( coeff % vaecrb )
      
      deallocate( coeff%dna0 )
      deallocate( coeff%dpa0 )
      deallocate( coeff%vla0 )
      deallocate( coeff%vsa0 )
      deallocate( coeff%hci0 )
      deallocate( coeff%hcib )
      deallocate( coeff%hce0 )
      deallocate( coeff%sig0 )
      deallocate( coeff%alf0 )
      
      end subroutine destroyB2Coeff
      
      
      
      subroutine createB2Derivatives( nCv, nFc, ns, dv )
      integer, intent(in) :: nCv, nFc, ns
      type(B2Derivatives), intent(inout) :: dv
      
! fluxes
      allocate( dv % fch( nFc, 0:1, 0:1 ) )
      allocate( dv % fch_32( nFc, 0:1, 0:1 ) )
      allocate( dv % fch_52( nFc, 0:1, 0:1 ) )
      allocate( dv % fch_p( nFc, 0:1, 0:1 ) )
      
      allocate( dv % fne( nFc, 0:1, 0:1 ) )
      allocate( dv % fne_32( nFc, 0:1, 0:1 ) )
      allocate( dv % fne_52( nFc, 0:1, 0:1 ) )
      
      allocate( dv % fni( nFc, 0:1, 0:1 ) )    
      allocate( dv % fni_32( nFc, 0:1, 0:1 ) )
      allocate( dv % fni_52( nFc, 0:1, 0:1 ) )
      
      allocate( dv % fna( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( dv % fna_32( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( dv % fna_52( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( dv % fna_nodrift( nFc, 0:1, 0:1, 0:ns-1 ) )
      
      allocate( dv % fhe( nFc, 0:1, 0:1 ) )    
      allocate( dv % fhet( nFc, 0:1, 0:1 ) )
      allocate( dv % fhi( nFc, 0:1, 0:1 ) )    
      allocate( dv % fhit( nFc, 0:1, 0:1 ) )
      
! kinetic energy
      allocate( dv % kinrgy( nCv, 0:ns-1 ) )    
      
! derived state quantities
      allocate( dv % ne( nCv ) )    
      allocate( dv % ne2( nCv ) )    
      allocate( dv % ni( nCv, 0:1 ) )    
      allocate( dv % ue( nCv ) )
      
! numerical coefficients
      allocate( dv % conc( nFc, 0:1, 0:1 ) )
      allocate( dv % fllime( nFc ) )    
      allocate( dv % floe( nFc, 0:1, 0:1 ) )
      allocate( dv % cone( nFc, 0:1, 0:1, 0:2) )
      allocate( dv % fllimi( nFc ) )    
      allocate( dv % floi( nFc, 0:1, 0:1 ) )
      allocate( dv % coni( nFc, 0:1, 0:1, 0:2 ) )
      
! residuals
      allocate( dv % resmo( nCv, 0:ns-1 ) )
      allocate( dv % resco( nCv, 0:ns-1 ) )
      allocate( dv % respo( nCv ) )
      allocate( dv % reshe( nCv ) )
      allocate( dv % reshi( nCv ) )
      allocate( dv % resht( nCv ) )
      allocate( dv % resmt( nCv ) )
      
! corrections
      allocate( dv % corua( nCv, 0:ns-1 ) )
      allocate( dv % corpa( nCv, 0:ns-1 ) )
      allocate( dv % corut( nCv ) )
      allocate( dv % corpo( nCv ) )
      allocate( dv % cortt( nCv ) )
      allocate( dv % corte( nCv ) )
      allocate( dv % corti( nCv ) )
      
! coupling coefficients
      allocate( dv % pcca( nCv, 0:1, 0:1, 0:ns-1 ) )
      allocate( dv % pccm( nCv, 0:1, 0:1 ) )
      
      end subroutine createB2Derivatives
      
      
      subroutine destroyB2Derivatives( dv )
      type(B2Derivatives), intent(inout) :: dv
      
      deallocate( dv % fne )
      deallocate( dv % fni )
      deallocate( dv % fne_32 )
      deallocate( dv % fne_52 )
      deallocate( dv % fna_32 )
      deallocate( dv % fna_52 )
      deallocate( dv % fna_nodrift )
      deallocate( dv % fni_32 )
      deallocate( dv % fni_52 )
      deallocate( dv % fna )
      deallocate( dv % fch )
      deallocate( dv % fch_32 )
      deallocate( dv % fch_52 )
      deallocate( dv % fch_p )
      deallocate( dv % fna )
      deallocate( dv % fna_32 )
      deallocate( dv % fna_52 )
      deallocate( dv % fna_nodrift )
      
      deallocate( dv % ne )
      deallocate( dv % ne2 )
      deallocate( dv % ni )
      deallocate( dv % ue )
      
! kinetic energy
      deallocate( dv % kinrgy )
      
      deallocate( dv % fhi )
      deallocate( dv % fhit )
      deallocate( dv % fhe )
      deallocate( dv % fhet )

      deallocate( dv % conc )
      deallocate( dv % fllime )
      deallocate( dv % floe )
      deallocate( dv % cone )
      deallocate( dv % fllimi )
      deallocate( dv % floi )
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
           
      end subroutine destroyB2Derivatives
      
      
      subroutine createB2Source( nCv, ns, sr )
      integer, intent(in) :: nCv, ns
      type(B2Source), intent(inout) :: sr
      
      allocate( sr % sch( nCv, 0:3 ) )
      allocate( sr % she( nCv, 0:3 ) )
      allocate( sr % shi( nCv, 0:3 ) )
      allocate( sr % smo( nCv, 0:3, 0:ns-1 ) )
      allocate( sr % smq( nCv, 0:3, 0:ns-1 ) )
      allocate( sr % sna( nCv, 0:1, 0:ns-1 ) )
      
      end subroutine createB2Source
      
      subroutine destroyB2Source( sr )
      type(B2Source), intent(inout) :: sr
      
      deallocate( sr % sch )
      deallocate( sr % she )
      deallocate( sr % shi )
      deallocate( sr % smo )
      deallocate( sr % smq )
      deallocate( sr % sna )
      end subroutine destroyB2Source
      
      subroutine createB2SourceWork( nCv, ns, srw )
      integer, intent(in) :: nCv, ns
      type(B2SourceWork), intent(inout) :: srw
      
      allocate( srw % sch0( nCv, 0:3 ) )
      allocate( srw % she0( nCv, 0:3 ) )
      allocate( srw % shi0( nCv, 0:3 ) )
      allocate( srw % smo0( nCv, 0:3, 0:ns-1 ) )
      allocate( srw % smq0( nCv, 0:3, 0:ns-1 ) )
      allocate( srw % sna0( nCv, 0:1, 0:ns-1 ) )
      
      end subroutine createB2SourceWork

      subroutine destroyB2SourceWork( srw )
      type(B2SourceWork), intent(inout) :: srw
      
      deallocate( srw % sch0 )
      deallocate( srw % she0 )
      deallocate( srw % shi0 )
      deallocate( srw % smo0 )
      deallocate( srw % smq0 )
      deallocate( srw % sna0 )
      end subroutine destroyB2SourceWork
      
      
      subroutine createB2Rates( nCv, ns, rt )
      use b2mod_rates
      
      integer, intent(in) :: nCv, ns
      type(B2Rates), intent(inout) :: rt
      
      allocate( rt%rlcx( nCv, 0:1, 0:ns-1, 0:nscxmax-1 ) ) ! nscxmax is from b2mod_rates
      
      allocate( rt%rlqa( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlrd( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlbr( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlra( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlsa( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlza( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlz2( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlpt( nCv, 0:1, 0:ns-1 ) )
      allocate( rt%rlpi( nCv, 0:1, 0:ns-1 ) )
      
      allocate( rt%rza( nCv, 0:ns-1 ) )
      allocate( rt%rz2( nCv, 0:ns-1 ) )
      allocate( rt%rpt( nCv, 0:ns-1 ) )
      allocate( rt%rpi( nCv, 0:ns-1 ) )
      end subroutine createB2Rates
      
      
      subroutine destroyB2Rates( rt )
      type(B2Rates), intent(inout) :: rt
      
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
      
      end subroutine destroyB2Rates
      
      
      subroutine createB2RatesWork( nCv, ns, rtw )
      integer, intent(in) :: nCv, ns
      type(B2RatesWork), intent(inout) :: rtw
      
      allocate( rtw%rsa( nCv, 0:ns-1 ) )
      allocate( rtw%rra( nCv, 0:ns-1 ) )
      allocate( rtw%rqa( nCv, 0:ns-1 ) )
      allocate( rtw%rrd( nCv, 0:ns-1 ) )
      allocate( rtw%rbr( nCv, 0:ns-1 ) )
      allocate( rtw%rcx( nCv, 0:ns-1 ) )
      end subroutine createB2RatesWork
      
      
      subroutine destroyB2RatesWork( rtw )
      type(B2RatesWork), intent(inout) :: rtw
      
      if ( allocated( rtw%rsa ) ) deallocate( rtw%rsa )
      if ( allocated( rtw%rra ) ) deallocate( rtw%rra )
      if ( allocated( rtw%rqa ) ) deallocate( rtw%rqa )
      if ( allocated( rtw%rrd ) ) deallocate( rtw%rrd )
      if ( allocated( rtw%rbr ) ) deallocate( rtw%rbr )
      if ( allocated( rtw%rcx ) ) deallocate( rtw%rcx )
      end subroutine destroyB2RatesWork
      
      subroutine createB2PlasmaSnapshot( nCv, nFc, ns, snap )
      integer, intent(in) :: nCv, nFc, ns
      type(B2PlasmaSnapshot), intent(inout) :: snap
      
      allocate( snap % na( nCv, 0:ns-1 ) )
      allocate( snap % ua( nCv, 0:ns-1 ) )
      allocate( snap % po( nCv ) )
      allocate( snap % te( nCv ) )
      allocate( snap % ti( nCv ) )
      
      allocate( snap % ne( nCv ) ) 
      allocate( snap % ni( nCv, 0:1 ) )    
      
      allocate( snap % fch( nFc, 0:1, 0:1 ) )
      allocate( snap % fna( nFc, 0:1, 0:1, 0:ns-1 ) )
      allocate( snap % fhe( nFc, 0:1, 0:1 ) )    
      allocate( snap % fhi( nFc, 0:1, 0:1 ) )    
      
      allocate( snap % kinrgy( nCv, 0:ns-1 ) )    
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
      end subroutine destroyB2PlasmaSnapshot
      
      subroutine createB2Diagnostic( nCv, ns, diag )
      use b2mod_indirect , only: nnreg
      
      integer, intent(in) :: nCv, ns
      type(B2Diagnostic), intent(inout) :: diag
      
      allocate( diag % srcna( nCv, 0:ns-1 ) )
      allocate( diag % srcmo( nCv, 0:ns-1 ) )
      allocate( diag % srcpo( nCv ) )
      allocate( diag % srche( nCv ) )
      allocate( diag % srchi( nCv ) )
      allocate( diag % srcmt( nCv ) )
      
      allocate( diag % aresco( 0:ns-1 ) )
      allocate( diag % aresmo( 0:ns-1 ) )
      
      allocate( diag % acorpa( 0:ns-1 ) )
      allocate( diag % acorua( 0:ns-1 ) )
      
      allocate( diag % rescoreg(0:nnreg(0), 0:ns-1 ) )
      allocate( diag % resmoreg(0:nnreg(0), 0:ns-1 ) )
      allocate( diag % reshereg(0:nnreg(0)) )
      allocate( diag % reshireg(0:nnreg(0)) )
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
      end subroutine destroyB2Diagnostic
      
      
      subroutine createB2Update( nCv, ns, update )
      integer, intent(in) :: nCv, ns
      type(B2Update), intent(inout) :: update
      
! corrections
      allocate( update % pa( nCv, 0:ns-1 ) )
      allocate( update % ua( nCv, 0:ns-1 ) )
      allocate( update % na( nCv, 0:ns-1 ) )
      allocate( update % po( nCv ) )
      allocate( update % te( nCv ) )
      allocate( update % ti( nCv ) )
      
      end subroutine createB2Update
      
      
      subroutine destroyB2Update( update )
      type(B2Update), intent(inout) :: update
      
      deallocate( update % pa )
      deallocate( update % ua )
      deallocate( update % na )
      deallocate( update % po )
      deallocate( update % te )
      deallocate( update % ti )
      end subroutine destroyB2Update
         
           
      end module b2us_plasma
      
