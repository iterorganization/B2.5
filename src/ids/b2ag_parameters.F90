module b2ag_parameters

  use b2mod_types
  use tradui_constants
  use b2mod_version
  use b2ag_ghostcells

  implicit none
  
  ! parg - (0:99) real*8 array.
  ! parg(0:*) contains a set of numerical parameters that specify
  ! details of the geometry. The interpretation of these parameters
  ! is problem-dependent; see the code.
  real (kind=R8) :: parg(0:99) 

  ! isolating cuts
  integer :: niso, nxiso(1:NISOMX)

  logical :: doGhostCells

contains

  !> Read and initialize parameters for b2ag
  !> nx, ny: dimensions for the B2 data structures (which are then
  !>         of size (-1:nx, -1:ny)
  !> nx1, ny1: size of input grid which is to be downscaled to nx, ny  
  subroutine b2ag_read_parameters_st(ninp, nout, nx, ny, nx1, ny1)
    integer, intent(in) :: ninp(0:1), nout(0:1)
    integer, intent(out) :: nx, ny, nx1, ny1

    ! internal
    character :: id*8, cnamip*80, cvalip*80
    integer :: nnx, nny, lun=96
    character*256 local_sonnet
    integer :: istyle

    external xertst, ipsetc, xerrab, streql, b2agx0_st, ipgetc, cfverr, &
         & open_file
    logical :: streql, open_file

    call b2agx0_st (ninp(0), nx, ny, nx1, ny1)

    !   ..read geometry parameters parg
    read (ninp(0),'(a8)',err=91) id
    call xertst (streql(id,'*param'),&
         & 'faulty input instead of *param, id='//id)
    read (ninp(0),*,err=92) parg
    ! extract istyle
    istyle = nint(parg(0))

    !   ..echo parg
    write (nout(0),'(/2x,a)') 'The input parameters are:'
    write (nout(0),*) parg
    !   ..read and echo code internal parameters
    write (nout(0),'(/2x,a)') 'non-default internal parameters:'
1   continue
    read (ninp(0),*,end=2,err=93) cnamip, cvalip
    if (cnamip(1:1).ne.'*') then
      write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
      call ipsetc (cnamip, cvalip)
    endif
    goto 1
2   continue
    write (nout(0),'(2x,a)') '(end of list of internal parameters)'

    ! open geometry file and get nx, ny,... from there
    nnx = 0
    nny = 0
    doGhostCells = .false.

    call ipgetc ('b2agfs_geometry', local_sonnet)
    if(open_file(lun,local_sonnet, 'data.local/meshes|data/meshes')) then
        call cfverr(lun,grid_version)
        if(istyle.eq.-1) then 
            ! Carre grid input
            if (grid_version.lt."01.001.028") then
                read(lun,*) nnx,nny
                niso = 0
            else
                read(lun,*) nnx,nny,niso,nxiso(1:nisomx)
                doGhostCells = .true.
            end if
        end if
        close(lun)
    end if

    ! check that b2ag and geometry file sizes are compatible
    if(nnx+2*niso.ne.nx1.or.nny.ne.ny1) then        
        write(*,'(a,4i4)') 'b2ag_read_parameters_st: nx,ny code and data disagree',nx1,ny1,nnx+2*niso,nny
        stop
    endif

    
    ! fix nx1, ny1 in case ghost cells are not provided
    if (doGhostCells) then
        call computeGridSizeWithGhostCells_st(nnx, nny, niso, nx1, ny1)
    end if

    return
    !-----------------------------------------------------------------------
    !.error conditions

91  call xerrab ('error trying to read heading *param')
92  call xerrab ('error trying to read parg')
93  call xerrab ('error trying to read non-default internal parameters')

  end subroutine b2ag_read_parameters_st

  !> Read and initialize parameters for b2ag
  !> nCv, nFc, nVx: dimensions for the B2 data structures  
  subroutine b2ag_read_parameters(ninp,nout,m)
    use b2us_map
    integer, intent(in) :: ninp(0:1), nout(0:1)
    type(mapping) :: m 
    !internal
    character :: id*8, cnamip*80, cvalip*80
    integer :: nCv0, nFc0, nVx0, nCi, nCg, lun=96, idum(0:3), idum2(0:4), &
        & nCmxFc0, nCmxVx0, nFmxCv0, nVmxCv0, nVmxFc0, &
        & nCmxVx, nCmxFc, nFmxCv, nVmxCv, nVmxFc, nCv, nFc, nVx, nncut
    character*256 local_sonnet
    integer :: istyle
    logical :: streql, open_file

    external xertst, ipsetc, xerrab, streql, b2agx0, ipgetc, cfverr, &
         & open_file, cfruin

     call b2agx0 (ninp(0), nCv, nFc, nVx) ! nog issue met deze functie


     ! istyle == -1
     ! means Carre grid input
     !   ..read geometry parameters parg
     read (ninp(0),'(a8)',err=93) id
     call xertst (streql(id,'*param'),&
    & 'faulty input instead of *param, id='//id)
     read (ninp(0),*,err=94) parg
     ! extract istyle (what does it mean?, still appl. for us)
     istyle = nint(parg(0))

     !   ..echo parg
     write (nout(0),'(/2x,a)') 'The input parameters are:'
     write (nout(0),*) parg
     !   ..read and echo code internal parameters
     write (nout(0),'(/2x,a)') 'non-default internal parameters:'
1    continue
     read (ninp(0),*,end=2,err=95) cnamip, cvalip
     if (cnamip(1:1).ne.'*') then
       write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
       call ipsetc (cnamip, cvalip)
     endif
     goto 1
2    continue
     write (nout(0),'(2x,a)') '(end of list of internal parameters)'

     ! open geometry file and get nx, ny,... from there
     nCv0 = 0
     nFc0 = 0
     nVx0 = 0
     doGhostCells = .false.

     call ipgetc ('b2agfs_geometry', local_sonnet)
     if(open_file(lun,local_sonnet, 'data.local/meshes|data/meshes')) then
         call cfverr(lun,grid_version)
         if(istyle.eq.-1) then 
            ! Carre grid input
            if (grid_version.eq."03.002.000") then
		! obtain nCv, nFc, nVx from geometry file
      	        call cfruin (lun, 4, idum, 'nCv,nFc,nVx,nCg')
                nCv0 = idum(0)
                nFc0 = idum(1)
	            nVx0 = idum(2)
                m%nCg = idum(3)
                doGhostCells = .true.
        ! obtain nCmxFc and nCmxVx from geometry file
                call cfruin (lun, 5, idum2, 'nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc')
                nCmxVx0 = idum2(0)
                nCmxFc0 = idum2(1)
                nFmxCv0 = idum2(2)
                nVmxCv0 = idum2(3)
                nVmxFc0 = idum2(4)
                call cfruin (lun, 2, idum, 'nx,ny,nncut') ! here still add nncut?? 
                m%nx = idum(0)
                m%ny = idum(1)
                nncut = idum(2)
            !else
                !read(lun,*) nnx,nny,niso,nxiso(1:nisomx)
                !doGhostCells = .true.
            end if
        end if
        
     end if

     ! check that b2ag and geometry file sizes are compatible
     if(nCv.ne.nCv0.or.nFc.ne.nFc0.or.nVx.ne.nVx0) then        
        write(*,'(a,4i4)') 'b2ag_read_parameters: nCv,nFc, nVx code and data disagree',nCv,nFc,nVx,nCv0,nFc0,nVx0
        stop
     endif
 
     ! fix nCv, in case ghost cells are not provided
     if (doGhostCells) then
        rewind lun
        nCi = nCv
        call computeGridSizeWithGhostCells(lun,nCv,m%nCg,nCmxVx0,nCmxFc0, &
              &  nFmxCv0, nVmxCv0, nVmxFc0,nCv,nCmxVx,nCmxFc, &
              &  nFmxCv, nVmxCv, nVmxFc)
        m%nCv = nCv 
        m%nFc = nFc0
        m%nVx = nVx0
        m%nCi = m%nCv - m%nCg
        m%nCmxVx = nCmxVx
        m%nCmxFc = nCmxFc
        !m%nFmxCv = nFmxCv
        m%nVmxCv = nVmxCv
        m%nVmxFc = nVmxFc

	 ! subroutine of b2ag_ghostcells
     end if
     close(lun)

      !introduce the flux tubes etc.
      !rough over-estimation
      m%nFs = max(1,(1 + 2*nncut)*m%ny)
      m%nFt = m%nFs
      !m%nFs = 0
      !m%nFt = 0



     return

    !-----------------------------------------------------------------------
    !.error conditions

91 call xerrab ('error trying to read header *dimens')
92 call xerrab ('error trying to read after '//id)
93  call xerrab ('error trying to read heading *param')
94  call xerrab ('error trying to read parg')
95  call xerrab ('error trying to read non-default internal parameters')
    

  end subroutine b2ag_read_parameters

    !-----------------------------------------------------------------------



end module b2ag_parameters
