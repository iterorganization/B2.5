      module b2mod_mdsplus
      use b2mod_types
      use b2mod_constants
*   ..database quantities
      integer :: dbsolpsindex=0, dbshot=0, dbnCv=0, dbnFc=0, dbns=0,
     &     dbversion=0, dbsvn=0
      real (SELECTED_REAL_KIND (14)) :: dbtime=0.0_R8, 
     &     dbimpsep_te=0.0_R8, dbimpsep_ti=0.0_R8, dbimpsep_ne=0.0_R8,
     &     dbimpsep_d=0.0_R8, dbimpsep_kye=0.0_R8, dbimpsep_kyi=0.0_R8,
     &     dbimpsep_kyi0=0.0_R8, dbimpsep_cr=0.0_R8, dbimpsep_cz=0.0_R8,
     &     dbompsep_te=0.0_R8, dbompsep_ti=0.0_R8, dbompsep_ne=0.0_R8,
     &     dbompsep_d=0.0_R8, dbompsep_kye=0.0_R8, dbompsep_kyi=0.0_R8,
     &     dbompsep_kyi0=0.0_R8, dbompsep_cr=0.0_R8, dbompsep_cz=0.0_R8,
     &     dbsep_fhe=0.0_R8, dbsep_fhi=0.0_R8, 
     &     dbsep_fne=0.0_R8, dbsep_fni=0.0_R8
      character*(32) :: dbexp='NOT_SET', dbuser='NOT_SET',
     &     dbhost='NOT_SET', dbdate='NOT_SET', dbobjectcode='NOT_SET',
     &     dbb2fstatiid='NOT_SET', dbb2fstateid='NOT_SET',
     &     dbb2fgmtryid='NOT_SET', dbb2fratesid='NOT_SET',
     &     dbb2fplasmaid='NOT_SET',
     &     dbsolpsversion='NOT_SET',
     &     db_version_b2fstati='NOT_SET',
     &     db_version_b2fstate='NOT_SET',
     &     db_version_b2frates='NOT_SET',
     &     db_version_b2fmovie='NOT_SET',
     &     db_version_b2ftrack='NOT_SET',
     &     db_version_plasmastate='NOT_SET',
     &     db_version_b2fplasma='NOT_SET',
     &     dbgit='None'
      character*(511) :: dbspecies='NOT_SET'
      character*(255) :: dbdirectory='NOT_SET',
     &     dbcomments='NOT_SET', dbgeometry='NOT_SET', 
     &     dbphysics='NOT_SET', dbmain='NOT_SET', dbsolpstop='NOT_SET'
      character*8192 dbline

      contains

#ifdef MDSPLUS
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      integer function descr2_str0(idtype,iflag,ilen)
#     include "mdslib.inc"
      integer, intent(in) :: idtype
      integer, intent(in) :: iflag
      integer, intent(in) :: ilen
      descr2_str0 = descr2(idtype,iflag,ilen)
      end function descr2_str0

      integer function descr2_str1(idtype,il1,iflag,ilen)
#     include "mdslib.inc"
      integer, intent(in) :: idtype
      integer, intent(in) :: iflag
      integer, intent(in) :: il1, ilen
      descr2_str1 = descr2(idtype,il1,iflag,ilen)
      end function descr2_str1

      integer function descr2_i0(idtype,iflag)
#     include "mdslib.inc"
      integer, intent(in) :: idtype
      integer, intent(in) :: iflag
      descr2_i0 = descr2(idtype,iflag)
      end function descr2_i0

      integer function descr2_i1(idtype,il1,iflag)
#     include "mdslib.inc"
      integer, intent(in) :: idtype
      integer, intent(in) :: iflag
      integer, intent(in) :: il1
      descr2_i1 = descr2(idtype,il1,iflag)
      end function descr2_i1

      integer function descr2_i2(idtype,il1,il2,iflag)
#     include "mdslib.inc"
      integer, intent(in) :: idtype
      integer, intent(in) :: iflag
      integer, intent(in) :: il1, il2
      descr2_i2 = descr2(idtype,il1,il2,iflag)
      end function descr2_i2

      integer function descr2_i3(idtype,il1,il2,il3,iflag)
#     include "mdslib.inc"
      integer, intent(in) :: idtype
      integer, intent(in) :: iflag
      integer, intent(in) :: il1, il2, il3
      descr2_i3 = descr2(idtype,il1,il2,il3,iflag)
      end function descr2_i3

      integer function descr2_i4(idtype,il1,il2,il3,il4,iflag)
#     include "mdslib.inc"
      integer, intent(in) :: idtype
      integer, intent(in) :: iflag
      integer, intent(in) :: il1, il2, il3, il4
      descr2_i4 = descr2(idtype,il1,il2,il3,il4,iflag)
      end function descr2_i4

      integer function mdsput2_i0(name,units,description,i0,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      integer, intent(in) :: i0
      integer, intent(in) :: iflg
      mdsput2_i0 = mdsput2(name,units,description,i0,iflg)
      end function mdsput2_i0

      integer function mdsput2_i1d(name,units,description,i1d,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      integer, allocatable, intent(in) :: i1d(:)
      integer, intent(in) :: iflg
      mdsput2_i1d = mdsput2(name,units,description,i1d,iflg)
      end function mdsput2_i1d

      integer function mdsput2_i2d(name,units,description,i2d,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      integer, allocatable, intent(in) :: i2d(:,:)
      integer, intent(in) :: iflg
      mdsput2_i2d = mdsput2(name,units,description,i2d,iflg)
      end function mdsput2_i2d

      integer function mdsput2_d0(name,units,description,d0,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      real (kind=R8), intent(in) :: d0
      integer, intent(in) :: iflg
      mdsput2_d0 = mdsput2(name,units,description,d0,iflg)
      end function mdsput2_d0

      integer function mdsput2_d1d(name,units,description,d1d,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      real (kind=R8), allocatable, intent(in) :: d1d(:)
      integer, intent(in) :: iflg
      mdsput2_d1d = mdsput2(name,units,description,d1d,iflg)
      end function mdsput2_d1d

      integer function mdsput2_d2d(name,units,description,d2d,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      real (kind=R8), allocatable, intent(in) :: d2d(:,:)
      integer, intent(in) :: iflg
      mdsput2_d2d = mdsput2(name,units,description,d2d,iflg)
      end function mdsput2_d2d

      integer function mdsput2_d3d(name,units,description,d3d,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      real (kind=R8), allocatable, intent(in) :: d3d(:,:,:)
      integer, intent(in) :: iflg
      mdsput2_d3d = mdsput2(name,units,description,d3d,iflg)
      end function mdsput2_d3d

      integer function mdsput2_d4d(name,units,description,d4d,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      real (kind=R8), allocatable, intent(in) :: d4d(:,:,:,:)
      integer, intent(in) :: iflg
      mdsput2_d4d = mdsput2(name,units,description,d4d,iflg)
      end function mdsput2_d4d

      integer function mdsput2_str(name,units,description,string,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
      character*(*), intent(in) :: string
      integer, intent(in) :: iflg
      mdsput2_str = mdsput2(name,units,description,string,iflg)
      end function mdsput2_str

      integer function mdsput2_str1(name,units,description,string,iflg)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      character*(*), intent(in) :: units
      integer, intent(in) :: description
#ifndef LEGACYCOMP
      character(len=:), allocatable, intent(in) :: string(:)
#else
      integer, parameter :: strmax = 1024
      character(len=strmax+1), allocatable, intent(in) :: string(:)
#endif
      integer, intent(in) :: iflg
      mdsput2_str1 = mdsput2(name,units,description,string,iflg)
      end function mdsput2_str1

      integer function mdsvalue2_i0(name,description,i0,iflg,retlen)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      integer, intent(in) :: description
      integer, intent(out) :: i0
      integer, intent(in) :: iflg
      integer, intent(out) :: retlen
      mdsvalue2_i0 = mdsvalue2(name,description,i0,iflg,retlen)
      end function mdsvalue2_i0

      integer function mdsvalue2_i2(name,description,i2d,iflg,retlen)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      integer, intent(in) :: description
      integer, intent(out), allocatable :: i2d(:,:)
      integer, intent(in) :: iflg
      integer, intent(out) :: retlen
      mdsvalue2_i2 = mdsvalue2(name,description,i2d,iflg,retlen)
      end function mdsvalue2_i2

      integer function mdsvalue2_r1(name,description,r1d,iflg,retlen)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      integer, intent(in) :: description
      real(kind=R8), intent(out), allocatable :: r1d(:)
      integer, intent(in) :: iflg
      integer, intent(out) :: retlen
      mdsvalue2_r1 = mdsvalue2(name,description,r1d,iflg,retlen)
      end function mdsvalue2_r1

      integer function mdsvalue2_r2(name,description,r2d,iflg,retlen)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      integer, intent(in) :: description
      real(kind=R8), intent(out), allocatable :: r2d(:,:)
      integer, intent(in) :: iflg
      integer, intent(out) :: retlen
      mdsvalue2_r2 = mdsvalue2(name,description,r2d,iflg,retlen)
      end function mdsvalue2_r2

      integer function mdsvalue2_r3(name,description,r3d,iflg,retlen)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      integer, intent(in) :: description
      real(kind=R8), intent(out), allocatable :: r3d(:,:,:)
      integer, intent(in) :: iflg
      integer, intent(out) :: retlen
      mdsvalue2_r3 = mdsvalue2(name,description,r3d,iflg,retlen)
      end function mdsvalue2_r3

      integer function mdsvalue2_r4(name,description,r4d,iflg,retlen)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      integer, intent(in) :: description
      real(kind=R8), intent(out), allocatable :: r4d(:,:,:,:)
      integer, intent(in) :: iflg
      integer, intent(out) :: retlen
      mdsvalue2_r4 = mdsvalue2(name,description,r4d,iflg,retlen)
      end function mdsvalue2_r4

      integer function mdsvalue2_str(name,description,str0,iflg,retlen)
#     include "mdslib.inc"
      character*(*), intent(in) :: name
      integer, intent(in) :: description
      character*(*), intent(out) :: str0
      integer, intent(in) :: iflg
      integer, intent(out) :: retlen
      mdsvalue2_str = mdsvalue2(name,description,str0,iflg,retlen)
      end function mdsvalue2_str
#endif
#endif

      subroutine b2mod_mdsplus_write_to_db(nwrite)
      use b2mod_subsys
      implicit none
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer status, retlen
#endif
      integer nwrite
      character*16 shotnumber

      call subini ('b2mod_mdsplus_write_to_db')
      write(dbline,
     &     '(a,a,i10,a,3a,a,i10,a,a,f20.12,a,9a,3(a,i6,a),3a,'//
     &     '22(a,e30.20,a),19(a,a,a),2(a,i10,a),6a,a,a)')
     &    'sqlshot("insert into id set ',
     &    'solpsindex=',dbsolpsindex,',',
     &    'exp=''',trim(dbexp),''',',
     &    'shot=',dbshot,',',
     &    'time=',dbtime,',',
     &    'user=''',trim(dbuser),''',',
     &    'host=''',trim(dbhost),''',',
     &    'date=''',trim(dbdate),''',',
     &    'nCv=',dbnCv,',',
     &    'nFc=',dbnFc,',',
     &    'ns=',dbns,',',
     &    'species=''',trim(dbspecies),''',',
     &    'impsep_te=',dbimpsep_te,',',
     &    'impsep_ti=',dbimpsep_ti,',',
     &    'impsep_ne=',dbimpsep_ne,',',
     &    'impsep_d=',dbimpsep_d,',',
     &    'impsep_kye=',dbimpsep_kye,',',
     &    'impsep_kyi=',dbimpsep_kyi,',',
     &    'impsep_kyi0=',dbimpsep_kyi0,',',
     &    'impsep_cr=',dbimpsep_cr,',',
     &    'impsep_cz=',dbimpsep_cz,',',
     &    'ompsep_te=',dbompsep_te,',',
     &    'ompsep_ti=',dbompsep_ti,',',
     &    'ompsep_ne=',dbompsep_ne,',',
     &    'ompsep_d=',dbompsep_d,',',
     &    'ompsep_kye=',dbompsep_kye,',',
     &    'ompsep_kyi=',dbompsep_kyi,',',
     &    'ompsep_kyi0=',dbompsep_kyi0,',',
     &    'ompsep_cr=',dbompsep_cr,',',
     &    'ompsep_cz=',dbompsep_cz,',',
     &    'sep_fhe=',dbsep_fhe,',',
     &    'sep_fhi=',dbsep_fhi,',',
     &    'sep_fne=',dbsep_fne,',',
     &    'sep_fni=',dbsep_fni,',',
     &    'directory=''',trim(dbdirectory),''',',
     &    'comments=''',trim(dbcomments),''',',
     &    'geometry=''',trim(dbgeometry),''',',
     &    'physics=''',trim(dbphysics),''',',
     &    'main=''',trim(dbmain),''',',
     &    'solpstop=''',trim(dbsolpstop),''',',
     &    'objectcode=''',trim(dbobjectcode),''',',
     &    'b2fstatiid=''',trim(dbb2fstatiid),''',',
     &    'b2fstateid=''',trim(dbb2fstateid),''',',
     &    'b2fgmtryid=''',trim(dbb2fgmtryid),''',',
     &    'b2fratesid=''',trim(dbb2fratesid),''',',
     &    'b2fplasmaid=''',trim(dbb2fplasmaid),''',',
     &    'b2fmovie_version=''',trim(db_version_b2fmovie),''',',
     &    'b2fplasma_version=''',trim(db_version_b2fplasma),''',',
     &    'b2frates_version=''',trim(db_version_b2frates),''',',
     &    'b2fstate_version=''',trim(db_version_b2fstate),''',',
     &    'b2fstati_version=''',trim(db_version_b2fstati),''',',
     &    'b2ftrack_version=''',trim(db_version_b2ftrack),''',',
     &    'plasmastate_version=''',trim(db_version_plasmastate),''',',
     &    'svn_version=',dbsvn,',',
     &    'version=',dbversion,',',
     &    'git_version=''',trim(dbgit),''',',
     &    'solpsversion=''',trim(dbsolpsversion),'''',
     &    '","select last_insert_id()",',
     &    '"solps")'
      write(nwrite,'(a,a)') '    DBLINE : ',trim(dbline)
      
#ifdef MDSPLUS
      status=mdsvalue2(
     *     trim(dbline) // char(0),
     *     descr2(idtype_cstring,0,16),shotnumber,
     *     0,retlen)
      if(mod(status,2).eq.1) then
         write(*,*) 'shotnumber: ',shotnumber
      else
         call mds_err_msg(status)
         call xerrab ('Invalid status value')
      endif
#else
      shotnumber = '0'
#endif

      if(shotnumber.ne.'0') then
        write(dbline,
     &   '(a,3a,a,i10,a,a,f20.12,a,9a,3(a,i6,a),3a,'//
     &    '22(a,e30.20,a),19(a,a,a),2(a,i10,a),6a,a,i10,a,a)')
     &   'sqlshot("update id set ',
     &   'exp=''',trim(dbexp),''',',
     &   'shot=',dbshot,',',
     &   'time=',dbtime,',',
     &   'user=''',trim(dbuser),''',',
     &   'host=''',trim(dbhost),''',',
     &   'date=''',trim(dbdate),''',',
     &   'nCv=',dbnCv,',',
     &   'nFc=',dbnFc,',',
     &   'ns=',dbns,',',
     &   'species=''',trim(dbspecies),''',',
     &   'impsep_te=',dbimpsep_te,',',
     &   'impsep_ti=',dbimpsep_ti,',',
     &   'impsep_ne=',dbimpsep_ne,',',
     &   'impsep_d=',dbimpsep_d,',',
     &   'impsep_kye=',dbimpsep_kye,',',
     &   'impsep_kyi=',dbimpsep_kyi,',',
     &   'impsep_kyi0=',dbimpsep_kyi0,',',
     &   'impsep_cr=',dbimpsep_cr,',',
     &   'impsep_cz=',dbimpsep_cz,',',
     &   'ompsep_te=',dbompsep_te,',',
     &   'ompsep_ti=',dbompsep_ti,',',
     &   'ompsep_ne=',dbompsep_ne,',',
     &   'ompsep_d=',dbompsep_d,',',
     &   'ompsep_kye=',dbompsep_kye,',',
     &   'ompsep_kyi=',dbompsep_kyi,',',
     &   'ompsep_kyi0=',dbompsep_kyi0,',',
     &   'ompsep_cr=',dbompsep_cr,',',
     &   'ompsep_cz=',dbompsep_cz,',',
     &   'sep_fhe=',dbsep_fhe,',',
     &   'sep_fhi=',dbsep_fhi,',',
     &   'sep_fne=',dbsep_fne,',',
     &   'sep_fni=',dbsep_fni,',',
     &   'directory=''',trim(dbdirectory),''',',
     &   'comments=''',trim(dbcomments),''',',
     &   'geometry=''',trim(dbgeometry),''',',
     &   'physics=''',trim(dbphysics),''',',
     &   'main=''',trim(dbmain),''',',
     &   'solpstop=''',trim(dbsolpstop),''',',
     &   'objectcode=''',trim(dbobjectcode),''',',
     &   'b2fstatiid=''',trim(dbb2fstatiid),''',',
     &   'b2fstateid=''',trim(dbb2fstateid),''',',
     &   'b2fgmtryid=''',trim(dbb2fgmtryid),''',',
     &   'b2fratesid=''',trim(dbb2fratesid),''',',
     &   'b2fplasmaid=''',trim(dbb2fplasmaid),''',',
     &   'b2fmovie_version=''',trim(db_version_b2fmovie),''',',
     &   'b2fplasma_version=''',trim(db_version_b2fplasma),''',',
     &   'b2frates_version=''',trim(db_version_b2frates),''',',
     &   'b2fstate_version=''',trim(db_version_b2fstate),''',',
     &   'b2fstati_version=''',trim(db_version_b2fstati),''',',
     &   'b2ftrack_version=''',trim(db_version_b2ftrack),''',',
     &   'plasmastate_version=''',trim(db_version_plasmastate),''',',
     &   'svn_version=',dbsvn,',',
     &   'version=',dbversion,',',
     &   'git_version=''',trim(dbgit),''',',
     &   'solpsversion=''',trim(dbsolpsversion),'''',
     &   ' where solpsindex=',dbsolpsindex,
     &   '","select last_insert_id()",',
     &   '"solps")'
        write(nwrite,'(a,a)') '    DBLINE : ',trim(dbline)

#ifdef MDSPLUS
        status=mdsvalue2(
     *   trim(dbline) // char(0),
     *   descr2(idtype_cstring,0,16),shotnumber,
     *   0,retlen)
        if(mod(status,2).eq.1) then
          write(*,*) 'shotnumber: ',shotnumber
        else
          call mds_err_msg(status)
          call xerrab ('Invalid status value')
        endif
#endif

      endif
      call subend ()
      return
      end subroutine b2mod_mdsplus_write_to_db

      end module b2mod_mdsplus
