      module b2mod_mdsplus
      use b2mod_types
      use b2mod_constants
*   ..database quantities
      integer :: dbsolpsindex=0, dbshot=0, dbnx=0, dbny=0, dbns=0,
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

      subroutine b2mod_mdsplus_write_to_db(nwrite)
      implicit none
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer status, retlen
#endif
      integer nwrite
      character*16 shotnumber
      external subini, subend, xerrab

      call subini ('b2mod_mdsplus_write_to_db')
      write(dbline,
     &     '(a,a,i10,a,3a,a,i10,a,a,f20.12,a,9a,3(a,i6,a),3a,
     &     22(a,e30.20,a),19(a,a,a),2(a,i10,a),6a,a,a)')
     &    'sqlshot("insert into id set ',
     &    'solpsindex=',dbsolpsindex,',',
     &    'exp=''',trim(dbexp),''',',
     &    'shot=',dbshot,',',
     &    'time=',dbtime,',',
     &    'user=''',trim(dbuser),''',',
     &    'host=''',trim(dbhost),''',',
     &    'date=''',trim(dbdate),''',',
     &    'nx=',dbnx,',',
     &    'ny=',dbny,',',
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
     &   '(a,3a,a,i10,a,a,f20.12,a,9a,3(a,i6,a),3a,
     &   22(a,e30.20,a),19(a,a,a),2(a,i10,a),6a,a,i10,a,a)')
     &   'sqlshot("update id set ',
     &   'exp=''',trim(dbexp),''',',
     &   'shot=',dbshot,',',
     &   'time=',dbtime,',',
     &   'user=''',trim(dbuser),''',',
     &   'host=''',trim(dbhost),''',',
     &   'date=''',trim(dbdate),''',',
     &   'nx=',dbnx,',',
     &   'ny=',dbny,',',
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
