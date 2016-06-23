      module b2mod_tallies
      use b2mod_types
      implicit none
      real (kind=R8), allocatable, save ::
     *  rsanareg(:,:),
     *  rsahireg(:,:),
     *  rsamoreg(:,:),
     *  rranareg(:,:),
     *  rrahireg(:,:),
     *  rramoreg(:,:),
     *  rqahereg(:,:),
     *  rqradreg(:,:),
     *  rqbrmreg(:,:),
     *  rcxnareg(:,:),
     *  rcxhireg(:,:),
     *  rcxmoreg(:,:)
      real (kind=R8), allocatable, save ::
     *  fnaxreg(:,:),
     *  fnayreg(:,:),
     *  fhmxreg(:,:),
     *  fhmyreg(:,:),
     *  fhpxreg(:,:),
     *  fhpyreg(:,:),
     *  fhixreg(:),
     *  fhiyreg(:),
     *  fhexreg(:),
     *  fheyreg(:),
     *  fchxreg(:),
     *  fchyreg(:),
     *  fhtxreg(:),                                               !cjk
     *  fhtyreg(:),
     *  fhjxreg(:),
     *  fhjyreg(:)
      real (kind=R8), allocatable, save ::
     *  nareg(:,:),
     *  tereg(:),
     *  nereg(:),
     *  ne2reg(:),
     *  tireg(:),
     *  nireg(:),
     *  poreg(:),
     *  volreg(:)
      real (kind=R8), allocatable, save ::                        !cjk
     *  qconvixreg(:,:),
     *  qconviyreg(:,:),
     *  qconvexreg(:),
     *  qconveyreg(:)
      real (kind=R8), allocatable, save ::
     *  b2stbr_sna_reg(:,:),
     *  b2stbr_smo_reg(:,:),
     *  b2stbr_sne_reg(:),
     *  b2stbr_she_reg(:),
     *  b2stbr_shi_reg(:),
     *  b2stbr_sch_reg(:),
     *  b2stbc_sna_reg(:,:),
     *  b2stbc_smo_reg(:,:),
     *  b2stbc_sne_reg(:),
     *  b2stbc_she_reg(:),
     *  b2stbc_shi_reg(:),
     *  b2stbc_sch_reg(:),
     *  b2stbm_sna_reg(:,:),
     *  b2stbm_smo_reg(:,:),
     *  b2stbm_sne_reg(:),
     *  b2stbm_she_reg(:),
     *  b2stbm_shi_reg(:),
     *  b2stbm_sch_reg(:),
     *  b2sext_sna_reg(:,:),
     *  b2sext_smo_reg(:,:),
     *  b2sext_sne_reg(:),
     *  b2sext_she_reg(:),
     *  b2sext_shi_reg(:),
     *  b2sext_sch_reg(:)
      real (kind=R8), allocatable, save ::
     *  b2qie(:),
     *  b2divue(:),
     *  b2divua(:),
     *  b2exbe(:),
     *  b2exba(:),
     *  b2visa(:),
     *  b2joule(:),
     *  b2fraa(:),
     *  b2str(:),
     *  b2wrong1(:),
     *  b2wrong2(:),
     *  b2wrong3(:),
     *  b2she(:),
     *  b2shi(:),
     *  b2she0(:),
     *  b2shi0(:),
     *  b2radreg(:),
     *  b2bremreg(:)
      real (kind=R8), allocatable, save ::
     *  resco_reg(:,:),
     *  resmo_reg(:,:),
     *  reshe_reg(:),
     *  reshi_reg(:),
     *  resmt_reg(:),
     *  respo_reg(:)

      contains

      subroutine alloc_b2mod_tallies(nnreg,nsd)
      implicit none
      integer nnreg(0:2), nsd

      allocate(rsanareg(0:nnreg(0),0:nsd-1))
      allocate(rsahireg(0:nnreg(0),0:nsd-1))
      allocate(rsamoreg(0:nnreg(0),0:nsd-1))
      allocate(rranareg(0:nnreg(0),0:nsd-1))
      allocate(rrahireg(0:nnreg(0),0:nsd-1))
      allocate(rramoreg(0:nnreg(0),0:nsd-1))
      allocate(rqahereg(0:nnreg(0),0:nsd-1))
      allocate(rqradreg(0:nnreg(0),0:nsd-1))
      allocate(rqbrmreg(0:nnreg(0),0:nsd-1))
      allocate(rcxnareg(0:nnreg(0),0:nsd-1))
      allocate(rcxhireg(0:nnreg(0),0:nsd-1))
      allocate(rcxmoreg(0:nnreg(0),0:nsd-1))
      rsanareg=0.0_R8
      rsahireg=0.0_R8
      rsamoreg=0.0_R8
      rranareg=0.0_R8
      rrahireg=0.0_R8
      rramoreg=0.0_R8
      rqahereg=0.0_R8
      rqradreg=0.0_R8
      rqbrmreg=0.0_R8
      rcxnareg=0.0_R8
      rcxhireg=0.0_R8
      rcxmoreg=0.0_R8

      allocate(fnaxreg(0:nnreg(1),0:nsd-1))
      allocate(fnayreg(0:nnreg(2),0:nsd-1))
      allocate(fhmxreg(0:nnreg(1),0:nsd-1))
      allocate(fhmyreg(0:nnreg(2),0:nsd-1))
      allocate(fhpxreg(0:nnreg(1),0:nsd-1))
      allocate(fhpyreg(0:nnreg(2),0:nsd-1))
      allocate(fhixreg(0:nnreg(1)))
      allocate(fhiyreg(0:nnreg(2)))
      allocate(fhexreg(0:nnreg(1)))
      allocate(fheyreg(0:nnreg(2)))
      allocate(fchxreg(0:nnreg(1)))
      allocate(fchyreg(0:nnreg(2)))
      allocate(fhtxreg(0:nnreg(1)))                                 !cjk
      allocate(fhtyreg(0:nnreg(2)))
      allocate(fhjxreg(0:nnreg(1)))
      allocate(fhjyreg(0:nnreg(2)))
      allocate(nareg(0:nnreg(0),0:nsd-1))
      allocate(tereg(0:nnreg(0)))
      allocate(nereg(0:nnreg(0)))
      allocate(ne2reg(0:nnreg(0)))
      allocate(tireg(0:nnreg(0)))
      allocate(nireg(0:nnreg(0)))
      allocate(poreg(0:nnreg(0)))
      allocate(volreg(0:nnreg(0)))
      allocate(qconvixreg(0:nnreg(1),0:nsd-1))                      !cjk
      allocate(qconviyreg(0:nnreg(2),0:nsd-1))
      allocate(qconvexreg(0:nnreg(1)))
      allocate(qconveyreg(0:nnreg(2)))
      allocate(b2stbr_sna_reg(0:nnreg(0),0:nsd-1))
      allocate(b2stbr_smo_reg(0:nnreg(0),0:nsd-1))
      allocate(b2stbr_sne_reg(0:nnreg(0)))
      allocate(b2stbr_she_reg(0:nnreg(0)))
      allocate(b2stbr_shi_reg(0:nnreg(0)))
      allocate(b2stbr_sch_reg(0:nnreg(0)))
      allocate(b2stbc_sna_reg(0:nnreg(0),0:nsd-1))
      allocate(b2stbc_smo_reg(0:nnreg(0),0:nsd-1))
      allocate(b2stbc_sne_reg(0:nnreg(0)))
      allocate(b2stbc_she_reg(0:nnreg(0)))
      allocate(b2stbc_shi_reg(0:nnreg(0)))
      allocate(b2stbc_sch_reg(0:nnreg(0)))
      allocate(b2stbm_sna_reg(0:nnreg(0),0:nsd-1))
      allocate(b2stbm_smo_reg(0:nnreg(0),0:nsd-1))
      allocate(b2stbm_sne_reg(0:nnreg(0)))
      allocate(b2stbm_she_reg(0:nnreg(0)))
      allocate(b2stbm_shi_reg(0:nnreg(0)))
      allocate(b2stbm_sch_reg(0:nnreg(0)))
      allocate(b2sext_sna_reg(0:nnreg(0),0:nsd-1))
      allocate(b2sext_smo_reg(0:nnreg(0),0:nsd-1))
      allocate(b2sext_sne_reg(0:nnreg(0)))
      allocate(b2sext_she_reg(0:nnreg(0)))
      allocate(b2sext_shi_reg(0:nnreg(0)))
      allocate(b2sext_sch_reg(0:nnreg(0)))
      allocate(b2qie(0:nnreg(0)))
      b2qie=0.0_R8
      allocate(b2divue(0:nnreg(0)))
      allocate(b2divua(0:nnreg(0)))
      allocate(b2exbe(0:nnreg(0)))
      allocate(b2exba(0:nnreg(0)))
      allocate(b2visa(0:nnreg(0)))
      allocate(b2joule(0:nnreg(0)))
      allocate(b2fraa(0:nnreg(0)))
      allocate(b2str(0:nnreg(0)))
      allocate(b2wrong1(0:nnreg(0)))
      b2wrong1=0.0_R8
      allocate(b2wrong2(0:nnreg(0)))
      b2wrong2=0.0_R8
      allocate(b2wrong3(0:nnreg(0)))
      b2wrong3=0.0_R8
      allocate(b2she(0:nnreg(0)))
      allocate(b2shi(0:nnreg(0)))
      allocate(b2she0(0:nnreg(0)))
      allocate(b2shi0(0:nnreg(0)))
      allocate(b2radreg(0:nnreg(0)))
      b2radreg=0.0_R8
      allocate(b2bremreg(0:nnreg(0)))
      b2bremreg=0.0_R8
      allocate(resco_reg(0:nnreg(0),0:nsd-1))
      allocate(resmo_reg(0:nnreg(0),0:nsd-1))
      allocate(reshe_reg(0:nnreg(0)))
      allocate(reshi_reg(0:nnreg(0)))
      allocate(resmt_reg(0:nnreg(0)))
      allocate(respo_reg(0:nnreg(0)))

      return
      end subroutine alloc_b2mod_tallies

      subroutine dealloc_b2mod_tallies
      implicit none

      deallocate(rsanareg)
      deallocate(rsahireg)
      deallocate(rsamoreg)
      deallocate(rranareg)
      deallocate(rrahireg)
      deallocate(rramoreg)
      deallocate(rqahereg)
      deallocate(rqradreg)
      deallocate(rqbrmreg)
      deallocate(rcxnareg)
      deallocate(rcxhireg)
      deallocate(rcxmoreg)
      deallocate(fnaxreg)
      deallocate(fnayreg)
      deallocate(fhmxreg)
      deallocate(fhmyreg)
      deallocate(fhpxreg)
      deallocate(fhpyreg)
      deallocate(fhixreg)
      deallocate(fhiyreg)
      deallocate(fhexreg)
      deallocate(fheyreg)
      deallocate(fchxreg)
      deallocate(fchyreg)
      deallocate(fhtxreg)
      deallocate(fhtyreg)
      deallocate(fhjxreg)
      deallocate(fhjyreg)
      deallocate(nareg)
      deallocate(tereg)
      deallocate(nereg)
      deallocate(ne2reg)
      deallocate(tireg)
      deallocate(nireg)
      deallocate(poreg)
      deallocate(volreg)
      deallocate(qconvixreg)
      deallocate(qconviyreg)
      deallocate(qconvexreg)
      deallocate(qconveyreg)
      deallocate(b2stbr_sna_reg)
      deallocate(b2stbr_smo_reg)
      deallocate(b2stbr_sne_reg)
      deallocate(b2stbr_she_reg)
      deallocate(b2stbr_shi_reg)
      deallocate(b2stbr_sch_reg)
      deallocate(b2stbc_sna_reg)
      deallocate(b2stbc_smo_reg)
      deallocate(b2stbc_sne_reg)
      deallocate(b2stbc_she_reg)
      deallocate(b2stbc_shi_reg)
      deallocate(b2stbc_sch_reg)
      deallocate(b2stbm_sna_reg)
      deallocate(b2stbm_smo_reg)
      deallocate(b2stbm_sne_reg)
      deallocate(b2stbm_she_reg)
      deallocate(b2stbm_shi_reg)
      deallocate(b2stbm_sch_reg)
      deallocate(b2sext_sna_reg)
      deallocate(b2sext_smo_reg)
      deallocate(b2sext_sne_reg)
      deallocate(b2sext_she_reg)
      deallocate(b2sext_shi_reg)
      deallocate(b2sext_sch_reg)
      deallocate(b2qie)
      deallocate(b2divue)
      deallocate(b2divua)
      deallocate(b2exbe)
      deallocate(b2exba)
      deallocate(b2visa)
      deallocate(b2joule)
      deallocate(b2fraa)
      deallocate(b2str)
      deallocate(b2wrong1)
      deallocate(b2wrong2)
      deallocate(b2wrong3)
      deallocate(b2she)
      deallocate(b2shi)
      deallocate(b2she0)
      deallocate(b2shi0)
      deallocate(b2radreg)
      deallocate(b2bremreg)
      deallocate(resco_reg)
      deallocate(resmo_reg)
      deallocate(reshe_reg)
      deallocate(reshi_reg)
      deallocate(resmt_reg)
      deallocate(respo_reg)

      return
      end subroutine dealloc_b2mod_tallies

      subroutine read_b2mod_tallies(nread)
      implicit none
      integer nread
      call xertst(0.lt.nread, 'faulty argument nread')

      return
      end subroutine read_b2mod_tallies

      subroutine write_b2mod_tallies(nwrite)
      implicit none
      integer nwrite
      call xertst(0.lt.nwrite, 'faulty argument nwrite')

      return
      end subroutine write_b2mod_tallies

#ifndef NO_CDF
      subroutine write_cdf_variable(ncid,ntstep,varname,var)
      implicit none
      integer :: ncid, ntstep
      character*(*) :: varname
      real (kind=R8) :: var(*)
#include <netcdf.inc>
      integer i
      integer :: status, varid, ndims, dimids(NF_MAX_VAR_DIMS),
     1 start(NF_MAX_VAR_DIMS), count(NF_MAX_VAR_DIMS),
     2 time_dep

      if(ntstep.ge.0) then
        time_dep=1
      else
        time_dep=0
      endif
      status=nf_inq_varid(ncid,varname,varid)
      call check_cdf_status(status)
      status=nf_inq_varndims(ncid,varid,ndims)
      call check_cdf_status(status)
      status=nf_inq_vardimid(ncid,varid,dimids)
      call check_cdf_status(status)
      do i=1,ndims-time_dep
        start(i)=1
        status=nf_inq_dimlen(ncid,dimids(i),count(i))
      enddo
      if(time_dep.eq.1) then
        start(ndims)=ntstep
        count(ndims)=1
      endif
      status=nf_put_vara_double(ncid,varid,start,count,var)
      call check_cdf_status(status)

      return
      end subroutine write_cdf_variable

      subroutine write_cdf_char_variable(ncid,ntstep,varname,var)
      implicit none
      integer :: ncid, ntstep
      character*(*) :: varname
      character*(*) :: var(*)
#include <netcdf.inc>
      integer i
      integer :: status, varid, ndims, dimids(NF_MAX_VAR_DIMS),
     1 start(NF_MAX_VAR_DIMS), count(NF_MAX_VAR_DIMS),
     2 time_dep

      if(ntstep.ge.0) then
        time_dep=1
      else
        time_dep=0
      endif
      status=nf_inq_varid(ncid,varname,varid)
      call check_cdf_status(status)
      status=nf_inq_varndims(ncid,varid,ndims)
      call check_cdf_status(status)
      status=nf_inq_vardimid(ncid,varid,dimids)
      call check_cdf_status(status)
      do i=1,ndims-time_dep
        start(i)=1
        status=nf_inq_dimlen(ncid,dimids(i),count(i))
      enddo
      if(time_dep.eq.1) then
        start(ndims)=ntstep
        count(ndims)=1
      endif
      status=nf_put_vara_text(ncid,varid,start,count,var)
      call check_cdf_status(status)

      return
      end subroutine write_cdf_char_variable

      subroutine create_cdf_var_with_attrs(
     1 ncid,varname,vartype,ndims,vdims,units,long_name)
      integer :: ncid, vartype, ndims, vdims(ndims)
      character*(*) :: varname, units, long_name
#include <netcdf.inc>
      integer :: status, varid

      status=nf_def_var(ncid,varname,vartype,ndims,vdims,varid)
      call check_cdf_status(status)
      status=nf_put_att_text(ncid,varid,'units',len(units),units)
      call check_cdf_status(status)
      status=nf_put_att_text(ncid,varid,'long_name',
     1 len(long_name),long_name)
      call check_cdf_status(status)

      return
      end subroutine create_cdf_var_with_attrs

      subroutine define_cdf_b2mod_tallies(nnreg,nsd)
      integer nnreg(0:2), nsd
#include <netcdf.inc>
      integer :: status, 
     1 ncid, vregid, xregid, yregid, nsid, tid, varid, specid
      integer :: vdims(3)

      status=nf_create('b2tallies.nc',NF_CLOBBER,ncid)
      call check_cdf_status(status)

      status=nf_def_dim(ncid,'vreg',nnreg(0)+1,vregid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'xreg',nnreg(1)+1,xregid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'yreg',nnreg(2)+1,yregid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'ns',nsd,nsid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'speclng',10,specid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'time',NF_UNLIMITED,tid)
      call check_cdf_status(status)

      vdims(1)=tid
      status=nf_def_var(ncid,'times',NF_DOUBLE,1,vdims,varid)
      call check_cdf_status(status)

      vdims(1)=specid
      vdims(2)=nsid
      call create_cdf_var_with_attrs(
     *  ncid,'species',NF_CHAR,
     *  2,vdims,'','Species')

      vdims(1)=vregid
      vdims(2)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'tereg',NF_DOUBLE,
     *  2,vdims,'eV','Average electron temperature')
      call create_cdf_var_with_attrs(
     *  ncid,'nereg',NF_DOUBLE,
     *  2,vdims,'m^-3','Average electron density')
      call create_cdf_var_with_attrs(
     *  ncid,'ne2reg',NF_DOUBLE,
     *  2,vdims,'m^-3','ne2reg')
      call create_cdf_var_with_attrs(
     *  ncid,'tireg',NF_DOUBLE,
     *  2,vdims,'eV','Average ion temperature')
      call create_cdf_var_with_attrs(
     *  ncid,'nireg',NF_DOUBLE,
     *  2,vdims,'m^-3','Average ion density (ni)')
      call create_cdf_var_with_attrs(
     *  ncid,'poreg',NF_DOUBLE,
     *  2,vdims,'eV','Average potential')
      call create_cdf_var_with_attrs(
     *  ncid,'volreg',NF_DOUBLE,
     *  2,vdims,'m^3','Volume')                     !!! change ???
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbr_sne_reg',NF_DOUBLE,
     *  2,vdims,'e','Recycling electron particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbr_she_reg',NF_DOUBLE,
     *  2,vdims,'W','Recycling electron energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbr_shi_reg',NF_DOUBLE,
     *  2,vdims,'W','Recycling ion energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbr_sch_reg',NF_DOUBLE,
     *  2,vdims,'A','Recycling charge source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbc_sne_reg',NF_DOUBLE,
     *  2,vdims,'e','Boundary electron particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbc_she_reg',NF_DOUBLE,
     *  2,vdims,'W','Boundary electron energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbc_shi_reg',NF_DOUBLE,
     *  2,vdims,'W','Boundary ion energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbc_sch_reg',NF_DOUBLE,
     *  2,vdims,'A','Boundary charge source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbm_sne_reg',NF_DOUBLE,
     *  2,vdims,'e','Additional electron particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbm_she_reg',NF_DOUBLE,
     *  2,vdims,'W','Additional electron energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbm_shi_reg',NF_DOUBLE,
     *  2,vdims,'W','Additional ion energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbm_sch_reg',NF_DOUBLE,
     *  2,vdims,'A','Additional charge source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2sext_sne_reg',NF_DOUBLE,
     *  2,vdims,'e','External electron particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2sext_she_reg',NF_DOUBLE,
     *  2,vdims,'W','External electron energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2sext_shi_reg',NF_DOUBLE,
     *  2,vdims,'W','External ion energy source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2sext_sch_reg',NF_DOUBLE,
     *  2,vdims,'A','External charge source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2qie',NF_DOUBLE,
     *  2,vdims,'W','Electron/ion exchange term')
      call create_cdf_var_with_attrs(
     *  ncid,'b2divue',NF_DOUBLE,
     *  2,vdims,'W','b2divue')
      call create_cdf_var_with_attrs(
     *  ncid,'b2divua',NF_DOUBLE,
     *  2,vdims,'W','b2divua')
      call create_cdf_var_with_attrs(
     *  ncid,'b2exbe',NF_DOUBLE,
     *  2,vdims,'W','b2exbe')
      call create_cdf_var_with_attrs(
     *  ncid,'b2exba',NF_DOUBLE,
     *  2,vdims,'W','b2exba')
      call create_cdf_var_with_attrs(
     *  ncid,'b2visa',NF_DOUBLE,
     *  2,vdims,'W','b2visa')
      call create_cdf_var_with_attrs(
     *  ncid,'b2joule',NF_DOUBLE,
     *  2,vdims,'W','b2joule')
      call create_cdf_var_with_attrs(
     *  ncid,'b2fraa',NF_DOUBLE,
     *  2,vdims,'W','b2fraa')
      call create_cdf_var_with_attrs(
     *  ncid,'b2str',NF_DOUBLE,
     *  2,vdims,'W','b2str')
      call create_cdf_var_with_attrs(
     *  ncid,'b2wrong1',NF_DOUBLE,
     *  2,vdims,'W','b2wrong1')
      call create_cdf_var_with_attrs(
     *  ncid,'b2wrong2',NF_DOUBLE,
     *  2,vdims,'W','b2wrong2')
      call create_cdf_var_with_attrs(
     *  ncid,'b2wrong3',NF_DOUBLE,
     *  2,vdims,'W','b2wrong3')
      call create_cdf_var_with_attrs(
     *  ncid,'b2she',NF_DOUBLE,
     *  2,vdims,'W','b2she')
      call create_cdf_var_with_attrs(
     *  ncid,'b2shi',NF_DOUBLE,
     *  2,vdims,'W','b2shi')
      call create_cdf_var_with_attrs(
     *  ncid,'b2she0',NF_DOUBLE,
     *  2,vdims,'W','b2she0')
      call create_cdf_var_with_attrs(
     *  ncid,'b2shi0',NF_DOUBLE,
     *  2,vdims,'W','b2shi0')
      call create_cdf_var_with_attrs(
     *  ncid,'reshe_reg',NF_DOUBLE,
     *  2,vdims,'W','reshe_reg')
      call create_cdf_var_with_attrs(
     *  ncid,'reshi_reg',NF_DOUBLE,
     *  2,vdims,'W','reshi_reg')
      call create_cdf_var_with_attrs(
     *  ncid,'resmt_reg',NF_DOUBLE,
     *  2,vdims,'kg.m.s^-1','resmt_reg')
      call create_cdf_var_with_attrs(
     *  ncid,'respo_reg',NF_DOUBLE,
     *  2,vdims,'eV','respo_reg')
      call create_cdf_var_with_attrs(
     *  ncid,'b2bremreg',NF_DOUBLE,
     *  2,vdims,'W','b2bremreg')
      call create_cdf_var_with_attrs(
     *  ncid,'b2radreg',NF_DOUBLE,
     *  2,vdims,'W','b2radreg')

      vdims(1)=vregid
      vdims(2)=nsid
      vdims(3)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'rsanareg',NF_DOUBLE,
     *  3,vdims,'s^-1','rsanareg')
      call create_cdf_var_with_attrs(
     *  ncid,'rsahireg',NF_DOUBLE,
     *  3,vdims,'W','rsahireg')
      call create_cdf_var_with_attrs(
     *  ncid,'rsamoreg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','rsamoreg')
      call create_cdf_var_with_attrs(
     *  ncid,'rranareg',NF_DOUBLE,
     *  3,vdims,'s^-1','rranareg')
      call create_cdf_var_with_attrs(
     *  ncid,'rrahireg',NF_DOUBLE,
     *  3,vdims,'W','rrahireg')
      call create_cdf_var_with_attrs(
     *  ncid,'rramoreg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','rramoreg')
      call create_cdf_var_with_attrs(
     *  ncid,'rqahereg',NF_DOUBLE,
     *  3,vdims,'W','rqahereg')
      call create_cdf_var_with_attrs(
     *  ncid,'rqradreg',NF_DOUBLE,
     *  3,vdims,'W','rqradreg')
      call create_cdf_var_with_attrs(
     *  ncid,'rqbrmreg',NF_DOUBLE,
     *  3,vdims,'W','rqbrmreg')
      call create_cdf_var_with_attrs(
     *  ncid,'rcxnareg',NF_DOUBLE,
     *  3,vdims,'s^-1','rcxnareg')
      call create_cdf_var_with_attrs(
     *  ncid,'rcxhireg',NF_DOUBLE,
     *  3,vdims,'W','rcxhireg')
      call create_cdf_var_with_attrs(
     *  ncid,'rcxmoreg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','rcxmoreg')
      call create_cdf_var_with_attrs(
     *  ncid,'nareg',NF_DOUBLE,
     *  3,vdims,'m^-3','Average ion density')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbr_sna_reg',NF_DOUBLE,
     *  3,vdims,'s^-1','Recycling particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbr_smo_reg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','Recycling parallel momentum source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbc_sna_reg',NF_DOUBLE,
     *  3,vdims,'s^-1','Boundary particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbc_smo_reg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','Boundary parallel momentum source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbm_sna_reg',NF_DOUBLE,
     *  3,vdims,'s^-1','Additional particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2stbm_smo_reg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','Additional parallel momentum source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2sext_sna_reg',NF_DOUBLE,
     *  3,vdims,'s^-1','External particle source')
      call create_cdf_var_with_attrs(
     *  ncid,'b2sext_smo_reg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','External parallel momentum source')
      call create_cdf_var_with_attrs(
     *  ncid,'resco_reg',NF_DOUBLE,
     *  3,vdims,'#','resco_reg')
      call create_cdf_var_with_attrs(
     *  ncid,'resmo_reg',NF_DOUBLE,
     *  3,vdims,'kg.m.s^-1','resmo_reg')

      vdims(1)=xregid
      vdims(2)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'fhixreg',NF_DOUBLE,
     *  2,vdims,'W','fhixreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fhexreg',NF_DOUBLE,
     *  2,vdims,'W','fhexreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fchxreg',NF_DOUBLE,
     *  2,vdims,'A','fchxreg')
      call create_cdf_var_with_attrs(                              !cjk
     *  ncid,'fhtxreg',NF_DOUBLE,
     *  2,vdims,'W','fhtxreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fhjxreg',NF_DOUBLE,
     *  2,vdims,'W','fhjxreg')

      vdims(1)=xregid                                              !cjk
      vdims(2)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'qconvexreg',NF_DOUBLE,
     *  2,vdims,'W','qconvexreg')

      vdims(1)=xregid
      vdims(2)=nsid
      vdims(3)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'fnaxreg',NF_DOUBLE,
     *  3,vdims,'s^-1','fnaxreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fhmxreg',NF_DOUBLE,
     *  3,vdims,'W','fhmxreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fhpxreg',NF_DOUBLE,
     *  3,vdims,'W','fhpxreg')
      call create_cdf_var_with_attrs(                              !cjk
     *  ncid,'qconvixreg',NF_DOUBLE,
     *  3,vdims,'W','qconvixreg')

      vdims(1)=yregid
      vdims(2)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'fhiyreg',NF_DOUBLE,
     *  2,vdims,'W','fhiyreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fheyreg',NF_DOUBLE,
     *  2,vdims,'W','fheyreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fchyreg',NF_DOUBLE,
     *  2,vdims,'A','fchyreg')
      call create_cdf_var_with_attrs(                              !cjk
     *  ncid,'fhtyreg',NF_DOUBLE,
     *  2,vdims,'W','fhtyreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fhjyreg',NF_DOUBLE,
     *  2,vdims,'W','fhjyreg')

      vdims(1)=yregid                                              !cjk
      vdims(2)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'qconveyreg',NF_DOUBLE,
     *  2,vdims,'W','qconveyreg')

      vdims(1)=yregid
      vdims(2)=nsid
      vdims(3)=tid
      call create_cdf_var_with_attrs(
     *  ncid,'fnayreg',NF_DOUBLE,
     *  3,vdims,'s^-1','fnayreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fhmyreg',NF_DOUBLE,
     *  3,vdims,'W','fhmyreg')
      call create_cdf_var_with_attrs(
     *  ncid,'fhpyreg',NF_DOUBLE,
     *  3,vdims,'W','fhpyreg')
      call create_cdf_var_with_attrs(                               !cjk
     *  ncid,'qconviyreg',NF_DOUBLE,
     *  3,vdims,'W','qconviyreg')

      status=nf_close(ncid)
      call check_cdf_status(status)

      return
      end subroutine define_cdf_b2mod_tallies
#endif

      end module b2mod_tallies
