c  version : 25.06.97 21:29
c
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_indirect

      integer wklng
      parameter (wklng=10)

      integer nlabels,nfns,nparameters,nsparameters,nconstants,npart
c typeofplot=1 physical domain only, =2 computational domain only, =3 both
      parameter(nlabels=112,nfns=130,nparameters=23,nsparameters=16,
     *     nconstants=7,npart=100000)

      real (kind=R8) ::
     *   values(nparameters),constant_value(nconstants),
     *   arrowl,earrowl,arroww,sepwidth,lwidth,wwidth,owidth,
     *   ewidth,vwidth,zwidth
      integer dimensions(nlabels),nidim,linlog,wlinlog
      integer ndim(2,wklng),typeofplot,iwrk,movie_skip,average
cank 960511
      integer nlimi,nstsi,natmi,nmoli,nioni
cank
      integer mcol,scol,vcol,ocol,ecol,acol,zcol,y_transform_l
      integer dcol,wcol,tcol,lfsz

      logical writing_to_file,color_set,use_compile,
     * set_resignore_regions_to_zero,plotting_wall

      character*4 labels(nlabels),fns(nfns),
     *	   parameters(nparameters),sparameters(nsparameters),
     *     constant_name(nconstants)
      character*64 descriptions(nlabels)
      character*64 graphlabel(wklng)
      character*60 runlabel
      character*256 solpstop
      character*2 elements(92)
      character*80 b2version
      character*256 moviename
      character*256 expfile

      save
     r values,constant_value,
     r arroww,arrowl,earrowl,sepwidth,lwidth,wwidth,owidth,
     r ewidth,vwidth,zwidth,
     r dimensions,ndim,typeofplot,nidim,linlog,wlinlog,iwrk,
     i nlimi,nstsi,natmi,nmoli,nioni,movie_skip,average,
     i mcol,scol,vcol,ocol,ecol,acol,zcol,y_transform_l,
     i dcol,wcol,tcol,lfsz,
     l writing_to_file,color_set,use_compile,
     l set_resignore_regions_to_zero,plotting_wall,
     c labels,descriptions,fns,parameters,sparameters,
     c constant_name,
     c graphlabel,runlabel,solpstop,elements,
     c b2version,moviename,expfile

