      module b2mod_graphics
      use b2mod_types
      implicit none

      real (kind=R4) :: xlo,xhi,ylo,yhi,lwidth,xsize,ysize,xpos,ypos,
     . wwidth,xwidth,zwidth,yheight,glsz,lbsz,
     . wallbar_left,wallbar_right,colorbar_left,colorbar_right,
     . wmin,wmax,trgmin,trgmax,targetbar_left,targetbar_right
      real clval(1000),eclval(1000),wclval(1000),tclval(1000),
     . runtime,partition_colors

      integer nx_page,ny_page,ix_page,iy_page,
     . ncol,zcol,hcol,gcol,bcol,xcol,fcol,lfsz,
     . nofcls,nofecls,nofwcls,noftcls,ncont,
     . fill_option,same_page,ndec,ngraph,wlinlog,nwork,noax,
     . trglinlog

      logical dirty_graph,dirty_page,new_page,vessel,outline,
     . colprt,lalong,labels,contours,grid,trig,wall,intpol,sptrix,
     . contplot,initialised_plotting,logx,logy,logz,target,
     . want_aspect,lzero,label_same_lines,color_map_loaded,
     . nice_axes,bold

      character*256 extralabel,globalheader,device,wallbarlabel,
     . targetbarlabel,color_file,xlabel,ylabel

c     the zoom common
      integer winid
      integer wsnshwin,wsnX11

      end module b2mod_graphics
