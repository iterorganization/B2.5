      real xlo,xhi,ylo,yhi,lwidth,xsize,ysize,xpos,ypos,
     . wwidth,xwidth,zwidth,yheight,lbsz,
     . wallbar_left,wallbar_right,colorbar_left,colorbar_right,
     . wmin,wmax,trgmin,trgmax,targetbar_left,targetbar_right
      real clval(1000),wclval(1000),tclval(1000),
     . runtime,partition_colors

      integer nx_page,ny_page,ix_page,iy_page,
     . ncol,zcol,hcol,gcol,bcol,xcol,lfsz,
     . nofcls,nofwcls,noftcls,ncont,
     . fill_option,same_page,ndec,ngraph,wlinlog,nwork,noax,
     . trglinlog

      logical dirty_graph,dirty_page,new_page,vessel,outline,
     . colprt,lalong,labels,contours,grid,wall,intpol,sptrix,
     . contplot,initialised_plotting,logx,logy,logz,target,
     . want_aspect,lzero,label_same_lines,color_map_loaded

      character*256 extralabel,globalheader,device,wallbarlabel,
     . targetbarlabel,color_file

      common /graphics_vars/
     r xlo,xhi,ylo,yhi,xsize,ysize,xpos,ypos,wwidth,zwidth,
     . xwidth,yheight,clval,wclval,tclval,runtime,lwidth,lbsz,
     . wallbar_left,wallbar_right,colorbar_left,colorbar_right,
     . wmin,wmax,trgmin,trgmax,partition_colors,
     . targetbar_left,targetbar_right,
     i nx_page,ny_page,ix_page,iy_page,ncont,fill_option,same_page,
     . ncol,zcol,hcol,gcol,bcol,xcol,nofcls,nofwcls,noftcls,lfsz,
     . ndec,ngraph,wlinlog,nwork,noax,trglinlog,
     l dirty_graph,dirty_page,new_page,labels,contours,grid,wall,
     . intpol,sptrix,vessel,outline,contplot,initialised_plotting,
     . colprt,lalong,logx,logy,logz,want_aspect,lzero,label_same_lines,
     . target,color_map_loaded,
     c extralabel,globalheader,device,wallbarlabel,targetbarlabel,
     . color_file

      save /graphics_vars/
