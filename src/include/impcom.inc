      parameter (nform=4, maxkey=8, nnum=17,
     ,   lterm=5, lenw=255, maxlin=255)
      logical dumanl,immed,nonam,ybrl,ybrr,yastr,
     ,                    exls,erenc,done,fixkey,yline
      character blnc,blnam*(maxkey),term,aform,number,
     ,                    ln,wr,fd*12,fe*12,fi*10,fr*8
      character fin*8,nam*8,promp*16
      real*8 ww(4),oldnam
      parameter (blnc=' ', blnam=' ')
      common/impcon/fin,promp,jmess,mimp,mess,keylen,
     ,                     lfin,lena,mfr,immed,fixkey
      common/impint/ww,nam,oldnam,len(nform),
     ,       lbeg,ka,kf,kw,krpb,krpf,nrpb,nrps,leng,
     ,       nin,ninp,kline,lim,dumanl,nonam,ybrl,
     ,                ybrr,yastr,exls,erenc,done,yline
      common/impchr/aform(nform),number(nnum),
     ,     ln(maxlin),term(lterm),wr(lenw),fd,fe,fi,fr
