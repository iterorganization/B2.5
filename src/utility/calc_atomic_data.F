      program calc_atomic_data

c  version : 08.08.2000 11:35

      use b2mod_types
      use b2mod_constants
      implicit none
#     include <eir_reac.F>
      integer iarg,iostat,nden,iden
      integer ir,i_k,i_l,i_m,i_n,i,j
      real (kind=R8) :: D11(0:8,0:8),density(101),csx(101)
      real (kind=R8) :: t1, t2, te, tei, tef, de, dej, def, tfac, dnnn
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      character arg*256
#ifdef NAGFOR
      integer lenval, ierror
#else
      integer iargc
#endif

#ifndef NAGFOR
      iarg=iargc()
#else
      iarg=command_argument_count()
#endif
      if(iarg.lt.6) then
        write(*,*) 'At least 6 arguments are necessary'
      endif

#ifndef NAGFOR
      call getarg(1,filnam)
      call getarg(2,h123)
      call getarg(3,reac)
      call getarg(4,crc)
      call getarg(5,arg)
#else
      call get_command_argument(1,filnam,lenval,ierror)
      call get_command_argument(2,h123,lenval,ierror)
      call get_command_argument(3,reac,lenval,ierror)
      call get_command_argument(4,crc,lenval,ierror)
      call get_command_argument(5,arg,lenval,ierror)
#endif
      read(arg,*,iostat=iostat) t1
      if(iostat.ne.0) t1=0.0_R8
#ifndef NAGFOR
      call getarg(6,arg)
#else
      call get_command_argument(6,arg,lenval,ierror)
#endif
      read(arg,*,iostat=iostat) t2
      if(iostat.ne.0) t2=0.0_R8
      if(t1.le.0.0_R8) t1=0.1_R8
      if(t2.le.t1) t2=10100.0_R8
      tfac=exp(log(t2/t1)/100.0_R8)
      nden=iarg-6
      if(nden.gt.101) nden=101
      do iden=1,nden
#ifndef NAGFOR
        call getarg(iden+6,arg)
#else
        call get_command_argument(iden+6,arg,lenval,ierror)
#endif
        read(arg,*) density(iden)
      enddo
      if(nden.le.0) then
        nden=1
        density(1)=1.0e14_R8
      endif

      if(iarg.ge.3) then
        IR=1
        CALL SLREAC_b2plot(IR,FILNAM,H123,REAC,CRC)
        write(*,'(''# '',i4)') modclf(ir)
        i_k=mod(modclf(ir)/1,10)
        i_l=mod(modclf(ir)/10,10)
        i_m=mod(modclf(ir)/100,10)
        i_n=mod(modclf(ir)/1000,10)
        DO J=1,9
          DO I=1,9
            if(i_k.eq.1) then
              D11(J-1,I-1)=CREAC(J,0,IR)
            elseif(i_l.eq.1 .or. i_n.eq.1) then
              D11(J-1,I-1)=CREAC(J,1,IR)
            elseif(i_l.ge.2.or.i_n.ge.2) then
              D11(J-1,I-1)=CREAC(J,I,IR)
            else
              write(*,*) 'Do not yet know about ',i_n,i_m,i_l,i_k
              D11(J-1,I-1)=CREAC(J,I,IR)
            endif
          enddo
        enddo
      endif
      
      write(*,'(''# start of data'')') 
      te=t1
      do while (te.lt.t2*sqrt(tfac))
        TEF=LOG(TE)
        do iden=1,nden
          de=density(iden)
          if(i_l.eq.2.or.i_n.eq.2) then
            DEF=LOG(DE)
          else
            DEF=LOG(DE*1.D-8)
          endif
          dnnn=0.0_R8
          DO I=0,8
            TEI=TEF**I
            if(i_k.eq.1.or.i_l.eq.1.or.i_n.eq.1) then
              DNNN=DNNN+D11(I,0)*TEI
            elseif(i_l.ge.2.or.i_n.ge.2.or.modclf(ir).eq.0) then
              DO J=0,8
                DEJ=DEF**J
                DNNN=DNNN+D11(I,J)*TEI*DEJ
              enddo
            else
c              write(*,*) 'Do not yet know about ',i_n,i_m,i_l,i_k
            endif
          enddo
          if(dnnn.gt.log(1.0e30_R8)) then
            dnnn=1.0e30_R8
          elseif(dnnn.lt.log(1.0e-30_R8)) then
            dnnn=1.0e-30_R8
          else
            DNNN=EXP(DNNN)
          endif
          csx(iden)=dnnn
        enddo
        write(*,'(1p,102e12.5)') te,(csx(iden),iden=1,nden)
        te=te*tfac
      enddo
      write(*,'(''# end of data'')') 

      stop

      end
