      include "other_code.f"
      include "mod_alloys.f"            
      include "mod_stress.f"            
      include "mod_wkcoup.f"            

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Uexternaldb...                                                        +
c +       lop = 0 beginning of analysis                                     +
c +             1 start of increment                                        +
c +             2 end of increment                                          +
c +             3 end of analysis                                           +
c +             4 beginning of restart                                      +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      use mod_wkcoup
      implicit none
      integer lop,lrestart,kstep,kinc
      real(8) time(2),dtime
      real(8) x1,x2,x3
      real(8) v3_1(3),v3_2(3),v3_3(3)
      integer i,j,ix,ip
      integer ix1,ix2,ix3,ix4,ix5
      integer iv5_1(5)
      integer np_x,np_y,np_z
      character(len=255) :: Fpathx
      character(len=255) :: PathFilex
      integer get_thread_id,getnumthreads,I_thread,N_thread
c
      I_thread=0
      N_thread=1
c      I_thread=get_thread_id()
c      N_thread=getnumthreads()
c      
      if(lop==0)then   !ini
         print '("First call",2I3,I5,10E15.5)' 
c     &         ,I_thread,N_thread,kinc,dtime,time
         call sub_wkcoup_ini
      endif
      if(lop==1)then
         print '("dt-start",  2I3,I5,10E15.5)' 
c     &         ,I_thread,N_thread,kinc,dtime,time
         call sub_wkcoup_evl
      endif
      if(lop==2)then
         print '("dt-end",    2I3,I5,10E15.5)' 
c     &         ,I_thread,N_thread,kinc,dtime,time
      endif
      

      return
      end



c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +  Abaqus Umat, Abaqus Umat, Abaqus Umat, Abaqus Umat, Abaqus Umat        +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine umat(stress,statev,ddsdde,
     &                sse,spd,scd,rpl,
     &                ddsddt,drplde,drpldt,stran,dstran, 
     &                time,dtime,temp,dtemp,
     &                predef,dpred,cmname,ndi,nshr,ntens, 
     &                nstatv,props,nprops,
     &                coords,drot,pnewdt,celent, 
     &                dfgrd0,dfgrd1,noel,npt,layer,
     &                kspt,kstep,kinc)
c                       
c      include "aba_param.inc"
c
      use mod_wkcoup
      implicit none
      character*8 cmname
      integer noel,npt,layer,kspt,kstep,kinc
      integer ntens,nstatv,nprops,ndi,nshr
      real(8) stress(ntens), statev(nstatv), 
     &        ddsdde(ntens, ntens),
     &        ddsddt(ntens), drplde(ntens), 
     &        stran(ntens), dstran(ntens),
     &        predef(1), dpred(1), props(nprops), 
     &        coords(3), drot(3,3),
     &        dfgrd0(3,3), dfgrd1(3,3)
      real(8) dtime,time(2),sse,spd,scd,rpl
      real(8) pnewdt,celent,temp,dtemp,drpldt
c
      integer ising,Ialloy,ie,ig,Nslp
      integer get_thread_id,getnumthreads,I_thread,N_thread
      real(8) Fg0(3,3),Fg(3,3),eang00(4),eang(4),QM(3,3)
      real(8) cs(6),Fp(3,3),Fe(3,3),MatJacb(6,6)
      real(8) IVB(48),IVB_gnd(48),IVB_trp(48),gam(48),bstrs(48)
      real(8) pk2i(6),pk2i_gnd(6),pk2i_msf(48,6),gamvf(48)
      real(8) Vtrp(48),IFtrp(3,3),Lprs(3,3)
      real(8) dt1,dgmdt(48)
c
      integer i,j,k
      real(8) x1,x2,x3
      integer IB1(9),IB2(9)
      real(8) XI33(3,3),XI66(6,6),XI99(9,9),XInn(48,48)
      real(8) XI333(3,3,3)
c--------------------------------------------------------------------
      I_thread=0
      N_thread=1
c      I_thread=get_thread_id()
c      N_thread=getnumthreads()
c      if(I_thread==0)then
c         print '(I3,2I5,10e15.3,10e15.5)', 
c     &   I_thread,noel,npt,dtime,time(1),time(2)
c      endif
c--------------------------------------------------------------------

      call constantsIDT(IB1,IB2,XI33,XI333,XI66,XI99,XInn)

      ie           = noel
      ig           = npt
      dt1          = dtime
      Fg0          = dfgrd0
      Fg           = dfgrd1
      ialloy       = int(props(1))
      eang00(1:3)  = props(2:4)
      cs    (1:6)  = stress( 1: 6)
      eang  (1:4)  = statev( 1: 4)                 
      pk2i  (1:6)  = statev( 5:10)                 
      Fp    (1,:)  = statev(11:13)            
      Fp    (2,:)  = statev(14:16)              
      Fp    (3,:)  = statev(17:19)              
      Fe    (1,:)  = statev(20:22)              
      Fe    (2,:)  = statev(23:25)              
      Fe    (3,:)  = statev(26:28)              
      IVB   (1:48) = statev(28+1:28+48)              
      gam   (1:48) = statev(28+48+1:28+48*2)              

      IVB_gnd      = 0
      pk2i_gnd     = 0
      Lprs         = fem_Lprs(ie,ig,:,:)
      bstrs        = fem_bastress(ie,ig,:)  

c--------------------------------------------------------------------
      call calculate_stress_IVB_matJacb( 
     &     ie,ig,ialloy,eang00,Fg0,Fg,dt1,time(1),time(2),
     &     IVB_gnd,pk2i_gnd,Lprs,bstrs,
     &     pk2i,Fe,Fp,eang,IVB,dgmdt,gam,Nslp,
     &     cs,MatJacb,ising)      
c--------------------------------------------------------------------
c      print*,ie,ig,ising
c      call pm(MatJacb,6,6)

      if(ising/=0)then
         pnewdt=0.5d0
c         if(I_thread==0)then
         print*,'umat did not converge, err=',ising 
         print '(10I5)',I_thread,ie,ig,Ialloy,Nslp
         call pm(Fg0,3,3)
         call pm(Fg,3,3)
         call pm(Fe,3,3)
         call pm(Fp,3,3)
         call pv(IVB(1:12),12)
c         endif
         return
      else         
         pnewdt=1.25
         if(ntens==4)then
            stress(1:4)=cs(1:4)
            ddsdde(1:4,1:4)=MatJacb(1:4,1:4)
         elseif(ntens==6)then
            stress=cs
            ddsdde=MatJacb
         else
            print*,'-----------------------------------------'
            print*
            print*,'this code can not work under plane stress!'
            print*
            print*,'-----------------------------------------'
            call xit
         endif
         statev( 1: 4)=eang                 
         statev( 5:10)=pk2i                 
         statev(11:13)=Fp(1,:)              
         statev(14:16)=Fp(2,:)              
         statev(17:19)=Fp(3,:)              
         statev(20:22)=Fe(1,:)              
         statev(23:25)=Fe(2,:)              
         statev(26:28)=Fe(3,:)              
         statev(28+1:28+Nslp)=IVB(1:Nslp)              
         statev(28+48+1:28+48+Nslp)=gam(1:Nslp)              

         fem_dgmdt(ie,ig,1:Nslp)=dgmdt(1:Nslp)
         fem_pk2i(ie,ig,:)=pk2i
         zeit1=time(2)
         dzeit1=dt1

      endif         
c
      return
      end
      
      
