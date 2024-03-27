c================================================================
c
c    weak coupling module: mod_FEMFFT_data
c    set mapping between FEM data and FFT data
c
c================================================================
      module mod_FEMFFT_data
         implicit none
         integer, parameter :: Tnel=10000     ! FEM element number 
         integer, parameter :: Tngp=8         ! FEM integration point number per element
         integer, parameter :: Tnfx=128       ! FFT grid number along X 
         integer, parameter :: Tnfy=128       ! FFT grid number along Y
         integer, parameter :: Tnfz=2         ! FFT grid number along Z
         real(8) detx,dety,detz               ! FFT unit lenghth
         real(8) dzeit1,zeit1                 ! time-step and time
         !---FEM-FFT communication
         real(8) fem_xyz   (Tnel,Tngp,3)
         integer fem_Inf   (Tnel,Tngp,4)
         integer fft_Inf   (Tnfx,Tnfy,Tnfz,3  )
         real(8) fft_Fqc   (Tnfx,Tnfy,Tnfz,3  )
         real(8) fft_xyz   (Tnfx,Tnfy,Tnfz,3  )
        !---partitial diffusion equations
         real(8) fft_P0    (Tnfx,Tnfy,Tnfz    )
         real(8) fft_A     (Tnfx,Tnfy,Tnfz    )
         real(8) fft_B     (Tnfx,Tnfy,Tnfz,3  )
         real(8) fft_C     (Tnfx,Tnfy,Tnfz    )
         real(8) fft_S0    (Tnfx,Tnfy,Tnfz)
         real(8) fft_V0    (Tnfx,Tnfy,Tnfz,3  )
         real(8) fft_M0    (Tnfx,Tnfy,Tnfz,3,3)
         real(8) fft_S1    (Tnfx,Tnfy,Tnfz    )
         real(8) fft_V1    (Tnfx,Tnfy,Tnfz,3  )
         real(8) fft_M1    (Tnfx,Tnfy,Tnfz,3,3)
         !---crystal plasticity data
         real(8) fem_Fp00    (Tnel,Tngp,3,3)
         real(8) fem_Fp      (Tnel,Tngp,3,3)
         real(8) fem_Fe      (Tnel,Tngp,3,3)
         real(8) fem_pk2i    (Tnel,Tngp,6)
         real(8) fem_cs      (Tnel,Tngp,6)
         real(8) fem_dgmdt   (Tnel,Tngp,48)
         real(8) fem_IVB     (Tnel,Tngp,48)
         real(8) fem_gam     (Tnel,Tngp,48)    
      contains
         !c========================================
         !c initialization of mod_FEMFFT_data  
         !c========================================
         subroutine sub_FEMFFT_ini
         implicit none
         integer i,j,k,iex,igx
         real(8) x1,x2,x3
         zeit1=0
         dzeit1=0
         fem_xyz=0
         fem_Inf=0
         fft_Inf=0
         fft_Fqc=0
         fft_xyz=0
         fft_P0=0
         fft_A=0
         fft_B=0
         fft_C=0
         fft_S0=0
         fft_V0=0
         fft_M0=0
         fft_S1=0
         fft_V1=0
         fft_M1=0
         fem_Fe=0
         fem_Fp=0
         fem_Fp00=0
         fem_cs=0
         fem_pk2i=0
         fem_dgmdt=0
         fem_gam=0
         fem_IVB=0
         return
         endsubroutine
         !c========================================
         !c evolution of mod_FEMFFT_data  
         !c========================================
         subroutine sub_FEMFFT_evl
         implicit none
         integer i,j,k
         real(8) x1,x2,x3
         call sub_FEMFFT_map(Tnel,Tngp,
     &                       Tnfx,Tnfy,Tnfz,
     &                       fem_xyz,
     &                       fft_xyz,
     &                       fft_Fqc,
     &                       fft_Inf,
     &                       fem_Inf,
     &                       detx,dety,detz)
         return
         endsubroutine
      endmodule mod_FEMFFT_data

c================================================================
c
c    weak coupling module: gnd
c    output for plastic deformation calculation: IVB_gnd, pk2i_gnd
c
c      abaqus default length unit mm 
c      in nonlocal calculate we use length unit m 
c      C_unit is used to rescale the mesh size
c        use the mm mesh C_unit=1.d-3
c        for example, reduce 100  times need: C_unit=1.d-3/100
c        for example, reduce 1000 times need: C_unit=1.d-3/1000
c
c================================================================
      module mod_wkcoup_gnd
         use mod_FEMFFT_data
         implicit none
         integer, parameter :: I_gnd_iso = 1       !-->with(1),without(0)
         integer, parameter :: I_gnd_kin = 1       !-->with(1),without(0)  
         real(8), parameter :: C_taui    = 0.01    !-->Taylor hardening parameter
         real(8), parameter :: L_size    = 1.0d-9  !-->unit m, length scale
         real(8), parameter :: B_lattice = 2.5d-10 !-->unit m, Burgers vector 
         real(8), parameter :: CD_smooth = 3.d-15  !-->unit m, Burgers vector 
         real(8), parameter :: C_unit=1.d-3/1000   !-->micrometer  
         real(8) fem_dFp     (Tnel,Tngp,3,3)
         real(8) fem_dFp_1gd (Tnel,Tngp,3,3,3)
         real(8) fem_dFp_2gd (Tnel,Tngp,3,3,3,3)
         real(8) fem_IVB_gnd (Tnel,Tngp,48)
         real(8) fem_pk2i_gnd(Tnel,Tngp,6)
      contains
         !c========================================
         !c initialization of mod_wkcoup_gnd  
         !c========================================
         subroutine sub_gnd_ini
         implicit none
         integer i,j,k
         real(8) x1,x2,x3
         fem_dFp=0
         fem_dFp_1gd=0
         fem_dFp_2gd=0
         fem_IVB_gnd=0
         fem_pk2i_gnd=0
         return
         endsubroutine
         !c========================================
         !c evolution of mod_wkcoup_gnd  
         !c========================================
         subroutine sub_gnd_evl
         implicit none
         integer i,j,k,l,m,n,ix1,ix2,ix3,ii,np,ising,ising_1
         integer ip1,ip2,id,Ialloy
         integer ix,iy,iz,idx
         integer Nslp,iex,igx,is
         integer i1,i2,i3,i4,i5,i6
         integer j1,j2,j3,j4,j5,j6
         real(8) x1,x2,x3,dydx,ddyddx
         real(8) mx33_1(3,3),mx33_2(3,3),mx33_3(3,3)
         real(8) dFp(3,3),dFp_1gd(3,3,3),dFp_2gd(3,3,3,3)
         real(8) Fp0(3,3),Fp(3,3)
c
         do iex=1,Tnel
         do igx=1,Tngp
         if(fem_Inf(iex,igx,1)/=0)then
            mx33_1=fem_Fp00(iex,igx,:,:)
            mx33_2=fem_Fp(iex,igx,:,:)
            fem_dFp(iex,igx,:,:)=matmul(mx33_2,transpose(mx33_1))
         endif
         enddo
         enddo
         do ip1=1,3
         do ip2=1,3
            do ix=1,Tnfx
            do iy=1,Tnfy
            do iz=1,Tnfz
               iex=fft_Inf(ix,iy,iz,1)
               igx=fft_Inf(ix,iy,iz,2)
               fft_S0(ix,iy,iz)=fem_dFp(iex,igx,ip1,ip2)
               fft_P0(ix,iy,iz)=fft_S0(ix,iy,iz)
               fft_A (ix,iy,iz)=-1
               fft_B (ix,iy,iz,:)=0
               fft_C (ix,iy,iz)=CD_smooth
            enddo
            enddo
            enddo
            call sub_PDEsolve_FFT(Tnfx,Tnfy,Tnfz,fft_Fqc,
     &                     fft_P0,fft_A,fft_B,fft_C,
     &                     fft_S0,
     &                     fft_S1,fft_V1,fft_M1)
            do iex=1,Tnel
            do igx=1,Tngp
            if(fem_Inf(iex,igx,1)/=0)then
               ix=fem_Inf(iex,igx,2)
               iy=fem_Inf(iex,igx,3)
               iz=fem_Inf(iex,igx,4)
               fem_dFp(iex,igx,ip1,ip2)=fft_S1(ix,iy,iz)
               fem_dFp_1gd(iex,igx,ip1,ip2,:)=fft_V1(ix,iy,iz,:)
               fem_dFp_2gd(iex,igx,ip1,ip2,:,:)=fft_M1(ix,iy,iz,:,:)
            endif
            enddo
            enddo
         enddo
         enddo
         return
         endsubroutine
      endmodule mod_wkcoup_gnd

c================================================================
c
c    weak coupling module: reslips of twin
c    output for plastic deformation calculation: Fprs, IFprs, IVB_rs
c
c================================================================
      module mod_wkcoup_reslip
         use mod_FEMFFT_data
         implicit none
         real(8) fem_Lprs  (Tnel,Tngp,3,3)
         real(8) fem_IVBrs (Tnel,Tngp,6,48)
         real(8) smdMi_rs  (6,18, 3,3)       ! 6 twin systems, 18 slip systems
         real(8) smdSMi_rs (6,18, 3,3)
         real(8) smdVi1_rs (6,18, 6  ) 
         real(8) smdVi2_rs (6,18, 6  )
         real(8) hdM_rs(18,18)
      contains
         !c========================================
         !c initialization of mod_wkcoup_reslip  
         !c========================================
         subroutine sub_rs_ini
         implicit none
         real(8), parameter:: crssba0  = 1.75d0  !* tau_0 
         real(8), parameter:: crsspr0  = 25.d0
         real(8), parameter:: crsspy0  = 40.d0
         real(8), parameter:: c_cpl    = 1.d0
         real(8), parameter:: c_oth    = 1.d0
         integer iex,igx
         integer i,j,is,js 
         real(8) vd(24,3),vn(24,3)
         real(8) Qm(3,3),mx33_1(3,3),vx3_1(3),vx3_2(3) 
         real(8) XI33(3,3),XI66(6,6),XI99(9,9)
         real(8) XInn(48,48),XI333(3,3,3) 
         integer IB1(9),IB2(9)
         real(8) x1,x2,x3 
c
         call constantsIDT(IB1,IB2,XI33,XI333,XI66,XI99,XInn)
c
         vd( 1,:)=[ 3.0,  0., 0.                ] 
         vd( 2,:)=[ -1.5,  2.598076211353316, 0.] 
         vd( 3,:)=[ -1.5, -2.598076211353316, 0.] 
         vd( 4,:)=[-1.5,  2.598076211353316, 0.] 
         vd( 5,:)=[ 3.0,  0., 0.               ] 
         vd( 6,:)=[ -1.5, -2.598076211353316,0.] 
         vd(7,:)=[-1.5,  -2.598076211353316, 4.8990] 
         vd(8,:)=[-3.0,  0.,  4.8990               ] 
         vd(9,:)=[ 1.5, 2.598076211353316,   4.8990]
         vd(10,:)=[ -1.5, 2.598076211353316, 4.8990] 
         vd(11,:)=[ 3.0,  0., 4.8990               ] 
         vd(12,:)=[ 1.5, -2.598076211353316, 4.8990] 
         vd(13,:)=[ 3.0,  0., 4.8990               ] 
         vd(14,:)=[1.5, 2.598076211353316,  4.8990 ] 
         vd(15,:)=[ -1.5, -2.598076211353316,4.8990] 
         vd(16,:)=[ -1.5, 2.598076211353316, 4.8990] 
         vd(17,:)=[ -3.0,  0., 4.8990              ] 
         vd(18,:)=[ -1.5, 2.598076211353316, 4.8990] 
         vd(19,:)=[ 1.5, -0.866025403784438, 1.6330] 
         vd(20,:)=[-1.5, 0.866025403784438,  1.6330] 
         vd(21,:)=[ -1.5, -0.866025403784438,1.6330]
         vd(22,:)=[ 1.50, 0.866025403784438, 1.6330] 
         vd(23,:)=[ 0., 1.732050807568877,   1.6330]
         vd(24,:)=[ 0., -1.732050807568877,  1.6330] 
c     
         vn( 1,:)=[ 0.,  0.,  1.6330]
         vn( 2,:)=[ 0.,  0.,  1.6330]
         vn( 3,:)=[ 0.,  0.,  1.6330]
         vn( 4,:)=[ 1.5, 0.866025403784438, 0.]
         vn( 5,:)=[ 0., -1.732050807568877, 0.]
         vn( 6,:)=[ -1.5,0.866025403784438, 0.]
         vn(7,:)=[ 1.5,  0.866025403784438, 1.6330]
         vn(8,:)=[ 1.5,  0.866025403784438, 1.6330]
         vn(9,:)=[0.,  -1.732050807568877,   1.6330]
         vn(10,:)=[0.,  -1.732050807568877,  1.6330]
         vn(11,:)=[ -1.5,  0.866025403784438,1.6330]
         vn(12,:)=[ -1.5,  0.866025403784438,1.6330]
         vn(13,:)=[ -1.5,-0.866025403784438, 1.6330]
         vn(14,:)=[-1.5,  -0.866025403784438,1.6330]
         vn(15,:)=[ 0., 1.732050807568877,   1.6330]
         vn(16,:)=[ 0., 1.732050807568877,   1.6330]
         vn(17,:)=[ 1.5,  -0.866025403784438,1.6330]
         vn(18,:)=[1.5,  -0.866025403784438, 1.6330]
         vn(19,:)=[ -1.5,0.866025403784438, 3.2660]
         vn(20,:)=[1.5,  -0.866025403784438,3.2660]
         vn(21,:)=[1.5,  0.866025403784438, 3.2660]
         vn(22,:)=[-1.5, -0.866025403784438,3.2660]
         vn(23,:)=[ 0., -1.732050807568877, 3.2660]
         vn(24,:)=[0., 1.732050807568877,   3.2660]
         do is=1,24
            x1=dsqrt(sum(vd(is,:)**2)); vd(is,:)=vd(is,:)/x1
            x1=dsqrt(sum(vn(is,:)**2)); vn(is,:)=vn(is,:)/x1
         enddo
c
         do is=1,18
         do js=1,18
            x1=sum(dabs(vn(is,:)-vn(js,:)))
            if(x1<1.d-10)then
               hdM_rs(is,js)=c_cpl
            else                             
               hdM_rs(is,js)=c_oth
            endif
         enddo
         enddo
c
         do is=1,6
            do i=1,3
            do j=1,3
               Qm(i,j)=2*vn(is+18,i)*vn(is+18,j)-XI33(i,j)
            enddo
            enddo
            do js=1,18
               vx3_1=matmul(Qm,vd(js,:))
               vx3_2=matmul(Qm,vn(js,:))
               do i=1,3
               do j=1,3
                  mx33_1(i,j)=vx3_1(i)*vx3_2(j)
               enddo
               enddo
               smdMi_rs(is,js,:,:)=mx33_1
               smdSMi_rs(is,js,:,:)=(mx33_1+transpose(mx33_1))/2
               do i=1,6
                  smdVi1_rs(is,js,i)=smdSMi_rs(is,js,ib1(i),ib2(i))
               enddo
               smdVi2_rs(is,js,1:3)=smdVi1_rs(is,js,1:3)
               smdVi2_rs(is,js,4:6)=smdVi1_rs(is,js,4:6)*2
            enddo
         enddo

         do iex=1,Tnel
         do igx=1,Tngp
            do is=1,6
               fem_IVBrs(iex,igx,is, 1: 3)=crssba0
               fem_IVBrs(iex,igx,is, 4: 6)=crsspr0
               fem_IVBrs(iex,igx,is, 7:18)=crsspy0
            enddo
         enddo
         enddo
         fem_Lprs=0
c
         return
         endsubroutine
         !c========================================
         !c evolution of mod_wkcoup_reslip  
         !c========================================
         subroutine sub_rs_evl
         implicit none
         real(8),parameter :: hdrtba0  = 20.d0          !* initial hardeing of basal slip
         real(8),parameter :: hdrtpr0  = 1500.d0        !* initial hardeing of prismatic slip
         real(8),parameter :: hdrtpy0  = 3000.d0        !* initial hardeing of pyramidal slip
         real(8),parameter :: shrts0   = 1.d-5          !* reference rate   
         real(8),parameter :: pwfl     = 10.d0          !* power of slip
         real(8),parameter :: pwhd     = 0.6d0          !* power of strain hardening
         real(8),parameter :: crssbas  = 40.d0          !* initial resistance of basal slip
         real(8),parameter :: crssprs  = 85.d0          !* initial resistance of prismatic slip
         real(8),parameter :: crsspys  = 150.d0         !* initial resistance of pyramidal slip
         real(8),parameter :: gammatw  = 0.129          !* characteristic shear of twinning mode
         real(8) hdrt0(18)        
         real(8) crsss(18)
         real(8) hdm(18,18)
         real(8) dgmdtrs(6,18)
         real(8) dIVBdtrs(6,18)
         real(8) Lprs(3,3)    
         integer i,j,is,js,ks
         integer iex,igx      
         real(8) x1,x2,x3
c
         crsss(1:3)  = crssbas 
         crsss(4:6)  = crssprs
         crsss(7:18) = crsspys
         hdrt0(1:3)  = hdrtba0
         hdrt0(4:6)  = hdrtpr0 
         hdrt0(7:18) = hdrtpy0
c
         do iex=1,Tnel
         do igx=1,Tngp

            dgmdtrs=0
            do is=1,6
            do js=1,18
               x1=dot_product(fem_pk2i(iex,igx,:),smdVi2_rs(is,js,:))
               x2=x1/fem_IVBrs(iex,igx,is,js)
               dgmdtrs(is,js)=shrts0*dabs(x2)**pwfl*dsign(1.d0,x2) 
            enddo
            enddo

            dIVBdtrs=0
            do is=1,6
            do js=1,18
               do ks=1,18
                  x1=1-fem_IVBrs(iex,igx,is,ks)/crsss(ks)
                  x1=max(1.d-50,x1)
                  dIVBdtrs(is,js)=dIVBdtrs(is,js)
     &              +hdM_rs(js,ks)*hdrt0(ks)*x1**pwhd
     &              *dabs(dgmdtrs(is,ks))
               enddo
            enddo
            enddo

            Lprs=0
            do is=1,6
            do js=1,18
               Lprs=Lprs+fem_gam(iex,igx,is+18)/gammatw
     &             *dgmdtrs(is,js)*smdMi_rs(is,js,:,:)
            enddo
            enddo

            fem_IVBrs(iex,igx,1:6,1:18)= 
     &      fem_IVBrs(iex,igx,1:6,1:18)+
     &      dIVBdtrs(1:6,1:18)*dzeit1

            fem_Lprs(iex,igx,:,:)=Lprs

         enddo
         enddo
         return
         endsubroutine
      endmodule mod_wkcoup_reslip
      
c================================================================
c
c    weak coupling module: Back stress of slip system
c
c================================================================     
      module mod_wkcoup_backStress
         use mod_FEMFFT_data
         implicit none
         real(8) fem_bastress(Tnel,Tngp,48)
      contains
         !c========================================
         !c initialization of mod_wkcoup_backStress  
         !c========================================
         subroutine sub_bs_ini
         implicit none
         integer i,j,k
         real(8) x1,x2,x3
         fem_bastress=0
         return
         endsubroutine
         !c========================================
         !c evolution of mod_wkcoup_backStress  
         !c========================================
         subroutine sub_bs_evl
         implicit none
         integer iex,igx,is
         real(8) x1,x2,x3 
         real(8) pcALL(24),pdALL(24),pkALL(24) 

         real(8),parameter:: fitParaC1=1.d0
         real(8),parameter:: fitParaC2=1.d0
         real(8),parameter:: fitParaC3=1.d0

         real(8),parameter :: pc = 4200.      
         real(8),parameter :: pd = 50. 
         real(8),parameter :: pk = 7.5  

c         real(8),parameter :: pc = fitParaC1      
c         real(8),parameter :: pd = fitParaC2 
c         real(8),parameter :: pk = fitParaC3  

         do iex=1,Tnel
         do igx=1,Tngp

         do is=1,18
            x1=fem_dgmdt(iex,igx,is)
            x2=fem_bastress(iex,igx,is)
            x3=pc*x1-pd*x2*dabs(x1)*(dabs(x2)/(pc/pd))**pk  
            fem_bastress(iex,igx,is)=fem_bastress(iex,igx,is)+x3*dzeit1
         enddo

         do is=19,24
            x1=fem_dgmdt(iex,igx,is) *1 
            x2=fem_bastress(iex,igx,is)
            x3=pc*x1-pd*x2*dabs(x1)*(dabs(x2)/(pc/pd))**pk  
            fem_bastress(iex,igx,is)=fem_bastress(iex,igx,is)+x3*dzeit1
         enddo

         enddo
         enddo
         return
         endsubroutine
      endmodule mod_wkcoup_backStress
c================================================================
c
c    weak coupling module
c
c================================================================
      module mod_wkcoup
         use mod_wkcoup_gnd
         use mod_wkcoup_reslip
         use mod_wkcoup_backStress
         implicit none
         integer, parameter :: Iwkcoup_gnd=0  !-->with(1), without(0) 
         integer, parameter :: Iwkcoup_rs=0   !-->with(1), without(0) 
         integer, parameter :: Iwkcoup_bs=1   !-->with(1), without(0)
      contains
c================================================================
c        state variable ini
c================================================================
         subroutine sub_wkcoup_ini
            implicit none
            print*,'1==================================1'
            print*,'1   wkcoupMod was initionilized.   1'
            print*,'1==================================1'
            call sub_FEMFFT_ini
            if(Iwkcoup_gnd/=0) call sub_gnd_ini
            if(Iwkcoup_rs /=0) call sub_rs_ini
            if(Iwkcoup_bs /=0) call sub_bs_ini
            return
         endsubroutine
c================================================================
c        state variable evolution
c================================================================
         subroutine sub_wkcoup_evl
            implicit none
c            call sub_FEMFFT_evl
            if(Iwkcoup_gnd/=0) call sub_gnd_evl
            if(Iwkcoup_rs /=0) call sub_rs_evl
            if(Iwkcoup_bs /=0) call sub_bs_evl
            return
         endsubroutine
c
      endmodule mod_wkcoup


