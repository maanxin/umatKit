! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +   subroutine calculate stress, IVB, stiffness                 +
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine calculate_stress_IVB_matJacb( 
     &   ie,ig,ialloy,eang00,Fg0,Fg,dt1,zeit0,zeit,
     &   IVB_gnd,pk2i_gnd,Lprs,bstrs,
     &   pk2i,Fe,Fp,eang,IVB,dgmdt,gam,Nslp,
     &   cs,MatJacb,ising)      

         implicit none
         integer Icurrent_dt,iNRloop
         integer :: Nnr_max=200                            
         real(8) :: toler_NRloop=1.d-5         
         integer :: Iexp_abq=0 
         integer :: Iexp_loc=0 
         integer :: Imth_add=0 
         integer :: first_call=1
         integer Ialloy,Nslp
         real(8) dt0,dt1,zeit0,zeit
!---------------------------------------------------------------c
!        general continuum mechanics variables                  c
!---------------------------------------------------------------c
         integer ie,ig
         real(8) det_Fg,det_Fe,det_Fp
         real(8) QM(3,3)
         real(8) Fg00(3,3),Fg0(3,3),TFg0(3,3),Fg(3,3),TFg(3,3)
         real(8) Fe00(3,3),Fe0(3,3),TFe0(3,3),Fe(3,3),TFe(3,3)
         real(8) Fp00(3,3),Fp0(3,3),TFp0(3,3),Fp(3,3),TFp(3,3)
         real(8) IFg0(3,3),TIFg0(3,3),IFg(3,3),TIFg(3,3)
         real(8) IFe0(3,3),TIFe0(3,3),IFe(3,3),TIFe(3,3)
         real(8) IFp0(3,3),TIFp0(3,3),IFp(3,3),TIFp(3,3)
         real(8) pk2i00(6),pk2i0(6),pk2i(6),pk2i_max(6)
         real(8) pk2r(6),cs(6),dpk2i(6)
         real(8) pk2r_M(3,3),pk2i_M(3,3),cs_M(3,3)
         real(8) Lp(3,3),Lptw(3,3), cs0(6)  !* Lptw is the contribution of the twin part
         real(8) DPv(6),pegrd_1st(6,3),pegrd_2st(6,3,3)
         real(8) dFpdt(3,3)
         real(8) MatJacb0(6,6),MatJacb(6,6)
         real(8) Cstr_max(3,3),Fem(3,3),Femt(3,3)
         real(8) Estr_max(3,3),VEstr_max(6)
         real(8) dpk2i_dt(6),Rvpk2i(6)
         real(8) dRvpk2i_dpk2i(6,6),IdRvpk2i_dpk2i(6,6)
         real(8) STFjc_66(6,6),STFtk_66(6,6)
         real(8) CGE(3,3),CGEe_max(3,3)
         real(8) eang00(4),eang0(4),eang(4),Lg_glb(3,3)
         real(8) Dmc(3,3),Dvc(6),Wmc(3,3),Lg(3,3)
         real(8) dDmc(3,3),dDvc(6),dWmc(3,3),dLg(3,3)
         real(8) smdMc(48,3,3),smdMi(48,3,3)
         real(8) smdSMc(48,3,3),smdAMc(48,3,3)
         real(8) smdSMi(48,3,3),smdAMi(48,3,3)
         real(8) smdVc1(48,6),smdVc2(48,6)
         real(8) smdVi1(48,6),smdVi2(48,6)
         real(8) vd_slp(48,3)
         real(8) vl_slp(48,3)
         real(8) vn_slp(48,3)
         real(8) STFrlx6(48,6),STFrlx33(3,3,48)
         real(8) dSTFrlx6_dE(48,6,6)
         real(8) STFec26(6,6),STFec29(9,9),STFec43(3,3,3,3)
         real(8) STFei26(6,6),STFei29(9,9),STFei43(3,3,3,3)
         real(8) dcsdt(6),dcsdt_max(6),dcs_max(6)
         real(8) csM0(3,3),csM(3,3),dcs(6)
!---------------------------------------------------------------c
!        microstructure(slip system based) relative variables   c
!---------------------------------------------------------------c
         real(8) shgrd00(48,13),shgrd0(48,13)
         real(8) shgrd(48,13),dIVB(48)
         real(8) gam00(48),gam0(48),gam(48),bstrs(48),gamvf(48)   
         real(8) IVB00(48),IVB0(48),IVB(48),hdM(48,48)
         real(8) dIVBdt(48)
         real(8) shrt_gnd(48)
         real(8) tau(48)
         real(8) dgmdt(48),detGam(48)
         real(8) vd_ltc(48,3)
         real(8) vl_ltc(48,3)
         real(8) vn_ltc(48,3)
         real(8) IVB_ini(48),refv_IVB,refv_pk2i
         real(8) pk2i_eq(6)
         real(8) tau_eff(48)
         real(8) IVB_eff(48)
!---------------------------------------------------------------c
!        strain gradient effect                                 c
!---------------------------------------------------------------c
         real(8) Rho_gnd(48)
         real(8) IVB_gnd(48)
         real(8) pk2i_gnd(6),tau_gnd(48)

!---------------------------------------------------------------c
!        misfit stress effect                                   c
!---------------------------------------------------------------c
         real(8) pk2i_msf(48,6),tau_msf(48)

!---------------------------------------------------------------c
!        phase transformation effect                            c
!---------------------------------------------------------------c
         real(8) Ftrp(3,3),IFtrp(3,3)
         real(8) IVB_trp(48)
         real(8) Vtrp(12),Vtrp_e(12),Vtrp_s(12)
!---------------------------------------------------------------c
!        non-Schmid effect                                      c
!---------------------------------------------------------------c
         real(8) tau_nsd(48)
!---------------------------------------------------------------c
!        reslips effect                                      c
!---------------------------------------------------------------c         
         real(8) Lprs(3,3)

!---------------------------------------------------------------c
!        vaules for newton-raphson algoriths                    c
!---------------------------------------------------------------c
         real(8) ddgmdt_dtau(48)
         real(8) ddgmdt_dIVB(48)
         real(8) ddIVBdt_ddgmdt(48,48)
         real(8) ddIVBdt_dIVB(48,48)
         real(8) ddgmdt_dpk2i(48,6)
         real(8) ddIVBdt_dpk2i(48,6)
         real(8) GV1(6),GV2(48)
         real(8) dGv1_dpk2i(6,6),dGv1_dIVB( 6,48)
         real(8) dGv2_dpk2i(48,6),dGv2_dIVB(48,48)
         real(8) IdGv1_dpk2i(6,6),IdGv2_dIVB(48,48)
         real(8) eqM66Gv1(6,6),eqMnnGv2(48,48)
         real(8) IeqM66Gv1(6,6),IeqMnnGv2(48,48)
         real(8) eqM6nGv1(6,48),eqMn6Gv2(48,6)
!---------------------------------------------------------------c
!        vaules for material tangent calcuation                 c
!---------------------------------------------------------------c
         real(8) dIVB_dpk2i(48,6)
         real(8) dTdgmdt_dpk2i(48,9)
         real(8) dGv1_dE(6,6)
         real(8) dpk2i_dE(6,6),dPk2r_dE(6,6)
         real(8) dIFp_dE(9,6),dTIFp_dE(9,6)
         real(8) dSTFsh_dE(48,6,6)
         real(8) dCGEe_mx_dE(6,6)
         real(8) dpk2i_mx_dE(6,6)
         real(8) dFp_dE(9,6),ddgmdt_dE(48,6)
         integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
         integer ising_1,ising
         real(8) x1,x2,x3,x4,y1,z1
         real(8) M1_66(6,6),M1_96(9,6),M1_99(9,9)
         real(8) M2_66(6,6),M2_96(9,6),M2_99(9,9),M3_99(9,9)
         real(8) M1_3333(3,3,3,3),M2_3333(3,3,3,3)
         real(8) MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
         integer IB1(9),IB2(9)
         real(8) XI33(3,3),XI66(6,6),XI99(9,9),XInn(48,48)
         real(8) XI333(3,3,3)
         !1-----------------------------1
         !1   Initialization            1
         !1-----------------------------1
         ising=0
         call constantsIDT(IB1,IB2,XI33,XI333,XI66,XI99,XInn)

         if(Ialloy==1)then
            call sub_Aluminum_ini(STFei26,Nslp,
     &           smdMi,smdSMi,smdAMi,smdVi1,smdVi2,hdM,
     &           refv_pk2i,refv_IVB,IVB_ini)
         elseif(Ialloy==2)then
            call sub_Ferrite_ini(STFei26,Nslp,
     &           smdMi,smdSMi,smdAMi,smdVi1,smdVi2,hdM,
     &           refv_pk2i,refv_IVB,IVB_ini)
         elseif(Ialloy==3)then
            call sub_Magnesium_ini(STFei26,Nslp,
     &           smdMi,smdSMi,smdAMi,smdVi1,smdVi2,hdM,
     &           refv_pk2i,refv_IVB,IVB_ini)
         else
            print*,'no this material:', Ialloy
            call xit
         endif
         
         if(zeit==0)then
            call icams_eang2Q(eang00(1),eang00(2),eang00(3),QM)
            Fe   = transpose(QM)
            Fp   = QM
            IVB  = IVB_ini              
            gam  = 0.d0              
            eang = eang00
         endif
         pk2i0 = pk2i
         cs0   = cs
         Fp0   = Fp
         Fe0   = Fe
         IVB0  = IVB
         eang0 = eang
         gam0  = gam
         do i=1,6; csM0(ib1(i),ib2(i))=cs0(i  ); enddo
         do i=7,9; csM0(ib1(i),ib2(i))=cs0(i-3); enddo
c
         IFtrp=XI33
c
         !1-----------------------------1
         !1   check strong distortion   1
         !1-----------------------------1
         call icams_determ(Fg,det_Fg)
         if(det_Fg < 1.d-10)then
            ising=1 
            pk2i=pk2i0
            cs=cs0
            Fp=Fp0
            Fe=Fe0
            IVB=IVB0
            eang=eang0
            gam=gam0
            write(6,*) 'Strongly distorted', ie,ig,dt1
            write(*,*) 'Strongly distorted', ie,ig,dt1
            return
         endif

         ising_1=0
         call gaussj(Fg,3,IFg,ising_1)
         if(ising_1/=0)then
            ising=2
            pk2i=pk2i0
            cs=cs0
            Fp=Fp0
            Fe=Fe0
            IVB=IVB0
            eang=eang0
            gam=gam0
            write(*,*) 'Fg is non-ivertable', ie,ig,dt1
            write(6,*) 'Fg is non-ivertable', ie,ig,dt1
            return
         endif

         det_Fe=det_Fg
         TFg=transpose(Fg)
         CGE=matmul(TFg,Fg)
         !1--------------------------------1
         !1   Whether Fp0 is ivertible     1
         !1--------------------------------1
         ising_1=0
         call gaussj(Fp0,3,IFp0,ising_1)
         if(ising_1/=0)then
            ising=3
            pk2i=pk2i0
            cs=cs0
            Fp=Fp0
            Fe=Fe0
            IVB=IVB0
            eang=eang0
            gam=gam0
            write(6,*) 'non ivertable for Fp0',ie,ig,dt1
            write(*,*) 'non ivertable for Fp0',ie,ig,dt1
            return
         endif
         TIFp0=transpose(IFp0)
         !1-------------------------------------------------1
         !1   Calculate pk2i_max(6)=C*(CGE_max(6)-I(6))/2   1
         !1-------------------------------------------------1
         MX1=matmul( matmul(TIFp0,CGE),IFp0 )                    
         CGEe_max=matmul( matmul(transpose(IFtrp),MX1), IFtrp )  !==> add trps effect
!         CGEe_max=matmul( matmul(transpose(IFprs),MX1), IFprs ) 

         !** the following changed by meijuan 10/25/2018         !==> add reslips effect
         MX1=(CGEe_max-XI33)/2
         MX1 = MX1 - matmul(CGEe_max,Lprs)*dt1/2
     &       - matmul(transpose(Lprs),CGEe_max)*dt1/2
     
         pk2i_max=0
         do i=1,6
         do j=1,6
            pk2i_max(i)=pk2i_max(i)
     &      +STFei26(i,j)*MX1(ib1(j),ib2(j))
         enddo
         enddo

         !1------------------------------------------------------------------1
         !1   Calculate STFrlx6=C*( CGE_max*smdMi + (CGE_max*smdMi)^T )/2    1
         !1------------------------------------------------------------------1
         do is=1,Nslp
            MX1=matmul(CGEe_max,smdMi(is,:,:))
            MX2=matmul( matmul(transpose(IFtrp), !==> add trp effect
     &            (MX1+transpose(MX1))/2), IFtrp )
!            MX2=matmul( matmul(transpose(IFprs), !==> add reslips effect
!     &            (MX1+transpose(MX1))/2), IFprs )
            do i=1,6
               STFrlx6(is,i)=0
               do j=1,6
                  STFrlx6(is,i)=STFrlx6(is,i)
     &            +STFei26(i,j)*MX2(ib1(j),ib2(j))
               enddo
            enddo
         enddo
!======================================================!
!                                                      !
!    Begin newton raphson method to solve pk2i, IVB    !
!                                                      !
!======================================================!
         do iNRloop=1,Nnr_max
            do is=1,Nslp
               tau(is)=dot_product(pk2i,smdVi2(is,:))
            enddo               
            do is=1,Nslp
               tau_gnd(is)=dot_product(pk2i_gnd,smdVi2(is,:))
            enddo               

            tau_eff = tau + tau_gnd - bstrs
            IVB_eff = IVB + IVB_gnd             

            !1----------------------------------------------1
            !1   get material flow and hardening equations  1
            !1----------------------------------------------1
            ising_1=0

            if(Ialloy==1)then 
               call sub_flhd_Aluminum(Iexp_loc,gam,
     &              tau,hdM,IVB,tau_eff,IVB_eff,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,
     &              ising_1)
            elseif(Ialloy==2)then 
               call sub_flhd_Ferrite(Iexp_loc,gam,
     &              tau,hdM,IVB,tau_eff,IVB_eff,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,
     &              ising_1)
            elseif(Ialloy==3)then 
               call sub_flhd_Magnesium(Iexp_loc,gam,
     &              tau,hdM,IVB,tau_eff,IVB_eff,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,
     &              ising_1)
            else
               print*,'no this material:', Ialloy
               call xit
            endif
 
            if(ising_1/=0)then
               ising=4
               pk2i=pk2i0
               cs=cs0
               Fp=Fp0
               Fe=Fe0
               IVB=IVB0
               eang=eang0
               gam=gam0
               write(6,*) 'Flow-Harden has problem',ie,ig,dt1
               write(*,*) 'Flow-Harden has problem',ie,ig,dt1,Nslp
               return
            endif

            do is=1,Nslp
               ddgmdt_dpk2i(is,:)=ddgmdt_dtau(is)*smdVi2(is,:)
            enddo
            ddIVBdt_dpk2i(1:Nslp,:)=
     &      matmul(ddIVBdt_ddgmdt(1:Nslp,1:Nslp), 
     &      ddgmdt_dpk2i(1:Nslp,:))

            Gv1=+pk2i-pk2i_max
     &      +matmul(transpose(STFrlx6(1:Nslp,:)),dgmdt(1:Nslp))*dt1     
            Gv2(1:Nslp)=+IVB(1:Nslp)-IVB0(1:Nslp)-dIVBdt(1:Nslp)*dt1   
            x1=sum(dabs(Gv1))/refv_pk2i
            x2=sum(dabs(Gv2(1:Nslp)))/refv_IVB
                                                                      
       if(ie==-1)then
          call pv(pk2i,6)
          call pv(IVB(1:12),12)
          print*,iNRloop,x1+x2,toler_NRloop
       endif
         
            if((x1+x2<toler_NRloop.and.iNRloop>1).or.Iexp_loc==1)then
               pk2i=+pk2i_max-matmul(transpose(STFrlx6(1:Nslp,:)),
     &              dgmdt(1:Nslp))*dt1
               IVB(1:Nslp)=IVB0(1:Nslp)+dIVBdt(1:Nslp)*dt1 
               goto 101 !Converge jump out NR loop and continue
            endif

            !-------------------------------------
            dGv1_dpk2i=XI66+matmul(transpose(STFrlx6(1:Nslp,:)),
     &                 ddgmdt_dpk2i(1:Nslp,:))*dt1
            ising_1=0
            call gaussj(dGv1_dpk2i,6,idGv1_dpk2i,ising_1)
            if(ising_1/=0)then
               ising=5
               pk2i=pk2i0
               cs=cs0
               Fp=Fp0
               Fe=Fe0
               IVB=IVB0
               eang=eang0
               gam=gam0
               write(6,*) 'dGv1_dpk2i is ivertable',ie,ig,dt1
               write(*,*) 'dGv1_dpk2i is ivertable',ie,ig,dt1
               return
            endif

            do is=1,Nslp
               dGv1_dIVB(:,is)=STFrlx6(is,:)*ddgmdt_dIVB(is)*dt1
            enddo
            dGv2_dpk2i(1:Nslp,:)=-ddIVBdt_dpk2i(1:Nslp,:)*dt1

            !-------------------------------------
            dGv2_dIVB(1:Nslp,1:Nslp) =+XInn(1:Nslp,1:Nslp)
     &                  -ddIVBdt_dIVB(1:Nslp,1:Nslp)*dt1

            ising_1=0
            call gaussj( dGv2_dIVB(1:Nslp,1:Nslp),Nslp,
     &                  idGv2_dIVB(1:Nslp,1:Nslp),ising_1)
            if(ising_1/=0)then

               call pv(pk2i,6)
               call pv(tau(1:Nslp),Nslp)
               call pv(IVB0(1:Nslp),Nslp)
               call pv(IVB(1:Nslp),Nslp)
               call pm(XInn(1:Nslp,1:Nslp),Nslp,Nslp)
               call pm(ddIVBdt_dIVB(1:Nslp,1:Nslp),Nslp,Nslp)
               call pm(dGv2_dIVB(1:Nslp,1:Nslp),Nslp,Nslp)

               ising=6
               pk2i=pk2i0
               cs=cs0
               Fp=Fp0
               Fe=Fe0
               IVB=IVB0
               eang=eang0
               gam=gam0
               write(6,*) 'dGv2_dIVB is ivertable',ie,ig,dt1
               write(*,*) 'dGv2_dIVB is ivertable',ie,ig,dt1
               return
            endif
            !-------------------------------------
            eqM6nGv1=matmul(dGv1_dIVB(:,1:Nslp),
     &                     IdGv2_dIVB(1:Nslp,1:Nslp))
            eqM66Gv1=-matmul(matmul(dGv1_dIVB(:,1:Nslp),
     &         IdGv2_dIVB(1:Nslp,1:Nslp)),dGv2_dpk2i(1:Nslp,:))
     &         +dGv1_dpk2i
            ising_1=0
            call gaussj(eqM66Gv1,6,IeqM66Gv1,ising_1)
            if(ising_1/=0)then
               ising=7
               pk2i=pk2i0
               cs=cs0
               Fp=Fp0
               Fe=Fe0
               IVB=IVB0
               eang=eang0
               gam=gam0
               write(6,*) 'eqM66Gv1 is ivertable',ie,ig,dt1
               write(*,*) 'eqM66Gv1 is ivertable',ie,ig,dt1
               return
            endif
            !-------------------------------------
            eqMn6Gv2(1:Nslp,:)=matmul(dGv2_dpk2i(1:Nslp,:),IdGv1_dpk2i)
            eqMnnGv2(1:Nslp,1:Nslp)=-matmul(matmul(dGv2_dpk2i(1:Nslp,:),
     &                          IdGv1_dpk2i),dGv1_dIVB(:,1:Nslp))
     &                         +dGv2_dIVB(1:Nslp,1:Nslp)
            ising_1=0
            call gaussj( eqMnnGv2(1:Nslp,1:Nslp),Nslp,
     &                  IeqMnnGv2(1:Nslp,1:Nslp),ising_1)
            if(ising_1/=0)then
               ising=8
               pk2i=pk2i0
               cs=cs0
               Fp=Fp0
               Fe=Fe0
               IVB=IVB0
               eang=eang0
               gam=gam0
               write(6,*) 'eqMnnGv2 is ivertable',ie,ig,dt1
               write(*,*) 'eqMnnGv2 is ivertable',ie,ig,dt1
               return
            endif
                                              
            !-------------------------------------
            dpk2i=matmul(IeqM66Gv1,
     &         -(Gv1-matmul(eqM6nGv1(:,1:Nslp),Gv2(1:Nslp))))
            dIVB(1:Nslp)=matmul( IeqMnnGv2(1:Nslp,1:Nslp),
     &         -(Gv2(1:Nslp)-matmul(eqMn6Gv2(1:Nslp,:),Gv1)) )
            pk2i=pk2i+dpk2i
            IVB(1:Nslp) =IVB(1:Nslp) + dIVB(1:Nslp) 
         enddo

         ising=9
         pk2i=pk2i0
         cs=cs0
         Fp=Fp0
         Fe=Fe0
         IVB=IVB0
         eang=eang0
         gam=gam0
         write(6,*) 'NR-algorithm is not converge',ie,ig,dt1
         write(*,*) 'NR-algorithm is not converge',ie,ig,dt1
         return

101      continue
!======================================================!
!                                                      !
!    End of Newton Raphson method to solve pk2i, IVB   !
!                                                      !
!======================================================!
         call icams_conv6to33(pk2i,ib1,ib2,pk2i_M)
         do i=1,3
         do j=1,3
            Lp(i,j)=dot_product(dgmdt(1:Nslp),smdMi(1:Nslp,i,j))
         enddo
         enddo
         Lp = Lp + Lprs     

         Fp=matmul(XI33+Lp*dt1,Fp0)                    
         call icams_determ(Fp,det_Fp)
         if(det_Fp<1.d-10)then
            ising=10
            pk2i=pk2i0
            cs=cs0
            Fp=Fp0
            Fe=Fe0
            IVB=IVB0
            eang=eang0
            gam=gam0
            write(6,*) 'det of Fp is almost zero',ie,ig,dt1
            write(*,*) 'det of Fp is almost zero',ie,ig,dt1
            call pm(Lp,3,3)
            return
         endif
         Fp=Fp/det_Fp**(1/3.0)
         det_Fp=1    ! unnecessary

         dFpdt=(Fp-Fp0)/dt1

         ising_1=0
         call gaussj(Fp,3,iFp,ising_1)
         if(ising_1/=0)then
            ising=11
            pk2i=pk2i0
            cs=cs0
            Fp=Fp0
            Fe=Fe0
            IVB=IVB0
            eang=eang0
            gam=gam0
            write(6,*) 'Fp is non ivertable', ie,ig,dt1
            write(*,*) 'Fp is non ivertable', ie,ig,dt1
            return
         endif                                          

         TIFp=transpose(iFp)
         Fe=matmul(Fg,iFp)
         call icams_determ(Fe,det_Fe)
         csM=matmul(matmul(Fe,pk2i_M),transpose(Fe))/det_Fe
         call icams_conv33to6(csM,ib1,ib2,cs)

         pk2r_M=matmul(matmul(IFg,csM),transpose(IFg))*det_Fe
         call icams_conv33to6(pk2r_M,ib1,ib2,pk2r)                      

         detGam=dgmdt*dt1            
         if(ialloy==3)then
            gam( 1:18)=gam( 1:18)+dabs(detGam(1:18))               
            gam(19:24)=gam(19:24)+detGam(19:24)               
         else
            gam=gam+dabs(detGam)               
         endif
         
         call caleulang(Fe,eang(1:3),ising_1)
         if(ising_1/=0) eang(1:3)=eang0(1:3)
         call icams_misori(eang00(1:3),eang(1:3),eang(4))

c=============================================================================
c
c   calculate material stiffness for user material and user element  
c
c           (STF_JC_3333)_ijkl = + (dpk2r_dE)_mnop*F_im*F_jn*F_ko*F_lp/det(F)
c                                +  I_ik * CS_lj
c                                + CS_ik *  I_lj
c                                - CS_ij *  I_kl
c        MatJacb=(STFjc_66+transpose(STFjc_66))/2
c
c   calculate material stiffness for user material and user element  
c
c           (STF_TK_3333)_ijkl = + (dpk2r_dE)_mnop*F_im*F_jn*F_ko*F_lp
c           MatJacb=(STF_TK_66+transpose(STF_TK_66))/2/det(F)
c
c=============================================================================
         dIVB_dpk2i=-matmul(IdGv2_dIVB,dGv2_dpk2i)
         do i=1,6
         do j=1,6
            if(j<=3)then
               dCGEe_mx_dE(i,j)=2*TIFp0(ib1(i),ib1(j))
     &                         *IFp0(ib2(j),ib2(i))
            else
               dCGEe_mx_dE(i,j)=
     &         (+2*TIFp0(ib1(i),ib1(j  ))*IFp0(ib2(j  ),ib2(i))
     &          +2*TIFp0(ib1(i),ib1(j+3))*IFp0(ib2(j+3),ib2(i)))/2
            endif
         enddo
         enddo
         !---------------------------------------------
         dpk2i_mx_dE=matmul( STFei26 , dCGEe_mx_dE )/2

         do is=1,Nslp
            MX1=matmul(IFp0,smdMi(is,:,:))
            MX2=transpose(MX1)
            do i=1,6
            do j=1,6
               M1_66(i,j)=
     &           +2*TIFp0(ib1(i),ib1(j))*MX1 (ib2(j),ib2(j))
     &           +2*MX2  (ib1(i),ib1(j))*IFp0(ib2(j),ib2(j))
            enddo
            enddo
            dSTFrlx6_dE(is,:,:)=matmul(STFei26 , M1_66)/2
         enddo
         !---------------------------------------------
         do i=1,6
         do j=1,6
            dGv1_dE(i,j)=
     &      -dpk2i_mx_dE(i,j)
     &      +dot_product(dgmdt(1:Nslp),dSTFrlx6_dE(1:Nslp,i,j))*dt1
         enddo
         enddo

         !---------------------------------------------
         dpk2i_dE=-matmul(IeqM66Gv1,dGv1_dE)
         !---------------------------------------------

c       if(ie==1)then
c          call pm(IeqM66Gv1,6,6)
c          call pm(dGv1_dE,6,6)
c          call pm(dpk2i_dE,6,6)
c       endif


         M1_99=0
         M2_99=0
         M3_99=0
         do is=1,Nslp
            MX1=matmul(IFp0,smdMi(is,:,:))
            MX2=transpose(MX1)
            MX4=matmul(smdMi(is,:,:),Fp0)
            do i=1,9
               if(i<=6) i1=i
               if(i >6) i1=i-3
               MX3(ib1(i),ib2(i))=
     &         +ddgmdt_dpk2i(is,i1)
     &         +ddgmdt_dIVB(is)*dIVB_dpk2i(is,i1)
               dTdgmdt_dpk2i(is,i)=+ddgmdt_dpk2i(is,i1)
     &         +ddgmdt_dIVB(is)*dIVB_dpk2i(is,i1)
            enddo
            do i=1,9
            do j=1,9
            M1_99(i,j)=M1_99(i,j)+MX1(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
            M2_99(i,j)=M2_99(i,j)+MX2(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
            M3_99(i,j)=M3_99(i,j)+MX4(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
            enddo
            enddo
         enddo
         dIFp_dE=0
         dTIFp_dE=0
         M1_96(1:6,:)=dpk2i_dE
         M1_96(7,:)=M1_96(4,:)
         M1_96(8,:)=M1_96(5,:)
         M1_96(9,:)=M1_96(6,:)
          dIFp_dE=-matmul(M1_99,M1_96)*dt1
         dTIFp_dE=-matmul(M2_99,M1_96)*dt1     
           dFp_dE=-matmul(M3_99,M1_96)*dt1
         ddgmdt_dE=matmul(dTdgmdt_dpk2i,M1_96)*dt1
         dFp_dE(:,4:6)=dFp_dE(:,4:6)*2
         ddgmdt_dE(:,4:6)=ddgmdt_dE(:,4:6)*2

         !---------------------------------------------
         dpk2r_dE=0
         MX1=matmul(pk2i_M,TIFp)
         MX2=matmul(IFp,pk2i_M)
         do i=1,6
         do k=1,6
            do m=1,9
               if(m<=6)m1=m
               if(m >6)m1=m-3
               dpk2r_dE(i,k)=dpk2r_dE(i,k)
     &        +XI33(ib1(i),ib1(m))* MX1(ib2(m),ib2(i))* dIFp_dE(m ,k)
     &        + IFp(ib1(i),ib1(m))*TIFp(ib2(m),ib2(i))*dpk2i_dE(m1,k)
     &        + MX2(ib1(i),ib1(m))*XI33(ib2(i),ib2(m))*dTIFp_dE(m ,k)
            enddo
         enddo
         enddo

         
         !-------------------------------------------------------------------------------------------
         ! Stiffness for user material:  use Jaummann rate of Cauchy stress STFjc_66 
         !-------------------------------------------------------------------------------------------
         do i=1,9
         do j=1,9
            if(i<=6) i1=i
            if(i >6) i1=i-3
            if(j<=6) j1=j
            if(j >6) j1=j-3
            M1_3333(ib1(i),ib2(i),ib1(j),ib2(j))=dpk2r_dE(i1,j1)
         enddo
         enddo
         M2_3333=0
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
            x1=0
            do i1=1,3
            do j1=1,3
            do k1=1,3
            do l1=1,3
               x1=x1+M1_3333(i1,j1,k1,l1)
     &        *Fg(i,i1)*Fg(j,j1)*Fg(k,k1)*Fg(l,l1)
            enddo
            enddo
            enddo
            enddo
            M2_3333(i,j,k,l)=x1/det_Fg
     &                   +XI33(i,k)*csM0(l,j)
     &                   +csM0(i,k)*XI33(l,j)
     &                   +csM0(i,j)*XI33(k,l)
         enddo
         enddo
         enddo
         enddo
         do i=1,6
         do j=1,6
            STFjc_66(I,J)=M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
         enddo
         enddo
         MatJacb=STFjc_66

c         !c-----------------------------------------------------------------------------------------
c         !c Stiffness for for user element: use Trusdell rate of Kirchhoff stress STFtk_66 
c         !c-----------------------------------------------------------------------------------------
c         do i=1,9
c         do j=1,9
c            if(i<=6) i1=i
c            if(i >6) i1=i-3
c            if(j<=6) j1=j
c            if(j >6) j1=j-3
c            M1_3333(ib1(i),ib2(i),ib1(j),ib2(j))=dpk2r_dE(i1,j1)
c         enddo
c         enddo
c         M2_3333=0
c         do i=1,3
c         do j=1,3
c         do k=1,3
c         do l=1,3
c            do k1=1,3
c            do l1=1,3
c               M2_3333(i,j,k,l)=M2_3333(i,j,k,l)
c     &        +M1_3333(i,j,k1,l1)*Fg(k,k1)*Fg(l,l1)
c            enddo
c            enddo
c         enddo
c         enddo
c         enddo
c         enddo
c         do i=1,6
c         do j=1,6
c            STFtk_66(I,J)=M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
c         enddo
c         enddo
c
         return
         end

