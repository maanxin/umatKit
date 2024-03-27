
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #1==>Aluminum                                       +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine sub_Aluminum_ini(STFei26,Nslp,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,hdM,
     &              refv_pk2i,refv_IVB,IVB_ini)
            implicit none
            real(8),parameter:: c11    = 107.3d3  !
            real(8),parameter:: c12    =  60.9d3  ! Hosford Book,P16
c            real(8),parameter:: c44    =  28.3d3  !
            real(8),parameter:: c44    =  (c11-c12)/2
            real(8),parameter:: c_cpl  = 1.d0
            real(8),parameter:: c_oth  = 1.4d0
            real(8),parameter:: crss0  = 16.d0
            real(8) hdM(48,48)
            integer Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(48,3,3)
            real(8) smdSMi(48,3,3)
            real(8) smdAMi(48,3,3)
            real(8) smdVi1(48,6)
            real(8) smdVi2(48,6)
            real(8) vd(48,3)
            real(8) vl(48,3)
            real(8) vn(48,3)
            real(8) refv_pk2i
            real(8) refv_IVB
            real(8) IVB_ini(48)
            integer i,j,is,js
            real(8) x1,x2,x3
            integer IB1(9),IB2(9)
            data IB1/1,2,3,1,1,2,2,3,3/
            data IB2/1,2,3,2,3,3,1,1,2/
c
            STFei26(:,:)=0
            STFei26(1,1)=c11
            STFei26(2,2)=c11
            STFei26(3,3)=c11
            STFei26(4,4)=c44*2
            STFei26(5,5)=c44*2
            STFei26(6,6)=c44*2
            STFei26(2,3)=c12
            STFei26(3,2)=c12
            STFei26(1,3)=c12
            STFei26(3,1)=c12
            STFei26(1,2)=c12
            STFei26(2,1)=c12
c
            Nslp=12
            vd( 1,:)=[ 0,  1, -1] ; vn( 1,:)=[  1,  1,  1]
            vd( 2,:)=[ 1,  0, -1] ; vn( 2,:)=[  1,  1,  1]
            vd( 3,:)=[ 1, -1,  0] ; vn( 3,:)=[  1,  1,  1]
            vd( 4,:)=[ 0,  1, -1] ; vn( 4,:)=[ -1,  1,  1]
            vd( 5,:)=[ 1,  0,  1] ; vn( 5,:)=[ -1,  1,  1]
            vd( 6,:)=[ 1,  1,  0] ; vn( 6,:)=[ -1,  1,  1]
            vd( 7,:)=[ 0,  1,  1] ; vn( 7,:)=[  1, -1,  1]
            vd( 8,:)=[ 1,  0, -1] ; vn( 8,:)=[  1, -1,  1]
            vd( 9,:)=[ 1,  1,  0] ; vn( 9,:)=[  1, -1,  1]
            vd(10,:)=[ 0,  1,  1] ; vn(10,:)=[  1,  1, -1]
            vd(11,:)=[ 1,  0,  1] ; vn(11,:)=[  1,  1, -1]
            vd(12,:)=[ 1, -1,  0] ; vn(12,:)=[  1,  1, -1]
            do is=1,12
               vl(is,1)=vn(is,2)*vd(is,3)-vn(is,3)*vd(is,2)
               vl(is,2)=vn(is,3)*vd(is,1)-vn(is,1)*vd(is,3)
               vl(is,3)=vn(is,1)*vd(is,2)-vn(is,2)*vd(is,1)
               x1=dsqrt(sum(vn(is,:)**2))
               x2=dsqrt(sum(vd(is,:)**2))
               x3=dsqrt(sum(vl(is,:)**2))
               vn(is,:)=vn(is,:)/x1
               vd(is,:)=vd(is,:)/x2
               vl(is,:)=vl(is,:)/x3
               do i=1,3
               do j=1,3
                  smdMi(is,i,j)=vd(is,i)*vn(is,j)
               enddo
               enddo
               smdSMi(is,:,:)=( smdMi(is,:,:)
     &                       +transpose(smdMi(is,:,:)) )/2
               smdAMi(is,:,:)=( smdMi(is,:,:)
     &                       -transpose(smdMi(is,:,:)) )/2
               do i=1,6
                  smdVi1(is,i)=smdSMi(is,ib1(i),ib2(i))
               enddo
               smdVi2(is,1:3)=smdVi1(is,1:3)
               smdVi2(is,4:6)=smdVi1(is,4:6)*2
            enddo
c
            do is=1,12
            do js=1,12
               x1=sum(dabs(vn(is,:)-vn(js,:)))
               if(x1<1.d-10)then
                  hdM(is,js)=c_cpl
               else
                  hdM(is,js)=c_oth
               endif
            enddo
            enddo
c
            refv_pk2i = c44*1.d-6
            refv_IVB  = c44*1.d-6
            IVB_ini   = crss0 
c
            return
         endsubroutine
c=================================================================
         subroutine sub_flhd_Aluminum(Iexp_loc,gam,
     &                 tau,hdM,IVB,tau_eff,IVB_eff,
     &                 dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &                 dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,
     &                 ising)
            implicit none

c            real(8),parameter:: shrt0  = 0.001    !
c            real(8),parameter:: pwfl   = 83.33    ! Kalidindi
c            real(8),parameter:: pwhd   = 2.25     ! JMPS1992
c            real(8),parameter:: crsss  = 148.0    ! Vol40, No3
c            real(8),parameter:: hdrt0  = 180.d0   ! pp 537-569

            real(8),parameter:: shrt0  = 0.001    !
c            real(8),parameter:: pwfl   = 20       ! Kalidindi
            real(8),parameter:: pwfl   = 10       ! Kalidindi
            real(8),parameter:: pwhd   = 2.25     ! JMPS1992
            real(8),parameter:: crsss  = 148.0    ! Vol40, No3
            real(8),parameter:: hdrt0  = 180.d0   ! pp 537-569

            real(8) hdM(48,48)
            integer i,j,is,js,ising
            integer Iexp_loc               
            real(8) gam(48)                        
            real(8) tau(48),IVB(48)
            real(8) tau_eff(48),IVB_eff(48)
            real(8) dgmdt(48)
            real(8) ddgmdt_dtau(48)
            real(8) ddgmdt_dIVB(48)
            real(8) dIVBdt(48)
            real(8) ddIVBdt_ddgmdt(48,48)
            real(8) ddIVBdt_dIVB(48,48)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,12
               x1=crsss*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
c               print*,Is, IVB(is), x1,x2
c               read*
                  return
               endif
            enddo

c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB
            do is=1,12
               dgmdt(is)=shrt0*(dabs(tau_eff(is))/IVB_eff(is))
     &                          **pwfl*dsign(1.d0,tau_eff(is))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0
     &            *(dabs(tau_eff(is))/IVB_eff(is))**(pwfl-1)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo
c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do is=1,12
            do js=1,12
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + hdM(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=hdM(is,js)*hdrt0
     &            *dsign(1.d0,tau_eff(js))*x1**pwhd 

                  ddIVBdt_dIVB(is,js)=hdM(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau_eff(js))*x1**pwhd 
     &            -hdM(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss

               endif
            enddo
            enddo

            return
         end


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #2==>Ferrite                                        +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine sub_Ferrite_ini(STFei26,Nslp,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,hdM,
     &              refv_pk2i,refv_IVB,IVB_ini)
            implicit none
            real(8),parameter:: c11    = 231.0d3
            real(8),parameter:: c12    = 134.7d3
            real(8),parameter:: c44    = 116.4d3 ! c44=2*(c11-c12) iso
            real(8),parameter:: crss0  = 40.d0
            real(8),parameter:: c_cpl  = 1.d0
            real(8),parameter:: c_oth  = 1.4d0
            real(8) hdm(48,48)
            integer Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(48,3,3)
            real(8) smdSMi(48,3,3)
            real(8) smdAMi(48,3,3)
            real(8) smdVi1(48,6)
            real(8) smdVi2(48,6)
            real(8) vd(48,3)
            real(8) vl(48,3)
            real(8) vn(48,3)
            real(8) refv_pk2i
            real(8) refv_IVB
            real(8) IVB_ini(48)
            integer i,j,is,js
            real(8) x1,x2,x3
            integer IB1(9),IB2(9)
            data IB1/1,2,3,1,1,2,2,3,3/
            data IB2/1,2,3,2,3,3,1,1,2/
c
            STFei26(:,:)=0
            STFei26(1,1)=c11
            STFei26(2,2)=c11
            STFei26(3,3)=c11
            STFei26(4,4)=c44*2
            STFei26(5,5)=c44*2
            STFei26(6,6)=c44*2
            STFei26(2,3)=c12
            STFei26(3,2)=c12
            STFei26(1,3)=c12
            STFei26(3,1)=c12
            STFei26(1,2)=c12
            STFei26(2,1)=c12
c
            Nslp=12
            vn( 1,:)=[ 0,  1, -1] ; vd( 1,:)=[  1,  1,  1]
            vn( 2,:)=[ 1,  0, -1] ; vd( 2,:)=[  1,  1,  1]
            vn( 3,:)=[ 1, -1,  0] ; vd( 3,:)=[  1,  1,  1]
            vn( 4,:)=[ 0,  1, -1] ; vd( 4,:)=[ -1,  1,  1]
            vn( 5,:)=[ 1,  0,  1] ; vd( 5,:)=[ -1,  1,  1]
            vn( 6,:)=[ 1,  1,  0] ; vd( 6,:)=[ -1,  1,  1]
            vn( 7,:)=[ 0,  1,  1] ; vd( 7,:)=[  1, -1,  1]
            vn( 8,:)=[ 1,  0, -1] ; vd( 8,:)=[  1, -1,  1]
            vn( 9,:)=[ 1,  1,  0] ; vd( 9,:)=[  1, -1,  1]
            vn(10,:)=[ 0,  1,  1] ; vd(10,:)=[  1,  1, -1]
            vn(11,:)=[ 1,  0,  1] ; vd(11,:)=[  1,  1, -1]
            vn(12,:)=[ 1, -1,  0] ; vd(12,:)=[  1,  1, -1]
            do is=1,12
               vl(is,1)=vn(is,2)*vd(is,3)-vn(is,3)*vd(is,2)
               vl(is,2)=vn(is,3)*vd(is,1)-vn(is,1)*vd(is,3)
               vl(is,3)=vn(is,1)*vd(is,2)-vn(is,2)*vd(is,1)
               x1=dsqrt(sum(vn(is,:)**2))
               x2=dsqrt(sum(vd(is,:)**2))
               x3=dsqrt(sum(vl(is,:)**2))
               vn(is,:)=vn(is,:)/x1
               vd(is,:)=vd(is,:)/x2
               vl(is,:)=vl(is,:)/x3
               do i=1,3
               do j=1,3
                  smdMi(is,i,j)=vd(is,i)*vn(is,j)
               enddo
               enddo
               smdSMi(is,:,:)=( smdMi(is,:,:)
     &                       +transpose(smdMi(is,:,:)) )/2
               smdAMi(is,:,:)=( smdMi(is,:,:)
     &                       -transpose(smdMi(is,:,:)) )/2
               do i=1,6
                  smdVi1(is,i)=smdSMi(is,ib1(i),ib2(i))
               enddo
               smdVi2(is,1:3)=smdVi1(is,1:3)
               smdVi2(is,4:6)=smdVi1(is,4:6)*2
            enddo
c
            do is=1,12
            do js=1,12
               x1=sum(dabs(vn(is,:)-vn(js,:)))
               if(x1<1.d-10)then
                  hdm(is,js)=c_cpl
               else
                  hdm(is,js)=c_oth
               endif
            enddo
            enddo

            refv_pk2i = c44*1.d-6
            refv_IVB  = c44*1.d-6
            IVB_ini   = crss0 
c
            return
         end
c=================================================================
         subroutine sub_flhd_Ferrite(Iexp_loc,gam,
     &                 tau,hdM,IVB,tau_eff,IVB_eff,
     &                 dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &                 dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,
     &                 ising)
            implicit none
            real(8),parameter:: shrt0  = 1.d-3
            real(8),parameter:: pwfl   = 45.5d0
            real(8),parameter:: pwhd   = 5.d0
            real(8),parameter:: crsss  = 400.d0
            real(8),parameter:: hdrt0  = 500.d0
            real(8) hdm(48,48)

            integer i,j,is,js,ising
            integer Iexp_loc
            real(8) gam(48)
            real(8) tau(48),IVB(48)
            real(8) tau_eff(48),IVB_eff(48)
            real(8) dgmdt(48)
            real(8) ddgmdt_dtau(48)
            real(8) ddgmdt_dIVB(48)
            real(8) dIVBdt(48)
            real(8) ddIVBdt_ddgmdt(48,48)
            real(8) ddIVBdt_dIVB(48,48)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,12
               x1=crsss*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
c               print*,Is, IVB(is), x1,x2
c               read*
                  return
               endif
            enddo

c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB

            do is=1,12
c            do is=2,2

               dgmdt(is)=shrt0*(dabs(tau_eff(is))/IVB_eff(is))
     &                          **pwfl*dsign(1.d0,tau_eff(is))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0
     &            *(dabs(tau_eff(is))/IVB_eff(is))**(pwfl-1)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo
c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do is=1,12
            do js=1,12
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + hdm(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=hdm(is,js)*hdrt0
     &            *dsign(1.d0,tau_eff(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=hdm(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau_eff(js))*x1**pwhd 
     &            -hdm(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo

            return
         endsubroutine
      

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #3==>Mg                                             +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine sub_Magnesium_ini(STFei26,Nslpt,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,hdM,
     &              refv_pk2i,refv_IVB,IVB_ini)
            implicit none
c
c
            real(8),parameter:: c11    = 59.4d3
            real(8),parameter:: c12    = 25.6d3
            real(8),parameter:: c44    = 16.4d3 
            real(8),parameter:: c33    = 61.6d3 
            real(8),parameter:: c13    = 21.4d3  
            real(8) c66     

            real(8),parameter:: fitParaA1=1.d0
            real(8),parameter:: fitParaA2=1.d0
            real(8),parameter:: fitParaA3=1.d0
            real(8),parameter:: fitParaA4=1.d0
            
c            real(8),parameter:: crssba0  = fitParaA1        
c            real(8),parameter:: crsspr0  = fitParaA2
c            real(8),parameter:: crsspy0  = fitParaA3                               
c            real(8),parameter:: crsstw0  = fitParaA4  

c            real(8),parameter:: crssba0  = 9.6        
c            real(8),parameter:: crsspr0  = 192.
c            real(8),parameter:: crsspy0  = 192.                               
c            real(8),parameter:: crsstw0  = 37.8 

            real(8),parameter:: crssba0  = 12.    !10.d0  *2
            real(8),parameter:: crsspr0  = 240.   !80.d0  *2.5
            real(8),parameter:: crsspy0  = 240.   !88.d0  *2.5                               
            real(8),parameter:: crsstw0  = 63.    !35        *1  


            real(8),parameter:: c_cpl  = 1.d0
            real(8),parameter:: c_oth  = 1.d0
            real(8),parameter:: c_t2s  = 2.d0   
            real(8),parameter:: c_s2t  = 0.d0
            real(8),parameter:: c_t2t  = 1.d0
 
c            real(8),parameter:: c_cpl  = 1.d0
c            real(8),parameter:: c_oth  = 1.d0
c            real(8),parameter:: c_t2s  = 0.d0   
c            real(8),parameter:: c_s2t  = 0.d0
c            real(8),parameter:: c_t2t  = 1.d0
 
            real(8) hdM(48,48)
            integer Nslp,Nslpt,Ntwn
            real(8) STFei26(6,6)
            real(8) smdMi(48,3,3)
            real(8) smdSMi(48,3,3)
            real(8) smdAMi(48,3,3)
            real(8) smdVi1(48,6)
            real(8) smdVi2(48,6)
            real(8) vd(48,3)
            real(8) vl(48,3)
            real(8) vn(48,3)
            real(8) crss0(48)
            real(8) refv_pk2i
            real(8) refv_IVB
            real(8) IVB_ini(48)
            integer i,j,is,js
            real(8) x1,x2,x3
            integer IB1(9),IB2(9)
            data IB1/1,2,3,1,1,2,2,3,3/
            data IB2/1,2,3,2,3,3,1,1,2/
c
            STFei26(:,:)=0
            c66=0.5*(c11-c12)
            STFei26(1,1)=c11
            STFei26(2,2)=c11
            STFei26(3,3)=c33
            STFei26(4,4)=c44*2
            STFei26(5,5)=c44*2
            STFei26(6,6)=c66*2
            STFei26(2,3)=c13
            STFei26(3,2)=c13
            STFei26(1,3)=c13
            STFei26(3,1)=c13
            STFei26(1,2)=c12
            STFei26(2,1)=c12
c       
            Nslpt = 24     ! total number 24  
            Nslp  = 18     ! 6 2nd order pyramidal plane are not considered here 
            Ntwn  = 6      ! 6 tensile twining 

            !* vd is the slip direction                               
            !* 3 basal
            vd( 1,:)=[ 3.0,  0., 0.                ] 
            vd( 2,:)=[ -1.5,  2.598076211353316, 0.] 
            vd( 3,:)=[ -1.5, -2.598076211353316, 0.] 
            !* 3 prismatic
            vd( 4,:)=[-1.5,  2.598076211353316, 0.] 
            vd( 5,:)=[ 3.0,  0., 0.               ] 
            vd( 6,:)=[ -1.5, -2.598076211353316,0.] 
            !* 12 1st order pyramidal
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
            !* 6 tensile twin
            vd(19,:)=[ 1.5, -0.866025403784438, 1.6330] 
            vd(20,:)=[-1.5, 0.866025403784438,  1.6330] 
            vd(21,:)=[ -1.5, -0.866025403784438,1.6330]
            vd(22,:)=[ 1.50, 0.866025403784438, 1.6330] 
            vd(23,:)=[ 0., 1.732050807568877,   1.6330]
            vd(24,:)=[ 0., -1.732050807568877,  1.6330] 
            ! vn is the normal direction of planes
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
            !* twining
            vn(19,:)=[ -1.5,0.866025403784438, 3.2660]
            vn(20,:)=[1.5,  -0.866025403784438,3.2660]
            vn(21,:)=[1.5,  0.866025403784438, 3.2660]
            vn(22,:)=[-1.5, -0.866025403784438,3.2660]
            vn(23,:)=[ 0., -1.732050807568877, 3.2660]
            vn(24,:)=[0., 1.732050807568877,   3.2660]
            !
            do is=1,24
               vl(is,1)=vn(is,2)*vd(is,3)-vn(is,3)*vd(is,2)
               vl(is,2)=vn(is,3)*vd(is,1)-vn(is,1)*vd(is,3)
               vl(is,3)=vn(is,1)*vd(is,2)-vn(is,2)*vd(is,1)
               x1=dsqrt(sum(vn(is,:)**2))
               x2=dsqrt(sum(vd(is,:)**2))
               x3=dsqrt(sum(vl(is,:)**2))
               vn(is,:)=vn(is,:)/x1
               vd(is,:)=vd(is,:)/x2
               vl(is,:)=vl(is,:)/x3
               do i=1,3
               do j=1,3
                  smdMi(is,i,j)=vd(is,i)*vn(is,j)
               enddo
               enddo
               smdSMi(is,:,:)=( smdMi(is,:,:)
     &                       +transpose(smdMi(is,:,:)) )/2
               smdAMi(is,:,:)=( smdMi(is,:,:)
     &                       -transpose(smdMi(is,:,:)) )/2
               do i=1,6
                  smdVi1(is,i)=smdSMi(is,ib1(i),ib2(i))
               enddo
               smdVi2(is,1:3)=smdVi1(is,1:3)
               smdVi2(is,4:6)=smdVi1(is,4:6)*2
            enddo
            !** slip-slip
            do is=1,18
            do js=1,18
               x1=sum(dabs(vn(is,:)-vn(js,:)))
               if(x1<1.d-10)then
                  hdM(is,js)=c_cpl
               else                             
                  hdM(is,js)=c_oth
               endif
            enddo
            enddo
            !** twin-slip  
            do is=19,24
            do js= 1,18
               hdM(js,is) = c_t2s
            enddo
            enddo    
            !** twin-twin
            do is=19,24
            do js=19,24
               hdM(js,is) = c_t2t
            enddo
            enddo
            !** slip-twin
            do is= 1,18
            do js=19,24
               hdM(js,is) = c_s2t !=0
            enddo
            enddo  
            !
            refv_pk2i      = c44*1.d-6  
            refv_IVB       = c44*1.d-6
            IVB_ini( 1: 3) = crssba0   !* sl
            IVB_ini( 4: 6) = crsspr0   !* sl
            IVB_ini( 7:18) = crsspy0   !* sl
            IVB_ini(19:24) = crsstw0   !* tw

            return
         endsubroutine
         
c=================================================================
         subroutine sub_flhd_Magnesium(Iexp_loc,gam,
     &                 tau,hdM,IVB,tau_eff,IVB_eff,
     &                 dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &                 dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,
     &                 ising)     
            implicit none
c
            real(8),parameter:: fitParaB1=1.d0
            real(8),parameter:: fitParaB2=1.d0
            real(8),parameter:: fitParaB3=1.d0
            real(8),parameter:: fitParaB4=1.d0
c
c            real(8),parameter:: shrts0   = 0.001
c            real(8),parameter:: shrtt0   = 0.001
c            real(8),parameter:: pwfl     = 10.d0
c            real(8),parameter:: pwhd     = 2.2  
c            real(8),parameter:: atw      = 2.2  
c
c            real(8),parameter:: crssbas  =  24.  
c            real(8),parameter:: crssprs  = 270. 
c            real(8),parameter:: crsspys  = 270. 
c            real(8),parameter:: crsstws  =  65. 
cc            
c            real(8),parameter:: hdrtba0  = 20.       
c            real(8),parameter:: hdrtpr0  = 400.     
c            real(8),parameter:: hdrtpy0  = 400.     
c            real(8),parameter:: htw0     = 20.       
c------------

            real(8),parameter:: shrts0  = 0.01  !*0   !  1.d-2     !* reference shear rate
            real(8),parameter:: shrtt0  = 0.1   !*0.1   !   2.d-0    !* twin rate 0
            real(8),parameter:: pwfl    = 10.0     !  10.d0     !* 1/ m (m = 0.1)
            real(8),parameter:: pwhd    = 0.6      !  0.6d0     !* asl
            real(8),parameter:: atw     = 1.d0      !* atw
            
            real(8),parameter:: crssbas  = 250.d0 
            real(8),parameter:: crssprs  = 1880.d0 
            real(8),parameter:: crsspys  = 3580.d0
            real(8),parameter:: crsstws  = 100.d0
                                                
            real(8),parameter:: hdrtba0  = 20.d0    
            real(8),parameter:: hdrtpr0  = 2831.d0 /2 
            real(8),parameter:: hdrtpy0  = 2990.d0 /2 
            real(8),parameter:: htw0     = 24.d0          





c
            real(8),parameter:: gammatw  = 0.129    !* characteristic shear of twinning mode
c
            integer i,j,is,js,ising
            integer Iexp_loc
            real(8) gam(48)
            real(8) tau(48),IVB(48)
            real(8) tau_eff(48),IVB_eff(48)
            real(8) dgmdt(48)
            real(8) ddgmdt_dtau(48)
            real(8) ddgmdt_dIVB(48)
            real(8) dIVBdt(48)
            real(8) ddIVBdt_ddgmdt(48,48)
            real(8) ddIVBdt_dIVB(48,48)
            real(8) hdM(48,48)
            real(8) crsss(48)
            real(8) hdrt0(48)
c
            real(8) x1,x2,x3,x10,x11,x12
            real(8) VFslip,VFtwin
            real(8) dr,dd 
            real(8) fct0,fct1,fct2,fct3 
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0

            crsss(1:3)   = crssbas
            crsss(4:6)   = crssprs
            crsss(7:18)  = crsspys
            crsss(19:24) = crsstws

            hdrt0(1:3)   = hdrtba0
            hdrt0(4:6)   = hdrtpr0
            hdrt0(7:18)  = hdrtpy0
            hdrt0(19:24) = htw0

            VFtwin = sum(dabs(gam(19:24)))/gammatw
            VFslip = 1-VFtwin 
            
c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB, 
            !** slip rate:  
            do is=1,18
               fct1=dsign(1.d0,tau_eff(is))
               fct2=dabs(tau_eff(is)/IVB_eff(is))
               dgmdt(is)= shrts0*fct2**pwfl*fct1 
               if(Iexp_loc/=1)then               
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrts0*fct2**(pwfl-1)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is) 
               endif
            enddo
            
            !** twin rate: 
            do is=19,24
               fct0=dsign(1.d0,tau_eff(is)) 
               fct1=1.0 !  (dsign(1.d0, 0.8-VFtwin)+1.d0)/2

               if(tau_eff(is)>0)then
                  fct2=fct0
               else
                  fct2=fct0 *(dsign(1.d0, gam(is))+1.d0)/2
               endif

               fct3=tau_eff(is)/IVB_eff(is)
               dgmdt(is)=fct1*fct2*shrtt0*(fct3)**pwfl 
               if(Iexp_loc/=1)then               
                  ddgmdt_dtau(is)=fct1*fct2*shrtt0
     &            *pwfl*(fct3)**(pwfl-1)/IVB_eff(is)           
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)          
               endif
 
            enddo

c            do is=19,24
c               fct1=(dsign(1.d0,0.8d0-VFtwin)+1.d0)/2
c               fct2=(dsign(1.d0,tau_eff(is)) +1.d0)/2
c               fct3=tau_eff(is)/IVB_eff(is)
c               if(tau_eff(is)>0.d0) then
c                  dgmdt(is)=fct1*fct2*shrtt0*(fct3)**pwfl 
c                  if(Iexp_loc/=1)then               
c                     ddgmdt_dtau(is)=fct1*fct2*shrtt0
c     &               *pwfl*(fct3)**(pwfl-1)/IVB_eff(is)           
c                     ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)          
c                  endif
c               endif
c            enddo


            return

c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            !** slip part, dIVBdt 
            do is=1,18
               do js=1,18          !--hardening due to slip
                  x1=1.d0-IVB(js)/crsss(js)
                  fct1=dabs(x1)
                  fct2=dsign(1.d0,x1)
                  fct3=dsign(1.d0,dgmdt(js))
                  dIVBdt(is)=dIVBdt(is) + hdM(is,js)*hdrt0(js)
     &                      *fct1**pwhd*fct2*dgmdt(js)*fct3 
                  if(Iexp_loc/=1)then
                     ddIVBdt_ddgmdt(is,js)=hdM(is,js)*hdrt0(js)      
     &               *fct1**pwhd*fct2*fct3              
                     ddIVBdt_dIVB(is,js)=hdM(is,js)*hdrt0(js)
     &               *( fct1**pwhd*fct2*ddgmdt_dIVB(js)*fct3                                                                              
     &                 -pwhd*fct1**(pwhd-1)/crsss(js)*dgmdt(js)*fct3 )
                  endif
               enddo
               do js=19,24  !--hardening due to twin
                  x1=1.d0-IVB(js)/crsss(js)
                  fct1=dabs(x1)
                  fct2=dsign(1.d0,x1)
                  fct3=dsign(1.d0,dgmdt(js))
                  dIVBdt(is)=dIVBdt(is) + hdM(is,js)*hdrt0(js)
     &                      *fct1**atw*fct2*dgmdt(js)*gammatw*fct3       !* gammatw is included in dgmdt(is) here . 0.8 is included in dgmdt
                  if(Iexp_loc/=1)then
                     ddIVBdt_ddgmdt(is,js)=hdM(is,js)*hdrt0(js)      
     &               *fct1**atw*fct2*gammatw*fct3              
                     ddIVBdt_dIVB(is,js)=hdM(is,js)*hdrt0(js)*gammatw
     &               *( fct1**atw*fct2*ddgmdt_dIVB(js)*fct3                                                                              
     &                 -atw*fct1**(atw-1)/crsss(js)*dgmdt(js)*fct3 )
                  endif
               enddo
            enddo
            !     
            !** twin part (Nslp+1:Nslp+Ntwn,Nslp+1:Nslp+Ntwn) dIVBdt and derivatives !** the following is checked
            !* dIVBdt first 
            do is=19,24
            do js=19,24
               x1=1.d0-IVB(js)/crsss(js)
               fct1=dabs(x1)
               fct2=dsign(1.d0,x1)
               fct3=dsign(1.d0,dgmdt(js))
               dIVBdt(is)=dIVBdt(is) + hdM(is,js)*hdrt0(js)     !* hdM(is,js) is qt2t here  
     &                   *fct1**atw*fct2*dgmdt(js)*gammatw*fct3 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=hdM(is,js)*hdrt0(js)      
     &            *fct1**atw*fct2*gammatw*fct3              
                  ddIVBdt_dIVB(is,js)=hdM(is,js)*hdrt0(js)*gammatw
     &            *( fct1**atw*fct2*ddgmdt_dIVB(js)*fct3                                                                              
     &              -atw*fct1**(atw-1)/crsss(js)*dgmdt(js)*fct3 )
               endif
            enddo           
            enddo

c            return

!------------------------------------------
!** add the constants multiplier
!* slips are multiplied with volume fraction for parent part, and twins are multiplied with 0.129
            dgmdt(1:18)        = dgmdt(1:18)*VFslip
            ddgmdt_dtau(1:18)  = ddgmdt_dtau(1:18)*VFslip
            ddgmdt_dIVB(1:18)  = ddgmdt_dIVB(1:18)*VFslip
            dgmdt(19:24)       = dgmdt(19:24)*gammatw
            ddgmdt_dtau(19:24) = ddgmdt_dtau(19:24)*gammatw
            ddgmdt_dIVB(19:24) = ddgmdt_dIVB(19:24)*gammatw
         ddIVBdt_ddgmdt(1:24,1:18)=ddIVBdt_ddgmdt(1:24,1:18)/VFslip 
         ddIVBdt_ddgmdt(1:24,19:24)=ddIVBdt_ddgmdt(1:24,19:24)/gammatw 
!------------------------------------------

            return
         end
