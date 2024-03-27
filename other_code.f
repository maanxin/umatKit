cccccccccccccccccccccccccccccccccccccccccccc
      subroutine caleulang(Mx,vv,ising)     
cccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer irdummy,ising,ising_1
      real(8) det
      real(8) vv(3)
      real(8) Mx(3,3),rm(3,3),um(3,3)
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
cccccccccccccccccccccccccccccccccccccccccccc
      ising_1=0
!      call PDECOMPOSITION(Mx,Um,Rm,ising_1)
      call polar_decomp(Mx,Um,Rm,ising_1)
      if(ising_1/=0)then
         write(6,*) 'pdecompsition of Mx is failure'
         call pm(Mx,3,3)
         ising=1
         return
      endif
      call icams_Q2Eang(Rm,vv(1),vv(2),vv(3))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_misori(v1,v2,ang)
ccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) v1(3),v2(3),ang,x1,x2
      real(8) QM1(3,3),QM2(3,3),dQM(3,3)
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
      call icams_Eang2Q(v1(1),v1(2),v1(3),QM1)
      call icams_Eang2Q(v2(1),v2(2),v2(3),QM2)
      dQM=matmul(QM2,transpose(QM1))
      x1=dQM(1,1)+dQM(2,2)+dQM(3,3)
      x2=(x1-1.d0)/2
      if(dabs(x2)>1.d0)x2=1.d0*sign(1.d0,x2)
      ang=dabs(pi/2-dasin(x2)) *r2g
      return
      end

C****************************************************************
      subroutine pdecomposition(Mx,UMx,RMx,ising)
C****************************************************************
      implicit none
      integer ising
      real(8) Mx(3,3),ce(3,3),RMx(3,3),UMx(3,3),IUMx(3,3)
      real(8) eb1(3,3),eb2(3,3),eb3(3,3)
      real(8) ev1,ev2,ev3,det
      ising=0
      ce=matmul(transpose(Mx),Mx)
      call spectral(ce,ev1,ev2,ev3,eb1,eb2,eb3,ising)
      if(ev1<=0 .or. ev2<=0 .or. ev3<=0 .or. ising/=0)then
         write(6,*) 'eigen value of ce <0'
         print*, ev1,ev2,ev3
         ising=1
         return
      endif
      UMx=dsqrt(ev1)*eb1+dsqrt(ev2)*eb2+dsqrt(ev3)*eb3
      IUMx=1/dsqrt(ev1)*eb1+1/dsqrt(ev2)*eb2+1/dsqrt(ev3)*eb3
      RMx=matmul(Mx,IUMx)
      return 
      end

      subroutine icams_determ(a,det)
c********************************************************************
c     This routine calculates the determinant of a
c********************************************************************
      implicit none
      real(8) a(3,3),v1,v2,v3,det
      v1= a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      v2= a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      v3= a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      det= v1-v2+v3
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_Eang2Q(p1,p,p2,QM)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     rotate from X[100],Y[010],Z[001] to v1,v2,v3
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) QM(3,3)
      real(8) p1,p,p2,xp1,xp,xp2
      real(8) c1,c,c2,s1,s,s2
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
      xp1=p1*g2r
      xp =p *g2r
      xp2=p2*g2r
      c1=dcos(xp1)
      s1=dsin(xp1)
      c =dcos(xp )
      s =dsin(xp )
      s2=dsin(xp2)
      c2=dcos(xp2)
      QM(1,1)=+c1*c2-s1*s2*c
      QM(1,2)=+s1*c2+c1*s2*c
      QM(1,3)=+s2*s
      QM(2,1)=-c1*s2-s1*c2*c
      QM(2,2)=-s1*s2+c1*c2*c
      QM(2,3)=+c2*s
      QM(3,1)=+s1*s
      QM(3,2)=-c1*s
      QM(3,3)=+c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_Q2Eang(QM,phi1,PHI,phi2)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) QM(3,3)
      real(8) phi1,PHI,phi2
      real(8) sqhkl,squvw,sqhk,val
      real(8) pi,r2g,g2r,Tol
      Tol=1.d-15
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
c---------------------------------------------
c
c             v1   v2    v3
c
c           | u   v2_1    h  |
c     QM =  | v   v2_2    k  |  with v2 = v3 x v1   
c           | w   v2_3    l  | 
c
c             100   010   001
c          X|  u   v2_1    h  |
c     QM = Y|  v   v2_2    k  |    
c          Z|  w   v2_3    l  | 
c
c---------------------------------------------


c---------------------------------------------
c
c    if the roation tensor is defined as following
c    QM_ij:= (dx/dX)_ij = dx_i/dX_j
c
c               X     Y     X
c          v1|  u   v2_1    h  |
c     QM'= v2|  v   v2_2    k  |    
c          v3|  w   v2_3    l  | 
c    the ratation must be take 
c
c    QM=transpose(QM')
c
c---------------------------------------------
      squvw=dsqrt(QM(1,1)**2+QM(2,1)**2+QM(3,1)**2)
      sqhkl=dsqrt(QM(1,3)**2+QM(2,3)**2+QM(3,3)**2)
      sqhk =dsqrt(QM(1,3)**2+QM(2,3)**2           )

      val=QM(3,3)/sqhkl
      if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
      PHI=dacos(val)

      if(PHI < TOL) then
         phi2=0.0
         val=QM(1,1)/squvw
         if(QM(2,1) <= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi1=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi1=2*pi-dacos(val)
            phi1=-dacos(val)
         endif
      else
c        val=QM(2,3)/sqhk
         val=QM(2,3)/dsin(PHI)
         if(QM(1,3) >= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi2=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi2=2*pi-dacos(val)
            phi2=-dacos(val)
         endif
         val=-QM(3,2)/dsin(PHI)
         if(QM(3,1) >= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi1=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi1=2*pi-dacos(val)
            phi1=-dacos(val)
         endif
      endif
      phi1=phi1*r2g
      PHI=PHI*r2g
      phi2=phi2*r2g
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_Q2Eang_old(QM,phi1,PHI,phi2)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) QM(3,3)
      real(8) phi1,PHI,phi2,x1,x2,x3
      real(8) pi,r2g,g2r,Tol
      Tol=1.d-15
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180
c---------------------------------------------
c            v1   v2=v3xv1    v3
c
c           | u   v2_1    h  |
c     QM =  | v   v2_2    k  |    
c           | w   v2_3    l  | 
c
c---------------------------------------------

      x1=dsqrt(QM(1,1)**2+QM(2,1)**2+QM(3,1)**2)
      x2=dsqrt(QM(1,2)**2+QM(2,2)**2+QM(3,2)**2)
      x3=dsqrt(QM(1,3)**2+QM(2,3)**2+QM(3,3)**2)
      QM(:,1)=QM(:,1)/x1
      QM(:,2)=QM(:,2)/x1
      QM(:,3)=QM(:,3)/x1

      x1=QM(3,3)
      if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

      PHI=dacos(x2)

      if(PHI<Tol)then
         phi2=0.d0

         x1=QM(1,1)
         if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

         if(QM(1,2)>=0.d0)then
            phi1=dacos(x2)
         else
c           phi1=2*pi-dacos(val)
            phi1=-dacos(x2)
         endif
      else

         x1=QM(2,3)/dsin(PHI)
         if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

         if(QM(1,3)>=0.d0) then
            phi2=dacos(x2)
         else
c           phi2=2*pi-dacos(val)
            phi2=-dacos(x2)
         endif

         x1=-QM(3,2)/dsin(PHI)
         if(dabs(x1)>1.d0)x2=1.d0*sign(1.d0,x1)

         if(QM(3,1) >= 0.d0) then
            phi1=dacos(x2)
         else
c           phi1=2*pi-dacos(val)
            phi1=-dacos(x2)
         endif
      endif
      phi1=phi1*r2g
      PHI=PHI*r2g
      phi2=phi2*r2g
      return
      end




cccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_angax2QM(ang,u,v,w,QM)
cccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) QM(3,3)
      real(8) s,c,u2,v2,w2,ang,u,v,w,x1
      real(8) pi,r2g,g2r
      pi=dacos(-1.d0)
      r2g=180/pi
      g2r=pi/180

      x1=dsqrt(u**2+v**2+w**2)
      u2=u/x1
      v2=v/x1
      w2=w/x1
      s=dsin(ang)
      c=dcos(ang)

      QM(1,1)=(1-u2**2)*c+u2**2
      QM(2,2)=(1-v2**2)*c+v2**2
      QM(3,3)=(1-w2**2)*c+w2**2

      QM(1,2)=u2*v2*(1-c)+w2*s
      QM(2,1)=u2*v2*(1-c)-w2*s

      QM(1,3)=u2*w2*(1-c)-v2*s
      QM(3,1)=u2*w2*(1-c)+v2*s

      QM(2,3)=v2*w2*(1-c)+u2*s
      QM(3,2)=v2*w2*(1-c)-u2*s


      return
      end

c
c
C********************************************************************** 
      SUBROUTINE HI(M,HI1M,HI2M,HI3M)
C**** HAUPTINVARIANTEN HI1M, HI2M, HI3M DER 3X3 MATRIX M
      IMPLICIT NONE
      real(8) M(3,3),HI1M,HI2M,HI3M 
      HI1M = M(1,1)+M(2,2)+M(3,3)
      HI2M =(M(1,1)+M(2,2)+M(3,3))**2/2.d0
     &      -M(1,1)**2/2.d0
     &      -M(2,2)**2/2.d0
     &      -M(3,3)**2/2.d0
     &      -M(1,2)*M(2,1)
     &      -M(1,3)*M(3,1)
     &      -M(2,3)*M(3,2)
      HI3M =+M(1,1)*M(2,2)*M(3,3)
     &      +M(2,1)*M(1,3)*M(3,2)
     &      +M(3,1)*M(1,2)*M(2,3)
     &      -M(1,1)*M(2,3)*M(3,2)
     &      -M(2,1)*M(1,2)*M(3,3)
     &      -M(3,1)*M(1,3)*M(2,2)
      RETURN  
      END


C**********************************************************************
      SUBROUTINE NORM(M,ZM,SM,NO)
C**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,ZM,SM
      DOUBLE PRECISION M(ZM,SM),NO  
      NO=0.d0
      DO I=1,ZM
      DO J=1,SM 
         NO=NO+M(I,J)**2.d0 
      END DO 
      END DO
      NO=dsqrt(NO)
      RETURN 
      END  


C**********************************************************************
      subroutine wm(M,ni,nj)
      implicit none
      integer ni,nj,i,j
      real(8) M(ni,nj)
      write(6,*) 
      do i=1,ni
         write(6,100) (m(i,j),j=1,nj)          
      enddo
      write(6,*) 
100   format(48e20.8)       
      return
      end


C**********************************************************************
      subroutine pm(M,ni,nj)
      implicit none
      integer ni,nj,i,j
      real(8) M(ni,nj)
      write(*,*) 
      do i=1,ni
         !write(*,100) (m(i,j),j=1,nj)          
         print '(48e15.4)', (M(i,j),j=1,nj)          
      enddo
      write(*,*) 
c100   format(48e20.8)        
100   format(48e12.4)       
c100   format(48f6.3)       
      return
      end

C**********************************************************************
      subroutine pv(M,ni)
      implicit none
      integer ni,i
      real(8) M(ni)
      write(*,100) (M(i),i=1,ni)          
c      write(6,*) 
c100   format(48e20.8)       
100   format(48e12.4)       
      return
      end


C******************************************** 
      subroutine gaussj(A,n,B,ising)
C******************************************** 
      implicit none
      integer   n,ising
      integer   i,icol,irow,j,k,l,ll,i1,i2
      integer   indxc(n),indxr(n),ipiv(n)
      real(8)   A0(n,n),A(n,n),B(n,n),C(n,n),vx(n)
      real(8)   big,dum,pivinv,x1,x2,x3
c-------------------------------------------
c      write(6,*) 'coming to jd'
c      print*, 'a = ', a
c      call pm(A,3,3)
c      call flush(6)

      ising=0
      A0=A
      C=A
      ipiv=0
c-----loop for the pivot procedures, from 1 to n
      do i=1,n
         big=0.d0
         do j=1,n
         do k=1,n
            if(ipiv(j).ne.1 .and. ipiv(k)==0)then
               if(dabs(a(j,k)).ge.big)then    
                  big=dabs(a(j,k))         
                  irow=j
                  icol=k
               endif
            endif
            if(ipiv(j).ne.1 .and. ipiv(k).gt.1)then
c              print*,'sigular matrix in gauss_jordan 1'
c              write(6,*)'sigular matrix in gauss_jordan'
c              call flush(6)
               ising=1
               A=A0
               return
            endif
         enddo
         enddo
c--------check whether the pivot element is zero or not
         if(a(irow,icol)==0.d0)then
c	      print*,'sigular matrix in gauss_jordan 2'
c	      print*,'indices is:', irow,icol
c	      print*, 'a = ', a
c           write(6,*)'sigular matrix in gauss_jordan'
c           write(6,*)'indices is:', irow,icol
c           call flush(6)
            ising=1
            A=A0
            return
         endif
c-------------------------------------------------------
c        if one component is selected as pivot element
c        the second indice is important, so it is marked
c        from 0 to 1 in ipiv(:) array
c        after convert, it is the row number
c-------------------------------------------------------
         ipiv(icol)=ipiv(icol)+1
c--------record the row and collum number for ith pivot element
         indxr(i)=irow
         indxc(i)=icol
c--------change pivot element to diagonal position
         if(irow.ne.icol)then
            vx=a(irow,:)
            a(irow,:)=a(icol,:)
            a(icol,:)=vx
         endif
c--------eliminate the elements besides a(icol,icol)
         pivinv=1.d0/a(icol,icol)
         a(icol,icol)=1.d0
         a(icol,:)=a(icol,:)*pivinv 
         do i2=1,n
            if(i2.ne.icol)then
               dum=a(i2,icol)
               a(i2,icol)=0.d0
               a(i2,:)=a(i2,:)-a(icol,:)*dum
            endif
         enddo
      enddo

c-----after maximum pivot strategy elimination
c-----rearrage the left matrix
      do l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            vx=a(:,indxr(l))
            a(:,indxr(l))=a(:,indxc(l))
            a(:,indxc(l))=vx
         endif
      enddo
      B=A
      A=C

      return
      end


C***************************************************************C
C      SUBROUTINE SPOLAR (F,U,R)				        C
C							  	        C
C *performs the polar decomposition of the	               C
C   tensor F=RU using a Cayley-Hamilton theorem.		        C
C *Then, U=R^T F is found					        C
C***************************************************************C
      SUBROUTINE polar_decomp(F,U,R,ising)
      implicit none
      integer i,j,k,iflag1,iflag2,ising
      real(8) F(3,3),U(3,3),R(3,3),RT(3,3),C(3,3),CS(3,3),UINV(3,3)
      real(8) x1,x2,x3,f1,f2,c1,c2,c3,u1,u2,u3,b,b1,b2,b3,b4
C
      C=matmul(transpose(F),F)
      CS=matmul(C,C)

c-----1st, 2st, 3st invariant for tensor C
      C1=C(1,1)+C(2,2)+C(3,3)
      C2=(C1**2.d0-(CS(1,1)+CS(2,2)+CS(3,3)))/2
      C3=+C(1,1)*(C(2,2)*C(3,3)-C(2,3)*C(3,2))
     1	 -C(1,2)*(C(2,1)*C(3,3)-C(2,3)*C(3,1))
     2	 +C(1,3)*(C(2,1)*C(3,2)-C(2,2)*C(3,1))

c-----3st invariant for tensor U
      U3=dsqrt(C3)

      X1= 2.0**5.0 /27.0 * (2.0*C1**3.0-9.0*C1*C2+27.0*C3)
      X2= 2.**10.0 /27.0 * (4.0*C2**3.0 - C1**2.0*C2**2.0 +
     1	  4.0*C1**3.0*C3 - 18.0 * C1*C2*C3 + 27.0 * C3**2.0)

      IF(X2<0)X2=0
      F1=X1+dsqrt(X2)
      IFLAG1=0
      IF(F1<0)IFLAG1=1
      F2=X1-dsqrt(X2)
      IFLAG2=0
      IF(F2<0)IFLAG2=1

      IF(IFLAG1==1) F1=-F1
      IF(IFLAG2==1) F2=-F2

      X3= -2.0/3.0*C1 + F1**(1.0/3.0) + F2**(1.0/3.0)
      IF(IFLAG1==1) X3= -2.0/3.0*C1 - F1**(1.0/3.0) + F2**(1.0/3.0)
      IF(IFLAG2==1) X3= -2.0/3.0*C1 + F1**(1.0/3.0) - F2**(1.0/3.0)

c-----1st, 2st invariant for tensor U
      B=-2.0*C1
      if(X3==B)then
         U1=dsqrt(C1+2.0*dsqrt(C2))
      else

         x1=dsqrt(2.0*C1+X3)
         if(x1==0)then
            ising=1
            return
         endif      
         U1= 0.5 * ( x1 + dsqrt(2.0*C1 - X3 + 16.0*dsqrt(C3)/x1) )
      endif
      U2=dsqrt(C2+2.0*U3*U1)

      B1= U3**2.0 * (U3+U1*C1) + U1**2.0 * (U1*C3+U3*C2)
      if(B1==0)then
         ising=1
         return
      endif      

      B2= U1*(U1*U2-U3) / B1
      B3=-(U1*U2-U3) * (U3+U1*C1) / B1
      B4= (U2*U3*(U3+U1*C1) + U1**2.0 * (U2*C2+C3))/B1

      UINV=B2*CS + B3*C
      Uinv(1,1)=Uinv(1,1)+B4
      Uinv(2,2)=Uinv(2,2)+B4
      Uinv(3,3)=Uinv(3,3)+B4

      R=matmul(F,UINV)
      U=matmul(transpose(R),F)

      RETURN
      END


C**********************************************************************
      SUBROUTINE spectral(M,EW1,EW2,EW3,EB1,EB2,EB3,ising)
C**********************************************************************
      implicit none
      integer ising
      real(8) M(3,3),EB1(3,3),EB2(3,3),EB3(3,3),EW1,EW2,EW3,
     &        HI1M,HI2M,HI3M,TOL,R,S,T,P,Q,RHO,PHI,Y1,Y2,Y3,D1,D2,D3,
     &        E(3,3),M1(3,3),M2(3,3),M3(3,3)
      real(8) x1,x2,x3,PI
      ising=0
      TOL=1.d-15
      PI=dacos(-1.d0)
      CALL HI(M,HI1M,HI2M,HI3M)
      R=-HI1M
      S= HI2M
      T=-HI3M
      P=S-R**2.d0/3.d0
      Q=2.d0/27.d0*R**3.d0-R*S/3.d0+T
      RHO=dsqrt(-3.d0*P**3.d0)/9.d0
      if(dabs(RHO)<=Tol)then
         ising=1
         return         
      endif
      !======================================
      !  modification for overcome overflow
      !======================================
      if(dabs(RHO)<=Tol) RHO=sign(1.d0,RHO)*TOL
      x1=-Q/RHO/2.d0
      if(dabs(x1)>1.d0) x1=1.d0*sign(1.d0,x1)
      PHI=dacos(x1)
      Y1=2*RHO**(1.d0/3.d0)*dcos(PHI/3.d0)
      Y2=2*RHO**(1.d0/3.d0)*dcos(PHI/3.d0+2.d0/3.d0*PI)
      Y3=2*RHO**(1.d0/3.d0)*dcos(PHI/3.d0+4.d0/3.d0*PI)
      EW1=Y1-R/3.d0
      EW2=Y2-R/3.d0
      EW3=Y3-R/3.d0
      E(1,:)=[1,0,0]
      E(2,:)=[0,1,0]
      E(3,:)=[0,0,1]
      EB1=0.d0
      EB2=0.d0
      EB3=0.d0
      IF(dabs(ew1-ew2)<tol .and. dabs(ew1-ew3)<tol)THEN
         EB1=E
      elseif(dabs(ew2-ew3)<TOL)THEN
         EB1=MATMUL(M-EW2*E,M-EW3*E)/((ew1-ew2)*(ew1-ew3))
         EB2=E-EB1
      elseif(dabs(ew1-ew3)<TOL)THEN
         EB2=MATMUL(M-EW1*E,M-EW3*E)/((ew2-ew1)*(ew2-ew3))
         EB1=E-EB2
      elseif(dabs(ew1-ew2)<TOL)THEN
         EB3=MATMUL(M-EW1*E,M-EW2*E)/((ew3-ew1)*(ew3-ew2))
         EB1=E-EB3
      else
         EB1=MATMUL(M-EW2*E,M-EW3*E)/((ew1-ew2)*(ew1-ew3))
         EB2=MATMUL(M-EW1*E,M-EW3*E)/((ew2-ew1)*(ew2-ew3))
         EB3=MATMUL(M-EW1*E,M-EW2*E)/((ew3-ew1)*(ew3-ew2))
      endif

      RETURN
      END


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Fast Fourier Transformation Code                            +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fourn(Datx,nn,ndim,isg)
      implicit none
      integer isg,ndim,nn(ndim)
      integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1
      integer k2,n,nprev,nrem,ntot,ifp2,ip1,ip2,ip3,k1
      real(8) Datx(*)
      real(8) tempi,tempr
      real(8) theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then

            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=Datx(i3)
                tempi=Datx(i3+1)
                Datx(i3)=Datx(i3rev)
                Datx(i3+1)=Datx(i3rev+1)
                Datx(i3rev)=tempr
                Datx(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif

          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isg*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*Datx(k2)-sngl(wi)*Datx(k2+1)
                tempi=sngl(wr)*Datx(k2+1)+sngl(wi)*Datx(k2)
                Datx(k2)=Datx(k1)-tempr

                Datx(k2+1)=Datx(k1+1)-tempi
                Datx(k1)=Datx(k1)+tempr
                Datx(k1+1)=Datx(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      end

c**********************************************
      subroutine icams_conv33to6(Am,ibx1,ibx2,Av)
c**********************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) Am(3,3),Av(6)
      do i=1,6
         Av(i)=Am(ibx1(i),ibx2(i))
      enddo
      return
      end


c***********************************************
      subroutine icams_conv33to9(Am,ibx1,ibx2,Av)
c***********************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) Am(3,3),Av(9)
      do i=1,9
         Av(i)=Am(ibx1(i),ibx2(i))
      enddo
      return
      end


c**************************************************
      subroutine icams_conv6to33(Av,ibx1,ibx2,Am)
c**************************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) Am(3,3),Av(6)
      do i=1,6
         Am(ibx1(i),ibx2(i))=Av(i)
      enddo
      do i=7,9
         Am(ibx1(i),ibx2(i))=Av(i-3)
      enddo
      return
      end


c====================================================================================================
c     Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
c     by n, stored in a physical np by np array. On output, elements of a above the diagonal are
c     destroyed. d returns the eigenvalues of a in its first n elements. v is a matrix with the same
c     logical and physical dimensions as a, whose columns contain, on output, the normalized
c     eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
c====================================================================================================
      subroutine jacobi(a,n,np,d,v,nrot,ising)
      implicit none
      integer n,np,nrot,nmax
      real(8) a(np,np),d(np),v(np,np)
      parameter (nmax=500)
      integer i,ip,iq,j,ising
      real(8) c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)
      ising=0
      do ip=1,n  !initialize to the identity matrix.
         do iq=1,n
            v(ip,iq)=0.
         enddo
         v(ip,ip)=1.
      enddo 
      do ip=1,n
         b(ip)=a(ip,ip) !initialize b and d to the diagonal of a.
         d(ip)=b(ip)
         z(ip)=0.  !this vector will accumulate terms of the form tapq
      enddo        !as in equation (11.1.14).
      nrot=0
      do i=1,50
         sm=0.
         do ip=1,n-1  !sum off-diagonal elements.
         do iq=ip+1,n
            sm=sm+dabs(a(ip,iq))
         enddo 
         enddo 
         if(sm==0.)return   
         !===================!
         ! sucessful return  !
         !===================!
         if(i.lt.4)then       
            tresh=0.2*sm/n**2  !...on the first three sweeps.
         else
            tresh=0.           !...thereafter.
         endif
         do ip=1,n-1
         do iq=ip+1,n
            g=100.*dabs(a(ip,iq))
            !after four sweeps, skip the rotation if the off-diagonal element is small.
            if((i.gt.4).and.(dabs(d(ip))+g==dabs(d(ip)))
     *         .and.(dabs(d(iq))+g==dabs(d(iq))))then
               a(ip,iq)=0.
            else if(dabs(a(ip,iq)).gt.tresh)then
               h=d(iq)-d(ip)
               if(dabs(h)+g==dabs(h))then
                  t=a(ip,iq)/h          !t = 1/(2*theta)
               else
                  theta=0.5*h/a(ip,iq)  !equation (11.1.10).
                  t=1./(dabs(theta)+dsqrt(1.+theta**2))
                  if(theta.lt.0.)t=-t
               endif
               c=1./dsqrt(1+t**2)
               s=t*c
               tau=s/(1.+c)
               h=t*a(ip,iq)
               z(ip)=z(ip)-h
               z(iq)=z(iq)+h
               d(ip)=d(ip)-h
               d(iq)=d(iq)+h
               a(ip,iq)=0.
               do j=1,ip-1  !case of rotations 1 = j < p.
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo 
               do j=ip+1,iq-1 !case of rotations p < j < q.
                  g=a(ip,j)
                  h=a(j,iq)
                  a(ip,j)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo 
               do j=iq+1,n !case of rotations q < j = n.
                  g=a(ip,j)
                  h=a(iq,j)
                  a(ip,j)=g-s*(h+g*tau)
                  a(iq,j)=h+s*(g-h*tau)
               enddo
               do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*tau)
                  v(j,iq)=h+s*(g-h*tau)
               enddo 
               nrot=nrot+1
            endif
         enddo 
         enddo 
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip) !update d with the sum of tapq ,
            z(ip)=0.    !and reinitialize z.
         enddo
      enddo 

      write(*,*) 'too many iterations in jacobi, give up'
      ising=1

      return
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+                                                           +
c+     calculate stress due to dislocation line              +
c+                                                           +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine stress_dislocation(p1,p2,t,b,mu,nu,p,sm)
      implicit none
      integer i,j,k,l,m,n,ip
      real(8) p1(3),p2(3),p(3),t(3),b(3),mu,nu,sm(3,3),sv(6)
      real(8) r(2,3),rL(2,3),rP(2,3),rY(2,3)
      real(8) c_r(2),c_rL(2),c_rP(2),c_rY(2)
      real(8) stf(2,3,3,3)
      real(8) x1,x2,x3,mx1(3,3),mx2(3,3),vx(3)
      real(8) eij(3,3),eijk(3,3,3)
      real(8) c1,c2,c3,c4,c5
      real(8) z1,z2,z3,z4,z5,z6,z7

      eij(1,:)=[1,0,0]
      eij(2,:)=[0,1,0]
      eij(3,:)=[0,0,1]
      eijk=0
      eijk(1,2,3)= 1
      eijk(1,3,2)=-1
      eijk(2,1,3)=-1
      eijk(2,3,1)= 1
      eijk(3,1,2)= 1
      eijk(3,2,1)=-1

      r(1,:)=p-p1
      c_r(1)=dsqrt(sum(r(1,:)**2))
      c_rL(1)=dot_product(r(1,:),t)
      rP(1,:)=r(1,:)-c_rL(1)*t;          
      rY(1,:)=r(1,:)+c_r(1)*t
      c_rY(1)=dsqrt(sum(rY(1,:)**2))

      r(2,:)=p-p2
      c_r(2)=dsqrt(sum(r(2,:)**2))
      c_rL(2)=dot_product(r(2,:),t)
      rP(2,:)=r(2,:)-c_rL(2)*t;          
      rY(2,:)=r(2,:)+c_r(2)*t
      c_rY(2)=dsqrt(sum(rY(2,:)**2))

c      print*,c_r
c      print*,c_rL
c      print*,c_rY
c      read*

      stf=0
      do ip=1,2
         c1=mu/(3.14*c_rY(ip)**2)
         c2=1/(1-nu)
         vx=0
         do i=1,3
         do j=1,3
         do k=1,3
            vx(k)=vx(k)+eijk(k,i,j)*rY(ip,i)*t(j)
         enddo
         enddo
         enddo
         vx=vx/(2*(1-nu))
         c4=2/c_rY(ip)**2
         c5=c_rL(ip)/c_r(ip)
         do i=1,3
         do j=1,3
         do k=1,3
            z1=( dot_product(eijk(i,k,:),rY(ip,:))*t(j)
     &          +dot_product(eijk(j,k,:),rY(ip,:))*t(i))/2
            z2=( dot_product(eijk(i,k,:),t       )*rY(ip,j)
     &          +dot_product(eijk(j,k,:),t       )*rY(ip,i))/2
            z3=eij(i,j)+t(i)*t(j)
            z4=rP(ip,i)*rY(ip,j)+rP(ip,j)*rY(ip,i)
            z5=rY(ip,i)*rY(ip,j)
            stf(ip,i,j,k)=c1*(z1-c2*z2-vx(k)*(z3+c4*(z4+c5*z5)))         
         enddo
         enddo
         enddo
      enddo
      do i=1,3
      do j=1,3
         sm(i,j)=dot_product(stf(2,i,j,:)-stf(1,i,j,:),b)
      enddo
      enddo

      return
      end




c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Crystal Plasticity Fast Fourier Transformation Code         +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine smooth(nx,ny,nz,L_size,rve_X)
      implicit none
      integer i,j,k,k1,ising
      integer nx,ny,nz,ix,iy,iz,mx,my,mz
      integer, parameter :: ax=0
      integer, parameter :: ay=0
      integer, parameter :: az=0
      real(8) fct_smooth,L_size
      real(8) WSfft((nx+ax*2)*(ny+ay*2)*(nz+az*2)*2)
      real(8) rve_V((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_R((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_I((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_X(nx,ny,nz)
      real(8) x1,x2,x3,Fqc(3)

      mx=nx+ax*2
      my=ny+ay*2
      mz=nz+az*2

!      fct_smooth = 4*3.14**2*10.d0
!      fct_smooth = 4*3.14**2*20.d0
!      fct_smooth = 4*3.14**2*40.d0
!      fct_smooth = 4*3.14**2*100.d0
!      fct_smooth = 4*3.14**2*50.d0
!      fct_smooth = 4*3.14**2*30.d0

      fct_smooth = 4*3.14**2*60.d0
!      fct_smooth = 4*3.14**2*40.d0
!      fct_smooth = 4*3.14**2*40.d0*1.d-6

!      fct_smooth = 4*3.14**2*4.d19*L_size**2

c-----transfer stress to frequency space
      rve_V=0
      rve_V( ax+1:ax+nx, ay+1:ay+ny, az+1:az+nz ) = rve_X 

      WSfft=0               
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; WSfft(k1)=rve_V(ix,iy,iz)
         k1=k1+1; WSfft(k1)=0               
      enddo
      enddo
      enddo
      call fourn(WSfft,[mx,my,mz],3,1)
      WSfft=WSfft/(mx*my*mz)
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; rve_R(ix,iy,iz)=WSfft(k1) !--real part
         k1=k1+1; rve_I(ix,iy,iz)=WSfft(k1) !--image part
      enddo
      enddo
      enddo

c-----calculate equilibrium value in frequency space
      do ix=1,mx
      do iy=1,my
      do iz=1,mz
         if(ix<=mx/2+1)then
            Fqc(1) = (ix-1.)/mx
         else
            Fqc(1) = (ix-(mx+1.))/mx
         endif
         if(iy<=my/2+1)then
            Fqc(2) = (iy-1.)/my
         else
            Fqc(2) = (iy-(my+1.))/my
         endif
         if(iz<=mz/2+1)then
            Fqc(3) = (iz-1.)/mz
         else
            Fqc(3) = (iz-(mz+1.))/mz
         endif
         x1= 1 + fct_smooth*sum(Fqc**2)
         rve_R(ix,iy,iz) = rve_R(ix,iy,iz)/x1
         rve_I(ix,iy,iz) = rve_I(ix,iy,iz)/x1
      enddo
      enddo
      enddo

c-----transfer new value from frequency space to physical space
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; WSfft(k1)=rve_R(ix,iy,iz)
         k1=k1+1; WSfft(k1)=rve_I(ix,iy,iz)
      enddo
      enddo
      enddo
      call fourn(WSfft,[mx,my,mz], 3, -1)
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; rve_V(ix,iy,iz)=WSfft(k1) 
         k1=k1+1
      enddo
      enddo
      enddo

      rve_X = rve_V( ax+1:ax+nx, ay+1:ay+ny, az+1:az+nz ) 
c
      return
      end










ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine GAUSSJ_new(a,n,np,b,m,mp,ierr)

c     http://www.public.iastate.edu/~jba/Fortran/gaussj.fc  
c     Purpose: Solution of the system of linear equations AX = B by
c     Gauss-Jordan elimination, where A is a matrix of order N and B is
c     an N x M matrix.  On output A is replaced by its matrix inverse
c     and B is preplaced by the corresponding set of solution vectors.

c  Source: W.H. Press et al, "Numerical Recipes," 1989, p. 28.

c  Modifications: 
c     1. Double  precision.
c     2. Error parameter IERR included.  0 = no error. 1 = singular 
c        matrix encountered; no inverse is returned.

c  Prepared by J. Applequist, 8/17/91.

      implicit real*8(a-h,o-z)

c        Set largest anticipated value of N.

      parameter (nmax=500)
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
      ierr=0
      do 11 j=1,n
      ipiv(j)=0
 11   continue
      do 22 i=1,n
      big=0.d0
      do 13 j=1,n
      if (ipiv(j).ne.1) then
      do 12 k=1,n
      if (ipiv(k).eq.0) then
      if (dabs(a(j,k)).ge.big) then
      big=dabs(a(j,k))
      irow=j
      icol=k
      endif
      else if (ipiv(k).gt.1) then
      ierr=1
      return
      endif
 12   continue
      endif
 13   continue
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
      do 14 l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
 14   continue
      do 15 l=1,m
      dum=b(irow,l)
      b(irow,l)=b(icol,l)
      b(icol,l)=dum
 15   continue
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.d0) then
      ierr=1
      return
      endif
      pivinv=1.d0/a(icol,icol)
      a(icol,icol)=1.d0
      do 16 l=1,n
      a(icol,l)=a(icol,l)*pivinv
 16   continue
      do 17 l=1,m
      b(icol,l)=b(icol,l)*pivinv
 17   continue
      do 21 ll=1,n
      if (ll.ne.icol) then
      dum=a(ll,icol)
      a(ll,icol)=0.d0
      do 18 l=1,n
      a(ll,l)=a(ll,l)-a(icol,l)*dum
 18   continue
      do 19 l=1,m
      b(ll,l)=b(ll,l)-b(icol,l)*dum
 19   continue
      endif
 21   continue
 22   continue
      do 24 l=n,1,-1
      if (indxr(l).ne.indxc(l)) then
      do 23 k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
 23   continue
      endif
 24   continue
      return
      end



c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Set connection between FFT and FEM mesh                     +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine sub_FEMFFT_map(Tnel,Tngp,
     &                          Tnfx,Tnfy,Tnfz,
     &                          fem_xyz,
     &                          fft_xyz,
     &                          fft_Fqc,
     &                          fft_Inf,
     &                          fem_Inf,
     &                          detx,dety,detz)
      implicit none
c
      integer Tnel,Tngp,Tnfx,Tnfy,Tnfz
      real(8) detx,dety,detz
      integer fem_Inf (Tnel,Tngp,4)
      real(8) fem_xyz (Tnel,Tngp,3)
      integer fft_Inf (Tnfx,Tnfy,Tnfz,3)
      real(8) fft_Fqc (Tnfx,Tnfy,Tnfz,3)
      real(8) fft_xyz (Tnfx,Tnfy,Tnfz,3)
c
      integer i,j,k,l,m,n,ix1,ix2,ix3,ii,np,ising
      integer ix,iy,iz,i1,i2,j1,j2,iex,igx,idx,ip1,ip2
      real(8) mx33_1(3,3),mx33_2(3,3),mx33_3(3,3)
      real(8) BOX_x(4),BOX_y(4),BOX_z(4)
      real(8) Lenx,Leny,Lenz
      real(8) x1,x2,x3,dydx,ddyddx
      real(8) x_min,x_max,x_add
      real(8) y_min,y_max,y_add
      real(8) z_min,z_max,z_add
      real(8) x_tmp
c
      x_min =+1.d50  
      x_max =-1.d50
      y_min =+1.d50  
      y_max =-1.d50
      z_min =+1.d50  
      z_max =-1.d50
      do i=1,Tnel
      do j=1,Tngp
         if(x_min>fem_xyz(i,j,1)) x_min=fem_xyz(i,j,1)
         if(y_min>fem_xyz(i,j,2)) y_min=fem_xyz(i,j,2)
         if(z_min>fem_xyz(i,j,3)) z_min=fem_xyz(i,j,3)
         if(x_max<fem_xyz(i,j,1)) x_max=fem_xyz(i,j,1)
         if(y_max<fem_xyz(i,j,2)) y_max=fem_xyz(i,j,2)
         if(z_max<fem_xyz(i,j,3)) z_max=fem_xyz(i,j,3)
      enddo
      enddo
      if(x_max==x_min) x_max=x_min+1
      if(y_max==y_min) y_max=y_min+1
      if(z_max==z_min) z_max=z_min+1
      x_add=(x_max-x_min)*1.d-1
      y_add=(y_max-y_min)*1.d-1
      z_add=(z_max-z_min)*1.d-1
      BOX_x=[x_min-x_add,x_min,x_max,x_max+x_add]
      BOX_y=[y_min-y_add,y_min,y_max,y_max+y_add]
      BOX_z=[z_min-z_add,z_min,z_max,z_max+z_add]
      Lenx = BOX_x(4)-BOX_x(1)
      Leny = BOX_y(4)-BOX_y(1)
      Lenz = BOX_z(4)-BOX_z(1)
      detx = Lenx/Tnfx
      dety = Leny/Tnfy
      detz = Lenz/Tnfz

      !-----------------------------------------
      !Frequency and coordinite for FFT grids
      !-----------------------------------------
      fft_Fqc=0
      do ix=1,Tnfx
      do iy=1,Tnfy
      do iz=1,Tnfz
         if(ix<=Tnfx/2+1)then
            x1 = (ix-1.       )/Lenx
         else
            x1 = (ix-(Tnfx+1.))/Lenx
         endif
         if(iy<=Tnfy/2+1)then
            x2 = (iy-1.       )/Leny
         else
            x2 = (iy-(Tnfy+1.))/Leny
         endif
         if(iz<=Tnfz/2+1)then
            x3 = (iz-1.       )/Lenz
         else
            x3 = (iz-(Tnfz+1.))/Lenz
         endif
         fft_Fqc(ix,iy,iz,:) = [x1,x2,x3]
         fft_xyz(ix,iy,iz,1) = BOX_x(1) + detx*(ix-1)
         fft_xyz(ix,iy,iz,2) = BOX_y(1) + dety*(iy-1)
         fft_xyz(ix,iy,iz,3) = BOX_z(1) + detz*(iz-1)
      enddo
      enddo
      enddo

      !----------------------------------------------------
      ! Connection between FFT and FEM and BCs
      !----------------------------------------------------
      fft_Inf=0
      do ix=1,Tnfx
      do iy=1,Tnfy
      do iz=1,Tnfz
         x1=fft_xyz(ix,iy,iz,1)
         x2=fft_xyz(ix,iy,iz,2)
         x3=fft_xyz(ix,iy,iz,3)
         x_min=1.d50; iex=0; igx=0
         do i1=1,Tnel
         do i2=1,Tngp
         if(fem_Inf(i1,i2,1)/=0)then
            x_tmp=dsqrt(sum((fem_xyz(i1,i2,:)-[x1,x2,x3])**2))
            if(x_min>x_tmp)then
               x_min=x_tmp; iex=i1; igx=i2
            endif                
         endif
         enddo
         enddo
         idx=1
         if( x1<BOX_x(2) ) idx=2
         if( x1>BOX_x(3) ) idx=3
         if( x2<BOX_y(2) ) idx=4
         if( x2>BOX_y(3) ) idx=5
         if( x3<BOX_z(2) ) idx=6
         if( x3>BOX_z(3) ) idx=7
         fft_Inf(ix,iy,iz,1:3)=[iex,igx,idx]
      enddo         
      enddo         
      enddo

      do i1=1,Tnel
      do i2=1,Tngp
      if(fem_Inf(i1,i2,1)/=0)then
         x_min=1.d50; iex=0; igx=0
         do ix=1,Tnfx
         do iy=1,Tnfy
         do iz=1,Tnfz
            x1=fft_xyz(ix,iy,iz,1)
            x2=fft_xyz(ix,iy,iz,2)
            x3=fft_xyz(ix,iy,iz,3)
            if(fft_Inf(ix,iy,iz,3)==1)then
               x_tmp=dsqrt(sum((fem_xyz(i1,i2,:)-[x1,x2,x3])**2))
               if(x_min>x_tmp)then
                  x_min=x_tmp 
                  fem_Inf(i1,i2,2:4)=[ix,iy,iz]
               endif                
            endif
         enddo
         enddo
         enddo
      endif
      enddo
      enddo
c
      return
      end


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +  FFT Solver for P0+A*Y+B_i*dYdx_i+C*ddYddx_ii=0 with P.B.C    +
c +         input: P0, A, B, C, Y0                                +
c +         output: Y, dydx, ddYddx                               +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine sub_PDEsolve_FFT(Nftx,Nfty,Nftz,fft_Fqc,
     &                     fft_P0,fft_A,fft_B,fft_C,
     &                     fft_S0,
     &                     fft_S1,fft_V1,fft_M1)
      implicit none
      integer Nftx,Nfty,Nftz,Nele,Ngpt
      real(8) WSfft  (Nftx*Nfty*Nftz*2)
      real(8) fft_Fqc(Nftx,Nfty,Nftz,3)
      real(8) fft_R  (Nftx,Nfty,Nftz)
      real(8) fft_I  (Nftx,Nfty,Nftz)
      real(8) fft_Rp (Nftx,Nfty,Nftz)
      real(8) fft_Ip (Nftx,Nfty,Nftz)
      real(8) fft_P  (Nftx,Nfty,Nftz)
      real(8) fft_P0 (Nftx,Nfty,Nftz)
      real(8) fft_A  (Nftx,Nfty,Nftz)
      real(8) fft_B  (Nftx,Nfty,Nftz,3)
      real(8) fft_C  (Nftx,Nfty,Nftz)
      real(8) fft_Ap (Nftx,Nfty,Nftz)
      real(8) fft_Bp (Nftx,Nfty,Nftz,3)
      real(8) fft_Cp (Nftx,Nfty,Nftz)
      real(8) fft_S0 (Nftx,Nfty,Nftz)
      real(8) fft_S1 (Nftx,Nfty,Nftz)
      real(8) fft_V1 (Nftx,Nfty,Nftz,3)
      real(8) fft_M1 (Nftx,Nfty,Nftz,3,3)
      real(8) fft_S1_old (Nftx,Nfty,Nftz)

      real(8) A_avg,B_avg(3),C_avg
      real(8) rsd,Crsd,Fqc(3)
c
      integer i,j,k,k1,ising,i1,i2,i3,i4,j1,j2,j3,j4
      integer ix,iy,iz,Iloop,Nloop,ip,jp,TNP
      real(8) x1,x2,x3
c
      Nloop=100
      Crsd=1.d-8
      TNP=Nftx*Nfty*Nftz
c
c-----calculate the average parameters A, B, C
      A_avg=0
      B_avg=0
      C_avg=0
      do ix=1,Nftx
      do iy=1,Nfty
      do iz=1,Nftz
         A_avg = A_avg + fft_A(ix,iy,iz  )
         B_avg = B_avg + fft_B(ix,iy,iz,:)
         C_avg = C_avg + fft_C(ix,iy,iz  )
      enddo
      enddo
      enddo
      A_avg = A_avg/TNP !*50
      B_avg = B_avg/TNP !*50
      C_avg = C_avg/TNP !*50
c-----Parameters Ap, Bp, Cp for polarized field
      do ix=1,Nftx
      do iy=1,Nfty
      do iz=1,Nftz
         fft_Ap(ix,iy,iz  )=fft_A(ix,iy,iz  )-A_avg
         fft_Bp(ix,iy,iz,:)=fft_B(ix,iy,iz,:)-B_avg
         fft_Cp(ix,iy,iz  )=fft_C(ix,iy,iz  )-C_avg
      enddo
      enddo
      enddo


      Iloop=0
      fft_S1=fft_S0
      call Euc2Fqc(Nftx,Nfty,Nftz,fft_Fqc, fft_S1, fft_R,fft_I)
      do while(Iloop<Nloop) 
         Iloop=Iloop+1

c--------check convergence in hetergenious materials
         if(Iloop>1)then
            rsd=sum(dabs(fft_S1_old-fft_S1))/(Nftx*Nfty*Nftz)

            print '(2I5,10e15.5)',Iloop,Nloop,rsd,Crsd
c            read*

            if(rsd < Crsd) goto 102
         endif
         fft_S1_old=fft_S1

c--------1st, 2nd gradient of field in space
         do ip=1,3
            do ix=1,Nftx
            do iy=1,Nfty
            do iz=1,Nftz
               x1=2*3.14*fft_fqc(ix,iy,iz,ip)
               fft_Rp(ix,iy,iz) = x1*fft_I(ix,iy,iz)
               fft_Ip(ix,iy,iz) =-x1*fft_R(ix,iy,iz)
            enddo
            enddo
            enddo
            call Fqc2Euc(Nftx,Nfty,Nftz,fft_Fqc,fft_Rp,fft_Ip,
     &                   fft_V1(:,:,:,ip))
         enddo
         do ip=1,3
         do jp=1,3
            do ix=1,Nftx
            do iy=1,Nfty
            do iz=1,Nftz
               x1=4*9.86*fft_fqc(ix,iy,iz,ip)*fft_fqc(ix,iy,iz,jp)
               fft_Rp(ix,iy,iz) =-x1*fft_R(ix,iy,iz)
               fft_Ip(ix,iy,iz) =-x1*fft_I(ix,iy,iz)
            enddo
            enddo
            enddo
            call Fqc2Euc(Nftx,Nfty,Nftz,fft_Fqc,fft_Rp,fft_Ip,
     &                   fft_M1(:,:,:,ip,jp))
         enddo
         enddo

c--------polarized field
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            fft_P(ix,iy,iz) = fft_P0(ix,iy,iz) +
     &                        fft_Ap(ix,iy,iz)  *fft_S1(ix,iy,iz) + 
     &                        fft_Bp(ix,iy,iz,1)*fft_V1(ix,iy,iz,1) +
     &                        fft_Bp(ix,iy,iz,2)*fft_V1(ix,iy,iz,2) +
     &                        fft_Bp(ix,iy,iz,3)*fft_V1(ix,iy,iz,3) +
     &                        fft_Cp(ix,iy,iz)  *fft_M1(ix,iy,iz,1,1) +
     &                        fft_Cp(ix,iy,iz)  *fft_M1(ix,iy,iz,2,2) +
     &                        fft_Cp(ix,iy,iz)  *fft_M1(ix,iy,iz,3,3) 
         enddo
         enddo
         enddo
         call Euc2Fqc(Nftx,Nfty,Nftz,fft_Fqc, fft_P, fft_Rp,fft_Ip)

c--------calculate new field
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            Fqc = fft_fqc(ix,iy,iz,:)
            x1  = A_avg-4*9.86*C_avg*sum(Fqc**2)
            x2  = 2*3.14*dot_product(B_avg,Fqc)
            x3  = x1**2+x2**2
            x1  = x1/x3
            x2  = x2/x3
            fft_R(ix,iy,iz) = -fft_Rp(ix,iy,iz)*x1 + fft_Ip(ix,iy,iz)*x2
            fft_I(ix,iy,iz) = -fft_Rp(ix,iy,iz)*x2 - fft_Ip(ix,iy,iz)*x1
         enddo
         enddo
         enddo
         call Fqc2Euc(Nftx,Nfty,Nftz,fft_Fqc, fft_R,fft_I, fft_S1)

      enddo

102   continue
      return
      end




c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   FFT transformation: from space to frequency                 +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine Euc2Fqc(Nftx,Nfty,Nftz,fft_Fqc,fft_V,fft_R,fft_I)
      implicit none
      integer Nftx,Nfty,Nftz
      real(8) WSfft  (Nftx*Nfty*Nftz*2)
      real(8) fft_Fqc(Nftx,Nfty,Nftz,3)
      real(8) fft_R (Nftx,Nfty,Nftz)
      real(8) fft_I (Nftx,Nfty,Nftz)
      real(8) fft_V (Nftx,Nfty,Nftz)
      integer i,j,k,ix,iy,iz,ip
      real(8) x1,x2,x3
      WSfft=0               
      ip=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         ip=ip+1; WSfft(ip)=fft_V(ix,iy,iz)
         ip=ip+1; WSfft(ip)=0               
      enddo
      enddo
      enddo
      call fourn(WSfft,[Nftx,Nfty,Nftz],3,1)
      WSfft=WSfft/(Nftx*Nfty*Nftz)
      ip=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         ip=ip+1; fft_R(ix,iy,iz)=WSfft(ip) !--real part
         ip=ip+1; fft_I(ix,iy,iz)=WSfft(ip) !--image part
      enddo
      enddo
      enddo
      return
      end

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   FFT transformation: from frequency to space                 +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine Fqc2Euc(Nftx,Nfty,Nftz,fft_Fqc,fft_R,fft_I,fft_V)
      implicit none
      integer Nftx,Nfty,Nftz
      real(8) WSfft  (Nftx*Nfty*Nftz*2)
      real(8) fft_Fqc(Nftx,Nfty,Nftz,3)
      real(8) fft_R (Nftx,Nfty,Nftz)
      real(8) fft_I (Nftx,Nfty,Nftz)
      real(8) fft_V (Nftx,Nfty,Nftz)
      integer i,j,k,ix,iy,iz,ip
      real(8) x1,x2,x3
      WSfft=0               
      ip=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         ip=ip+1; WSfft(ip)=fft_R(ix,iy,iz)
         ip=ip+1; WSfft(ip)=fft_I(ix,iy,iz)
      enddo
      enddo
      enddo
      call fourn(WSfft,[Nftx,Nfty,Nftz], 3, -1)
      ip=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         ip=ip+1; fft_V(ix,iy,iz)=WSfft(ip) 
         ip=ip+1
      enddo
      enddo
      enddo
      return
      end


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Calculate first and second gradient by finite difference    +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine cal_grad_FD(fft_V1,
     &                      Nftx,Nfty,Nftz,
     &                      detx,dety,detz,
     &                      fft_dVdx,
     &                      fft_ddVddx)
      implicit none
      integer Nftx,Nfty,Nftz
      real(8) fft_V2    (Nftx,Nfty,Nftz)
      real(8) fft_V1    (Nftx,Nfty,Nftz)
      real(8) fft_dVdx  (Nftx,Nfty,Nftz,3)
      real(8) fft_ddVddx(Nftx,Nfty,Nftz,3,3)
      real(8) detx,dety,detz
c
      integer i,j,k
      integer ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2
      real(8) x1,x2,x3

c--------------------
      fft_dVdx=0
c--------------------
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         if(ix==1)then
            ix1=Nftx
            ix2=ix+1
         elseif(ix==Nftx)then
            ix1=Nftx-1
            ix2=1
         else
            ix1=ix-1
            ix2=ix+1
         endif
         if(iy==1)then
            iy1=Nfty
            iy2=iy+1
         elseif(iy==Nfty)then
            iy1=Nfty-1
            iy2=1
         else
            iy1=iy-1
            iy2=iy+1
         endif
         if(iz==1)then
            iz1=Nftz
            iz2=iz+1
         elseif(iz==Nftz)then
            iz1=Nftz-1
            iz2=1
         else
            iz1=iz-1
            iz2=iz+1
         endif

         if(Nftx>2)then
            x1=fft_V1(ix1, iy, iz  )
            x2=fft_V1(ix2, iy, iz  )
            fft_dVdx (ix,  iy, iz,1)=(x2-x1)/2/detx
         endif
         if(Nfty>2)then
            x1=fft_V1(ix, iy1, iz  )
            x2=fft_V1(ix, iy2, iz  )
            fft_dVdx (ix, iy,  iz,2)=(x2-x1)/2/dety
         endif
         if(Nftz>2)then
            x1=fft_V1(ix, iy, iz1  )
            x2=fft_V1(ix, iy, iz2  )
            fft_dVdx (ix, iy, iz, 3)=(x2-x1)/2/detz
         endif

      enddo
      enddo
      enddo

c-----------------------
      fft_ddVddx=0
c-----------------------
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         if(ix==1)then
            ix1=Nftx
            ix2=ix+1
         elseif(ix==Nftx)then
            ix1=Nftx-1
            ix2=1
         else
            ix1=ix-1
            ix2=ix+1
         endif
         if(iy==1)then
            iy1=Nfty
            iy2=iy+1
         elseif(iy==Nfty)then
            iy1=Nfty-1
            iy2=1
         else
            iy1=iy-1
            iy2=iy+1
         endif
         if(iz==1)then
            iz1=Nftz
            iz2=iz+1
         elseif(iz==Nftz)then
            iz1=Nftz-1
            iz2=1
         else
            iz1=iz-1
            iz2=iz+1
         endif

         do i=1,3
            if(Nftx>2)then
               x1=fft_dVdx(ix1, iy, iz, i  )
               x2=fft_dVdx(ix2, iy, iz, i  )
               fft_ddVddx (ix,  iy, iz, i,1)=(x2-x1)/2/detx
            endif
            if(Nfty>2)then
               x1=fft_dVdx(ix, iy1, iz, i  )
               x2=fft_dVdx(ix, iy2, iz, i  )
               fft_ddVddx (ix,  iy, iz, i,2)=(x2-x1)/2/dety
            endif
            if(Nftz>2)then
               x1=fft_dVdx(ix, iy, iz1, i  )
               x2=fft_dVdx(ix, iy, iz2, i  )
               fft_ddVddx (ix, iy, iz,  i,3)=(x2-x1)/2/detz
            endif
         enddo
      enddo
      enddo
      enddo
c
      return
      end


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Output date to ovito                                        +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine output_ovito(iframe,
     &                        fft_S0,
     &                        fft_V1,
     &                        fft_dVdx,
     &                        fft_ddVddx,
     &                        Nftx,Nfty,Nftz)
         implicit none
         integer Nftx,Nfty,Nftz,iframe
         real(8) fft_S0(Nftx,Nfty,Nftz)
         real(8) fft_V1(Nftx,Nfty,Nftz)
         real(8) fft_dVdx(Nftx,Nfty,Nftz,3)
         real(8) fft_ddVddx(Nftx,Nfty,Nftz,3,3)
c
         integer i,j,k
         integer ix,iy,iz
         real(8) x1,x2,x3,v3(3)
         character supp*4 
c      
      
         write(supp,'(i4.4)') iframe
         open(923,file='smth_'//supp//'.chkpt',status='unknown')
c         open(923,file='smth.chkpt')
c
         write(923,'("#F A 1 1 1 3 1 1")')
         write(923,'("#C Np Ng Nl x y z v0 v1 g1 g2 g3")')
         write(923,'("#X 1 0 0")')
         write(923,'("#Y 0 1 0")')
         write(923,'("#Z 0 0 1")')
         write(923,'("#E")')
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            v3=[ix,iy,iz]*2.d0  
            write(923,'(3I10,100e12.3)') 
     &         1,1,1,v3,
     &         fft_S0(ix,iy,iz),
     &         fft_V1(ix,iy,iz),
     &         fft_dVdx(ix,iy,iz,1:3)
         enddo
         enddo
         enddo
         close(923)
         return
      endsubroutine



c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   calculate Eschby matrix for 4vlm RVE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine Cal_Amatrix(nx,ny,nz,
     &                          pnx,pny,pnz,
     &                          c11_m,c12_m,c44_m,
     &                          c11_p,c12_p,c44_p,
     &                          IB1,IB2,
     &                          M_4vlm)
         implicit none
         character(len=200) :: Path=''         
         integer i1,j1,l1,m1,i,j,nx,ny,nz,pnx,pny,pnz
         integer IB1(9),IB2(9)
         real(8) c11(3),c12(3),c44(3)
         real(8) sv_4vlm(24),egv_4vlm(24),M_4vlm(24,24)
         real(8) c11_m,c12_m,c44_m
         real(8) c11_p,c12_p,c44_p
         real(8) det
c
         call getcwd(Path)        
         det=1.d-3
c
         print*
         print*, 'calcualte A-tensors for internal stresses:'
         print*
         do i=1,24
            print '("pertibation for componite:",I5)',i
            egv_4vlm=0.d0
            egv_4vlm(i)=det
            call Cal_INTstress_FFT(nx,ny,nz,
     &                      pnx,pny,pnz,
     &                      c11_m,c12_m,c44_m,
     &                      c11_p,c12_p,c44_p,
     &                      IB1,IB2,
     &                      egv_4vlm,
     &                      sv_4vlm)
            M_4vlm(:,i)=sv_4vlm/det
         enddo 

c
c         egv_4vlm=0.d0
c         egv_4vlm(19:21)=-0.005
c         call Cal_INTstress_FFT(nx,ny,nz,
c     &                      pnx,pny,pnz,
c     &                      c11_m,c12_m,c44_m,
c     &                      c11_p,c12_p,c44_p,
c     &                      IB1,IB2,
c     &                      egv_4vlm,
c     &                      sv_4vlm)
c         print*,'calclated by FFT'
c         call pv(sv_4vlm( 1: 6),6)
c         call pv(sv_4vlm( 7:12),6)
c         call pv(sv_4vlm(13:18),6)
c         call pv(sv_4vlm(19:24),6)
ccc
c         egv_4vlm=0.d0
c         egv_4vlm(19:21)=-0.005
c         sv_4vlm=matmul(M_4vlm,egv_4vlm)
c         print*,'calclated by Amatrix'
c         call pv(sv_4vlm( 1: 6),6)
c         call pv(sv_4vlm( 7:12),6)
c         call pv(sv_4vlm(13:18),6)
c         call pv(sv_4vlm(19:24),6)
c         read*
c

         return
         end

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   calculate internal stress of eigen strain for 4vlm RVE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine Cal_INTstress_FFT(nx,ny,nz,
     &                                pnx,pny,pnz,
     &                                c11_m,c12_m,c44_m,
     &                                c11_p,c12_p,c44_p,
     &                                IB1,IB2,
     &                                egv_4vlm,
     &                                sv_4vlm)
         implicit none
         integer nx,ny,nz,ix,iy,iz,ip,i,j,k,l,k1,m,npp,jj
         integer ii,i1,j1,l1,m1,ising,np_fft,IP_fft,pnx,pny,pnz
         real(8) epsfxR(6),epsfyR(6),epsfzR(6),epse2R(6),epse3R(6)
         real(8) eps0R(6),epscR(6),epse1R(6)
         real(8) c11(3),c12(3),c44(3),stf66_spc(6,6,nx,ny,nz)
         real(8) stf66(6,6,3),stf99(9,9,3),stf3333(3,3,3,3,3)
         real(8) cpc66_fqc(6,6,nx,ny,nz),freq_all(3,nx,ny,nz)
         real(8) epsR(6,nx,ny,nz),sigR(6,nx,ny,nz)         
         real(8) mx66(6,6),mx99(9,9),fqv3(3),v3(3),sigRtot(6)
         real(8) mx33_1(3,3),mx33_2(3,3),WSfft(nx*ny*nz*2)
         real(8) absvectorF,csR_avg(6),cpc33_fqc(3,3,nx,ny,nz)
         real(8) CPC6_Fqc(6,6),rve_vxnI(6,nx,ny,nz),rve_DvFR(6,nx,ny,nz)
         real(8) rve_DvFI(6,nx,ny,nz),rve_LDv(6,nx,ny,nz)               
         real(8) vx6_1(6),vx6_2(6),fqv3_all(3,nx,ny,nz),Fqc(3)
         real(8) CvgR_glb,rsdS,CrsdS,rve_vxnR(6,nx,ny,nz),x1,x2,x3,x4
         real(8) uni(3,3),avgsig(7,6),avgeps(7,6)
         real(8) sv_4vlm(24),egv_4vlm(24)
         real(8) px0,py0,pz0,px1,py1,pz1
         real(8) sigtot(7,6),epstot(7,6),vv3(3) 
         real(8) ev_x(6),ev_y(6),ev_z(6)
         real(8) ev_xy(6),ev_xz(6),ev_yz(6)
         real(8) ev_xyz(6),ev_p(6)
         real(8) sv_x(6),sv_y(6),sv_z(6)
         real(8) sv_xy(6),sv_xz(6),sv_yz(6)
         real(8) sv_xyz(6),sv_p(6)
         integer ii_p,ii_x,ii_y,ii_z,ii_xy,ii_xz,ii_yz,ii_xyz
         real(8) c11_m,c12_m,c44_m
         real(8) c11_p,c12_p,c44_p
         real(8) Astrain,Pstrain
         integer IB1(9),IB2(9),IJB(9,9)
c-----------------------------------------------------
         Astrain=sum(dabs(egv_4vlm))/24
c         ev_p   = egv_4vlm( 1: 6)
c         ev_x   = egv_4vlm( 7:12)
c         ev_y   = egv_4vlm(13:18)
c         ev_z   = egv_4vlm(19:24)

         ev_x   = egv_4vlm( 1: 6)
         ev_y   = egv_4vlm( 7:12)
         ev_z   = egv_4vlm(13:18)
         ev_p   = egv_4vlm(19:24)

         ev_xy  = (ev_x + ev_y)/2
         ev_xz  = (ev_x + ev_z)/2
         ev_yz  = (ev_y + ev_z)/2
         ev_xyz = (ev_x + ev_y + ev_z)/3
         c11=[ c11_m, c11_p, (c11_m+c11_p)/2 ]
         c12=[ c12_m, c12_p, (c12_m+c12_p)/2 ]
         c44=[ c44_m, c44_p, (c44_m+c44_p)/2 ]
         do i=1,9 
            IJB(IB1(i),IB2(i))=i
         enddo

c        Allocation
         stf66(:,:,:) = 0 
         stf99(:,:,:) = 0 
         stf3333(:,:,:,:,:) = 0 
         stf66_spc(:,:,:,:,:) = 0
         cpc33_fqc(:,:,:,:,:) = 0
         cpc66_fqc(:,:,:,:,:) = 0 
         fqv3_all(:,:,:,:) = 0
         epsR(:,:,:,:) = 0
         sigR(:,:,:,:) = 0
         fqv3(:)=0

         mx66(:,:)=0
         mx99(:,:)=0
         mx33_1(:,:)=0
         mx33_2(:,:)=0

         sigRtot=0
         sigtot=0
         epstot=0
         avgsig=0
         avgeps=0
         sv_4vlm=0

         do ip = 1,3
         stf99(1,1:6,ip) = [c11(ip),c12(ip),c12(ip),0.d0,0.d0,0.d0]
         stf99(2,1:6,ip) = [c12(ip),c11(ip),c12(ip),0.d0,0.d0,0.d0]
         stf99(3,1:6,ip) = [c12(ip),c12(ip),c11(ip),0.d0,0.d0,0.d0]
         stf99(4,1:6,ip) = [0.d0, 0.d0, 0.d0, c44(ip), 0.d0, 0.d0]
         stf99(5,1:6,ip) = [0.d0, 0.d0, 0.d0, 0.d0, c44(ip), 0.d0]
         stf99(6,1:6,ip) = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0, c44(ip)]
         stf99(1:6,7:9,ip) = stf99(1:6,4:6,ip)
         stf99(7:9,1:9,ip) = stf99(4:6,1:9,ip)
         stf66(1:6,1:3,ip) = stf99(1:6,1:3,ip)
         stf66(1:6,4:6,ip) = stf99(1:6,4:6,ip)*2
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
            stf3333(i,j,k,l,ip) = stf99(IJB(i,j),IJB(k,l),ip)
         enddo
         enddo
         enddo
         enddo 
         enddo
        
c ------------------------------------
c        microstructure
c ------------------------------------
         do ix=1,nx
         do iy=1,ny
         do iz=1,nz                                                   

c--Set up elastic constants
            if (ix<pnx .and. iy<pny .and. iz<pnz) then
               stf66_spc(1:6,1:6,ix,iy,iz) = stf66(1:6,1:6,2)
            else
               stf66_spc(1:6,1:6,ix,iy,iz) = stf66(1:6,1:6,1)                
            endif
               
c--Set up frequency space 
            if(ix<=nx/2+1)then
               fqv3_all(1,ix,iy,iz) = (ix-1.)/nx
            else
               fqv3_all(1,ix,iy,iz) = (ix-(nx+1.))/nx
            endif
            if(iy<=ny/2+1)then
               fqv3_all(2,ix,iy,iz) = (iy-1.)/ny
            else
               fqv3_all(2,ix,iy,iz) = (iy-(ny+1.))/ny
            endif
            if(iz<=nz/2+1)then
               fqv3_all(3,ix,iy,iz) = (iz-1.)/nz
            else
               fqv3_all(3,ix,iy,iz) = (iz-(nz+1.))/nz
            endif

c--Complience tensor in frequency space
            fqv3(:)=fqv3_all(:,ix,iy,iz)
            absvectorF = dsqrt(fqv3(1)**2+fqv3(2)**2+fqv3(3)**2)
            if(absvectorF/=0.0) then
               mx33_1(:,:)=0
               do i=1,3
               do j=1,3
               do k=1,3
               do l=1,3
                  mx33_1(i,j)=mx33_1(i,j)+stf3333(i,k,j,l,3)*
     &            fqv3(k)*fqv3(l)
               enddo
               enddo
               enddo
               enddo

               call gaussj(mx33_1,3,mx33_2,ising)
               if(ising/=0) stop    

               do i=1,3
               do j=1,3
               do k=1,3
               do l=1,3
               mx99(IJB(i,j),IJB(k,l))=(mx33_2(j,k)*fqv3(i)*fqv3(l)+
     &         mx33_2(i,k)*fqv3(j)*fqv3(l)+mx33_2(j,l)*fqv3(i)*fqv3(k)+
     &         mx33_2(i,l)*fqv3(j)*fqv3(k))/4
               enddo               
               enddo
               enddo
               enddo
            endif
            cpc33_fqc(:,:,ix,iy,iz)=mx33_2(:,:)
            cpc66_fqc(1:6,1:3,ix,iy,iz)=mx99(1:6,1:3);
            cpc66_fqc(1:6,4:6,ix,iy,iz)=mx99(1:6,4:6)*2;

c--Eigen strain
            if(ix<=pnx.and.iy<=pny.and.iz<=pnz) epsR(:,ix,iy,iz)=-ev_p
            if(ix<=pnx.and.iy<=pny.and.iz> pnz) epsR(:,ix,iy,iz)=-ev_z
            if(ix<=pnx.and.iy> pny.and.iz<=pnz) epsR(:,ix,iy,iz)=-ev_y
            if(ix<=pnx.and.iy> pny.and.iz> pnz) epsR(:,ix,iy,iz)=-ev_yz
            if(ix> pnx.and.iy<=pny.and.iz<=pnz) epsR(:,ix,iy,iz)=-ev_x
            if(ix> pnx.and.iy<=pny.and.iz> pnz) epsR(:,ix,iy,iz)=-ev_xz
            if(ix> pnx.and.iy> pny.and.iz<=pnz) epsR(:,ix,iy,iz)=-ev_xy
            if(ix> pnx.and.iy> pny.and.iz> pnz) epsR(:,ix,iy,iz)=-ev_xyz

         enddo
         enddo
         enddo

c -----------------------------------------
c        Start of FFT Iteration Procedure
c -----------------------------------------
         np_fft   = 100          ! max FFT loop number
         CvgR_glb = 1.d-03       ! cvg ratio for global FFT loop
         rve_vxnR = 0
         rve_vxnI = 0
         rve_DvFR = 0
         rve_DvFI = 0
         rve_LDv  = 0
         rsdS     = 1.d50
         CrsdS    = CvgR_glb*rsdS
         
c--------begin of fft loop       
 
         IP_fft=0
1110     IP_fft=IP_fft+1

c--------calcualte stress due to deformation
         do ix=1,nx
         do iy=1,ny
         do iz=1,nz
            sigR(:,ix,iy,iz)=matmul(stf66_spc(:,:,ix,iy,iz),
     &      epsR(:,ix,iy,iz))                   
         enddo
         enddo
         enddo

c--------transfer stress to frequency space
         do ii=1,6
            k1=0
            do iz=1,nz
            do iy=1,ny
            do ix=1,nx
               k1=k1+1; WSfft(k1)=sigR(ii,ix,iy,iz)      !--real part
               k1=k1+1; WSfft(k1)=0                      !--image part
            enddo
            enddo
            enddo
            call fourn(WSfft,[nx,ny,nz],3,1)

            WSfft=WSfft/(nx*ny*nz)

            k1=0
            do iz=1,nz
            do iy=1,ny
            do ix=1,nx
               k1=k1+1; rve_vxnR(ii,ix,iy,iz)=WSfft(k1) !--real part
               k1=k1+1; rve_vxnI(ii,ix,iy,iz)=WSfft(k1) !--image part
            enddo
            enddo
            enddo
         enddo
         
c--------check equilibrium in the frequency space
         rsdS=0
         do ix=1,nx
         do iy=1,ny
         do iz=1,nz
            Fqc=fqv3_all(:,ix,iy,iz)
            do i=1,9
               if(i<=6) j=i
               if(i>6 ) j=i-3
               mx33_1(IB1(i),IB2(i))=rve_vxnR(j,ix,iy,iz)

               mx33_2(IB1(i),IB2(i))=rve_vxnI(j,ix,iy,iz)
            enddo
            x1=sum(dabs(matmul(Mx33_1,Fqc)))/(nx*ny*nz)
            x2=sum(dabs(matmul(Mx33_2,Fqc)))/(nx*ny*nz)
            rsdS = rsdS + x1 + x2
         enddo
         enddo
         enddo

         CrsdS=CvgR_glb*1.d0
         

c--------calculate new polarization strain by green function in frequency space
         do ix=1,nx
         do iy=1,ny
         do iz=1,nz
            Fqc=fqv3_all(:,ix,iy,iz)
            if(Fqc(1)==0 .and. Fqc(2)==0 .and. Fqc(3)==0)then
               rve_DvFR(:,ix,iy,iz) = 0
               rve_DvFI(:,ix,iy,iz) = 0
            else
               CPC6_Fqc=cpc66_fqc(:,:,ix,iy,iz)
               vx6_1=rve_vxnR(:,ix,iy,iz)
               vx6_2=rve_vxnI(:,ix,iy,iz)
               rve_DvFR(:,ix,iy,iz) = -matmul(CPC6_Fqc,vx6_1)
               rve_DvFI(:,ix,iy,iz) = -matmul(CPC6_Fqc,vx6_2)
            endif
         enddo
         enddo
         enddo

c--------transfer new polarization strain from frequency space to physical space
         do m=1,6
            k1=0
            do iz=1,nz
            do iy=1,ny
            do ix=1,nx
               k1=k1+1; WSfft(k1)=rve_DvFR(m,ix,iy,iz)
               k1=k1+1; WSfft(k1)=rve_DvFI(m,ix,iy,iz)
            enddo
            enddo
            enddo
            call fourn(WSfft,[nx,ny,nz], 3, -1)
            k1=0
            do iz=1,nz
            do iy=1,ny
            do ix=1,nx
               
               k1=k1+1; rve_LDv(m,ix,iy,iz)=WSfft(k1) 
               k1=k1+1
            enddo
            enddo
            enddo
         enddo

c--------push forword to current configuration, add boundary condition
         do iz=1,nz
         do iy=1,ny
         do ix=1,nx
            epsR(:,ix,iy,iz)=epsR(:,ix,iy,iz)+rve_LDv(:,ix,iy,iz)
         enddo
         enddo
         enddo


c         print*,IP_fft,np_fft,rsdS,CrsdS,x1
c         if ( IP_fft<np_fft .and. rsdS>CrsdS ) go to 1110

         Pstrain=sum(dabs(rve_LDv))/(nx*ny*nz)
c         print*,IP_fft,np_fft,Pstrain,Astrain*CvgR_glb
c         read*
         if ( IP_fft<np_fft .and. Pstrain>Astrain*CvgR_glb ) go to 1110

c -----------------------------------------
c        End of FFT Iteration Procedure
c -----------------------------------------

         do ix=1,nx
         do iy=1,ny
         do iz=1,nz     
            sigRtot=sigRtot+sigR(:,ix,iy,iz)
         enddo
         enddo
         enddo 
         do ix=1,nx
         do iy=1,ny
         do iz=1,nz
            sigR(:,ix,iy,iz)=sigR(:,ix,iy,iz)-sigRtot(:)/(nx*ny*nz)     
         enddo
         enddo
         enddo 
          
         sv_p=0
         sv_x=0
         sv_y=0
         sv_z=0
         sv_xy=0
         sv_xz=0
         sv_yz=0
         sv_xyz=0

         ii_p=1
         ii_x=1
         ii_y=1
         ii_z=1
         ii_xy=1
         ii_xz=1
         ii_yz=1
         ii_xyz=1
         do ix=1,nx
         do iy=1,ny
         do iz=1,nz
            if(ix<=pnx.and.iy<=pny.and.iz<=pnz) then
               sv_p=sv_p+sigR(:,ix,iy,iz); ii_p=ii_p+1
            endif
            if(ix<=pnx.and.iy<=pny.and.iz> pnz)then
               sv_z=sv_z+sigR(:,ix,iy,iz); ii_z=ii_z+1
            endif
            if(ix<=pnx.and.iy> pny.and.iz<=pnz)then
               sv_y=sv_y+sigR(:,ix,iy,iz); ii_y=ii_y+1
            endif
            if(ix<=pnx.and.iy> pny.and.iz> pnz)then
               sv_yz=sv_yz+sigR(:,ix,iy,iz); ii_yz=ii_yz+1
            endif
            if(ix> pnx.and.iy<=pny.and.iz<=pnz)then
               sv_x=sv_x+sigR(:,ix,iy,iz); ii_x=ii_x+1
            endif
            if(ix> pnx.and.iy<=pny.and.iz> pnz)then
               sv_xz=sv_xz+sigR(:,ix,iy,iz); ii_xz=ii_xz+1
            endif
            if(ix> pnx.and.iy> pny.and.iz<=pnz)then
               sv_xy=sv_xy+sigR(:,ix,iy,iz); ii_xy=ii_xy+1
            endif
            if(ix> pnx.and.iy> pny.and.iz> pnz)then
               sv_xyz=sv_xyz+sigR(:,ix,iy,iz); ii_xyz=ii_xyz+1
            endif
        enddo
        enddo
        enddo

        sv_p   = sv_p   / ii_p
        sv_x   = sv_x   / ii_x
        sv_y   = sv_y   / ii_y
        sv_z   = sv_z   / ii_z
        sv_xy  = sv_xy  / ii_xy
        sv_xz  = sv_xz  / ii_xz
        sv_yz  = sv_yz  / ii_yz
        sv_xyz = sv_xyz / ii_xyz

        sv_4vlm( 1: 6)=sv_x  
        sv_4vlm( 7:12)=sv_y  
        sv_4vlm(13:18)=sv_z  
        sv_4vlm(19:24)=sv_p
       return
       end
c___________________________________       
       subroutine upsdv4_6t0(v,upsdv4_6t)
       
       implicit none
       real*8 upsdv4_6t(3,3,3,3),v(6,6)
       real*8, parameter:: r2 = 1.414213562373095d0
       upsdv4_6t(1,1,1,1) = v(1,1)
       upsdv4_6t(1,1,2,2) = v(1,2)
       upsdv4_6t(1,1,3,3) = v(1,3)
       upsdv4_6t(1,1,1,2) = v(1,4) / r2
       upsdv4_6t(1,1,2,3) = v(1,5) / r2
       upsdv4_6t(1,1,1,3) = v(1,6) / r2
       upsdv4_6t(2,2,1,1) = v(2,1)
       upsdv4_6t(2,2,2,2) = v(2,2)
       upsdv4_6t(2,2,3,3) = v(2,3)
       upsdv4_6t(2,2,1,2) = v(2,4) / r2
       upsdv4_6t(2,2,2,3) = v(2,5) / r2
       upsdv4_6t(2,2,1,3) = v(2,6) / r2
       upsdv4_6t(3,3,1,1) = v(3,1)
       upsdv4_6t(3,3,2,2) = v(3,2)
       upsdv4_6t(3,3,3,3) = v(3,3)
       upsdv4_6t(3,3,1,2) = v(3,4) / r2
       upsdv4_6t(3,3,2,3) = v(3,5) / r2
       upsdv4_6t(3,3,1,3) = v(3,6) / r2
       upsdv4_6t(1,2,1,1) = v(4,1)
       upsdv4_6t(1,2,2,2) = v(4,2)
       upsdv4_6t(1,2,3,3) = v(4,3)
       upsdv4_6t(1,2,1,2) = v(4,4) / 2.d0
       upsdv4_6t(1,2,2,3) = v(4,5) / 2.d0
       upsdv4_6t(1,2,1,3) = v(4,6) / 2.d0
       upsdv4_6t(2,3,1,1) = v(5,1)
       upsdv4_6t(2,3,2,2) = v(5,2)
       upsdv4_6t(2,3,3,3) = v(5,3)
       upsdv4_6t(2,3,1,2) = v(5,4) / 2.d0
       upsdv4_6t(2,3,2,3) = v(5,5) / 2.d0
       upsdv4_6t(2,3,1,3) = v(5,6) / 2.d0
       upsdv4_6t(1,3,1,1) = v(6,1)
       upsdv4_6t(1,3,2,2) = v(6,2)
       upsdv4_6t(1,3,3,3) = v(6,3)
       upsdv4_6t(1,3,1,2) = v(6,4) / 2.d0
       upsdv4_6t(1,3,2,3) = v(6,5) / 2.d0
       upsdv4_6t(1,3,1,3) = v(6,6) / 2.d0
       !
       upsdv4_6t(2,1,:,:) = upsdv4_6t(1,2,:,:) 
       upsdv4_6t(3,2,:,:) = upsdv4_6t(2,3,:,:) 
       upsdv4_6t(3,1,:,:) = upsdv4_6t(1,3,:,:) 
       !
       upsdv4_6t(:,:,2,1) = upsdv4_6t(:,:,1,2) 
       upsdv4_6t(:,:,3,2) = upsdv4_6t(:,:,2,3) 
       upsdv4_6t(:,:,3,1) = upsdv4_6t(:,:,1,3) 
       !
       return
       end
!_____________________________________________
       subroutine psdv4_6t0(v,psdv4_6t)
       
       implicit none
       real*8 psdv4_6t(6,6),v(3,3,3,3)
       real*8, parameter:: r2 = 1.414213562373095d0
       psdv4_6t(1,1) = v(1,1,1,1)
       psdv4_6t(1,2) = v(1,1,2,2)
       psdv4_6t(1,3) = v(1,1,3,3)
       psdv4_6t(1,4) = v(1,1,1,2) * r2
       psdv4_6t(1,5) = v(1,1,2,3) * r2
       psdv4_6t(1,6) = v(1,1,1,3) * r2
       psdv4_6t(2,1) = v(2,2,1,1)
       psdv4_6t(2,2) = v(2,2,2,2)
       psdv4_6t(2,3) = v(2,2,3,3)
       psdv4_6t(2,4) = v(2,2,1,2) * r2
       psdv4_6t(2,5) = v(2,2,2,3) * r2
       psdv4_6t(2,6) = v(2,2,1,3) * r2
       psdv4_6t(3,1) = v(3,3,1,1)
       psdv4_6t(3,2) = v(3,3,2,2)
       psdv4_6t(3,3) = v(3,3,3,3)
       psdv4_6t(3,4) = v(3,3,1,2) * r2
       psdv4_6t(3,5) = v(3,3,2,3) * r2
       psdv4_6t(3,6) = v(3,3,1,3) * r2
       psdv4_6t(4,1) = v(1,2,1,1)
       psdv4_6t(4,2) = v(1,2,2,2)
       psdv4_6t(4,3) = v(1,2,3,3)
       psdv4_6t(4,4) = v(1,2,1,2) * 2.d0
       psdv4_6t(4,5) = v(1,2,2,3) * 2.d0
       psdv4_6t(4,6) = v(1,2,1,3) * 2.d0
       psdv4_6t(5,1) = v(2,3,1,1)
       psdv4_6t(5,2) = v(2,3,2,2)
       psdv4_6t(5,3) = v(2,3,3,3)
       psdv4_6t(5,4) = v(2,3,1,2) * 2.d0
       psdv4_6t(5,5) = v(2,3,2,3) * 2.d0
       psdv4_6t(5,6) = v(2,3,1,3) * 2.d0
       psdv4_6t(6,1) = v(1,3,1,1)
       psdv4_6t(6,2) = v(1,3,2,2)
       psdv4_6t(6,3) = v(1,3,3,3)
       psdv4_6t(6,4) = v(1,3,1,2) * 2.d0
       psdv4_6t(6,5) = v(1,3,2,3) * 2.d0
       psdv4_6t(6,6) = v(1,3,1,3) * 2.d0
      !
       return 
       end     
!________________________________________

       subroutine rotation(r,r6,n,typ,i12,i23,i13)
!
!**** subroutine to create a rotation matrix for constitutive tensors
!     for mat is tensor format with dsqrt(2) in shear terms
!
!
!
       implicit none
       !
       integer n
       real*8  r(3,3), r6(n,n)         !* rotation tensors
       character typ                   !* flag for type of rotation matrix
                                       !  'p' = for pseusovectors (with root2 in shear)
                                       !  'e' = Voigt for strains 
                                       !  'v' = Voigt for stresses or constit. tensors (RCRt)
                                       !
       real*8, parameter:: root2 = 1.4142135623731d0
       real*8, parameter:: two = 2.d0
       !
       integer i12, i23, i13           !* stress orders
       !
       if (typ == 'p') then
          r6(1,1) = r(1,1)**2
          r6(1,2) = r(1,2)**2
          r6(1,3) = r(1,3)**2
          r6(2,1) = r(2,1)**2
          r6(2,2) = r(2,2)**2
          r6(2,3) = r(2,3)**2
          r6(3,1) = r(3,1)**2
          r6(3,2) = r(3,2)**2
          r6(3,3) = r(3,3)**2
          !
          r6(1,i12) = root2 * r(1,1) * r(1,2)
          r6(2,i12) = root2 * r(2,1) * r(2,2)
          r6(3,i12) = root2 * r(3,1) * r(3,2)
          !
          r6(i12,1) = root2 * r(1,1) * r(2,1) 
          r6(i12,2) = root2 * r(1,2) * r(2,2) 
          r6(i12,3) = root2 * r(1,3) * r(2,3) 
          !
          r6(i12,i12) = r(1,1) * r(2,2) + r(1,2) * r(2,1)
          !
        if (n > 4) then
        !
        r6(1,i23) = root2 * r(1,2) * r(1,3)
        r6(1,i13) = root2 * r(1,3) * r(1,1) 
        r6(2,i23) = root2 * r(2,2) * r(2,3)
        r6(2,i13) = root2 * r(2,3) * r(2,1)
        r6(3,i23) = root2 * r(3,2) * r(3,3)
        r6(3,i13) = root2 * r(3,3) * r(3,1)
        !
        r6(i23,1) = root2 * r(3,1) * r(2,1) 
        r6(i23,2) = root2 * r(3,2) * r(2,2) 
        r6(i23,3) = root2 * r(3,3) * r(2,3) 
        r6(i13,1) = root2 * r(3,1) * r(1,1) 
        r6(i13,2) = root2 * r(3,2) * r(1,2) 
        r6(i13,3) = root2 * r(3,3) * r(1,3) 
        !
        r6(i12,i23) = r(1,2) * r(2,3) + r(1,3) * r(2,2)
        r6(i12,i13) = r(1,3) * r(2,1) + r(1,1) * r(2,3)
        r6(i23,i12) = r(2,1) * r(3,2) + r(2,2) * r(3,1)
        r6(i23,i23) = r(2,2) * r(3,3) + r(2,3) * r(3,2)
        r6(i23,i13) = r(2,3) * r(3,1) + r(2,1) * r(3,3)
        r6(i13,i12) = r(3,1) * r(1,2) + r(3,2) * r(1,1)
        r6(i13,i23) = r(3,2) * r(1,3) + r(3,3) * r(1,2)
        r6(i13,i13) = r(3,3) * r(1,1) + r(3,1) * r(1,3)
        !
        endif
        !
       endif
        !
        !
        if (typ == 'e') then
          r6(1,1) = r(1,1)**2
          r6(1,2) = r(1,2)**2
          r6(1,3) = r(1,3)**2
          r6(2,1) = r(2,1)**2
          r6(2,2) = r(2,2)**2
          r6(2,3) = r(2,3)**2
          r6(3,1) = r(3,1)**2
          r6(3,2) = r(3,2)**2
          r6(3,3) = r(3,3)**2
        !
         r6(1,i12) = r(1,1) * r(1,2)
         r6(2,i12) = r(2,1) * r(2,2)
         r6(3,i12) = r(3,1) * r(3,2)
       !
         r6(i12,1) = 2.d0 * r(1,1) * r(2,1) 
         r6(i12,2) = 2.d0 * r(1,2) * r(2,2) 
         r6(i12,3) = 2.d0 * r(1,3) * r(2,3) 
       !
         r6(i12,i12) = r(1,1) * r(2,2) + r(1,2) * r(2,1)
       !
       if (n > 4) then
        !
        r6(1,i23) = r(1,2) * r(1,3)
        r6(1,i13) = r(1,3) * r(1,1) 
        r6(2,i23) = r(2,2) * r(2,3)
        r6(2,i13) = r(2,3) * r(2,1)
        r6(3,i23) = r(3,2) * r(3,3)
        r6(3,i13) = r(3,3) * r(3,1)
        !
        r6(i23,1) = 2.d0 * r(3,1) * r(2,1) 
        r6(i23,2) = 2.d0 * r(3,2) * r(2,2) 
        r6(i23,3) = 2.d0 * r(3,3) * r(2,3) 
        r6(i13,1) = 2.d0 * r(3,1) * r(1,1) 
        r6(i13,2) = 2.d0 * r(3,2) * r(1,2) 
        r6(i13,3) = 2.d0 * r(3,3) * r(1,3) 
        !
        r6(i12,i23) = r(1,2) * r(2,3) + r(1,3) * r(2,2)
        r6(i12,i13) = r(1,3) * r(2,1) + r(1,1) * r(2,3)
        r6(i23,i12) = r(2,1) * r(3,2) + r(2,2) * r(3,1)
        r6(i23,i23) = r(2,2) * r(3,3) + r(2,3) * r(3,2)
        r6(i23,i13) = r(2,3) * r(3,1) + r(2,1) * r(3,3)
        r6(i13,i12) = r(3,1) * r(1,2) + r(3,2) * r(1,1)
        r6(i13,i23) = r(3,2) * r(1,3) + r(3,3) * r(1,2)
        r6(i13,i13) = r(3,3) * r(1,1) + r(3,1) * r(1,3)
        !
       endif
    !
       endif
!
!
       if (typ == 'v' .or. typ == 's') then
         r6(1,1) = r(1,1)**2
         r6(1,2) = r(1,2)**2
         r6(1,3) = r(1,3)**2
         r6(2,1) = r(2,1)**2
         r6(2,2) = r(2,2)**2
         r6(2,3) = r(2,3)**2
         r6(3,1) = r(3,1)**2
         r6(3,2) = r(3,2)**2
         r6(3,3) = r(3,3)**2
        !
        r6(1,i12) = 2.d0 * r(1,1) * r(1,2)
        r6(2,i12) = 2.d0 * r(2,1) * r(2,2)
        r6(3,i12) = 2.d0 * r(3,1) * r(3,2)
        !
        r6(i12,1) = r(1,1) * r(2,1) 
        r6(i12,2) = r(1,2) * r(2,2) 
        r6(i12,3) = r(1,3) * r(2,3) 
        !
        r6(i12,i12) = r(1,1) * r(2,2) + r(1,2) * r(2,1)
        !
        if (n > 4) then
        !
        r6(1,i23) = 2.d0 * r(1,2) * r(1,3)
        r6(1,i13) = 2.d0 * r(1,3) * r(1,1) 
        r6(2,i23) = 2.d0 * r(2,2) * r(2,3)
        r6(2,i13) = 2.d0 * r(2,3) * r(2,1)
        r6(3,i23) = 2.d0 * r(3,2) * r(3,3)
        r6(3,i13) = 2.d0 * r(3,3) * r(3,1)
        !
        r6(i23,1) = r(3,1) * r(2,1) 
        r6(i23,2) = r(3,2) * r(2,2) 
        r6(i23,3) = r(3,3) * r(2,3) 
        r6(i13,1) = r(3,1) * r(1,1) 
        r6(i13,2) = r(3,2) * r(1,2) 
        r6(i13,3) = r(3,3) * r(1,3) 
        !
        r6(i12,i23) = r(1,2) * r(2,3) + r(1,3) * r(2,2)
        r6(i12,i13) = r(1,3) * r(2,1) + r(1,1) * r(2,3)
        r6(i23,i12) = r(2,1) * r(3,2) + r(2,2) * r(3,1)
        r6(i23,i23) = r(2,2) * r(3,3) + r(2,3) * r(3,2)
        r6(i23,i13) = r(2,3) * r(3,1) + r(2,1) * r(3,3)
        r6(i13,i12) = r(3,1) * r(1,2) + r(3,2) * r(1,1)
        r6(i13,i23) = r(3,2) * r(1,3) + r(3,3) * r(1,2)
        r6(i13,i13) = r(3,3) * r(1,1) + r(3,1) * r(1,3)
        !
       endif
       !
      endif
       !
      return
      end 
!___________________________________________________________
       
        
      subroutine constantsIDT(IB1,IB2,XI33,XI333,XI66,XI99,XInn)
      implicit none
      integer i,j,k
      integer IB1(9),IB2(9)
      real(8) XI33(3,3),XI66(6,6),XI99(9,9),XInn(48,48)
      real(8) XI333(3,3,3)
      XI333=0
      XI333(1,2,3)=+1
      XI333(1,3,2)=-1
      XI333(2,1,3)=-1
      XI333(2,3,1)=+1
      XI333(3,1,2)=+1
      XI333(3,2,1)=-1
      XI33=0; do i=1,3;  XI33(i,i)=1; enddo         
      XI66=0; do i=1,6;  XI66(i,i)=1; enddo         
      XI99=0; do i=1,9;  XI99(i,i)=1; enddo         
      XInn=0; do i=1,48; XInn(i,i)=1; enddo         
      IB1=[1,2,3,1,1,2,2,3,3]
      IB2=[1,2,3,2,3,3,1,1,2]
      return
      end




