#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "Damkoh_F.H"
#include "BC_TYPES.H"
#include "SPACE.H"
#include "EdgeSpeed_F.H"      
#define SDIM BL_SPACEDIM

      subroutine FORT_PUSHVTOG(lo, hi, dlo, dhi, U, DIMS(U), nc)
      implicit none
      integer nc
      integer  lo(SDIM),  hi(SDIM)
      integer dlo(SDIM), dhi(SDIM)
      integer DIMDEC(U)
      REAL_T U(DIMV(U),nc)
      integer i,j,k,n
      integer bc(SDIM,2)
      REAL_T xlo(SDIM), dx(SDIM)

      REAL_T x,y

c     Make up something for these that gets what we want
      do n=1,SDIM
         bc(n,1) = HOEXTRAP
         bc(n,2) = HOEXTRAP
         dx(n) = one
         xlo(n) = one
      end do

      do n=1,nc
         call filcc(U(D_DECL(ARG_L1(U),ARG_L2(U),ARG_L3(U)),n),DIMS(U),lo,hi,dx,xlo,bc)
      end do
      end

#if BL_SPACEDIM==2
      subroutine FORT_COMPDIFFTERMS(lo, hi, Y, DIMS(Y),
     &     RhoD, DIMS(RhoD),density,DIMS(density),vel,DIMS(vel), DT, DIMS(DT),
     &     dx, ns,si, nwhich) 
     
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
      integer si
      
      integer DIMDEC(Y)
      integer DIMDEC(density)
      integer DIMDEC(vel)
      integer DIMDEC(RhoD)
      
      integer DIMDEC(DT)
     
      REAL_T Y(DIMV(Y),ns)
      REAL_T density(DIMV(density))
      REAL_T vel(DIMV(vel),BL_SPACEDIM)
      REAL_T RhoD_xp(si), RhoD_xm(si)
      REAL_T RhoD_yp(si), RhoD_ym(si)     
      REAL_T RhoD(DIMV(RhoD),ns)
       
      REAL_T DT(DIMV(DT),(2*(BL_SPACEDIM)*si))
      REAL_T dx(SDIM)
      integer nwhich(si)
      integer i,j,n,p,s
      
      REAL_T dxInv2, dyInv2

      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      print*, si, nwhich(1)

       
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
            do p= 1,si      
              n= nwhich(si)+1
               RhoD_xp(p) = 0.5*( RhoD(i,j,n) + RhoD(i+1,j,n) )
               RhoD_xm(p) = 0.5*( RhoD(i,j,n) + RhoD(i-1,j,n) )
               RhoD_yp(p) = 0.5*( RhoD(i,j,n) + RhoD(i,j+1,n) )
               RhoD_ym(p) = 0.5*( RhoD(i,j,n) + RhoD(i,j-1,n) )
               
c               H_xm(n) = 0.5*( H(i,j,n) + H(i-1,j,n) )
c               H_xp(n) = 0.5*( H(i,j,n) + H(i+1,j,n) )
c               H_yp(n) = 0.5*( H(i,j,n) + H(i,j+1,n) )
c               H_ym(n) = 0.5*( H(i,j,n) + H(i,j-1,n) )
            
            
c            lambda_xm = 0.5*( lambda(i,j) + lambda(i-1,j) )
c            lambda_xp = 0.5*( lambda(i,j) + lambda(i+1,j) )
c            lambda_yp = 0.5*( lambda(i,j) + lambda(i,j+1) )
c            lambda_ym = 0.5*( lambda(i,j) + lambda(i,j-1) )
            
c            DT(i,j,ns+1) =
c    &           ( lambda_xp*(T(i+1,j)-T(i  ,j))
c     &           - lambda_xm*(T(i  ,j)-T(i-1,j)) ) *dxInv2 +
c     &           ( lambda_yp*(T(i,j+1)-T(i,j  ))
c     &           - lambda_ym*(T(i,j  )-T(i,j-1)) ) *dyInv2
             enddo
           do n=1,si
               s= nwhich(si)+1
        
          DT(i,j,(n-1)*BL_SPACEDIM+3) = 
     &              ( RhoD_xp(n)*(Y(i+1,j,s)-Y(i  ,j,s))
     &              - RhoD_xm(n)*(Y(i  ,j,s)-Y(i-1,j,s)) ) * dxInv2 
          DT(i,j,(n-1)*BL_SPACEDIM+4) =
     &              ( RhoD_yp(n)*(Y(i,j+1,s)-Y(i,j  ,s))
     &              - RhoD_ym(n)*(Y(i,j  ,s)-Y(i,j-1,s)) ) * dyInv2
               
c               DT(i,j,ns+1) = DT(i,j,ns+1) + 
c     &              ( RhoD_xp(n)*H_xp(n)*(Y(i+1,j,n)-Y(i  ,j,n))
c     &              - RhoD_xm(n)*H_xm(n)*(Y(i  ,j,n)-Y(i-1,j,n)) ) * dxInv2 +
c     &              ( RhoD_yp(n)*H_yp(n)*(Y(i,j+1,n)-Y(i,j  ,n))
c     &              - RhoD_ym(n)*H_ym(n)*(Y(i,j  ,n)-Y(i,j-1,n)) ) * dyInv2
            
             
         DT(i,j,2*(n-1)*BL_SPACEDIM+1) = (density(i+1,j)*vel(i+1,j,1)*Y(i+1,j,s)-
     &                    density(i-1,j)*vel(i-1,j,1)*Y(i-1,j,s))/(2.0*dx(1))
         DT(i,j,2*(n-1)*BL_SPACEDIM+2) = (density(i,j+1)*vel(i,j+1,2)*Y(i,j+1,s)-
     &                    density(i,j-1)*vel(i,j-1,2)*Y(i,j-1,s))/(2.0*dx(2))
c         DT(i,j,3) = (density(i+1,j)*vel(i+1,j,1)*Y(i+1,j,n)-
c     &                      density(i-1,j)*vel(i-1,j,1)*Y(i-1,j,n))
           enddo
         end do
      end do
      end
#else
      subroutine FORT_COMPDIFFTERMS(lo, hi, Y, DIMS(Y),
     &   RhoD, DIMS(RhoD),lambda,DIMS(lambda), Reac,DIMS(Reac),
     &   density,DIMS(density),Z,DIMS(Z), Speed, DIMS(Speed), dx,ws,ns)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(lambda)
      integer DIMDEC(Z)
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      integer DIMDEC(density)
      integer DIMDEC(Reac)
      integer DIMDEC(Speed)
      integer ws, ns,wsp
      REAL_T Y(DIMV(Y),ns), Speed(DIMV(Speed),4)
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T density(DIMV(density))
      REAL_T DT(DIMV(Speed),2)
      REAL_T lambda(DIMV(lambda))   
      REAL_T Reac(DIMV(Reac))
      REAL_T Z(DIMV(Z))
      REAL_T dx(SDIM)
      REAL_T dx2, dy2 , dz2,sumFlux
      REAL_T DOT(DIMV(Z))
      REAL_T Flux_Z(DIMV(Z),3), grad_Z(DIMV(Z),4)
      REAL_T Flux_Y(DIMV(Y),ns,3), grad_Y(DIMV(Y),4)
      REAL_T small, HAvg, YAvg
      integer i,j,k,n,p,s
      REAL_T const, tmp1, tmp2 
      REAL_T dxInv2, dyInv2, dzInv2
      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      dzInv2 = 1.d0 / (dx(3)*dx(3))
      dx2= 0.5/dx(1)
      dy2 = 0.5/dx(2)
      dz2 = 0.5/dx(3)
      small = 1e-2
      wsp = ws+1
       do k= lo(3), hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

           do n=1,ns
            HAvg = 1.0/(RhoD(i-1,j,k,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,1) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i-1  ,j,k,n))/dx(1)
            HAvg = 1.0/(lambda(i-1,j,k)) + 1.0/(lambda(i,j,k))
            HAvg = 2.0/HAvg
            Flux_Z(i,j,k,1) = HAvg*(Z(i-1,j,k)-Z(i,j,k))
           
           enddo
         end do
      end do
      enddo
    
       

       do k= lo(3), hi(3)
        do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

           do n=1,ns
            HAvg = 1.0/(RhoD(i,j-1,k,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,2) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i  ,j-1,k,n))/dx(2)
            HAvg = 1.0/(lambda(i,j-1,k)) + 1.0/(lambda(i,j,k))
            HAvg = 2.0/HAvg
            Flux_Z(i,j,k,2) = HAvg*(Z(i,j-1,k)-Z(i,j,k))
           
           enddo
         end do
      end do
      enddo
 
      do k= lo(3), hi(3)+1
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)

           do n=1,ns
            HAvg = 1.0/(RhoD(i,j,k-1,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,3) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i  ,j,k-1,n))/dx(3)
            HAvg = 1.0/(lambda(i,j,k)) + 1.0/(lambda(i,j,k-1))
            HAvg = 2.0/HAvg
            Flux_Z(i,j,k,3) = HAvg*(Z(i,j,k)-Z(i,j,k-1))
           
           enddo
         end do
      end do
      enddo
 

!Correction for species fluxes
   
        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux_Y(i,j,k,n,1)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i-1,j,k,n)+Y(i,j,k,n))
                    Flux_Y(i,j,k,n,1) = Flux_Y(i,j,k,n,1)-  sumFlux*YAvg
                  end do
               end do
            end do
         end do

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux_Y(i,j,k,n,2)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i,j-1,k,n)+Y(i,j,k,n))
                    Flux_Y(i,j,k,n,2) = Flux_Y(i,j,k,n,2)-  sumFlux*YAvg
                  end do
               end do
            end do
         end do


       do k = lo(3),hi(3)+1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux_Y(i,j,k,n,3)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i,j,k-1,n)+Y(i,j,k,n))
                    Flux_Y(i,j,k,n,3) = Flux_Y(i,j,k,n,3)-  sumFlux*YAvg
                  end do
               end do
            end do
         end do
  
       do k=lo(3), hi(3)
        do j=lo(2), hi(2)
         do i=lo(1), hi(1)
           grad_Z(i,j,k,1) = (Z(i+1,j,k)-Z(i-1,j,k))*dx2
           grad_Z(i,j,k,2) = (Z(i,j+1,k)-Z(i,j-1,k))*dy2
           grad_Z(i,j,k,3) = (Z(i,j,k+1)-Z(i,j,k-1))*dz2
           grad_Z(i,j,k,4) = sqrt(grad_Z(i,j,k,1)**2+grad_Z(i,j,k,2)**2+grad_Z(i,j,k,3)**2)
           
           grad_Y(i,j,k,1) = (Y(i+1,j,k,wsp)-Y(i-1,j,k,wsp))*dx2
           grad_Y(i,j,k,2) = (Y(i,j+1,k,wsp)-Y(i,j-1,k,wsp))*dy2
           grad_Y(i,j,k,3) = (Y(i,j,k+1,wsp)-Y(i,j,k-1,wsp))*dz2
           grad_Y(i,j,k,4) = sqrt(grad_Y(i,j,k,1)**2+grad_Y(i,j,k,2)**2+grad_Y(i,j,k,3)**2) 
         
         enddo
        enddo
       enddo

        do k=lo(3), hi(3)
         do j= lo(2), hi(2)
          do i=lo(1), hi(1)

           DT(i,j,k,1) = (Flux_Y(i+1,j,k,wsp,1)-Flux_Y(i,j,k,wsp,1))/(dx(1))+
     &               (Flux_Y(i,j+1,k,wsp,2)-Flux_Y(i,j,k,wsp,2))/dx(2) +
     &               (Flux_Y(i,j,k+1,wsp,3)-Flux_Y(i,j,k,wsp,3))/dx(3)
          
           DT(i,j,k,2) = (Flux_Z(i+1,j,k,1)-Flux_Z(i,j,k,1))/(dx(1))+
     &               (Flux_Z(i,j+1,k,2)-Flux_Z(i,j,k,2))/dx(2) +
     &               (Flux_Z(i,j,k+1,3)-Flux_Z(i,j,k,3))/dx(3)
           DOT(i,j,k) = (grad_Z(i,j,k,1)*grad_Y(i,j,k,1)+grad_Z(i,j,k,2)*grad_Y(i,j,k,2) +
     &                   grad_Z(i,j,k,3)*grad_Y(i,j,k,3))/(grad_Z(i,j,k,4)*grad_Y(i,j,k,4)+small)
           

       enddo
      enddo
      enddo

      do  i= lo(1), hi(1)
        do j= lo(2), hi(2)
         do k=lo(3), hi(3) 
!          if(Y(i,j,k,wsp).gt.1e-4)  then
         Speed(i,j,k,1) = -DT(i,j,k,2)/(density(i,j,k)*grad_Z(i,j,k,4)+small)

         Speed(i,j,k,2) = (-Reac(i,j,k)-DT(i,j,k,1))/(density(i,j,k)*grad_Y(i,j,k,4)+small)
!          else 
!          Speed(i,j,k,2)=0
!          Speed(i,j,k,1) =0 
!          endif
         Speed(i,j,k,3) = (Speed(i,j,k,2)-DOT(i,j,k)*Speed(i,j,k,1))/(sqrt(1-DOT(i,j,k)**2)+small)
         Speed(i,j,k,4) = DOT(i,j,k)
        enddo
       enddo
      enddo
    
            end
#endif
