#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "SPACE.H"
#include "compEnergyTerms_F.H"
      
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
c convection term calculation

      subroutine FORT_COMPCONVECTERMS(lo, hi, Rho, DIMS(Rho),
     &     Y, DIMS(Y),  H, DIMS(H),  Ht, DIMS(Ht),  Vx, DIMS(Vx), Vy, DIMS(Vy), CT, DIMS(CT), dx, ns)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
      integer DIMDEC(Rho)
      integer DIMDEC(Y)
      integer DIMDEC(H)
      integer DIMDEC(Vx)
      integer DIMDEC(Vy)
      integer DIMDEC(CT)
      integer DIMDEC(Ht)
      REAL_T Rho(DIMV(Rho))
      REAL_T Ht(DIMV(Ht))
      REAL_T Y(DIMV(Y),ns)
      REAL_T H(DIMV(H),ns)
      REAL_T Vx(DIMV(Vx))
      REAL_T Vy(DIMV(Vy))
      REAL_T CT(DIMV(CT), ns+1)
      REAL_T dx(SDIM)
     
      REAL_T h_xp
      REAL_T h_xm
      REAL_T h_yp
      REAL_T h_ym
      REAL_T rho_xp
      REAL_T rho_xm
      REAL_T rho_yp
      REAL_T rho_ym
      REAL_T Y_xp(maxspec)
      REAL_T Y_xm(maxspec)
      REAL_T Y_yp(maxspec)
      REAL_T Y_ym(maxspec)
      REAL_T U_xp
      REAL_T U_xm
      REAL_T V_yp
      REAL_T V_ym
      integer i,j, n
      
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            
            do n=1,ns
               Y_xm(n) = 0.5*( Y(i,j,n) + Y(i-1,j,n) )
               Y_xp(n) = 0.5*( Y(i,j,n) + Y(i+1,j,n) )
               Y_yp(n) = 0.5*( Y(i,j,n) + Y(i,j+1,n) )
               Y_ym(n) = 0.5*( Y(i,j,n) + Y(i,j-1,n) )
            enddo
            
            rho_xp = 0.5*( Rho(i,j) + Rho(i+1,j) )
            rho_xm = 0.5*( Rho(i,j) + Rho(i-1,j) )
            rho_yp = 0.5*( Rho(i,j) + Rho(i,j+1) )
            rho_ym = 0.5*( Rho(i,j) + Rho(i,j-1) )
            
            h_xp = 0.5*( Ht(i,j) + Ht(i+1,j) )
            h_xm = 0.5*( Ht(i,j) + Ht(i-1,j) )
            h_yp = 0.5*( Ht(i,j) + Ht(i,j+1) )
            h_ym = 0.5*( Ht(i,j) + Ht(i,j-1) )

            U_xp = 0.5*( Vx(i,j) + Vx(i+1,j) )
            U_xm = 0.5*( Vx(i,j) + Vx(i-1,j) )

            V_yp = 0.5*( Vy(i,j) + Vy(i,j+1) )
            V_ym = 0.5*( Vy(i,j) + Vy(i,j-1) )
          
            do n=1,ns
               CT(i,j,n) = 
     &              (rho_xp*Y_xp(n)*U_xp - rho_xm*Y_xm(n)*U_xm)/dx(1) + 
     &              (rho_yp*Y_yp(n)*V_yp - rho_ym*Y_ym(n)*V_ym)/dx(2)  
               
            end do
            
            CT(i,j, ns+1) =
     &           (rho_xp*h_xp*U_xp - rho_xm*h_xm*U_xm)/dx(1) + 
     &           (rho_yp*h_yp*V_yp - rho_ym*h_ym*V_ym)/dx(2)  
            
         enddo
      end do
      end
      

c diffusion term calculation
#if BL_SPACEDIM==2
      subroutine FORT_COMPDIFFTERMS(lo, hi, T, DIMS(T), Y, DIMS(Y),
     &     H, DIMS(H), RhoD, DIMS(RhoD), lambda, DIMS(lambda), DT, DIMS(DT),
     &     dx, ns)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(H)
      integer DIMDEC(RhoD)
      integer DIMDEC(lambda)
      integer DIMDEC(DT)
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),ns)
      REAL_T H(DIMV(H),ns)
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T lambda(DIMV(lambda))
      REAL_T DT(DIMV(DT),ns+2)
      REAL_T dx(SDIM)
      REAL_T H_xp(maxspec), RhoD_xp(maxspec), Y_xp(maxspec)
      REAL_T H_xm(maxspec), RhoD_xm(maxspec), Y_xm(maxspec)
      REAL_T H_yp(maxspec), RhoD_yp(maxspec), Y_yp(maxspec)
      REAL_T H_ym(maxspec), RhoD_ym(maxspec), Y_ym(maxspec)
     
      integer i,j,n
      REAL_T lambda_xp,lambda_xm,lambda_yp,lambda_ym
      REAL_T dxInv2, dyInv2

      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
           
            do n=1,ns
               RhoD_xp(n) = 0.5*( RhoD(i,j,n) + RhoD(i+1,j,n) )
               RhoD_xm(n) = 0.5*( RhoD(i,j,n) + RhoD(i-1,j,n) )
               RhoD_yp(n) = 0.5*( RhoD(i,j,n) + RhoD(i,j+1,n) )
               RhoD_ym(n) = 0.5*( RhoD(i,j,n) + RhoD(i,j-1,n) )
               
               H_xm(n) = 0.5*( H(i,j,n) + H(i-1,j,n) )
               H_xp(n) = 0.5*( H(i,j,n) + H(i+1,j,n) )
               H_yp(n) = 0.5*( H(i,j,n) + H(i,j+1,n) )
               H_ym(n) = 0.5*( H(i,j,n) + H(i,j-1,n) )

   
               Y_xm(n) = 0.5*( Y(i,j,n) + Y(i-1,j,n) )
               Y_xp(n) = 0.5*( Y(i,j,n) + Y(i+1,j,n) )
               Y_yp(n) = 0.5*( Y(i,j,n) + Y(i,j+1,n) )
               Y_ym(n) = 0.5*( Y(i,j,n) + Y(i,j-1,n) )
            enddo
            
            lambda_xm = 0.5*( lambda(i,j) + lambda(i-1,j) )
            lambda_xp = 0.5*( lambda(i,j) + lambda(i+1,j) )
            lambda_yp = 0.5*( lambda(i,j) + lambda(i,j+1) )
            lambda_ym = 0.5*( lambda(i,j) + lambda(i,j-1) )
            
            DT(i,j,ns+1) =
     &           ( lambda_xp*(T(i+1,j)-T(i  ,j))
     &           - lambda_xm*(T(i  ,j)-T(i-1,j)) ) *dxInv2 +
     &           ( lambda_yp*(T(i,j+1)-T(i,j  ))
     &           - lambda_ym*(T(i,j  )-T(i,j-1)) ) *dyInv2
            
            do n=1,ns
               DT(i,j,n) = 
     &              ( RhoD_xp(n)*(Y(i+1,j,n)-Y(i  ,j,n))
     &              - RhoD_xm(n)*(Y(i  ,j,n)-Y(i-1,j,n)) ) * dxInv2 +
     &              ( RhoD_yp(n)*(Y(i,j+1,n)-Y(i,j  ,n))
     &              - RhoD_ym(n)*(Y(i,j  ,n)-Y(i,j-1,n)) ) * dyInv2
               
               DT(i,j,ns+2) = DT(i,j,ns+2) + 
     &               (RhoD_xp(n)*H_xp(n) - RhoD_xm(n)*H_xm(n))*(Y_xp(n)-Y_xm(n))* dxInv2
     &               + (RhoD_yp(n)*H_yp(n) - RhoD_ym(n)*H_ym(n))*(Y_yp(n)-Y_ym(n))* dyInv2         

               DT(i,j,ns+3) = DT(i,j,ns+3) + 
     &               RhoD(i,j,n)*H(i,j,n)*(((Y(i+1,j,n)-Y(i  ,j,n))-
     &               (Y(i,j,n)-Y(i-1,j,n)))* dxInv2 +
     &               ((Y(i,j+1,n)-Y(i ,j,n))-
     &               (Y(i,j,n)-Y(i  ,j-1,n)))* dyInv2)
   
            enddo
            
         end do
      end do
      end
#else
      subroutine FORT_COMPDIFFTERMS(lo, hi, T, DIMS(T), Y, DIMS(Y),
     &     H, DIMS(H), RhoD, DIMS(RhoD), lambda, DIMS(lambda), DT, DIMS(DT),
     &     dx, ns)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(H)
      integer DIMDEC(RhoD)
      integer DIMDEC(lambda)
      integer DIMDEC(DT)
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),ns)
      REAL_T H(DIMV(H),ns)
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T lambda(DIMV(lambda))
      REAL_T DT(DIMV(DT),ns+2)
      REAL_T dx(SDIM)
      REAL_T H_xp(maxspec), RhoD_xp(maxspec)
      REAL_T H_xm(maxspec), RhoD_xm(maxspec)
      REAL_T H_yp(maxspec), RhoD_yp(maxspec)
      REAL_T H_ym(maxspec), RhoD_ym(maxspec)
      REAL_T H_zp(maxspec), RhoD_zp(maxspec)
      REAL_T H_zm(maxspec), RhoD_zm(maxspec)
      integer i,j,k,n
      REAL_T lambda_xp,lambda_xm,lambda_yp,lambda_ym,lambda_zp,lambda_zm
      REAL_T dxInv2, dyInv2, dzInv2

      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      dzInv2 = 1.d0 / (dx(3)*dx(3))

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,ns
                  RhoD_xp(n) = 0.5*( RhoD(i,j,k,n) + RhoD(i+1,j,k,n) )
                  RhoD_xm(n) = 0.5*( RhoD(i,j,k,n) + RhoD(i-1,j,k,n) )
                  RhoD_yp(n) = 0.5*( RhoD(i,j,k,n) + RhoD(i,j+1,k,n) )
                  RhoD_ym(n) = 0.5*( RhoD(i,j,k,n) + RhoD(i,j-1,k,n) )
                  RhoD_zp(n) = 0.5*( RhoD(i,j,k,n) + RhoD(i,j,k+1,n) )
                  RhoD_zm(n) = 0.5*( RhoD(i,j,k,n) + RhoD(i,j,k-1,n) )

                  H_xp(n) = 0.5*( H(i,j,k,n) + H(i+1,j,k,n) )
                  H_xm(n) = 0.5*( H(i,j,k,n) + H(i-1,j,k,n) )
                  H_yp(n) = 0.5*( H(i,j,k,n) + H(i,j+1,k,n) )
                  H_ym(n) = 0.5*( H(i,j,k,n) + H(i,j-1,k,n) )
                  H_zp(n) = 0.5*( H(i,j,k,n) + H(i,j,k+1,n) )
                  H_zm(n) = 0.5*( H(i,j,k,n) + H(i,j,k-1,n) )
               enddo

               lambda_xp = 0.5*( lambda(i,j,k) + lambda(i+1,j,k) )
               lambda_xm = 0.5*( lambda(i,j,k) + lambda(i-1,j,k) )
               lambda_yp = 0.5*( lambda(i,j,k) + lambda(i,j+1,k) )
               lambda_ym = 0.5*( lambda(i,j,k) + lambda(i,j-1,k) )
               lambda_zp = 0.5*( lambda(i,j,k) + lambda(i,j,k+1) )
               lambda_zm = 0.5*( lambda(i,j,k) + lambda(i,j,k-1) )

               DT(i,j,k,ns+1) =
     &              ( lambda_xp*(T(i+1,j,k)-T(i  ,j,k))
     &              - lambda_xm*(T(i  ,j,k)-T(i-1,j,k)) ) *dxInv2 +
     &              ( lambda_yp*(T(i,j+1,k)-T(i,j  ,k))
     &              - lambda_ym*(T(i,j  ,k)-T(i,j-1,k)) ) *dyInv2 +
     &              ( lambda_zp*(T(i,j,k+1)-T(i,j,k  ))
     &              - lambda_zm*(T(i,j,k  )-T(i,j,k-1)) ) *dzInv2

               do n=1,ns
                  DT(i,j,k,n) = 
     &                 ( RhoD_xp(n)*(Y(i+1,j,k,n)-Y(i  ,j,k,n))
     &                 - RhoD_xm(n)*(Y(i  ,j,k,n)-Y(i-1,j,k,n)) ) * dxInv2 +
     &                 ( RhoD_yp(n)*(Y(i,j+1,k,n)-Y(i,j  ,k,n))
     &                 - RhoD_ym(n)*(Y(i,j  ,k,n)-Y(i,j-1,k,n)) ) * dyInv2 +
     &                 ( RhoD_zp(n)*(Y(i,j,k+1,n)-Y(i,j,k  ,n))
     &                 - RhoD_zm(n)*(Y(i,j,k  ,n)-Y(i,j,k-1,n)) ) * dzInv2

                  DT(i,j,k,ns+2) = DT(i,j,k,ns+2) + 
     &                 ( RhoD_xp(n)*H_xp(n)*(Y(i+1,j,k,n)-Y(i  ,j,k,n))
     &                 - RhoD_xm(n)*H_xm(n)*(Y(i  ,j,k,n)-Y(i-1,j,k,n)) ) * dxInv2 +
     &                 ( RhoD_yp(n)*H_yp(n)*(Y(i,j+1,k,n)-Y(i,j  ,k,n))
     &                 - RhoD_ym(n)*H_ym(n)*(Y(i,j  ,k,n)-Y(i,j-1,k,n)) ) * dyInv2 +
     &                 ( RhoD_zp(n)*H_zp(n)*(Y(i,j,k+1,n)-Y(i,j,k  ,n))
     &                 - RhoD_zm(n)*H_zm(n)*(Y(i,j,k  ,n)-Y(i,j,k-1,n)) ) * dzInv2
               enddo

            end do
         end do
      end do
      end
#endif
c convection term calculation for energy equation written as sensible enthalpy
      subroutine FORT_COMPCDDRTERMS(lo, hi, Rho, DIMS(Rho),
     &     Y, DIMS(Y),  T, DIMS(T),  Cpmix, DIMS(Cpmix), H, DIMS(H),  
     &     Ht, DIMS(Ht),  Vx, DIMS(Vx), Vy, DIMS(Vy), RhoD, DIMS(RhoD), lambda, DIMS(lambda),
     &     CT, DIMS(CT), dx, ns)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
      integer DIMDEC(Rho)
      integer DIMDEC(Cpmix)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(H)
      integer DIMDEC(Vx)
      integer DIMDEC(Vy)
      integer DIMDEC(RhoD)
      integer DIMDEC(lambda)
      integer DIMDEC(CT)
      integer DIMDEC(Ht)
      REAL_T Rho(DIMV(Rho))
      REAL_T Cpmix(DIMV(Cpmix))
      REAL_T T(DIMV(T))
      REAL_T Ht(DIMV(Ht))
      REAL_T Y(DIMV(Y),ns)
      REAL_T H(DIMV(H),ns)
      REAL_T Vx(DIMV(Vx))
      REAL_T Vy(DIMV(Vy))
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T lambda(DIMV(lambda))
      REAL_T CT(DIMV(CT), ns+1)
      REAL_T dx(SDIM)
     
      REAL_T ht_xp
      REAL_T ht_xm
      REAL_T ht_yp
      REAL_T ht_ym
      REAL_T H_xp(maxspec)
      REAL_T H_xm(maxspec)
      REAL_T H_yp(maxspec)
      REAL_T H_ym(maxspec)
      REAL_T T_xp
      REAL_T T_xm
      REAL_T T_yp
      REAL_T T_ym
      REAL_T rho_xp
      REAL_T rho_xm
      REAL_T rho_yp
      REAL_T rho_ym
      REAL_T Cpmix_xp
      REAL_T Cpmix_xm
      REAL_T Cpmix_yp
      REAL_T Cpmix_ym
      REAL_T Y_xp(maxspec)
      REAL_T Y_xm(maxspec)
      REAL_T Y_yp(maxspec)
      REAL_T Y_ym(maxspec)
      REAL_T U_xp
      REAL_T U_xm
      REAL_T V_yp
      REAL_T V_ym
      REAL_T RhoD_xp(maxspec)
      REAL_T RhoD_xm(maxspec)
      REAL_T RhoD_yp(maxspec)
      REAL_T RhoD_ym(maxspec)
      REAL_T lambda_xp,lambda_xm,lambda_yp,lambda_ym
      REAL_T dxInv2, dyInv2
      integer i,j, n

      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            
            do n=1,ns
            
               RhoD_xp(n) = 0.5*( RhoD(i,j,n) + RhoD(i+1,j,n) )
               RhoD_xm(n) = 0.5*( RhoD(i,j,n) + RhoD(i-1,j,n) )
               RhoD_yp(n) = 0.5*( RhoD(i,j,n) + RhoD(i,j+1,n) )
               RhoD_ym(n) = 0.5*( RhoD(i,j,n) + RhoD(i,j-1,n) )
               
               H_xm(n) = 0.5*( H(i,j,n) + H(i-1,j,n) )
               H_xp(n) = 0.5*( H(i,j,n) + H(i+1,j,n) )
               H_yp(n) = 0.5*( H(i,j,n) + H(i,j+1,n) )
               H_ym(n) = 0.5*( H(i,j,n) + H(i,j-1,n) )

               Y_xm(n) = 0.5*( Y(i,j,n) + Y(i-1,j,n) )
               Y_xp(n) = 0.5*( Y(i,j,n) + Y(i+1,j,n) )
               Y_yp(n) = 0.5*( Y(i,j,n) + Y(i,j+1,n) )
               Y_ym(n) = 0.5*( Y(i,j,n) + Y(i,j-1,n) )
            enddo
            
            rho_xp = 0.5*( Rho(i,j) + Rho(i+1,j) )
            rho_xm = 0.5*( Rho(i,j) + Rho(i-1,j) )
            rho_yp = 0.5*( Rho(i,j) + Rho(i,j+1) )
            rho_ym = 0.5*( Rho(i,j) + Rho(i,j-1) )
           
            T_xp = 0.5*( T(i,j) + T(i+1,j) )
            T_xm = 0.5*( T(i,j) + T(i-1,j) )
            T_yp = 0.5*( T(i,j) + T(i,j+1) )
            T_ym = 0.5*( T(i,j) + T(i,j-1) )
           
            Cpmix_xp = 0.5*( Cpmix(i,j) + Cpmix(i+1,j) )
            Cpmix_xm = 0.5*( Cpmix(i,j) + Cpmix(i-1,j) )
            Cpmix_yp = 0.5*( Cpmix(i,j) + Cpmix(i,j+1) )
            Cpmix_ym = 0.5*( Cpmix(i,j) + Cpmix(i,j-1) )
          
            lambda_xm = 0.5*( lambda(i,j) + lambda(i-1,j) )
            lambda_xp = 0.5*( lambda(i,j) + lambda(i+1,j) )
            lambda_yp = 0.5*( lambda(i,j) + lambda(i,j+1) )
            lambda_ym = 0.5*( lambda(i,j) + lambda(i,j-1) )
           
            ht_xp = 0.5*( Ht(i,j) + Ht(i+1,j) )
            ht_xm = 0.5*( Ht(i,j) + Ht(i-1,j) )
            ht_yp = 0.5*( Ht(i,j) + Ht(i,j+1) )
            ht_ym = 0.5*( Ht(i,j) + Ht(i,j-1) )

            U_xp = 0.5*( Vx(i,j) + Vx(i+1,j) )
            U_xm = 0.5*( Vx(i,j) + Vx(i-1,j) )

            V_yp = 0.5*( Vy(i,j) + Vy(i,j+1) )
            V_ym = 0.5*( Vy(i,j) + Vy(i,j-1) )
            
            do n=1,ns
               CT(i,j,n) = 
     &              ( RhoD_xp(n)*(Y(i+1,j,n)-Y(i  ,j,n))
     &              - RhoD_xm(n)*(Y(i  ,j,n)-Y(i-1,j,n)) ) * dxInv2 +
     &              ( RhoD_yp(n)*(Y(i,j+1,n)-Y(i,j  ,n))
     &              - RhoD_ym(n)*(Y(i,j  ,n)-Y(i,j-1,n)) ) * dyInv2
               
               CT(i,j,ns+n) = 
     &              (rho_xp*Y_xp(n)*U_xp - rho_xm*Y_xm(n)*U_xm)/dx(1) + 
     &              (rho_yp*Y_yp(n)*V_yp - rho_ym*Y_ym(n)*V_ym)/dx(2)  
               
               CT(i,j, 2*ns+3) = CT(i,j,2*ns+3) + 
     &              RhoD(i,j,n)*((H(i+1,j,n)-H(i-1,j,n))*(Y(i+1,j,n)-Y(i-1,j,n))*0.25*dxInv2 +
     &              (H(i,j+1,n)-H(i,j-1,n))*(Y(i,j+1,n)-Y(i,j-1,n))*0.25*dyInv2)
            end do
            
            CT(i,j, 2*ns+1) =
     &           (rho_xp*Cpmix_xp*T_xp*U_xp - rho_xm*Cpmix_xm*T_xm*U_xm)/dx(1) + 
     &           (rho_yp*Cpmix_yp*T_yp*V_yp - rho_ym*Cpmix_ym*T_ym*V_ym)/dx(2) 
            
            CT(i,j,2*ns+2) =
     &           ( lambda_xp*(T(i+1,j)-T(i  ,j))
     &           - lambda_xm*(T(i  ,j)-T(i-1,j)) ) *dxInv2 +
     &           ( lambda_yp*(T(i,j+1)-T(i,j  ))
     &           - lambda_ym*(T(i,j  )-T(i,j-1)) ) *dyInv2
           
            CT(i,j,2*ns+4) = 0.25*(rho_xp+rho_xm+rho_yp+rho_yp)*0.25*(Cpmix_xp+ Cpmix_xm+
     &           Cpmix_yp +  Cpmix_ym)


         enddo
      end do
      end
      
