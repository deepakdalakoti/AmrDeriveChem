#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "compDiffTerms_F.H"
#include "BC_TYPES.H"
#include "SPACE.H"
      
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
      REAL_T DT(DIMV(DT),ns+1)
      REAL_T dx(SDIM)
      REAL_T H_xp(maxspec), RhoD_xp(maxspec)
      REAL_T H_xm(maxspec), RhoD_xm(maxspec)
      REAL_T H_yp(maxspec), RhoD_yp(maxspec)
      REAL_T H_ym(maxspec), RhoD_ym(maxspec)
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
               
               DT(i,j,ns+1) = DT(i,j,ns+1) + 
     &              ( RhoD_xp(n)*H_xp(n)*(Y(i+1,j,n)-Y(i  ,j,n))
     &              - RhoD_xm(n)*H_xm(n)*(Y(i  ,j,n)-Y(i-1,j,n)) ) * dxInv2 +
     &              ( RhoD_yp(n)*H_yp(n)*(Y(i,j+1,n)-Y(i,j  ,n))
     &              - RhoD_ym(n)*H_ym(n)*(Y(i,j  ,n)-Y(i,j-1,n)) ) * dyInv2
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
      REAL_T DT(DIMV(DT),ns+1)
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

                  DT(i,j,k,ns+1) = DT(i,j,k,ns+1) + 
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
