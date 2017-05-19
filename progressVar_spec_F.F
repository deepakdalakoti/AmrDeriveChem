#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "progressVar_spec_F.H"
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
      subroutine FORT_COMPDIFFTERMS(lo, hi, Y, DIMS(Y),
     &     RhoD, DIMS(RhoD), DT, DIMS(DT),
     &     dx, ns,ind) 
     
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
  
      
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      
      integer DIMDEC(DT)
     
      REAL_T Y(DIMV(Y),ns)
      REAL_T RhoD_xp, RhoD_xm
      REAL_T RhoD_yp, RhoD_ym     
      REAL_T RhoD(DIMV(RhoD),ns)
       
      REAL_T DT(DIMV(DT))
      REAL_T dx(SDIM)
      integer ind
      integer i,j,n,p,s
      REAL_T small = 1e-7;
      REAL_T dxInv2, dyInv2

      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
   
      n= ind+1
       
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
               RhoD_xp = 0.5*( RhoD(i,j,n) + RhoD(i+1,j,n) )
               RhoD_xm = 0.5*( RhoD(i,j,n) + RhoD(i-1,j,n) )
               RhoD_yp = 0.5*( RhoD(i,j,n) + RhoD(i,j+1,n) )
               RhoD_ym = 0.5*( RhoD(i,j,n) + RhoD(i,j-1,n) )
               
            enddo
          
           DT(i,j) = 
     &              ( RhoD_xp*(Y(i+1,j,k,n)-Y(i  ,j,k,n))
     &              - RhoD_xm*(Y(i  ,j,k,n)-Y(i-1,j,k,n)) ) * dxInv2 
     &             +  ( RhoD_yp*(Y(i,j+1,k,n)-Y(i,j  ,k,n))
     &              - RhoD_ym*(Y(i,j  ,k,n)-Y(i,j-1,k,n)) ) * dyInv2
           DT(i,j) =-1.0/(small+DT(i,j))
    
         end do
      end do
      end
#else
      subroutine FORT_COMPDIFFTERMS(lo, hi, T, DIMS(T),
     &     outs, maxim)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer sz
      integer DIMDEC(T)
      
      REAL_T T(DIMV(T))
      REAL_T outs(DIMV(T))
      integer i,j,k,l
      REAL_T maxim
           print*, maxim
           do k=lo(3), hi(3)
             do j= lo(2), hi(2)
              do i= lo(1),hi(1)

                outs(i,j,k) = T(i,j,k)/maxim
               
              enddo
             enddo
           enddo









            end
#endif
