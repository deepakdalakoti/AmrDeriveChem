#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "gradT_F.H"
#include "BC_TYPES.H"
#include "SPACE.H"
      
#define SDIM BL_SPACEDIM

      subroutine FORT_HACKVAL(lo, hi, U, DIMS(U), dx, plo)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(U)
      REAL_T U(DIMV(U))
      REAL_T dx(SDIM),plo(SDIM)
      integer i,j
      REAL_T x,y,x0,y0
#if BL_SPACEDIM==2
      x0 = 0.d0
      y0 = 0.d0
c      do j=lo(2),hi(2)
      do j=ARG_L2(U),ARG_H2(U)
         y = (j+0.5d0)*dx(2) + plo(2)
c         do i=lo(1),hi(1)
         do i=ARG_L1(U),ARG_H1(U)
            x = (i+0.5d0)*dx(1) + plo(1)
            U(i,j) = y
c            U(i,j) = SQRT( (x-x0)**2 + (y-y0)**2 )
         end do
      end do
#else
      integer k
      REAL_T z,z0
      x0 = 0.07374d0
      y0 = 0.0d0
      z0 = 0.0255d0
      do k=lo(3),hi(3)
         z = (k+0.5d0)*dx(3) + plo(3)
         do j=lo(2),hi(2)
            y = (j+0.5d0)*dx(2) + plo(2)
            do i=lo(1),hi(1)
               x = (i+0.5d0)*dx(1) + plo(1)
               U(i,j,k) = x
c               U(i,j,k) = SQRT( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
            end do
         end do
      end do
#endif
      end

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

      subroutine FORT_GRAD(lo, hi, T, DIMS(T), gr, DIMS(gr), delta)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(gr)
      REAL_T    T(DIMV(T))
      REAL_T   gr(DIMV(gr),SDIM+1) 
      REAL_T delta(SDIM)
      integer i,j,k,n

      REAL_T tdxI(3)
      
      do i=1,BL_SPACEDIM
         tdxI(i) =  one / (two*delta(i))
      enddo

#if BL_SPACEDIM==2
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            
            gr(i,j,1) = tdxI(1)*(T(i+1,j) - T(i-1,j))
            gr(i,j,2) = tdxI(2)*(T(i,j+1) - T(i,j-1))
            gr(i,j,3) = SQRT(gr(i,j,1)**2 + gr(i,j,2)**2)

         end do
      end do
#else
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               gr(i,j,k,1) = tdxI(1)*(T(i+1,j,k) - T(i-1,j,k))
               gr(i,j,k,2) = tdxI(2)*(T(i,j+1,k) - T(i,j-1,k))
               gr(i,j,k,3) = tdxI(3)*(T(i,j,k+1) - T(i,j,k-1))
               gr(i,j,k,4) = SQRT(gr(i,j,k,1)**2 + gr(i,j,k,2)**2 + gr(i,j,k,3)**2)

            end do
         end do
      end do
#endif
      end

