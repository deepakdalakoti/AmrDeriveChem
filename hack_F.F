#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "hack_F.H"
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
      x0 = 0.0d0
      y0 = 0.0d0
      z0 = 0.0d0
      do k=lo(3),hi(3)
         z = (k+0.5d0)*dx(3) + plo(3)
         do j=lo(2),hi(2)
            y = (j+0.5d0)*dx(2) + plo(2)
            do i=lo(1),hi(1)
               x = (i+0.5d0)*dx(1) + plo(1)
c               U(i,j,k) = x
               U(i,j,k) = SQRT( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
            end do
         end do
      end do
#endif
      end

