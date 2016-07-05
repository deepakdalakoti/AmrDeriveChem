#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "summassCyl_F.H"
#include "SPACE.H"
      
#define SDIM BL_SPACEDIM

      double precision function FORT_SUM_IN_RZ(lo, hi, U, DIMS(U), dx, plo, rMax, zMax)
      implicit none
      integer  lo(SDIM),  hi(SDIM)
      integer DIMDEC(U)
      REAL_T U(DIMV(U))
      REAL_T dx(SDIM), plo(SDIM), rMax, zMax
      integer i,j,k
      REAL_T x,y,z,r
      
      FORT_SUM_IN_RZ = 0.d0

#if BL_SPACEDIM==2
      do j=lo(2),hi(2)
         y = (j + 0.5)*dx(2) + plo(2)
         do i=lo(1),hi(1)
            x = (i + 0.5)*dx(1) + plo(1)
            r = SQRT(x*x + y*y)
            if (r.le.rMax) then
               FORT_SUM_IN_RZ = FORT_SUM_IN_RZ + U(i,j)
            endif
         enddo
      enddo
#else
      do k=lo(3),hi(3)
         z = (k + 0.5)*dx(3) + plo(3)
         if (z.le.zMax) then
            do j=lo(2),hi(2)
               y = (j + 0.5)*dx(2) + plo(2)
               do i=lo(1),hi(1)
                  x = (i + 0.5)*dx(1) + plo(1)
                  r = SQRT(x*x + y*y)
                  if (r.le.rMax) then
                     FORT_SUM_IN_RZ = FORT_SUM_IN_RZ + U(i,j,k)
                  endif
               enddo
            enddo
         endif
      enddo
#endif

      end

