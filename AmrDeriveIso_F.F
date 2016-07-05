#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "AmrDeriveIso_F.H"
#include "SPACE.H"
      
#define SDIM BL_SPACEDIM

      subroutine FORT_SETLOC(lo, hi, U, DIMS(U), dx, domnlo)
      implicit none
      integer  lo(SDIM),  hi(SDIM)
      integer DIMDEC(U)
      REAL_T U(DIMV(U),SDIM)
      REAL_T dx(SDIM), domnlo(SDIM)
      integer i,j,k
      REAL_T y,z
#if BL_SPACEDIM==2
      do j=lo(2),hi(2)
         y = (j + 0.5d0)*dx(2) + domnlo(2)
         do i=lo(1),hi(1)
            U(i,j,1) = (i + 0.5d0)*dx(1) + domnlo(1)
            U(i,j,2) = y
         enddo
      enddo
#else
      do k=lo(3),hi(3)
         z = (k + 0.5d0)*dx(3) + domnlo(3)
         do j=lo(2),hi(2)
            y = (j + 0.5d0)*dx(2) + domnlo(2)
            do i=lo(1),hi(1)
               U(i,j,k,1) = (i + 0.5d0)*dx(1) + domnlo(1)
               U(i,j,k,2) = y
               U(i,j,k,3) = z
            enddo
         enddo
      enddo
#endif

      end

      subroutine FORT_SETCLOC(lo, hi, U, DIMS(U), fdx, domnlo, ratio)
      implicit none
      integer  lo(SDIM),  hi(SDIM)
      integer DIMDEC(U)
      REAL_T U(DIMV(U),SDIM)
      REAL_T fdx(SDIM), domnlo(SDIM)
      integer ratio
      integer i,j,k,ci,cj,ck
      REAL_T y,z

#if BL_SPACEDIM==2
      do j=lo(2),hi(2)
         cj = j/ratio
         if (j.lt.0) cj = cj-1
         y = (cj + 0.5d0)*fdx(2)*ratio + domnlo(2)
         do i=lo(1),hi(1)
            ci = i/ratio
            if (i.lt.0) ci = ci-1
            U(i,j,1) = (ci + 0.5d0)*fdx(1)*ratio + domnlo(1)
            U(i,j,2) = y
         enddo
      enddo
#else
      do k=lo(3),hi(3)
         ck = k/ratio
         if (k.lt.0) ck = ck-1
         z = (ck + 0.5d0)*fdx(3)*ratio + domnlo(3)
         do j=lo(2),hi(2)
            cj = j/ratio
            if (j.lt.0) cj = cj-1
            y = (cj + 0.5d0)*fdx(2)*ratio + domnlo(2)
            do i=lo(1),hi(1)
               ci = i/ratio
               if (i.lt.0) ci = ci-1
               U(i,j,k,1) = (ci + 0.5d0)*fdx(1)*ratio + domnlo(1)
               U(i,j,k,2) = y
               U(i,j,k,3) = z
            enddo
         enddo
      enddo
#endif
      end

