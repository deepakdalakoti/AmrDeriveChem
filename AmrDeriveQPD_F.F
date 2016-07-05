#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "SPACE.H"
#include "AmrDeriveQPD_F.H"

#define SDIM BL_SPACEDIM
 
   
      subroutine FORT_CONDITIONEDRR(lo, hi, X, DIMS(X), F, DIMS(F), B, DIMS(B), lowval, upval, Nreacs)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(X)
      REAL_T   X(DIMV(X))
      integer DIMDEC(F)
      REAL_T   F(DIMV(F), Nreacs)
      integer DIMDEC(B)
      REAL_T   B(DIMV(B), Nreacs)
      REAL_T   lowval, upval
      integer i, j, Nreacs, n
      
      do i=lo(1),hi(1)
         do j=lo(2),hi(2)
            if ((X(i,j) .GE. lowval) .AND. (X(i,j) .LT. upval)) then
c     retain original values                 
            else
c     reset the FWR and RVS reaction rates to be zero
               do n=1, Nreacs
                  F(i, j, n) = 0.0
                  B(i, j, n) = 0.0
               enddo
            endif
         enddo
      enddo
      end
      
    
