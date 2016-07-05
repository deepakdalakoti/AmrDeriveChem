#include "REAL.H"
#include "ArrayLim.H"
#include "getReactionRates_F.H"

#define SDIM 2
 
      
      subroutine FORT_GETSUMRR(lo,hi, FWD,DIMS(FWD), 
     &     REV,DIMS(REV),SUM,DIMS(SUM), Nreacs)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(FWD)
      integer DIMDEC(REV)
      integer DIMDEC(SUM)
      REAL_T   FWD(DIMV(FWD),Nreacs)
      REAL_T   REV(DIMV(REV),Nreacs)
      REAL_T   SUM(DIMV(SUM),Nreacs)
      integer Nreacs
      integer i, j, n
      
      do i=lo(1),hi(1)
         do j=lo(2),hi(2)
            do n=1, Nreacs
               SUM(i,j, n) = FWD(i,j,n) - REV(i,j,n)
            enddo
         enddo
      enddo
      end
      
      subroutine FORT_GETSUMHEATRELEASE(lo,hi, HR,DIMS(HR), 
     &     T,DIMS(T), minT, maxT, deltaT, Ns, Sum)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(HR)
      integer DIMDEC(T)

      REAL_T   T(DIMV(T))
      REAL_T   HR(DIMV(HR))
      REAL_T   minT, maxT
      REAL_T   deltaT
      integer Ns
      integer i, j, n, ntot
      REAL_T  integr(Ns)
      REAL_T  Sum
      integer icount(Ns)
      Sum = 0.d0
      ntot = 0
      do n=1, Ns
         integr(n) = 0.d0
         icount(n) = 0
      enddo
      
      do i=lo(1),hi(1)
         do j=lo(2),hi(2)
            do n=1, Ns-1
               if ((T(i,j) .GE. (minT+(n-1)*deltaT)) .AND. (T(i,j)
     &          .LT. (minT+n*deltaT))) then
                  integr(n) = integr(n) + HR(i,j)
                  icount(n) = icount(n) +1
               endif
            enddo
         enddo
      enddo
      do n=1, Ns-1
         if (icount(n) .GT. 0) then
         integr(n) = integr(n)/icount(n)
         Sum = Sum + integr(n)*deltaT
         ntot = ntot +icount(n)
         endif
      enddo
c      print *, 'integral heat release in temperature space ', Sum
      end

   

      subroutine FORT_GETPDF(lo,hi, X,DIMS(X), min, max, delta, Nt)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(X)
      REAL_T   X(DIMV(X))
      REAL_T   min, max, delta
      integer Nt
      integer i, j, n, ntot
      integer icount(Nt)
  
      ntot = 0
      do n=1, Nt
         icount(n) = 0
      enddo
      
      do i=lo(1),hi(1)
         do j=lo(2),hi(2)
            do n=1, Nt-1
               if ((X(i,j) .GE. (min+(n-1)*delta)) .AND. (X(i,j)
     &          .LT. (min+n*delta))) then
                  icount(n) = icount(n) +1
               endif
            enddo
         enddo
      enddo
      do n=1, Nt-1
         ntot = ntot +icount(n)
      enddo
      print *, 'VARIABLES ="T","pdf"'
      print *, 'ZONE I= ', Nt, 'DATAPACKING=POINT'
      do n=1, Nt
         print*, min+(n-1)*delta, icount(n)/float(ntot)
      enddo
c      print *, 'ntot=',  ntot
      end
