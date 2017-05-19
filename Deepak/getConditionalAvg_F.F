#include <REAL.H>
#include <ArrayLim.H>
#include <CONSTANTS.H>
#include <SPACE.H>
#include <getJPDF_F.H>

#define SDIM BL_SPACEDIM
 
#if SDIM==2   
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
      do n=1, Nt
         ntot = ntot +icount(n)
      enddo
      print *, 'VARIABLES ="T","pdf"'
      print *, 'ZONE I= ', Nt, 'DATAPACKING=POINT'
      do n=1, Nt
         print*, min+(n-1)*delta, icount(n)/float(ntot)
      enddo
c      print *, 'ntot=',  ntot
      end

      subroutine FORT_GETJPDF(lo,hi, X,DIMS(X),  Y,DIMS(Y), Bins,DIMS(Bins), minX, maxX, 
     &     deltaX, NX, minY, maxY, deltaY, NY, vol,mask,DIMS(mask))
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(X)
      integer DIMDEC(Y)
      integer DIMDEC(Bins)
      integer DIMDEC(mask)
      REAL_T   X(DIMV(X))
      REAL_T   Y(DIMV(Y))
      REAL_T   mask(DIMV(mask))
      REAL_T   Bins(DIMV(Bins))
      REAL_T   minX, maxX, deltaX
      REAL_T   minY, maxY, deltaY
      REAL_T   vol
      integer NX, NY
      integer i, j, n, m
      print*, "DIMV" , Bins_l1     
      do i=lo(1),hi(1)
         do j=lo(2),hi(2)
            do n=1, NX
               if ((X(i,j) .GE. (minX+(n-1)*deltaX)) .AND. (X(i,j)
     &              .LT. (minX+n*deltaX))) then
                    if(mask(i,j).gt.0) then
                      Bins(n-1,0) = Bins(n-1,0)+ vol
                      Bins(n-1,1) = Bins(n-1,1) +Y(i,j)*vol
                    endif
c                      Bins(ARG_L1(Bins)+n,ARG_L2(Bins)+m,ARG_L3(Bins)) = 
c     &                       Bins(ARG_L1(Bins)+n,ARG_L2(Bins)+m,ARG_L3(Bins)) + vol
                        
                   
               endif
               if(X(i,j).eq.(minX+nX*deltaX)) then
                  Bins(NX-1,0) = Bins(NX-1,0)+vol
                  Bins(NX-1,1) = Bins(NX-1,1) +Y(i,j)*vol
               endif
            enddo
         enddo
      enddo
      end
#else
      subroutine FORT_GETPDF(lo,hi, X,DIMS(X), min, max, delta, Nt)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(X)
      REAL_T   X(DIMV(X))
      REAL_T   min, max, delta
      integer Nt
      integer i, j, k, n, ntot
      integer icount(Nt)
  
      ntot = 0
      do n=1, Nt
         icount(n) = 0
      enddo
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1, Nt-1
                  if ((X(i,j,k) .GE. (min+(n-1)*delta))
     &                 .AND. (X(i,j,k)  .LT. (min+n*delta))) then
                     icount(n) = icount(n) +1
                  endif
               enddo
            enddo
         enddo
      enddo
      do n=1, Nt
         ntot = ntot +icount(n)
      enddo
      print *, 'VARIABLES ="T","pdf"'
      print *, 'ZONE I= ', Nt, 'DATAPACKING=POINT'
      do n=1, Nt
         print*, min+(n-1)*delta, icount(n)/float(ntot)
      enddo
c      print *, 'ntot=',  ntot
      end

      subroutine FORT_GETJPDF(lo,hi, X,DIMS(X),  Y,DIMS(Y), Bins,DIMS(Bins), minX, maxX, 
     &     deltaX, NX, minY, maxY, deltaY, NY, vol,mask,DIMS(mask))
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(X)
      integer DIMDEC(Y)
      integer DIMDEC(Bins)
      integer DIMDEC(mask)
      REAL_T   X(DIMV(X))
      REAL_T   Y(DIMV(Y))
      REAL_T   BinAr(nX)
      REAL_T   mask(DIMV(mask))
      REAL_T   Bins(DIMV(Bins))
      REAL_T   minX, maxX, deltaX
      REAL_T   minY, maxY, deltaY
      REAL_T   vol
      integer NX, NY
      integer i, j, k, n, m
   
        
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1, NX
                  if ((X(i,j,k) .GE. (minX+(n-1)*deltaX)) .AND.
     &                 (X(i,j,k) .LT. (minX+n*deltaX))) then
                         if(mask(i,j,k).gt.0) then 
c                         if(X(i,j,k) > 930)  then
c                    print*, X(i,j,k), minX+(n-1)*deltaX,minX+n*deltaX,n
c                         endif
                          Bins(n-1,0,:) = Bins(n-1,0,:)+ vol
                          Bins(n-1,1,:) = Bins(n-1,1,:)+ vol*Y(i,j,k)
c                          if (n>53) print*, Bins(n,1,ARG_L3(Bins)), n
                          
                         endif
c                             Bins(ARG_L1(Bins)+n,ARG_L2(Bins)+m,ARG_L3(Bins)) = 
c     &                       Bins(ARG_L1(Bins)+n,ARG_L2(Bins)+m,ARG_L3(Bins)) + vol
                       
                  endif
                  if(X(i,j,k).eq.maxX) then
                    Bins(NX-1,0,:) = Bins(NX-1,0,:)+vol
                    Bins(NX-1,1,:) = Bins(NX-1,1,:) +vol*Y(i,j,k)
                  endif
               enddo
            enddo
         enddo
      enddo
      end
      
#endif      
