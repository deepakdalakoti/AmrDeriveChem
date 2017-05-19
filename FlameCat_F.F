#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "FlameCat_F.H"
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

      subroutine FORT_GRAD(lo, hi, dat, DIMS(dat), delta,
     &                OutBox,DIMS(OutBox),lo_thresh,hi_thresh,Z_st,PerHRR)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(dat)
      integer DIMDEC(OutBox)
      REAL_T    dat(DIMV(dat),5)
      REAL_T    delta(SDIM)
      REAL_T    OutBox(DIMV(OutBox)), PerHRR(5)
      integer i,j,k,n
      REAL_T small, lo_thresh,hi_thresh,Z_st 
      REAL_T tdxI(3)
      
      do i=1,BL_SPACEDIM
         tdxI(i) =  one / (two*delta(i))
      enddo
      small = 1e-7
#if BL_SPACEDIM==2
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            
            gr(i,j,1) = tdxI(1)*(Fuel(i+1,j) -Fuel(i-1,j))
            gr(i,j,2) = tdxI(2)*(Fuel(i,j+1) - Fuel(i,j-1))
            gr(i,j,3) = SQRT(gr(i,j,1)**2 + gr(i,j,2)**2)

            gr_ox(i,j,1) = tdxI(1)*(Ox(i+1,j) -Ox(i-1,j))
            gr_ox(i,j,2) = tdxI(2)*(Ox(i,j+1) - Ox(i,j-1))
            gr_ox(i,j,3) = SQRT(gr_ox(i,j,1)**2 + gr_ox(i,j,2)**2)
            
            OutBox(i,j) = (gr(i,j,1)*gr_ox(i,j,1)+gr(i,j,2)
     &               *gr_ox(i,j,2))/(gr(i,j,3)*gr_ox(i,j,3)+small)
         end do
      end do
#else
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
            
              if(dat(i,j,k,4).gt.1e9) then
                  PerHRR(1) = PerHRR(1)+dat(i,j,k,4)/1e9
                if(dat(i,j,k,1).gt.lo_thresh.and.dat(i,j,k,5).lt.1200) then
                  PerHRR(2) = PerHRR(2)+dat(i,j,k,4)/1e9  
                 OutBox(i,j,k)= 1
                elseif(dat(i,j,k,1).lt.lo_thresh.and.dat(i,j,k,2).lt.
     &                    hi_thresh.and.dat(i,j,k,3).gt.Z_st) then
                 OutBox(i,j,k)=2
                 PerHRR(3) = PerHRR(3)+dat(i,j,k,4)/1e9
                elseif(dat(i,j,k,1).lt.lo_thresh.and.dat(i,j,k,2).lt.
     &             hi_thresh.and.dat(i,j,k,3).lt.Z_st) then
                 PerHRR(4) = PerHRR(4)+dat(i,j,k,4)/1e9
                 OutBox(i,j,k)=3
                elseif(dat(i,j,k,1).lt.lo_thresh.and.dat(i,j,k,2).gt.hi_thresh
     &            ) then
                PerHRR(5)= dat(i,j,k,4)/1e9+PerHRR(5)
                OutBox(i,j,k)=4
                endif
               else
                OutBox(i,j,k)=0
            
             endif 
           end do
         end do
      end do


c      do k = lo(3),hi(3)
c         do j = lo(2),hi(2)
c            do i = lo(1),hi(1)
c            
c              if(dat(i,j,k,4).gt.1e8) then
c                  PerHRR(1) = PerHRR(1)+dat(i,j,k,4)/1e9
c                if(dat(i,j,k,5).lt.1200) then
c                  PerHRR(2) = PerHRR(2)+dat(i,j,k,4)/1e9  
c                 OutBox(i,j,k)= 1
c                elseif(dat(i,j,k,5).gt.
c     &                    1200.and.dat(i,j,k,3).gt.(Z_st+0.01)) then
c                 OutBox(i,j,k)=2
c                 PerHRR(3) = PerHRR(3)+dat(i,j,k,4)/1e9
c                elseif(dat(i,j,k,5).gt.
c     &             1200.and.dat(i,j,k,3).lt.(Z_st-0.01)) then
c                 PerHRR(4) = PerHRR(4)+dat(i,j,k,4)/1e9
c                 OutBox(i,j,k)=3
c                else
c                 PerHRR(5)= dat(i,j,k,4)/1e9+PerHRR(5)
c                 OutBox(i,j,k)=4
c                endif
c               else
c               OutBox(i,j,k)=-1
            
c             endif 
c           end do
c         end do
c      end do



#endif
      end

