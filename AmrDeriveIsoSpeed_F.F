#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "AmrDeriveIsoSpeed_F.H"
#include "SPACE.H"
#include "BC_TYPES.H"
      
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

          subroutine FORT_COMPDIFFTERMS(lo, hi, dx,indices,GridPts,nNodes,datain,DIMS(datain),auxcomps,nComps,sz)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(datain)
      integer nNodes,nComps,sz,auxcomps
      integer flo(SDIM),fhi(SDIM)
      REAL_T datain(DIMV(datain),auxcomps), dx(3)
      REAL_T indices(sz),GridPts(nNodes,nComps)
      REAL_T dxInv2,dyInv2, dzInv2,dx2,dy2,dz2,small
      integer i,j,k
      logical ok
      REAL_T grad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4)
      REAL_T vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      REAL_T pos(3),OutP(2),plo(3)
      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      dzInv2 = 1.d0 / (dx(3)*dx(3))
      dx2= 0.5/dx(1)
      dy2 = 0.5/dx(2)
      dz2 = 0.5/dx(3)
      small = 1e-3
      plo(1) = -0.00256
      plo(2) = -0.00256
      plo(3) = 0
      print*, 'size ' , sz, indices(1)
      if (sz.eq.0) return
      do k = lo(3),hi(3)
        do j= lo(2), hi(2)
         do i= lo(1), hi(1)
          
            grad(i,j,k,1)= (datain(i+1,j,k,4)- datain(i-1,j,k,4))*dx2
            grad(i,j,k,2)= (datain(i,j+1,k,4)- datain(i,j-1,k,4))*dy2
            grad(i,j,k,3)= (datain(i,j,k+1,4)- datain(i,j,k-1,4))*dz2
            grad(i,j,k,4) = sqrt(grad(i,j,k,1)**2 + grad(i,j,k,2)**2 + grad(i,j,k,3)**2)
            grad(i,j,k,1) = grad(i,j,k,1)/(grad(i,j,k,4)+small)
            grad(i,j,k,1) = grad(i,j,k,1)/(grad(i,j,k,4)+small)
            grad(i,j,k,1) = grad(i,j,k,1)/(grad(i,j,k,4)+small)
            vel(i,j,k) = datain(i,j,k,1)*grad(i,j,k,1)+datain(i,j,k,2)*grad(i,j,k,2)+
     &                    datain(i,j,k,3)*grad(i,j,k,3)
          enddo
         enddo
        enddo
      
      do i= 1, sz
       print*, indices(i),'in fortran'
!       pos(1) = GridPts(indices(i),1)
!       pos(2) = GridPts(indices(i),2)
!       pos(3) = GridPts(indices(i),3)
!
!       call ntrpv(pos,vel,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),dx,plo,OutP,2,ok)
!       GridPts(indices(i),nComps-auxcomps+1 +1) = OutP(1)
!       GridPts(indices(i),nComps-auxcomps+1 +2) = OutP(2)


       enddo
!      if(abs(Speed_Out(2)) >40 ) then
!         print*, Speed_Out(2), " HIGH"
!       endif
!       call ntrpv(pos,Y(:,:,:,4),DIMS(Y),dx,plo,Speed_Out,1,ok)
!        if(.not.ok) then
!         print*, 'not ok'
!        endif
!      print*, Speed_Out(1), Speed_Out(2), Speed_Out(3), Speed_Out(4)
!      print*, ok
      end

      subroutine ntrpv(x,g,DIMS(g),dx,plo,u,nc,ok)
      implicit none
      integer nc
      REAL_T x(SDIM),dx(SDIM),plo(SDIM),u(nc)
      integer DIMDEC(g)
      REAL_T g(DIMV(g),nc)
      logical ok
      integer b(SDIM), i
      REAL_T tmp
      REAL_T n(SDIM)
      do i=1,SDIM
         tmp = (x(i) - plo(i)) / dx(i) - 0.5d0
         b(i) = INT( tmp )
         n(i) = ( x(i) - ( (b(i)+0.5d0)*dx(i) + plo(i) ) )/dx(i)
         n(i) = MAX(0.d0,MIN(1.d0,n(i)))
c         if(n(i).lt.0.5) then
c           n(i)= 0.9*n(i)
c         endif
      enddo
      ok = .true.
      if (      b(1).lt.ARG_L1(g) .or. b(1).ge.ARG_H1(g)
     &     .or. b(2).lt.ARG_L2(g) .or. b(2).ge.ARG_H2(g)
#if BL_SPACEDIM==3
     &     .or. b(3).lt.ARG_L3(g) .or. b(3).ge.ARG_H3(g)
#endif
     &     ) then

         print *,'b:',b
         print *,'DIMS:',DIMS(g)

         ok = .false.
         return
      endif

      do i=1,nc
#if BL_SPACEDIM==2
         u(i) =
     &  +       n(1)  *    n(2)   * g(b(1)+1,b(2)+1,i)
     &  +       n(1)  *(1.d0-n(2))* g(b(1)+1,b(2)  ,i)
     &  +  (1.d0-n(1))*    n(2)   * g(b(1)  ,b(2)+1,i)
     &  +  (1.d0-n(1))*(1.d0-n(2))* g(b(1)  ,b(2)  ,i)
#else
         u(i) =
     &  +       n(1)  *    n(2)   *    n(3)    * g(b(1)+1,b(2)+1,b(3)+1,i)
     &  +       n(1)  *(1.d0-n(2))*    n(3)    * g(b(1)+1,b(2)  ,b(3)+1,i)
     &  +       n(1)  *    n(2)   *(1.d0-n(3)) * g(b(1)+1,b(2)+1,b(3)  ,i)
     &  +       n(1)  *(1.d0-n(2))*(1.d0-n(3)) * g(b(1)+1,b(2)  ,b(3)  ,i)
     &  +  (1.d0-n(1))*    n(2)   *    n(3)    * g(b(1)  ,b(2)+1,b(3)+1,i)
     &  +  (1.d0-n(1))*(1.d0-n(2))*    n(3)    * g(b(1)  ,b(2)  ,b(3)+1,i)
     &  +  (1.d0-n(1))*    n(2)   *(1.d0-n(3)) * g(b(1)  ,b(2)+1,b(3)  ,i)
     &  +  (1.d0-n(1))*(1.d0-n(2))*(1.d0-n(3)) * g(b(1)  ,b(2)  ,b(3)  ,i)

#endif
      enddo
      end
      subroutine ntrpv_cubic(x,g,DIMS(g),dx,plo,u,nc,ok)
      implicit none
      integer nc
      REAL_T x(SDIM),dx(SDIM),plo(SDIM),u(nc)
      integer DIMDEC(g)
      REAL_T g(DIMV(g),nc)
      logical ok
      integer b(SDIM), i
      REAL_T n(SDIM)
      REAL_T tester
      integer j,k
      REAL_T tmp

      REAL_T a(64),f(8),dfdx(8),dfdy(8),dfdz(8),d2fdxdy(8),d2fdxdz(8),d2fdydz(8),d3fdxdydz(8)
      logical, save :: first = .true.
      logical reuse_a
      REAL_T, allocatable, save :: a_SAVE(:,:)
      REAL_T, save :: b_SAVE(SDIM)
      integer, save :: a_nc_SAVE
      logical doTricubic, doTrilinear

      do i=1,SDIM         
         tmp = (x(i) - plo(i)) / dx(i) - 0.5d0
         b(i) = FLOOR( tmp )
         n(i) = ( x(i) - ( (b(i)+0.5d0)*dx(i) + plo(i) ) )/dx(i)
         n(i) = MAX(0.d0,MIN(1.d0,n(i)))
      enddo

      ok = .true.

      doTrilinear = .false.
      doTricubic = .true.

#if BL_SPACEDIM==2
      doTrilinear = .true.
      doTricubic = .false.
#endif
      
#if BL_SPACEDIM==3
c     Drop order on boundary.  
      if (doTricubic) then
         if (      b(1).eq.ARG_L1(g) .or. b(1).eq.ARG_H1(g)-1
     &        .or. b(2).eq.ARG_L2(g) .or. b(2).eq.ARG_H2(g)-1
     &        .or. b(3).eq.ARG_L3(g) .or. b(3).eq.ARG_H3(g)-1
     &        ) then

            doTricubic = .false.
            doTrilinear = .true.

         endif 

      endif
         
      if (doTricubic) then
         
         if (      b(1).lt.ARG_L1(g)+1 .or. b(1).gt.ARG_H1(g)-2
     &        .or. b(2).lt.ARG_L2(g)+1 .or. b(2).gt.ARG_H2(g)-2
     &        .or. b(3).lt.ARG_L3(g)+1 .or. b(3).gt.ARG_H3(g)-2
     &        ) then
            
            ok = .false.
            return
         endif

c     Have we done this before
c         reuse_a = (.not. first) .and. all(b(:).eq.b_SAVE(:))
         reuse_a = .false.
         
         if (first) first = .false.
         
         do i=1,nc
            
            if (reuse_a) then
               a(:) = a_SAVE(:,nc)
            else
               f(1) = g(b(1)  ,b(2)  ,b(3)  ,i)
               f(2) = g(b(1)+1,b(2)  ,b(3)  ,i)
               f(3) = g(b(1)  ,b(2)+1,b(3)  ,i)
               f(4) = g(b(1)+1,b(2)+1,b(3)  ,i)
               f(5) = g(b(1)  ,b(2)  ,b(3)+1,i)
               f(6) = g(b(1)+1,b(2)  ,b(3)+1,i)
               f(7) = g(b(1)  ,b(2)+1,b(3)+1,i)
               f(8) = g(b(1)+1,b(2)+1,b(3)+1,i)
               
               dfdx(1) = 0.5d0*( g(b(1)+1,b(2)  ,b(3)  ,i) - g(b(1)-1,b(2)  ,b(3)  ,i) )
               dfdx(2) = 0.5d0*( g(b(1)+2,b(2)  ,b(3)  ,i) - g(b(1)  ,b(2)  ,b(3)  ,i) )
               dfdx(3) = 0.5d0*( g(b(1)+1,b(2)+1,b(3)  ,i) - g(b(1)-1,b(2)+1,b(3)  ,i) )
               dfdx(4) = 0.5d0*( g(b(1)+2,b(2)+1,b(3)  ,i) - g(b(1)  ,b(2)+1,b(3)  ,i) )
               dfdx(5) = 0.5d0*( g(b(1)+1,b(2)  ,b(3)+1,i) - g(b(1)-1,b(2)  ,b(3)+1,i) )
               dfdx(6) = 0.5d0*( g(b(1)+2,b(2)  ,b(3)+1,i) - g(b(1)  ,b(2)  ,b(3)+1,i) )
               dfdx(7) = 0.5d0*( g(b(1)+1,b(2)+1,b(3)+1,i) - g(b(1)-1,b(2)+1,b(3)+1,i) )
               dfdx(8) = 0.5d0*( g(b(1)+2,b(2)+1,b(3)+1,i) - g(b(1)  ,b(2)+1,b(3)+1,i) )
               
               dfdy(1) = 0.5d0*( g(b(1)  ,b(2)+1,b(3)  ,i) - g(b(1)  ,b(2)-1,b(3)  ,i) )
               dfdy(2) = 0.5d0*( g(b(1)+1,b(2)+1,b(3)  ,i) - g(b(1)+1,b(2)-1,b(3)  ,i) )
               dfdy(3) = 0.5d0*( g(b(1)  ,b(2)+2,b(3)  ,i) - g(b(1)  ,b(2)  ,b(3)  ,i) )
               dfdy(4) = 0.5d0*( g(b(1)+1,b(2)+2,b(3)  ,i) - g(b(1)+1,b(2)  ,b(3)  ,i) )
               dfdy(5) = 0.5d0*( g(b(1)  ,b(2)+1,b(3)+1,i) - g(b(1)  ,b(2)-1,b(3)+1,i) )
               dfdy(6) = 0.5d0*( g(b(1)+1,b(2)+1,b(3)+1,i) - g(b(1)+1,b(2)-1,b(3)+1,i) )
               dfdy(7) = 0.5d0*( g(b(1)  ,b(2)+2,b(3)+1,i) - g(b(1)  ,b(2)  ,b(3)+1,i) )
               dfdy(8) = 0.5d0*( g(b(1)+1,b(2)+2,b(3)+1,i) - g(b(1)+1,b(2)  ,b(3)+1,i) )
               
               dfdz(1) = 0.5d0*( g(b(1)  ,b(2)  ,b(3)+1,i) - g(b(1)  ,b(2)  ,b(3)-1,i) )
               dfdz(2) = 0.5d0*( g(b(1)+1,b(2)  ,b(3)+1,i) - g(b(1)+1,b(2)  ,b(3)-1,i) )
               dfdz(3) = 0.5d0*( g(b(1)  ,b(2)+1,b(3)+1,i) - g(b(1)  ,b(2)+1,b(3)-1,i) )
               dfdz(4) = 0.5d0*( g(b(1)+1,b(2)+1,b(3)+1,i) - g(b(1)+1,b(2)+1,b(3)-1,i) )
               dfdz(5) = 0.5d0*( g(b(1)  ,b(2)  ,b(3)+2,i) - g(b(1)  ,b(2)  ,b(3)  ,i) )
               dfdz(6) = 0.5d0*( g(b(1)+1,b(2)  ,b(3)+2,i) - g(b(1)+1,b(2)  ,b(3)  ,i) )
               dfdz(7) = 0.5d0*( g(b(1)  ,b(2)+1,b(3)+2,i) - g(b(1)  ,b(2)+1,b(3)  ,i) )
               dfdz(8) = 0.5d0*( g(b(1)+1,b(2)+1,b(3)+2,i) - g(b(1)+1,b(2)+1,b(3)  ,i) )
               
               d2fdxdy(1) = 0.25d0*( g(b(1)+1,b(2)+1,b(3)  ,i) - g(b(1)-1,b(2)+1,b(3)  ,i)
     &              -                g(b(1)+1,b(2)-1,b(3)  ,i) + g(b(1)-1,b(2)-1,b(3)  ,i) )
               d2fdxdy(2) = 0.25d0*( g(b(1)+2,b(2)+1,b(3)  ,i) - g(b(1)  ,b(2)+1,b(3)  ,i)
     &              -                g(b(1)+2,b(2)-1,b(3)  ,i) + g(b(1)  ,b(2)-1,b(3)  ,i) )
               d2fdxdy(3) = 0.25d0*( g(b(1)+1,b(2)+2,b(3)  ,i) - g(b(1)-1,b(2)+2,b(3)  ,i)
     &              -                g(b(1)+1,b(2)  ,b(3)  ,i) + g(b(1)-1,b(2)  ,b(3)  ,i) )
               d2fdxdy(4) = 0.25d0*( g(b(1)+2,b(2)+2,b(3)  ,i) - g(b(1)  ,b(2)+2,b(3)  ,i)
     &              -                g(b(1)+2,b(2)  ,b(3)  ,i) + g(b(1)  ,b(2)  ,b(3)  ,i) )
               d2fdxdy(5) = 0.25d0*( g(b(1)+1,b(2)+1,b(3)+1,i) - g(b(1)-1,b(2)+1,b(3)+1,i)
     &              -                g(b(1)+1,b(2)-1,b(3)+1,i) + g(b(1)-1,b(2)-1,b(3)+1,i) )
               d2fdxdy(6) = 0.25d0*( g(b(1)+2,b(2)+1,b(3)+1,i) - g(b(1)  ,b(2)+1,b(3)+1,i)
     &              -                g(b(1)+2,b(2)-1,b(3)+1,i) + g(b(1)  ,b(2)-1,b(3)+1,i) )
               d2fdxdy(7) = 0.25d0*( g(b(1)+1,b(2)+2,b(3)+1,i) - g(b(1)-1,b(2)+2,b(3)+1,i)
     &              -                g(b(1)+1,b(2)  ,b(3)+1,i) + g(b(1)-1,b(2)  ,b(3)+1,i) )
               d2fdxdy(8) = 0.25d0*( g(b(1)+2,b(2)+2,b(3)+1,i) - g(b(1)  ,b(2)+2,b(3)+1,i)
     &              -                g(b(1)+2,b(2)  ,b(3)+1,i) + g(b(1)  ,b(2)  ,b(3)+1,i) )
               
               d2fdxdz(1) = 0.25d0*( g(b(1)+1,b(2)  ,b(3)+1,i) - g(b(1)-1,b(2)  ,b(3)+1,i)
     &              -                g(b(1)+1,b(2)  ,b(3)-1,i) + g(b(1)-1,b(2)  ,b(3)-1,i) )
               d2fdxdz(2) = 0.25d0*( g(b(1)+2,b(2)  ,b(3)+1,i) - g(b(1)  ,b(2)  ,b(3)+1,i)
     &              -                g(b(1)+2,b(2)  ,b(3)-1,i) + g(b(1)  ,b(2)  ,b(3)-1,i) )
               d2fdxdz(3) = 0.25d0*( g(b(1)+1,b(2)  ,b(3)+2,i) - g(b(1)-1,b(2)  ,b(3)+2,i)
     &              -                g(b(1)+1,b(2)  ,b(3)  ,i) + g(b(1)-1,b(2)  ,b(3)  ,i) )
               d2fdxdz(4) = 0.25d0*( g(b(1)+2,b(2)  ,b(3)+2,i) - g(b(1)  ,b(2)  ,b(3)+2,i)
     &              -                g(b(1)+2,b(2)  ,b(3)  ,i) + g(b(1)  ,b(2)  ,b(3)  ,i) )
               d2fdxdz(5) = 0.25d0*( g(b(1)+1,b(2)+1,b(3)+1,i) - g(b(1)-1,b(2)+1,b(3)+1,i)
     &              -                g(b(1)+1,b(2)+1,b(3)-1,i) + g(b(1)-1,b(2)+1,b(3)-1,i) )
               d2fdxdz(6) = 0.25d0*( g(b(1)+2,b(2)+1,b(3)+1,i) - g(b(1)  ,b(2)+1,b(3)+1,i)
     &              -                g(b(1)+2,b(2)+1,b(3)-1,i) + g(b(1)  ,b(2)+1,b(3)-1,i) )
               d2fdxdz(7) = 0.25d0*( g(b(1)+1,b(2)+1,b(3)+2,i) - g(b(1)-1,b(2)+1,b(3)+2,i)
     &              -                g(b(1)+1,b(2)+1,b(3)  ,i) + g(b(1)-1,b(2)+1,b(3)  ,i) )
               d2fdxdz(8) = 0.25d0*( g(b(1)+2,b(2)+1,b(3)+2,i) - g(b(1)  ,b(2)+1,b(3)+2,i)
     &              -                g(b(1)+2,b(2)+1,b(3)  ,i) + g(b(1)  ,b(2)+1,b(3)  ,i) )
               
               
               d2fdydz(1) = 0.25d0*( g(b(1)  ,b(2)+1,b(3)+1,i) - g(b(1)  ,b(2)-1,b(3)+1,i)
     &              -                g(b(1)  ,b(2)+1,b(3)-1,i) + g(b(1)  ,b(2)-1,b(3)-1,i) )
               d2fdydz(2) = 0.25d0*( g(b(1)+1,b(2)+1,b(3)+1,i) - g(b(1)+1,b(2)-1,b(3)+1,i)
     &              -                g(b(1)+1,b(2)+1,b(3)-1,i) + g(b(1)+1,b(2)-1,b(3)-1,i) )
               d2fdydz(3) = 0.25d0*( g(b(1)  ,b(2)+2,b(3)+1,i) - g(b(1)  ,b(2)  ,b(3)+1,i)
     &              -                g(b(1)  ,b(2)+2,b(3)-1,i) + g(b(1)  ,b(2)  ,b(3)-1,i) )
               d2fdydz(4) = 0.25d0*( g(b(1)+1,b(2)+2,b(3)+1,i) - g(b(1)+1,b(2)  ,b(3)+1,i)
     &              -                g(b(1)+1,b(2)+2,b(3)-1,i) + g(b(1)+1,b(2)  ,b(3)-1,i) )
               d2fdydz(5) = 0.25d0*( g(b(1)  ,b(2)+1,b(3)+2,i) - g(b(1)  ,b(2)-1,b(3)+2,i)
     &              -                g(b(1)  ,b(2)+1,b(3)  ,i) + g(b(1)  ,b(2)-1,b(3)  ,i) )
               d2fdydz(6) = 0.25d0*( g(b(1)+1,b(2)+1,b(3)+2,i) - g(b(1)+1,b(2)-1,b(3)+2,i)
     &              -                g(b(1)+1,b(2)+1,b(3)  ,i) + g(b(1)+1,b(2)-1,b(3)  ,i) )
               d2fdydz(7) = 0.25d0*( g(b(1)  ,b(2)+2,b(3)+2,i) - g(b(1)  ,b(2)  ,b(3)+2,i)
     &              -                g(b(1)  ,b(2)+2,b(3)  ,i) + g(b(1)  ,b(2)  ,b(3)  ,i) )
               d2fdydz(8) = 0.25d0*( g(b(1)+1,b(2)+2,b(3)+2,i) - g(b(1)+1,b(2)  ,b(3)+2,i)
     &              -                g(b(1)+1,b(2)+2,b(3)  ,i) + g(b(1)+1,b(2)  ,b(3)  ,i) )
               
               d3fdxdydz(1) = 0.125d0*( g(b(1)+1,b(2)+1,b(3)+1,i) - g(b(1)-1,b(2)+1,b(3)+1,i)
     &              -                   g(b(1)+1,b(2)-1,b(3)+1,i) + g(b(1)-1,b(2)-1,b(3)+1,i)
     &              -                   g(b(1)+1,b(2)+1,b(3)-1,i) + g(b(1)-1,b(2)+1,b(3)-1,i)
     &              +                   g(b(1)+1,b(2)-1,b(3)-1,i) - g(b(1)-1,b(2)-1,b(3)-1,i) )
               d3fdxdydz(2) = 0.125d0*( g(b(1)+2,b(2)+1,b(3)+1,i) - g(b(1)  ,b(2)+1,b(3)+1,i)
     &              -                   g(b(1)+2,b(2)-1,b(3)+1,i) + g(b(1)  ,b(2)-1,b(3)+1,i)
     &              -                   g(b(1)+2,b(2)+1,b(3)-1,i) + g(b(1)  ,b(2)+1,b(3)-1,i)
     &              +                   g(b(1)+2,b(2)-1,b(3)-1,i) - g(b(1)  ,b(2)-1,b(3)-1,i) )
               d3fdxdydz(3) = 0.125d0*( g(b(1)+1,b(2)+2,b(3)+1,i) - g(b(1)-1,b(2)+2,b(3)+1,i)
     &              -                   g(b(1)+1,b(2)  ,b(3)+1,i) + g(b(1)-1,b(2)  ,b(3)+1,i)
     &              -                   g(b(1)+1,b(2)+2,b(3)-1,i) + g(b(1)-1,b(2)+2,b(3)-1,i)
     &              +                   g(b(1)+1,b(2)  ,b(3)-1,i) - g(b(1)-1,b(2)  ,b(3)-1,i) )
               d3fdxdydz(4) = 0.125d0*( g(b(1)+2,b(2)+2,b(3)+1,i) - g(b(1)  ,b(2)+2,b(3)+1,i)
     &              -                   g(b(1)+2,b(2)  ,b(3)+1,i) + g(b(1)  ,b(2)  ,b(3)+1,i)
     &              -                   g(b(1)+2,b(2)+2,b(3)-1,i) + g(b(1)  ,b(2)+2,b(3)-1,i)
     &              +                   g(b(1)+2,b(2)  ,b(3)-1,i) - g(b(1)  ,b(2)  ,b(3)-1,i) )
               d3fdxdydz(5) = 0.125d0*( g(b(1)+1,b(2)+1,b(3)+2,i) - g(b(1)-1,b(2)+1,b(3)+2,i)
     &              -                   g(b(1)+1,b(2)-1,b(3)+2,i) + g(b(1)-1,b(2)-1,b(3)+2,i)
     &              -                   g(b(1)+1,b(2)+1,b(3)  ,i) + g(b(1)-1,b(2)+1,b(3)  ,i)
     &              +                   g(b(1)+1,b(2)-1,b(3)  ,i) - g(b(1)-1,b(2)-1,b(3)  ,i) )
               d3fdxdydz(6) = 0.125d0*( g(b(1)+2,b(2)+1,b(3)+2,i) - g(b(1)  ,b(2)+1,b(3)+2,i)
     &              -                   g(b(1)+2,b(2)-1,b(3)+2,i) + g(b(1)  ,b(2)-1,b(3)+2,i)
     &              -                   g(b(1)+2,b(2)+1,b(3)  ,i) + g(b(1)  ,b(2)+1,b(3)  ,i)
     &              +                   g(b(1)+2,b(2)-1,b(3)  ,i) - g(b(1)  ,b(2)-1,b(3)  ,i) )
               d3fdxdydz(7) = 0.125d0*( g(b(1)+1,b(2)+2,b(3)+2,i) - g(b(1)-1,b(2)+2,b(3)+2,i)
     &              -                   g(b(1)+1,b(2)  ,b(3)+2,i) + g(b(1)-1,b(2)  ,b(3)+2,i)
     &              -                   g(b(1)+1,b(2)+2,b(3)  ,i) + g(b(1)-1,b(2)+2,b(3)  ,i)
     &              +                   g(b(1)+1,b(2)  ,b(3)  ,i) - g(b(1)-1,b(2)  ,b(3)  ,i) )
               d3fdxdydz(8) = 0.125d0*( g(b(1)+2,b(2)+2,b(3)+2,i) - g(b(1)  ,b(2)+2,b(3)+2,i)
     &              -                   g(b(1)+2,b(2)  ,b(3)+2,i) + g(b(1)  ,b(2)  ,b(3)+2,i)
     &              -                   g(b(1)+2,b(2)+2,b(3)  ,i) + g(b(1)  ,b(2)+2,b(3)  ,i)
     &              +                   g(b(1)+2,b(2)  ,b(3)  ,i) - g(b(1)  ,b(2)  ,b(3)  ,i) )
               
               call FORT_TRICUBIC_GET_COEFF(a(1),f(1),dfdx(1),dfdy(1),dfdz(1),
     &              d2fdxdy(1),d2fdxdz(1),d2fdydz(1),d3fdxdydz(1))
               
               if (allocated(a_SAVE)) then
                  if (a_nc_SAVE.lt.nc) then
                     deallocate(a_SAVE)
                  endif
               endif
               
               if (.not. allocated(a_SAVE)) then
                  allocate(a_SAVE(64,nc))
               endif
               
               a_SAVE(:,i) = a(:)
               
            endif
            
            call FORT_TRICUBIC_EVAL(a(1),n(1),n(2),n(3),u(i))
            
         enddo
         
         b_SAVE(:) = b(:)

      else if (doTrilinear) then
#endif
         print*, 'in trilinear'  
         if (      b(1).lt.ARG_L1(g) .or. b(1).ge.ARG_H1(g)
     &        .or. b(2).lt.ARG_L2(g) .or. b(2).ge.ARG_H2(g)
#if BL_SPACEDIM==3
     &        .or. b(3).lt.ARG_L3(g) .or. b(3).ge.ARG_H3(g)
#endif
     &        ) then
         
            ok = .false.
            return
         endif

         do i=1,nc
#if BL_SPACEDIM==2
            u(i) =
     &           +       n(1)  *    n(2)   * g(b(1)+1,b(2)+1,i)
     &           +       n(1)  *(1.d0-n(2))* g(b(1)+1,b(2)  ,i)
     &           +  (1.d0-n(1))*    n(2)   * g(b(1)  ,b(2)+1,i)
     &           +  (1.d0-n(1))*(1.d0-n(2))* g(b(1)  ,b(2)  ,i)
#else 
            u(i) =
     &           +       n(1)  *    n(2)   *    n(3)    * g(b(1)+1,b(2)+1,b(3)+1,i)
     &           +       n(1)  *(1.d0-n(2))*    n(3)    * g(b(1)+1,b(2)  ,b(3)+1,i)
     &           +       n(1)  *    n(2)   *(1.d0-n(3)) * g(b(1)+1,b(2)+1,b(3)  ,i)
     &           +       n(1)  *(1.d0-n(2))*(1.d0-n(3)) * g(b(1)+1,b(2)  ,b(3)  ,i)
     &           +  (1.d0-n(1))*    n(2)   *    n(3)    * g(b(1)  ,b(2)+1,b(3)+1,i)
     &           +  (1.d0-n(1))*(1.d0-n(2))*    n(3)    * g(b(1)  ,b(2)  ,b(3)+1,i)
     &           +  (1.d0-n(1))*    n(2)   *(1.d0-n(3)) * g(b(1)  ,b(2)+1,b(3)  ,i)
     &           +  (1.d0-n(1))*(1.d0-n(2))*(1.d0-n(3)) * g(b(1)  ,b(2)  ,b(3)  ,i)
#endif
         enddo


#if BL_SPACEDIM==3
      endif
#endif

      end

