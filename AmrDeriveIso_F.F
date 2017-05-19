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

          subroutine FORT_COMPDIFFTERMS(lo, hi, Y, DIMS(Y),
     &   RhoD, DIMS(RhoD), Reac,DIMS(Reac),
     &   density,DIMS(density), Speed_Out, dx,ws,ns,velocity,pos,Z)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      integer DIMDEC(density)
      integer DIMDEC(Reac)
      integer ws, ns,wsp
      REAL_T Y(DIMV(Y),ns), Speed_Out(4), pos(3)
      REAL_T RhoD(DIMV(RhoD),ns),plo(SDIM)
      REAL_T Z(DIMV(Y))
      REAL_T velocity(DIMV(Reac),3)
      REAL_T density(DIMV(density))
      REAL_T DT(DIMV(Reac)), normal(DIMV(Reac),3)
      REAL_T Reac(DIMV(Reac)), Speed(DIMV(Reac),4)   
      REAL_T dx(SDIM)
      REAL_T dx2, dy2 , dz2,sumFlux
      REAL_T Flux_Y(DIMV(Y),ns,3), grad_Y(DIMV(Y),4)
      REAL_T small, HAvg, YAvg, DatOut(8)
      integer i,j,k,n,p,s,imid,jmid,kmid
      REAL_T const, tmp1, tmp2, denT 
      REAL_T dxInv2, dyInv2, dzInv2
      logical ok
      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      dzInv2 = 1.d0 / (dx(3)*dx(3))
      dx2= 0.5/dx(1)
      dy2 = 0.5/dx(2)
      dz2 = 0.5/dx(3)
      small = 1e-3
      wsp = ws+1
!     As thi box is not grown so I am making dimensions smaller here
      lo(3) = lo(3)+1
      hi(3) = hi(3)-1
      lo(2) = lo(2)+1
      hi(2) = hi(2)-1
      lo(1) = lo(1)+1
      hi(1) = hi(1)-1
      plo(1) = -0.00256
      plo(2) = -0.00256
      plo(3) = 0
       do k= lo(3), hi(3)
         do j=lo(2),hi(2)
           do i=lo(1),hi(1)+1

           do n=1,ns
            HAvg = 1.0/(RhoD(i-1,j,k,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,1) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i-1  ,j,k,n))/dx(1)
           
           enddo
         end do
      end do
      enddo
    
       

       do k= lo(3), hi(3)
        do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

           do n=1,ns
            HAvg = 1.0/(RhoD(i,j-1,k,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,2) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i  ,j-1,k,n))/dx(2)
           
           enddo
         end do
      end do
      enddo
 
      do k= lo(3), hi(3)+1
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)

           do n=1,ns
            HAvg = 1.0/(RhoD(i,j,k-1,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,3) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i  ,j,k-1,n))/dx(3)
           
           enddo
         end do
      end do
      enddo
 

!Correction for species fluxes
   
        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux_Y(i,j,k,n,1)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i-1,j,k,n)+Y(i,j,k,n))
                    Flux_Y(i,j,k,n,1) = Flux_Y(i,j,k,n,1)-  sumFlux*YAvg
                  end do
               end do
            end do
         end do

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux_Y(i,j,k,n,2)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i,j-1,k,n)+Y(i,j,k,n))
                    Flux_Y(i,j,k,n,2) = Flux_Y(i,j,k,n,2)-  sumFlux*YAvg
                  end do
               end do
            end do
         end do


       do k = lo(3),hi(3)+1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux_Y(i,j,k,n,3)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i,j,k-1,n)+Y(i,j,k,n))
                    Flux_Y(i,j,k,n,3) = Flux_Y(i,j,k,n,3)-  sumFlux*YAvg
                  end do
               end do
            end do
         end do
  
       do k=lo(3), hi(3)
        do j=lo(2), hi(2)
         do i=lo(1), hi(1)
           
           grad_Y(i,j,k,1) = (Y(i+1,j,k,wsp)-Y(i-1,j,k,wsp))*dx2
           grad_Y(i,j,k,2) = (Y(i,j+1,k,wsp)-Y(i,j-1,k,wsp))*dy2
           grad_Y(i,j,k,3) = (Y(i,j,k+1,wsp)-Y(i,j,k-1,wsp))*dz2
           grad_Y(i,j,k,4) = sqrt(grad_Y(i,j,k,1)**2+grad_Y(i,j,k,2)**2+grad_Y(i,j,k,3)**2) 
           normal(i,j,k,1) = -grad_Y(i,j,k,1)/(grad_Y(i,j,k,4)+small)
	   normal(i,j,k,2) = -grad_Y(i,j,k,2)/(grad_Y(i,j,k,4)+small)
           normal(i,j,k,3) = -grad_Y(i,j,k,3)/(grad_Y(i,j,k,4)+small)


         enddo
        enddo
       enddo

        do k=lo(3), hi(3)
         do j= lo(2), hi(2)
          do i=lo(1), hi(1)

           DT(i,j,k) = (Flux_Y(i+1,j,k,wsp,1)-Flux_Y(i,j,k,wsp,1))/(dx(1))+
     &               (Flux_Y(i,j+1,k,wsp,2)-Flux_Y(i,j,k,wsp,2))/dx(2) +
     &               (Flux_Y(i,j,k+1,wsp,3)-Flux_Y(i,j,k,wsp,3))/dx(3)
          
           

       enddo
      enddo
      enddo

       do i= lo(1), hi(1)
        do j= lo(2), hi(2)
         do k=lo(3), hi(3)
         denT = 23.2+((71.5-23.2)/0.45)*Z(i,j,k) 
         Speed(i,j,k,1) = DT(i,j,k)/(denT*grad_Y(i,j,k,4)+small)

         Speed(i,j,k,2) = (Reac(i,j,k))/(denT*grad_Y(i,j,k,4)+small)
         Speed(i,j,k,3) = Speed(i,j,k,1)+Speed(i,j,k,2)
         Speed(i,j,k,4) = velocity(i,j,k,1)*normal(i,j,k,1) + 
     &             velocity(i,j,k,2)*normal(i,j,k,2)+
     &             velocity(i,j,k,3)*normal(i,j,k,3)

         enddo
       enddo
      enddo

      imid = (lo(1)+hi(1))*0.5
      jmid = (lo(2)+hi(2))*0.5
      kmid = (lo(3)+hi(3))*0.5

      
!      print*, Speed(imid,jmid,kmid,1), Speed(imid,jmid,kmid,2), velocity(imid,jmid,kmid,1),Reac(imid,jmid,kmid), density(imid,jmid,kmid)
!      Speed_Out(1) = Speed(imid,jmid,kmid,1)
!      Speed_Out(2) = Speed(imid,jmid,kmid,2)
!      Speed_Out(3) = Speed(imid,jmid,kmid,3)
     
      call ntrpv(pos,Speed,DIMS(Reac),dx,plo,Speed_Out,4 ,ok)
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
         b(i) = FLOOR( tmp )
         n(i) = ( x(i) - ( (b(i)+0.5d0)*dx(i) + plo(i) ) )/dx(i)
         n(i) = MAX(0.d0,MIN(1.d0,n(i)))
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
