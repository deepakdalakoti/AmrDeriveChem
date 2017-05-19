#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "AmrDeriveIsoSec_F.H"
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

          subroutine FORT_COMPDIFFTERMS(lo, hi,flo, fhi, slo,shi,DIMS(fbox),DIMS(sbox),Y, DIMS(Y),
     &   RhoD, DIMS(RhoD), Reac,DIMS(Reac),
     &   density,DIMS(density), Speed_Out, dx,ws,ns,velocity,pos,Z,alpha)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Y), DIMDEC(fbox), DIMDEC(sbox)
      integer DIMDEC(RhoD)
      integer DIMDEC(density)
      integer DIMDEC(Reac)
      integer ws, ns,wsp
      integer flo(SDIM),fhi(SDIM),slo(SDIM),shi(SDIM)
      REAL_T Y(DIMV(Y),ns), Speed_Out(14), pos(3)
      REAL_T Z(DIMV(Y)), alpha(DIMV(Reac))
      REAL_T RhoD(DIMV(RhoD),ns),plo(SDIM)
      REAL_T velocity(DIMV(Reac),3)
      REAL_T density(DIMV(density))
      REAL_T DT(DIMV(sbox),4)
      REAL_T Reac(DIMV(Reac)), Speed(DIMV(sbox),14)   
      REAL_T dx(SDIM), NZ(DIMV(fbox),4), NYp(DIMV(fbox),4), T1(DIMV(fbox),4)
      REAL_T dx2, dy2 , dz2,sumFlux, T2(DIMV(fbox),4)
      REAL_T Flux_Y(DIMV(fbox),ns,3),Flux_Z(DIMV(fbox),4), dot(DIMV(fbox))
      REAL_T small, HAvg,alphaAvg, YAvg, DatOut(8)
      REAL_T FluxNZ(DIMV(fbox)), FluxT1(DIMV(fbox)), FluxNYp(DIMV(fbox))
      integer i,j,k,n,p,s,imid,jmid,kmid
      REAL_T const, tmp1, tmp2,denT 
      REAL_T dxInv2, dyInv2, dzInv2
      logical ok
      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      dzInv2 = 1.d0 / (dx(3)*dx(3))
      dx2= 0.5/dx(1)
      dy2 = 0.5/dx(2)
      dz2 = 0.5/dx(3)
      small = 1e-2
      wsp = ws+1
!     As thi box is not grown so I am making dimensions smaller here
      plo(1) = -0.00256
      plo(2) = -0.00256
      plo(3) = 0
!      print*, "inside" 

    
       do k= flo(3), fhi(3)
         do j=flo(2),fhi(2)
           do i=flo(1),fhi(1)+1

           do n=1,ns
            HAvg = 1.0/(RhoD(i-1,j,k,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,1) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i-1  ,j,k,n))/dx(1)
           
           alphaAvg = 1.0/alpha(i-1,j,k)+1.0/alpha(i,j,k)
           alphaAvg =2.0/alphaAvg
           Flux_Z(i,j,k,1) = alphaAvg*(Z(i,j,k)-Z(i-1,j,k))/dx(1)

           enddo
         end do
      end do
      enddo
    
       

       do k= flo(3), fhi(3)
        do j=flo(2),fhi(2)+1
          do i=flo(1),fhi(1)

           do n=1,ns
            HAvg = 1.0/(RhoD(i,j-1,k,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,2) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i  ,j-1,k,n))/dx(2)

           alphaAvg = 1.0/alpha(i,j-1,k)+1.0/alpha(i,j,k)
           alphaAvg =2.0/alphaAvg
           Flux_Z(i,j,k,2) = alphaAvg*(Z(i,j,k)-Z(i,j-1,k))/dx(2)

           enddo
         end do
      end do
      enddo
 
      do k= flo(3), fhi(3)+1
        do j=flo(2),fhi(2)
          do i=flo(1),fhi(1)

           do n=1,ns
            HAvg = 1.0/(RhoD(i,j,k-1,n))+1.0/(RhoD(i,j,k,n))
            HAvg = 2.0/HAvg
            Flux_Y(i,j,k,n,3) =
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i  ,j,k-1,n))/dx(3)
           alphaAvg = 1.0/alpha(i,j,k-1)+1.0/alpha(i,j,k)
           alphaAvg =2.0/alphaAvg
           Flux_Z(i,j,k,3) = alphaAvg*(Z(i,j,k)-Z(i,j,k-1))/dx(3)

           enddo
         end do
      end do
      enddo
 

!Correction for species fluxes
   
        do k = flo(3),fhi(3)
            do j = flo(2),fhi(2)
               do i = flo(1),fhi(1)+1
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

        do k = flo(3),fhi(3)
            do j = flo(2),fhi(2)+1
               do i = flo(1),fhi(1)
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


       do k = flo(3),fhi(3)+1
            do j = flo(2),fhi(2)
               do i = flo(1),fhi(1)
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
       call gradient(flo,fhi,Z,DIMS(Y),NZ,DIMS(fbox),dx,6)
       call gradient(flo,fhi,Y(:,:,:,wsp),DIMS(Y),NYp,DIMS(fbox),dx,6)
       
       do k=flo(3), fhi(3)
        do j=flo(2), fhi(2)
         do i=flo(1), fhi(1)
           
c              NZ(i,j,k,1) = -1.0*(Z(i+1,j,k) - Z(i-1,j,k))*dx2
c              NZ(i,j,k,2) = -1.0*(Z(i,j+1,k) - Z(i,j-1,k))*dy2
c              NZ(i,j,k,3) = -1.0*(Z(i,j,k+1) - Z(i,j,k-1))*dz2
          
              NZ(i,j,k,4) = sqrt( NZ(i,j,k,1)**2+NZ(i,j,k,2)**2+NZ(i,j,k,3)**2)
              NZ(i,j,k,1) = -1.0* NZ(i,j,k,1)/(NZ(i,j,k,4)+small)
              NZ(i,j,k,2) =  -1.0*NZ(i,j,k,2)/(NZ(i,j,k,4)+small)
              NZ(i,j,k,3) = -1.0* NZ(i,j,k,3)/(NZ(i,j,k,4)+small) 
!minus beacuse grdient just gives grsdirnt
                 
c              NYp(i,j,k,1) = -1.0*(Y(i+1,j,k,wsp) - Y(i-1,j,k,wsp))*dx2
c              NYp(i,j,k,2) = -1.0*(Y(i,j+1,k,wsp) - Y(i,j-1,k,wsp))*dy2
c              NYp(i,j,k,3) = -1.0*(Y(i,j,k+1,wsp) - Y(i,j,k-1,wsp))*dz2
              NYp(i,j,k,4) = sqrt(NYp(i,j,k,1)**2+NYp(i,j,k,2)**2+NYp(i,j,k,3)**2)
              NYp(i,j,k,1) =  -1.0*NYp(i,j,k,1)/(NYp(i,j,k,4)+small)
              NYp(i,j,k,2) =  -1.0*NYp(i,j,k,2)/(NYp(i,j,k,4)+small)
              NYp(i,j,k,3) =  -1.0*NYp(i,j,k,3)/(NYp(i,j,k,4)+small)
                   

             T1(i,j,k,1) = NYp(i,j,k,2)*NZ(i,j,k,3)-NZ(i,j,k,2)*NYp(i,j,k,3)
             T1(i,j,k,2)= NYp(i,j,k,3)*NZ(i,j,k,1)-NYp(i,j,k,1)*NZ(i,j,k,3)
             T1(i,j,k,3) =NYp(i,j,k,1)*NZ(i,j,k,2)-NYp(i,j,k,2)*NZ(i,j,k,1)
             T1(i,j,k,4) = sqrt(T1(i,j,k,1)**2+T1(i,j,k,2)**2+T1(i,j,k,3)**2)
             T1(i,j,k,1) = T1(i,j,k,1)/(T1(i,j,k,4)+small)
             T1(i,j,k,2) = T1(i,j,k,2)/(T1(i,j,k,4)+small)
             T1(i,j,k,3) = T1(i,j,k,3)/(T1(i,j,k,4)+small)

            T2(i,j,k,1) = NZ(i,j,k,2)*T1(i,j,k,3)-T1(i,j,k,2)*NZ(i,j,k,3)
            T2(i,j,k,2) = NZ(i,j,k,3)*T1(i,j,k,1)-T1(i,j,k,3)*NZ(i,j,k,1)
            T2(i,j,k,3) = NZ(i,j,k,1)*T1(i,j,k,2)-T1(i,j,k,1)*NZ(i,j,k,2)

            T2(i,j,k,4) = sqrt(T2(i,j,k,1)**2+T2(i,j,k,2)**2+T2(i,j,k,3)**2)
            T2(i,j,k,1)= T2(i,j,k,1)/(T2(i,j,k,4)+small) 
            T2(i,j,k,2) = T2(i,j,k,2)/(T2(i,j,k,4)+small)
            T2(i,j,k,3) = T2(i,j,k,3)/(T2(i,j,k,4)+small)
            dot(i,j,k) =NYp(i,j,k,1)*NZ(i,j,k,1)+NYp(i,j,k,2)*NZ(i,j,k,2)
     &        +NYp(i,j,k,3)*NZ(i,j,k,3)
             
         enddo
        enddo
       enddo
       n=wsp
       do k=flo(3), fhi(3)
       do j= flo(2), fhi(2)
        do i = flo(1), fhi(1)
c  calculating flux at cell centers at respective directions

           FluxNZ(i,j,k) = 0.5*(Flux_Z(i+1,j,k,1)+Flux_Z(i,j,k,1))*NZ(i,j,k,1) +
     &                       0.5*(Flux_Z(i,j+1,k,2)+Flux_Z(i,j,k,2))*NZ(i,j,k,2) +
     &                       0.5*(Flux_Z(i,j,k+1,3)+Flux_Z(i,j,k,3))*NZ(i,j,k,3)

           FluxNYp(i,j,k) = 0.5*(Flux_Y(i+1,j,k,n,1)+Flux_Y(i,j,k,n,1))*NYp(i,j,k,1) +
     &                       0.5*(Flux_Y(i,j+1,k,n,2)+Flux_Y(i,j,k,n,2))*NYp(i,j,k,2) +
     &                       0.5*(Flux_Y(i,j,k+1,n,3)+Flux_Y(i,j,k,n,3))*NYp(i,j,k,3)





       enddo
       enddo
      enddo

        do k=slo(3), shi(3)
         do j= slo(2), shi(2)
          do i=slo(1), shi(1)

           DT(i,j,k,1) = (Flux_Y(i+1,j,k,wsp,1)-Flux_Y(i,j,k,wsp,1))/(dx(1))+
     &               (Flux_Y(i,j+1,k,wsp,2)-Flux_Y(i,j,k,wsp,2))/dx(2) +
     &               (Flux_Y(i,j,k+1,wsp,3)-Flux_Y(i,j,k,wsp,3))/dx(3)
          
           DT(i,j,k,2) =   (Flux_Z(i+1,j,k,1)-Flux_Z(i,j,k,1))/(dx(1))+
     &               (Flux_Z(i,j+1,k,2)-Flux_Z(i,j,k,2))/dx(2) +
     &               (Flux_Z(i,j,k+1,3)-Flux_Z(i,j,k,3))/dx(3)
         

       enddo
      enddo
      enddo
       n=wsp
       do k= slo(3), shi(3)
         do j= slo(2), shi(2)
          do i=slo(1), shi(1)
              

          DT(i,j,k,3)=0.5*NYp(i,j,k,1)*(FluxNYp(i+1,j,k)- FluxNYp(i-1,j,k))/dx(1) + 
     &              0.5*NYp(i,j,k,2)*(FluxNYp(i,j+1,k)- FluxNYp(i,j-1,k))/dx(2) + 
     &              0.5*NYp(i,j,k,3)*(FluxNYp(i,j,k+1)- FluxNYp(i,j,k-1))/dx(3)  

          DT(i,j,k,4) = 0.5*NZ(i,j,k,1)*(FluxNZ(i+1,j,k)- FluxNZ(i-1,j,k))/dx(1) + 
     &               0.5*NZ(i,j,k,2)*(FluxNZ(i,j+1,k)- FluxNZ(i,j-1,k))/dx(2) + 
     &                0.5*NZ(i,j,k,3)*(FluxNZ(i,j,k+1)- FluxNZ(i,j,k-1))/dx(3)  

         enddo
        enddo
        enddo
       do k= slo(3), shi(3)
        do j= slo(2), shi(2)
         do i=slo(1), shi(1)
! harcoded min max of density along mixing line
!  
         denT = 23.2+((71.5-23.2)/0.45)*(Z(i,j,k))     
         Speed(i,j,k,1) = DT(i,j,k,1)/(density(i,j,k)*NYp(i,j,k,4)+small) 

         Speed(i,j,k,2) = (Reac(i,j,k))/(density(i,j,k)*NYp(i,j,k,4)+small)
         Speed(i,j,k,3) = Speed(i,j,k,1)+Speed(i,j,k,2)
         Speed(i,j,k,4) = DT(i,j,k,2)/(density(i,j,k)*NZ(i,j,k,4)+small)
         Speed(i,j,k,5) = (Speed(i,j,k,3)-dot(i,j,k)*Speed(i,j,k,4))/(sqrt(1.0-dot(i,j,k)**2)+small)
         Speed(i,j,k,6) = dot(i,j,k)
         Speed(i,j,k,7) = Speed(i,j,k,3)/(sqrt(1.0-dot(i,j,k)**2)+small)
         Speed(i,j,k,8) = -Speed(i,j,k,4)*dot(i,j,k)/(sqrt(1.0-dot(i,j,k)**2)+small)
         Speed(i,j,k,9) =  DT(i,j,k,3)/(density(i,j,k)*NYp(i,j,k,4)+small)
         Speed(i,j,k,10) = (DT(i,j,k,1)-DT(i,j,k,3))/(density(i,j,k)*NYp(i,j,k,4)+small)
         Speed(i,j,k,11) = DT(i,j,k,4)/(density(i,j,k)*NZ(i,j,k,4)+small)
         Speed(i,j,k,12) = (DT(i,j,k,2)-DT(i,j,k,4))/(density(i,j,k)*NZ(i,j,k,4)+small)
         Speed(i,j,k,13) = velocity(i,j,k,1)*T2(i,j,k,1)+velocity(i,j,k,2)*T2(i,j,k,2)+
     &                     velocity(i,j,k,3)*T2(i,j,k,3)
         Speed(i,j,k,14) = Reac(i,j,k)
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
     
      call ntrpv(pos,Speed,DIMS(sbox),dx,plo,Speed_Out,14,ok)
!      if(abs(Speed_Out(2)) >40 ) then
!         print*, Speed_Out(2), " HIGH"
!       endif
        if(.not.ok) then
         print*, 'not ok'
        endif
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

