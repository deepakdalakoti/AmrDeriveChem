#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "SPACE.H"
#include "budget2d_F.H"      
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

      subroutine FORT_COMPDIFFTERMS(flo,fhi,lo, hi, Y, DIMS(Y),
     &      RhoD, DIMS(RhoD), density,DIMS(density),vel,DIMS(vel), DT, DIMS(DT),
     &     dx, ns,si,nwhich1,Adv,total,cp, T,DiffT)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns,si,flo(SDIM),fhi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      integer DIMDEC(DT)
      integer DIMDEC(density)
      integer DIMDEC(vel)
      REAL_T Y(DIMV(Y),ns), T(DIMV(Y))
      REAL_T RhoD(DIMV(RhoD),ns+1)
      REAL_T density(DIMV(density))
      REAL_T DT(DIMV(DT),(SDIM*ns))
      REAL_T Flux(DIMV(Y),ns,SDIM)
      REAL_T vel(DIMV(vel),SDIM),FluxT(DIMV(Y),2)
      REAL_T dx(SDIM)
      REAL_T sumFlux, HAvg, YAvg, dot, lamAvg
      REAL_T Adv(DIMV(DT),si*SDIM)
      REAL_T total(DIMV(DT),si)
      REAL_T cp(DIMV(DT)), DiffT(DIMV(DT),4)
      REAL_T tot,x1,x2,y1,y2,z1,z2,tdx,tdy,tdz,dotp
      REAL_T small
      integer nwhich1(si), nwhich(si)
      integer i,j,k,n,p,s
      small = 1e-2
      tdx= 0.5/dx(1)
      tdy = 0.5/dx(2)
      do i=1,si 
        nwhich(i)=nwhich1(i)+1
 
      enddo 
       do j=flo(2),fhi(2)
        do i=flo(1),fhi(1)+1
           lamAvg = 1.0/(RhoD(i-1,j,ns+1)) + 1.0/(RhoD(i,j,ns+1))
           lamAvg = 2.0/lamAvg
           FluxT(i,j,1) = lamAvg*(T(i,j)-T(i-1,j))/dx(1)
           do n=1,ns
          HAvg = 1.0/(RhoD(i-1,j,n))+1.0/(RhoD(i,j,n))
          HAvg = 2.0/HAvg
            Flux(i,j,n,1) = 
     &             HAvg*
     &                    (Y(i,j,n)-Y(i-1  ,j,n))/dx(1)
                    
          
             
           enddo
         end do
      end do

       do j=flo(2),fhi(2)+1
        do i=flo(1),fhi(1)
           lamAvg = 1.0/(RhoD(i,j-1,ns+1)) + 1.0/(RhoD(i,j,ns+1))
           lamAvg = 2.0/lamAvg
           FluxT(i,j,2) = lamAvg*(T(i,j)-T(i,j-1))/dx(2)

         do n=1,ns
          HAvg = 1.0/(RhoD(i,j-1,n))+1.0/(RhoD(i,j,n))
          HAvg = 2.0/HAvg
          Flux(i,j,n,2) =
     &            HAvg*
     &                    (Y(i,j,n)-Y(i  ,j-1,n))/dx(2)



           enddo
         end do
      end do



            do j = flo(2),fhi(2)
               do i = flo(1),fhi(1)+1
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux(i,j,n,1)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i-1,j,n)+Y(i,j,n))
                    Flux(i,j,n,1) = Flux(i,j,n,1)-  sumFlux*YAvg
                  end do
               end do
         end do



            do j = flo(2),fhi(2)+1
               do i = flo(1),fhi(1)
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux(i,j,n,2)
c                 if(Y(i,j,k,n).gt.1e-5)  print*, Flux(i,j,k,n,2)
                  end do
            
                  do n=1,ns
                     YAvg = 0.5*(Y(i,j-1,n)+Y(i,j,n))
                     Flux(i,j,n,2) = Flux(i,j,n,2)- sumFlux*YAvg
                  end do
               end do
         end do


     
        do n=1,ns     
         do j= lo(2), hi(2)
          do i=lo(1), hi(1)
c  take care these fluxes are at cell faces

          DT(i,j,n) = (Flux(i+1,j,n,1)-Flux(i,j,n,1))/dx(1)
    

          DT(i,j,n+ns) = (Flux(i,j+1,n,2)-Flux(i,j,n,2))/dx(2)        
          
            
        enddo
       enddo
      enddo

      do j=lo(2), hi(2)
       do i = lo(1), hi(1)
      
        DiffT(i,j,1) = (FluxT(i+1,j,1)-FluxT(i,j,1))/dx(1)
        DiffT(i,j,2) = (FluxT(i,j+1,1)-FluxT(i,j,1))/dx(2)
c        DiffT(i,j,3) =0.5* density(i,j)*cp(i,j)*(vel(i+1,j,1)*T(i+1,j)-
c     &                 vel(i-1,j,1)*T(i-1,j))/dx(1)
c        DiffT(i,j,4) =0.5* density(i,j)*cp(i,j)*(vel(i,j+1,1)*T(i,j+1)-
c     &                 vel(i,j-1,1)*T(i,j-1))/dx(2)
       DiffT(i,j,3) =0.5* density(i,j)*cp(i,j)*vel(i,j,1)*(T(i+1,j)-
     &                 T(i-1,j))/dx(1)
        DiffT(i,j,4) =0.5* density(i,j)*cp(i,j)*vel(i,j+1,2)*(T(i,j+1)-
     &                 T(i,j-1))/dx(2)

 
       enddo
      enddo




  
       do n=1,si
         do j= lo(2), hi(2)
          do i = lo(1), hi(1)
    
            Adv(i,j,(n-1)*2+1) =0.5*(density(i+1,j)*vel(i+1,j,1)*Y(i+1,j,nwhich(n))
     &                          - density(i-1,j)*vel(i-1,j,1)*Y(i-1,j,nwhich(n)))/dx(1) 



            Adv(i,j,(n-1)*2+2) = 0.5*(density(i,j+1)*vel(i,j+1,2)*Y(i,j+1,nwhich(n))
     &                      - density(i,j-1)*vel(i,j-1,2)*Y(i,j-1,nwhich(n)))/dx(2) 

            Adv(i,j,(n-1)*2+1) = -Adv(i,j,(n-1)*2+1)
            Adv(i,j,(n-1)*2+2) = -Adv(i,j,(n-1)*2+2)

         
             total(i,j,n) = DT(i,j,nwhich(n))+DT(i,j,n+nwhich(n)) +
     &                     Adv(i,j,(n-1)*2+1)+Adv(i,j,(n-1)*2+2)

           
          enddo
         enddo
       enddo







            end
