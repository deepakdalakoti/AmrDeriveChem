#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "Damkoh_F.H"
#include "BC_TYPES.H"
#include "SPACE.H"
      
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

#if BL_SPACEDIM==2
      subroutine FORT_COMPDIFFTERMS(lo, hi, Y, DIMS(Y),
     &     RhoD, DIMS(RhoD),density,DIMS(density),vel,DIMS(vel), DT, DIMS(DT),
     &     dx, ns,si, nwhich) 
     
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
      integer si
      
      integer DIMDEC(Y)
      integer DIMDEC(density)
      integer DIMDEC(vel)
      integer DIMDEC(RhoD)
      
      integer DIMDEC(DT)
     
      REAL_T Y(DIMV(Y),ns)
      REAL_T density(DIMV(density))
      REAL_T vel(DIMV(vel),BL_SPACEDIM)
      REAL_T RhoD_xp(si), RhoD_xm(si)
      REAL_T RhoD_yp(si), RhoD_ym(si)     
      REAL_T RhoD(DIMV(RhoD),ns)
       
      REAL_T DT(DIMV(DT),(2*(BL_SPACEDIM)*si))
      REAL_T dx(SDIM)
      integer nwhich(si)
      integer i,j,n,p,s
      
      REAL_T dxInv2, dyInv2

      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
      print*, si, nwhich(1)

       
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
            do p= 1,si      
              n= nwhich(si)+1
               RhoD_xp(p) = 0.5*( RhoD(i,j,n) + RhoD(i+1,j,n) )
               RhoD_xm(p) = 0.5*( RhoD(i,j,n) + RhoD(i-1,j,n) )
               RhoD_yp(p) = 0.5*( RhoD(i,j,n) + RhoD(i,j+1,n) )
               RhoD_ym(p) = 0.5*( RhoD(i,j,n) + RhoD(i,j-1,n) )
               
            enddo
           do n=1,si
               s= nwhich(si)+1
        
          DT(i,j,(n-1)*BL_SPACEDIM+3) = 
     &              ( RhoD_xp(n)*(Y(i+1,j,s)-Y(i  ,j,s))
     &              - RhoD_xm(n)*(Y(i  ,j,s)-Y(i-1,j,s)) ) * dxInv2 
          DT(i,j,(n-1)*BL_SPACEDIM+4) =
     &              ( RhoD_yp(n)*(Y(i,j+1,s)-Y(i,j  ,s))
     &              - RhoD_ym(n)*(Y(i,j  ,s)-Y(i,j-1,s)) ) * dyInv2
               
            
             
         DT(i,j,2*(n-1)*BL_SPACEDIM+1) = (density(i+1,j)*vel(i+1,j,1)*Y(i+1,j,s)-
     &                    density(i-1,j)*vel(i-1,j,1)*Y(i-1,j,s))/(2.0*dx(1))
         DT(i,j,2*(n-1)*BL_SPACEDIM+2) = (density(i,j+1)*vel(i,j+1,2)*Y(i,j+1,s)-
     &                    density(i,j-1)*vel(i,j-1,2)*Y(i,j-1,s))/(2.0*dx(2))
           enddo
         end do
      end do
      end
#else
      subroutine FORT_COMPDIFFTERMS(lo, hi, Y, DIMS(Y),
     &      RhoD, DIMS(RhoD), density,DIMS(density),vel,DIMS(vel), DT, DIMS(DT),
     &     dx, ns,si,nwhich,Adv)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns,si
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      integer DIMDEC(DT)
      integer DIMDEC(density)
      integer DIMDEC(vel)
      REAL_T Y(DIMV(Y),ns)
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T density(DIMV(density))
      REAL_T DT(DIMV(DT),(3*ns))
      REAL_T Flux(DIMV(Y),ns,3)
      REAL_T vel(DIMV(vel),3)
      REAL_T dx(SDIM)
      REAL_T sumFlux, HAvg, YAvg
      REAL_T Adv(DIMV(DT),3)
      integer nwhich(si)
      integer i,j,k,n,p,s

      do k=lo(3), hi(3)
       do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
            
           do n=1,ns
          HAvg = 1.0/(RhoD(i-1,j,k,n))+1.0/(RhoD(i,j,k,n))
          HAvg = 2.0/HAvg
            Flux(i,j,k,n,1) = 
     &             HAvg*
     &                    (Y(i,j,k,n)-Y(i-1  ,j,k,n))/dx(1)
                    
          
             
           enddo
         end do
      end do
      enddo

      do k=lo(3), hi(3)
       do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)

         do n=1,ns
          HAvg = 1.0/(RhoD(i,j-1,k,n))+1.0/(RhoD(i,j,k,n))
          HAvg = 2.0/HAvg
          Flux(i,j,k,n,2) =
     &            HAvg*
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

           Flux(i,j,k,n,3) =
     &                  HAvg* (Y(i,j,k,n)-Y(i  ,j,k-1,n))/dx(3)

           enddo
         end do
      end do
      enddo

      do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux(i,j,k,n,1)
                  end do
                 
                  do n=1,ns
                    YAvg = 0.5*(Y(i-1,j,k,n)+Y(i,j,k,n))
                    Flux(i,j,k,n,1) = Flux(i,j,k,n,1)-  sumFlux*YAvg
                  end do
               end do
            end do
         end do



          do k=lo(3),hi(3)
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux(i,j,k,n,2)
c                 if(Y(i,j,k,n).gt.1e-5)  print*, Flux(i,j,k,n,2)
                  end do
            
                  do n=1,ns
                     YAvg = 0.5*(Y(i,j-1,k,n)+Y(i,j,k,n))
                     Flux(i,j,k,n,2) = Flux(i,j,k,n,2)- sumFlux*YAvg
                  end do
               end do
            end do
         end do

         do k = lo(3),hi(3)+1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  sumFlux = 0.d0
                  do n=1,Nspec
                     sumFlux = sumFlux + Flux(i,j,k,n,3)
                  end do
                  do n=1,Nspec
                     YAvg = 0.5*(Y(i,j,k-1,n)+Y(i,j,k,n))

                     Flux(i,j,k,n,3) =Flux(i,j,k,n,3)- sumFlux*YAvg
                  end do
               end do
            end do
         end do

       do k= lo(3), hi(3)
         do j= lo(2), hi(2)
          do i=lo(1), hi(1)
              
           do n=1,ns
           DT(i,j,k,n) = (Flux(i+1,j,k,n,1)-Flux(i,j,k,n,1))/(dx(1))
           DT(i,j,k,n+ns) = (Flux(i,j+1,k,n,2)-Flux(i,j,k,n,2))/(dx(2))
           DT(i,j,k,n+ns+ns) = (Flux(i,j,k+1,n,3)-Flux(i,j,k,n,3))/(dx(3))
c           DT(i,j,k,n) = DT(i,j,k,n)+DT(i,j,k,n+ns)+DT(i,j,k,n+ns+ns)
          enddo
       enddo
      enddo
      enddo
   
       do k=lo(3),hi(3)
        do j= lo(2), hi(2)
         do i=lo(1), hi(1)
          Adv(i,j,k,1) = -1.0*(density(i+1,j,k)*Y(i+1,j,k,nwhich(1))*vel(i+1,j,k,1)-
     &                   density(i-1,j,k)*Y(i-1,j,k,nwhich(1))*vel(i-1,j,k,1))/dx(1)
            
          Adv(i,j,k,2) = -1.0*(density(i,j+1,k)*Y(i,j+1,k,nwhich(1))*vel(i,j+1,k,1)-
     &                   density(i,j-1,k)*Y(i,j-1,k,nwhich(1))*vel(i,j-1,k,1))/dx(1)
 
          Adv(i,j,k,3) = -1.0*(density(i,j,k+1)*Y(i,j,k+1,nwhich(1))*vel(i,j,k+1,1)-
     &                   density(i,j,k-1)*Y(i,j,k-1,nwhich(1))*vel(i,j,k-1,1))/dx(1)

        enddo
      enddo
      enddo


            end
#endif
