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
     &     RhoD, DIMS(RhoD), DT, DIMS(DT),
     &     dx, ns,ind) 
     
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns
  
      
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      
      integer DIMDEC(DT)
     
      REAL_T Y(DIMV(Y),ns)
      REAL_T RhoD_xp, RhoD_xm
      REAL_T RhoD_yp, RhoD_ym     
      REAL_T RhoD(DIMV(RhoD),ns)
       
      REAL_T DT(DIMV(DT))
      REAL_T dx(SDIM)
      integer ind
      integer i,j,n,p,s
      REAL_T small = 1e-7;
      REAL_T dxInv2, dyInv2

      dxInv2 = 1.d0 / (dx(1)*dx(1))
      dyInv2 = 1.d0 / (dx(2)*dx(2))
   
      n= ind+1
       
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
               RhoD_xp = 0.5*( RhoD(i,j,n) + RhoD(i+1,j,n) )
               RhoD_xm = 0.5*( RhoD(i,j,n) + RhoD(i-1,j,n) )
               RhoD_yp = 0.5*( RhoD(i,j,n) + RhoD(i,j+1,n) )
               RhoD_ym = 0.5*( RhoD(i,j,n) + RhoD(i,j-1,n) )
               
            enddo
          
           DT(i,j) = 
     &              ( RhoD_xp*(Y(i+1,j,k,n)-Y(i  ,j,k,n))
     &              - RhoD_xm*(Y(i  ,j,k,n)-Y(i-1,j,k,n)) ) * dxInv2 
     &             +  ( RhoD_yp*(Y(i,j+1,k,n)-Y(i,j  ,k,n))
     &              - RhoD_ym*(Y(i,j  ,k,n)-Y(i,j-1,k,n)) ) * dyInv2
           DT(i,j) =-1.0/(small+DT(i,j))
    
         end do
      end do
      end
#else
      subroutine FORT_COMPDIFFTERMS(lo, hi, Y, DIMS(Y),
     &      RhoD, DIMS(RhoD), outs, DIMS(outs),
     &     dx, ns,ind1,Reac)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns,si
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      integer DIMDEC(outs)
      
      REAL_T Y(DIMV(Y),ns)
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T outs(DIMV(outs))
      REAL_T Reac(DIMV(outs))
      REAL_T dx(SDIM), Flux(DIMV(Y),ns,3), DT(DIMV(outs),3*ns)
      integer ind1,ind
      integer i,j,k,n,p,s
      REAL_T small,tmp
      REAL_T HAvg, YAvg, sumFlux
      ind = ind1+1
      small = 1e-6
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
          enddo
           outs(i,j,k)= Reac(i,j,k)/abs((DT(i,j,k,ind)+DT(i,j,k,ns+ind)
     &                     +DT(i,j,k,ind+2*ns)+small))
       enddo
      enddo
      enddo

 
            end
#endif
