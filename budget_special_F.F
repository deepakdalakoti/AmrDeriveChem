#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "SPACE.H"
#include "budget_special_F.H"      
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
      subroutine FORT_COMPDIFFTERMS(flo,fhi,lo, hi, Y, DIMS(Y),
     &      RhoD, DIMS(RhoD), density,DIMS(density),vel,DIMS(vel), DT, DIMS(DT),
     &     dx, ns,si,nwhich1,Adv,total,Ncomp, var, RR_OH, tot_in)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns,si,flo(SDIM),fhi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      integer DIMDEC(DT)
      integer DIMDEC(density)
      integer DIMDEC(vel), tot_in
      REAL_T Y(DIMV(Y),ns)
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T density(DIMV(density))
      REAL_T DT(DIMV(DT),(3*ns))
      REAL_T Flux(DIMV(Y),ns,3),FluxNZ(DIMV(Y),ns)
      REAL_T vel(DIMV(vel),3),FluxT1(DIMV(Y),ns)
      REAL_T dx(SDIM),FluxT2(DIMV(Y),ns)
      REAL_T sumFlux, HAvg, YAvg, dot
      REAL_T Adv(DIMV(DT),si)
      REAL_T total(DIMV(DT),si)
      REAL_T NZ(DIMV(Y),4),NYp(DIMV(Y),4),T1(DIMV(Y),4,ns),T2(DIMV(Y),4,ns)
      REAL_T var(DIMV(DT),tot_in)
      REAL_T Ncomp(DIMV(Y))
      REAL_T tot,x1,x2,y1,y2,z1,z2,tdx,tdy,tdz
      REAL_T small
      integer nwhich1(si), nwhich(si), RR_OH
      integer i,j,k,n,p,s
      small = 1e-2
      tdx= 0.5/dx(1)
      tdy = 0.5/dx(2)
      tdz = 0.5/dx(3)
      do i=1,si 
        nwhich(i)=nwhich1(i)+1
 
      enddo 
       print*, tot_in, RR_OH
      do k=flo(3), fhi(3)
       do j=flo(2),fhi(2)
        do i=flo(1),fhi(1)+1
         
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

      do k=flo(3), fhi(3)
       do j=flo(2),fhi(2)+1
        do i=flo(1),fhi(1)

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


      do k= flo(3), fhi(3)+1
       do j=flo(2),fhi(2)
        do i=flo(1),fhi(1)
   
           do n=1,ns
          HAvg = 1.0/(RhoD(i,j,k-1,n))+1.0/(RhoD(i,j,k,n))
          HAvg = 2.0/HAvg

           Flux(i,j,k,n,3) =
     &                  HAvg* (Y(i,j,k,n)-Y(i  ,j,k-1,n))/dx(3)

           enddo
         end do
      end do
      enddo

      do k = flo(3),fhi(3)
            do j = flo(2),fhi(2)
               do i = flo(1),fhi(1)+1
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



          do k=flo(3),fhi(3)
            do j = flo(2),fhi(2)+1
               do i = flo(1),fhi(1)
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

         do k = flo(3),fhi(3)+1
            do j = flo(2),fhi(2)
               do i = flo(1),fhi(1)
                  sumFlux = 0.d0
                  do n=1,ns
                     sumFlux = sumFlux + Flux(i,j,k,n,3)
                  end do
                  do n=1,ns
                     YAvg = 0.5*(Y(i,j,k-1,n)+Y(i,j,k,n))

                     Flux(i,j,k,n,3) =Flux(i,j,k,n,3)- sumFlux*YAvg
                  end do
               end do
            end do
         end do

      do n=1,ns
       do k=flo(3), fhi(3)
        do j= flo(2), fhi(2)
          do i= flo(1), fhi(1)
                  NZ(i,j,k,1) = -(Ncomp(i+1,j,k) - Ncomp(i-1,j,k))*tdx
                  NZ(i,j,k,2) = -(Ncomp(i,j+1,k) - Ncomp(i,j-1,k))*tdy
                  NZ(i,j,k,3) = -(Ncomp(i,j,k+1) - Ncomp(i,j,k-1))*tdz
                  NZ(i,j,k,4) = sqrt( NZ(i,j,k,1)**2+NZ(i,j,k,2)**2+NZ(i,j,k,3)**2)
                  NZ(i,j,k,1) =  NZ(i,j,k,1)/(NZ(i,j,k,4)+small)
                  NZ(i,j,k,2) =  NZ(i,j,k,2)/(NZ(i,j,k,4)+small)
                  NZ(i,j,k,3) =  NZ(i,j,k,3)/(NZ(i,j,k,4)+small)
                 
                  NYp(i,j,k,1) = -(Y(i+1,j,k,n) - Y(i-1,j,k,n))*tdx
                  NYp(i,j,k,2) = -(Y(i,j+1,k,n) - Y(i,j-1,k,n))*tdy
                  NYp(i,j,k,3) = -(Y(i,j,k+1,n) - Y(i,j,k-1,n))*tdz
                  NYp(i,j,k,4) = sqrt(NYp(i,j,k,1)**2+NYp(i,j,k,2)**2+NYp(i,j,k,3)**2)
                  NYp(i,j,k,1) =  NYp(i,j,k,1)/(NYp(i,j,k,4)+small)
                  NYp(i,j,k,2) =  NYp(i,j,k,2)/(NYp(i,j,k,4)+small)
                  NYp(i,j,k,3) =  NYp(i,j,k,3)/(NYp(i,j,k,4)+small)
                   

                 T1(i,j,k,1,n) = NYp(i,j,k,2)*NZ(i,j,k,3)-NZ(i,j,k,2)*NYp(i,j,k,3)
                 T1(i,j,k,2,n)= NYp(i,j,k,3)*NZ(i,j,k,1)-NYp(i,j,k,1)*NZ(i,j,k,3)
                 T1(i,j,k,3,n) =NYp(i,j,k,1)*NZ(i,j,k,2)-NYp(i,j,k,2)*NZ(i,j,k,1)
                 T1(i,j,k,4,n) = sqrt(T1(i,j,k,1,n)**2+T1(i,j,k,2,n)**2+T1(i,j,k,3,n)**2)
                 T1(i,j,k,1,n) = T1(i,j,k,1,n)/(T1(i,j,k,4,n)+small)
                 T1(i,j,k,2,n) = T1(i,j,k,2,n)/(T1(i,j,k,4,n)+small)
                 T1(i,j,k,3,n) = T1(i,j,k,3,n)/(T1(i,j,k,4,n)+small)
c               T2 is NZ cross T1 so
                 
            T2(i,j,k,1,n) = NZ(i,j,k,2)*T1(i,j,k,3,n)-T1(i,j,k,2,n)*NZ(i,j,k,3)
            T2(i,j,k,2,n) = NZ(i,j,k,3)*T1(i,j,k,1,n)-T1(i,j,k,3,n)*NZ(i,j,k,1)
            T2(i,j,k,3,n) = NZ(i,j,k,1)*T1(i,j,k,2,n)-T1(i,j,k,1,n)*NZ(i,j,k,2)

            T2(i,j,k,4,n) = sqrt(T2(i,j,k,1,n)**2+T2(i,j,k,2,n)**2+T2(i,j,k,3,n)**2)
            T2(i,j,k,1,n)= T2(i,j,k,1,n)/(T2(i,j,k,4,n)+small) 
            T2(i,j,k,2,n) = T2(i,j,k,2,n)/(T2(i,j,k,4,n)+small)
            T2(i,j,k,3,n) = T2(i,j,k,3,n)/(T2(i,j,k,4,n)+small)
            


           enddo
          enddo
         enddo
       enddo
      do n=1,ns       
      do k=flo(3), fhi(3)
       do j= flo(2), fhi(2)
        do i = flo(1), fhi(1)
c  calculating flux at cell centers at respective directions

           FluxNZ(i,j,k,n) = 0.5*(Flux(i+1,j,k,n,1)+Flux(i,j,k,n,1))*NZ(i,j,k,1) +
     &                       0.5*(Flux(i,j+1,k,n,2)+Flux(i,j,k,n,2))*NZ(i,j,k,2) +
     &                       0.5*(Flux(i,j,k+1,n,3)+Flux(i,j,k,n,3))*NZ(i,j,k,3)

           FluxT1(i,j,k,n) = 0.5*(Flux(i+1,j,k,n,1)+Flux(i,j,k,n,1))*T1(i,j,k,1,n) +
     &                       0.5*(Flux(i,j+1,k,n,2)+Flux(i,j,k,n,2))*T1(i,j,k,2,n) +
     &                       0.5*(Flux(i,j,k+1,n,3)+Flux(i,j,k,n,3))*T1(i,j,k,3,n)

           FluxT2(i,j,k,n) = 0.5*(Flux(i+1,j,k,n,1)+Flux(i,j,k,n,1))*T2(i,j,k,1,n) +
     &                       0.5*(Flux(i,j+1,k,n,2)+Flux(i,j,k,n,2))*T2(i,j,k,2,n) +
     &                       0.5*(Flux(i,j,k+1,n,3)+Flux(i,j,k,n,3))*T2(i,j,k,3,n)




         enddo
       enddo
       enddo
      enddo
     
        do n=1,ns     
       do k= lo(3), hi(3)
         do j= lo(2), hi(2)
          do i=lo(1), hi(1)
              

           DT(i,j,k,n) =0.5*NZ(i,j,k,1)*(FluxNZ(i+1,j,k,n)- FluxNZ(i-1,j,k,n))/dx(1) + 
     &              0.5*NZ(i,j,k,2)*(FluxNZ(i,j+1,k,n)- FluxNZ(i,j-1,k,n))/dx(2) + 
     &              0.5*NZ(i,j,k,3)*(FluxNZ(i,j,k+1,n)- FluxNZ(i,j,k-1,n))/dx(3)  

          DT(i,j,k,n+ns) = 0.5*T1(i,j,k,1,n)*(FluxT1(i+1,j,k,n)- FluxT1(i-1,j,k,n))/dx(1) + 
     &               0.5*T1(i,j,k,2,n)*(FluxT1(i,j+1,k,n)- FluxT1(i,j-1,k,n))/dx(2) + 
     &                0.5*T1(i,j,k,3,n)*(FluxT1(i,j,k+1,n)- FluxT1(i,j,k-1,n))/dx(3)  

          DT(i,j,k,n+ns+ns) = 0.5*T2(i,j,k,1,n)*(FluxT2(i+1,j,k,n)- FluxT2(i-1,j,k,n))/dx(1) + 
     &                   0.5*T2(i,j,k,2,n)*(FluxT2(i,j+1,k,n)- FluxT2(i,j-1,k,n))/dx(2) + 
     &                   0.5*T2(i,j,k,3,n)*(FluxT2(i,j,k+1,n)- FluxT2(i,j,k-1,n))/dx(3)  
c
c           DT(i,j,k,n) =NZ(i,j,k,1)*(FluxNZ(i+1,j,k,n)- FluxNZ(i,j,k,n))/dx(1) + 
c     &              NZ(i,j,k,2)*(FluxNZ(i,j+1,k,n)- FluxNZ(i,j,k,n))/dx(2) + 
c     &              NZ(i,j,k,3)*(FluxNZ(i,j,k+1,n)- FluxNZ(i,j,k,n))/dx(3)  
c
c          DT(i,j,k,n+ns) = T1(i,j,k,1,n)*(FluxT1(i+1,j,k,n)- FluxT1(i,j,k,n))/dx(1) + 
c     &               T1(i,j,k,2,n)*(FluxT1(i,j+1,k,n)- FluxT1(i,j,k,n))/dx(2) + 
c     &                T1(i,j,k,3,n)*(FluxT1(i,j,k+1,n)- FluxT1(i,j,k,n))/dx(3)  

c          DT(i,j,k,n+ns+ns) = T2(i,j,k,1,n)*(FluxT2(i+1,j,k,n)- FluxT2(i,j,k,n))/dx(1) + 
c     &                   T2(i,j,k,2,n)*(FluxT2(i,j+1,k,n)- FluxT2(i,j,k,n))/dx(2) + 
c     &                   T2(i,j,k,3,n)*(FluxT2(i,j,k+1,n)- FluxT2(i,j,k,n))/dx(3)  
            
        enddo
       enddo
      enddo
      enddo
  

      do n=1,si
       do k=lo(3),hi(3)
        do j= lo(2), hi(2)
         do i=lo(1), hi(1)
          
                      

       x1 = density(i+1,j,k)*Y(i+1,j,k,nwhich(n))*(vel(i+1,j,k,1)*NZ(i+1,j,k,1)
     &     +    vel(i+1,j,k,2)*NZ(i+1,j,k,2) + vel(i+1,j,k,3)*NZ(i+1,j,k,3))

      x2 =   density(i-1,j,k)*Y(i-1,j,k,nwhich(n))*(vel(i-1,j,k,1)*NZ(i-1,j,k,1)
     
     &   +    vel(i-1,j,k,2)*NZ(i-1,j,k,2) +vel(i-1,j,k,3)*NZ(i-1,j,k,3))

      y1 = density(i,j+1,k)*Y(i,j+1,k,nwhich(n))*(vel(i,j+1,k,1)*NZ(i,j+1,k,1)
     &     + vel(i,j+1,k,2)*NZ(i,j+1,k,2) + vel(i,j+1,k,3)*NZ(i,j+1,k,3))

      y2 =   density(i,j-1,k)*Y(i,j-1,k,nwhich(n))*(vel(i,j-1,k,1)*NZ(i,j-1,k,1)
     
     &   +    vel(i,j-1,k,2)*NZ(i,j-1,k,2) +vel(i,j-1,k,3)*NZ(i,j-1,k,3))

      z1 = density(i,j,k+1)*Y(i,j,k+1,nwhich(n))*(vel(i,j,k+1,1)*NZ(i,j,k+1,1)
     &     +  vel(i,j,k+1,2)*NZ(i,j,k+1,2) + vel(i,j,k+1,3)*NZ(i,j,k+1,3))

      z2 =   density(i,j,k-1)*Y(i,j,k-1,nwhich(n))*(vel(i,j,k-1,1)*NZ(i,j,k-1,1)
     
     &   +    vel(i,j,k-1,2)*NZ(i,j,k-1,2) +vel(i,j,k-1,3)*NZ(i,j,k-1,3))

c Negative sign in advection terms so that RHS = -grad(rhoUY) will be positive when advection bring species in
      Adv(i,j,k,(n-1)*3+1) = -0.5*NZ(i,j,k,1)*(x1-x2)/dx(1)-0.5*NZ(i,j,k,2)    
     &                        *(y1-y2)/dx(2) - 0.5*NZ(i,j,k,3)*(z1-z2)/dx(3)

      x1 = density(i+1,j,k)*Y(i+1,j,k,nwhich(n))*(vel(i+1,j,k,1)*T1(i+1,j,k,1,nwhich(n))
     &     +    vel(i+1,j,k,2)*T1(i+1,j,k,2,nwhich(n)) + vel(i+1,j,k,3)*T1(i+1,j,k,3,nwhich(n)))

      x2 =   density(i-1,j,k)*Y(i-1,j,k,nwhich(n))*(vel(i-1,j,k,1)*T1(i-1,j,k,1,nwhich(n))
     
     &   +    vel(i-1,j,k,2)*T1(i-1,j,k,2,nwhich(n)) +vel(i-1,j,k,3)*T1(i-1,j,k,3,nwhich(n)))

      y1 = density(i,j+1,k)*Y(i,j+1,k,nwhich(n))*(vel(i,j+1,k,1)*T1(i,j+1,k,1,nwhich(n))
     &     + vel(i,j+1,k,2)*T1(i,j+1,k,2,nwhich(n)) + vel(i,j+1,k,3)*T1(i,j+1,k,3,nwhich(n)))

      y2 =   density(i,j-1,k)*Y(i,j-1,k,nwhich(n))*(vel(i,j-1,k,1)*T1(i,j-1,k,1,nwhich(n))
     
     &   +    vel(i,j-1,k,2)*T1(i,j-1,k,2,nwhich(n)) +vel(i,j-1,k,3)*T1(i,j-1,k,3,nwhich(n)))

      z1 = density(i,j,k+1)*Y(i,j,k+1,nwhich(n))*(vel(i,j,k+1,1)*T1(i,j,k+1,1,nwhich(n))
     &     +  vel(i,j,k+1,2)*T1(i,j,k+1,2,nwhich(n)) + vel(i,j,k+1,3)*T1(i,j,k+1,3,nwhich(n)))

      z2 =   density(i,j,k-1)*Y(i,j,k-1,nwhich(n))*(vel(i,j,k-1,1)*T1(i,j,k-1,1,nwhich(n))
    
     &   +    vel(i,j,k-1,2)*T1(i,j,k-1,2,nwhich(n)) +vel(i,j,k-1,3)*T1(i,j,k-1,3,nwhich(n)))

c Negative sign in advection terms so that RHS = -grad(rhoUY) will be positive when advection bring species in
      Adv(i,j,k,(n-1)*3+2) = -0.5*T1(i,j,k,1,nwhich(n))*(x1-x2)/dx(1)-0.5*T1(i,j,k,2,nwhich(n))    
     &                        *(y1-y2)/dx(2) - 0.5*T1(i,j,k,3,nwhich(n))*(z1-z2)/dx(3)



      x1 = density(i+1,j,k)*Y(i+1,j,k,nwhich(n))*(vel(i+1,j,k,1)*T2(i+1,j,k,1,nwhich(n))
     &     +    vel(i+1,j,k,2)*T2(i+1,j,k,2,nwhich(n)) + vel(i+1,j,k,3)*T2(i+1,j,k,3,nwhich(n)))

      x2 =   density(i-1,j,k)*Y(i-1,j,k,nwhich(n))*(vel(i-1,j,k,1)*T2(i-1,j,k,1,nwhich(n))
     
     &   +    vel(i-1,j,k,2)*T2(i-1,j,k,2,nwhich(n)) +vel(i-1,j,k,3)*T2(i-1,j,k,3,nwhich(n)))
 
       y1 = density(i,j+1,k)*Y(i,j+1,k,nwhich(n))*(vel(i,j+1,k,1)*T2(i,j+1,k,1,nwhich(n))
     &     + vel(i,j+1,k,2)*T2(i,j+1,k,2,nwhich(n)) + vel(i,j+1,k,3)*T2(i,j+1,k,3,nwhich(n)))

      y2 =   density(i,j-1,k)*Y(i,j-1,k,nwhich(n))*(vel(i,j-1,k,1)*T2(i,j-1,k,1,nwhich(n))
  
     &   +    vel(i,j-1,k,2)*T2(i,j-1,k,2,nwhich(n)) +vel(i,j-1,k,3)*T2(i,j-1,k,3,nwhich(n)))

      z1 = density(i,j,k+1)*Y(i,j,k+1,nwhich(n))*(vel(i,j,k+1,1)*T2(i,j,k+1,1,nwhich(n))
     &     +  vel(i,j,k+1,2)*T2(i,j,k+1,2,nwhich(n)) + vel(i,j,k+1,3)*T2(i,j,k+1,3,nwhich(n)))

      z2 =   density(i,j,k-1)*Y(i,j,k-1,nwhich(n))*(vel(i,j,k-1,1)*T2(i,j,k-1,1,nwhich(n))
     
     &   +    vel(i,j,k-1,2)*T2(i,j,k-1,2,nwhich(n)) +vel(i,j,k-1,3)*T2(i,j,k-1,3,nwhich(n)))

cc Negative sign in advection terms so that RHS = -grad(rhoUY) will be positive when advection bring species in
      Adv(i,j,k,(n-1)*3+3) = -0.5*T2(i,j,k,1,nwhich(n))*(x1-x2)/dx(1)-0.5*T2(i,j,k,2,nwhich(n))    
     &                        *(y1-y2)/dx(2) - 0.5*T2(i,j,k,3,nwhich(n))*(z1-z2)/dx(3)

c      total(i,j,k,n) = Adv(i,j,k,(n-1)*3+1)+Adv(i,j,k,(n-1)*3+2)+Adv(i,j,k,(n-1)*3+3)
c     &                 +DT(i,j,k,nwhich(n)) + DT(i,j,k,nwhich(n)+ns)+DT(i,j,k,2*ns+nwhich(n))
       total(i,j,k,n) = abs(var(i,j,k,RR_OH))/(abs(DT(i,j,k,nwhich(n))+DT(i,j,k,nwhich(n)+ns)
     &                +                DT(i,j,k,nwhich(n)+2*ns))+small*5)
        enddo
      enddo
      enddo
      enddo

            end
#endif
