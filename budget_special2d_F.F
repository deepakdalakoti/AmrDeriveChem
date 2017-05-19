#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "SPACE.H"
#include "budget_special2d_F.H"      
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
     &     dx, ns,si,nwhich1,Adv,total,Ncomp)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer ns,si,flo(SDIM),fhi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(RhoD)
      integer DIMDEC(DT)
      integer DIMDEC(density)
      integer DIMDEC(vel)
      REAL_T Y(DIMV(Y),ns)
      REAL_T RhoD(DIMV(RhoD),ns)
      REAL_T density(DIMV(density))
      REAL_T DT(DIMV(DT),(SDIM*ns))
      REAL_T Flux(DIMV(Y),ns,SDIM),FluxNZ(DIMV(Y),ns)
      REAL_T vel(DIMV(vel),SDIM),FluxT1(DIMV(Y),ns)
      REAL_T dx(SDIM),FluxT2(DIMV(Y),ns)
      REAL_T sumFlux, HAvg, YAvg, dot
      REAL_T Adv(DIMV(DT),si)
      REAL_T total(DIMV(DT),si)
      REAL_T NZ(DIMV(Y),SDIM+1),NYp(DIMV(Y),SDIM+1),T1(DIMV(Y),SDIM+1,ns),T2(DIMV(Y),SDIM+1,ns)
      REAL_T Ncomp(DIMV(Y))
      REAL_T tot,x1,x2,y1,y2,z1,z2,tdx,tdy,tdz,dotp
      REAL_T small
      integer nwhich1(si), nwhich(si)
      integer i,j,k,n,p,s
      small = 1e-8
      tdx= 0.5/dx(1)
      tdy = 0.5/dx(2)
      do i=1,si 
        nwhich(i)=nwhich1(i)+1
 
      enddo 
       do j=flo(2),fhi(2)
        do i=flo(1),fhi(1)+1
         
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
        do j= flo(2), fhi(2)
          do i= flo(1), fhi(1)
                  NZ(i,j,1) = -(Ncomp(i+1,j) - Ncomp(i-1,j))*tdx
                  NZ(i,j,2) = -(Ncomp(i,j+1) - Ncomp(i,j-1))*tdy
                  NZ(i,j,3) = sqrt( NZ(i,j,1)**2+NZ(i,j,2)**2)
                  NZ(i,j,1) =  NZ(i,j,1)/(NZ(i,j,3)+small)
                  NZ(i,j,2) =  NZ(i,j,2)/(NZ(i,j,3)+small)
                 
                  NYp(i,j,1) = -(Y(i+1,j,n) - Y(i-1,j,n))*tdx
                  NYp(i,j,2) = -(Y(i,j+1,n) - Y(i,j-1,n))*tdy
                  NYp(i,j,3) = sqrt(NYp(i,j,1)**2+NYp(i,j,2)**2)
                  NYp(i,j,1) =  NYp(i,j,1)/(NYp(i,j,3)+small)
                  NYp(i,j,2) =  NYp(i,j,2)/(NYp(i,j,3)+small)
                  dotp = NYp(i,j,1)*NZ(i,j,1)+ NYp(i,j,2)*NZ(i,j,2)
                  T2(i,j,1,n) = (NYp(i,j,1)-NZ(i,j,1)*dotp)/(sqrt(1-dotp**2)+small)
                  T2(i,j,2,n) = (NYp(i,j,2)-NZ(i,j,2)*dotp)/(sqrt(1-dotp**2)+small)
                  T2(i,j,3,n)= sqrt(T2(i,j,1,n)**2 + T2(i,j,2,n)**2)
                  T2(i,j,1,n) = T2(i,j,1,n)/(T2(i,j,3,n)+small)
                  T2(i,j,2,n) = T2(i,j,2,n)/(T2(i,j,3,n)+small)

                          
                 


           enddo
          enddo
       enddo

      do n=1,ns       
       do j= flo(2), fhi(2)
        do i = flo(1), fhi(1)
c  calculating flux at cell centers at respective directions

           FluxNZ(i,j,n) = 0.5*(Flux(i+1,j,n,1)+Flux(i,j,n,1))*NZ(i,j,1) +
     &                       0.5*(Flux(i,j+1,n,2)+Flux(i,j,n,2))*NZ(i,j,2) 
     &                       


           FluxT2(i,j,n) = 0.5*(Flux(i+1,j,n,1)+Flux(i,j,n,1))*T2(i,j,1,n) +
     &                       0.5*(Flux(i,j+1,n,2)+Flux(i,j,n,2))*T2(i,j,2,n) 
     &                      




       enddo
       enddo
      enddo
     
        do n=1,ns     
         do j= lo(2), hi(2)
          do i=lo(1), hi(1)
              

           DT(i,j,n) =0.5*NZ(i,j,1)*(FluxNZ(i+1,j,n)- FluxNZ(i-1,j,n))/dx(1) + 
     &              0.5*NZ(i,j,2)*(FluxNZ(i,j+1,n)- FluxNZ(i,j-1,n))/dx(2)  
     &                


          DT(i,j,n+ns) = 0.5*T2(i,j,1,n)*(FluxT2(i+1,j,n)- FluxT2(i-1,j,n))/dx(1) + 
     &                   0.5*T2(i,j,2,n)*(FluxT2(i,j+1,n)- FluxT2(i,j-1,n))/dx(2)  
     &                     
          
            
        enddo
       enddo
      enddo
  

      do n=1,si
        do j= lo(2), hi(2)
         do i=lo(1), hi(1)
          
                      

       x1 = density(i+1,j)*Y(i+1,j,nwhich(n))*(vel(i+1,j,1)*NZ(i+1,j,1)
     &     +    vel(i+1,j,2)*NZ(i+1,j,2))

      x2 =   density(i-1,j)*Y(i-1,j,nwhich(n))*(vel(i-1,j,1)*NZ(i-1,j,1)
     
     &   +    vel(i-1,j,2)*NZ(i-1,j,2))

      y1 = density(i,j+1)*Y(i,j+1,nwhich(n))*(vel(i,j+1,1)*NZ(i,j+1,1)
     &     + vel(i,j+1,2)*NZ(i,j+1,2))

      y2 =   density(i,j-1)*Y(i,j-1,nwhich(n))*(vel(i,j-1,1)*NZ(i,j-1,1)
     
     &   +    vel(i,j-1,2)*NZ(i,j-1,2))


c Negative sign in advection terms so that RHS = -grad(rhoUY) will be positive when advection bring species in

	Adv(i,j,(n-1)*2+1) = -0.5*NZ(i,j,1)*(x1-x2)/dx(1)-0.5*NZ(i,j,2)    
     &                        *(y1-y2)/dx(2) 
        Adv(i,j,(n-1)*2+1) = vel(i,j,1)

      x1 = density(i+1,j)*Y(i+1,j,nwhich(n))*(vel(i+1,j,1)*T2(i+1,j,1,nwhich(n))
     &     +    vel(i+1,j,2)*T2(i+1,j,2,nwhich(n)) )

      x2 =   density(i-1,j)*Y(i-1,j,nwhich(n))*(vel(i-1,j,1)*T2(i-1,j,1,nwhich(n))
     
     &   +    vel(i-1,j,2)*T2(i-1,j,2,nwhich(n)) )
 
       y1 = density(i,j+1)*Y(i,j+1,nwhich(n))*(vel(i,j+1,1)*T2(i,j+1,1,nwhich(n))
     &     + vel(i,j+1,2)*T2(i,j+1,2,nwhich(n)) )

      y2 =   density(i,j-1)*Y(i,j-1,nwhich(n))*(vel(i,j-1,1)*T2(i,j-1,1,nwhich(n))
  
     &   +    vel(i,j-1,2)*T2(i,j-1,2,nwhich(n)) )

cc Negative sign in advection terms so that RHS = -grad(rhoUY) will be positive when advection bring species in
      Adv(i,j,(n-1)*2+2) = -0.5*T2(i,j,1,nwhich(n))*(x1-x2)/dx(1)-0.5*T2(i,j,2,nwhich(n))    
     &                        *(y1-y2)/dx(2)

      Adv(i,j,(n-1)*2+2) = vel(i,j,2)

!    Hack to get velocity along T2 by storing it in total 
       total(i,j,n) = T2(i,j,2,nwhich(n)) 
        enddo
      enddo
      enddo

            end
