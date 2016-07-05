#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "projectVintoIso_F.H"
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
         call filcc(U(D_DECL(ARG_L1(U),ARG_L2(U),ARG_L3(U)),n),DIMS(U),dlo,dhi,dx,xlo,bc)
      end do
      end

      subroutine FORT_PRJECTV(lo, hi, T, DIMS(T), v, DIMS(v), pv, DIMS(pv), delta)
      implicit none
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(T)
      integer DIMDEC(v)
      integer DIMDEC(pv)
      REAL_T    T(DIMV(T))
      REAL_T    v(DIMV(v),SDIM)
      REAL_T   pv(DIMV(pv),SDIM)
      REAL_T delta(SDIM)

      integer i,j,k,n
      REAL_T tdxI(3), norm(SDIM), mag, vdotn
      
      do i=1,BL_SPACEDIM
         tdxI(i) =  one / (two*delta(i))
      enddo

#if BL_SPACEDIM==2
      print *,'add code for 2D here'
#else
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               norm(1) = tdxI(1)*(T(i+1,j,k) - T(i-1,j,k))
               norm(2) = tdxI(2)*(T(i,j+1,k) - T(i,j-1,k))
               norm(3) = tdxI(3)*(T(i,j,k+1) - T(i,j,k-1))

               mag = SQRT( norm(1)*norm(1) + norm(2)*norm(2) + norm(3)*norm(3) )

               if (mag.gt.1.e-8) then
                  norm(1) = norm(1)/mag
                  norm(2) = norm(2)/mag
                  norm(3) = norm(3)/mag
               endif

               vdotn = v(i,j,k,1)*norm(1) + v(i,j,k,2)*norm(2) + v(i,j,k,3)*norm(3)

               pv(i,j,k,1) = v(i,j,k,1) - vdotn*norm(1)
               pv(i,j,k,2) = v(i,j,k,2) - vdotn*norm(2)
               pv(i,j,k,3) = v(i,j,k,3) - vdotn*norm(3)

            end do
         end do
      end do
#endif
      end

