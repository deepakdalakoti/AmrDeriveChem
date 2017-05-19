#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "SPACE.H"

	subroutine gradient(lo,hi,var,DIMS(var),res,DIMS(res),dx,order)
	implicit none
	integer  lo(3), hi(3), DIMDEC(var),DIMDEC(res)
	REAL_T var(DIMV(var))
	REAL_T res(DIMV(res),3),dx(3)
	integer order,i,j,k

	Select case(order)
	  case(2)
	   do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
		 res(i,j,k,1) = (var(i+1,j,k)-var(i-1,j,k))*0.5/dx(1)
                 res(i,j,k,2) = (var(i,j+1,k)-var(i,j-1,k))*0.5/dx(2)

		 res(i,j,k,3) = (var(i,j,k+1)-var(i,j,k-1))*0.5/dx(3)
	      enddo
	    enddo
	   enddo
 
	 
	 case(4)
	   do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
		 res(i,j,k,1) = (-1.0*var(i+2,j,k)+8.0*var(i+1,j,k)-8.0*var(i-1,j,k)+var(i-2,j,k))/(12.0*dx(1))

		 res(i,j,k,2) = (-1.0*var(i,j+2,k)+8.0*var(i,j+1,k)-8.0*var(i,j-1,k)+var(i,j-2,k))/(12.0*dx(2))

	         res(i,j,k,3) = (-1.0*var(i,j,k+2)+8.0*var(i,j,k+1)-8.0*var(i,j,k-1)+var(i,j,k-2))/(12.0*dx(3))
	      enddo
	    enddo
	   enddo
         case(6)
          do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
		 res(i,j,k,1) = (-1.0*var(i-3,j,k)+9.0*var(i-2,j,k)-45.0*var(i-1,j,k)+45.0*var(i+1,j,k)-9.0*var(i+2,j,k)+var(i+3,j,k))/(60.0*dx(1))
                 res(i,j,k,2) = (-1.0*var(i,j-3,k)+9.0*var(i,j-2,k)-45.0*var(i,j-1,k)+45.0*var(i,j+1,k)-9.0*var(i,j+2,k)+var(i,j+3,k))/(60.0*dx(2))
                 res(i,j,k,3) = (-1.0*var(i,j,k-3)+9.0*var(i,j,k-2)-45.0*var(i,j,k-1)+45.0*var(i,j,k+1)-9.0*var(i,j,k+2)+var(i,j,k+3))/(60.0*dx(3))

	      enddo
	    enddo
	   enddo


         



	   case default 
	    print*, " This orcder not implemented yet"
	   end select

	 end subroutine

