#include "REAL.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"
#include "SPACE.H"

	subroutine derivative(lo,hi,var,DIMS(var),res,DIMS(res),dx,order,dir)
	implicit none
	integer  lo(3), hi(3), DIMDEC(var),DIMDEC(res), dir
	REAL_T var(DIMV(var))
	REAL_T res(DIMV(res)),dx(3)
	integer order,i,j,k

	Select case(order)
	  case(2)
	    select case(dir)
	     case(1)
	   do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
		 res(i,j,k) = (var(i+1,j,k)-var(i-1,j,k))*0.5/dx(1)

	      enddo
	    enddo
	   enddo


	   case(2)
	  do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
		 res(i,j,k) = (var(i,j+1,k)-var(i,j-1,k))*0.5/dx(1)

	      enddo
	    enddo
	   enddo

	  case(3)
	  do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
		 res(i,j,k) = (var(i,j,k+1)-var(i,j,k-1))*0.5/dx(1)

	      enddo
	    enddo
	   enddo
	   case default 
	    print*, "wrong direction"
	   end select
	 case(4)
	    select case(dir)
	     case(1)
	   do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
		 res(i,j,k) = (-1.0*var(i+2,j,k)+8.0*var(i+1,j,k)-8.0*var(i-1,j,k)+var(i-2,j,k))/(12.0*dx(1))

	      enddo
	    enddo
	   enddo


	   case(2)
	  do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		  res(i,j,k) = (-1.0*var(i,j+2,k)+8.0*var(i,j+1,k)-8.0*var(i,j-1,k)+var(i,j-2,k))/(12.0*dx(1))

	      

	      enddo
	    enddo
	   enddo

	  case(3)
	  do k= lo(3),hi(3)
	    do j=lo(2), hi(2)
	      do i = lo(1),hi(1)
		     
	   res(i,j,k) = (-1.0*var(i,j,k+2)+8.0*var(i,j,k+1)-8.0*var(i,j,k-1)+var(i,j,k-2))/(12.0*dx(1))


	      enddo
	    enddo
	   enddo
          case default 
            print*, "Wrong Direction "
          end select

	   case default 
	    print*, " This orcder not implemented yet"
	   end select

	 end subroutine

