module prob_functions_mod

	public :: randn, std, idum, normal_pdf,uniform_pdf,normal_sample,nor2par,par2nor

	double precision :: idum

contains

!-----------------------------------------------------
!-  Normal PDF
!-----------------------------------------------------
double precision FUNCTION normal_pdf( x, a, b)
	!
	!! NORMAL_PDF evaluates the Normal PDF.
	!
	!  Discussion:
	!
	!    PDF(A,B;X)
	!      = exp ( - 0.5D+00 * ( ( X - A ) / B )^2 ) / ( B * sqrt ( 2 * PI ) )
	!
	!    The normal PDF is also known as the Gaussian PDF.
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    10 February 1999
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) X, the argument of the PDF.
	!
	!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
	!    0.0D+00 < B.
	!
	!    Output, real ( kind = 8 ) PDF, the value of the PDF.
	!

  	double precision ::  a
  	double precision :: b
  	double precision, parameter :: pi = 3.141592653589793D+00
  	double precision :: x
  	double precision :: y

  	y = ( x - a ) / b

  	normal_pdf = exp ( - 0.5D+00 * y * y )  / ( b * sqrt ( 2.0D+00 * pi ) )

end Function

!-----------------------------------------------------
!-  Uniform PDF
!-----------------------------------------------------
double precision function uniform_pdf ( x, a, b)
	!
	!UNIFORM_PDF evaluates the Uniform PDF.
	!
	!  Discussion:
	!
	!    The Uniform PDF is also known as the "Rectangular" or "de Moivre" PDF.
	!
	!    PDF(A,B;X) = 1 / ( B - A ) for A <= X <= B
	!               = 0 otherwise
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    01 February 1999
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) X, the argument of the PDF.
	!
	!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
	!    A < B.
	!
	!    Output, real ( kind = 8 ) PDF, the value of the PDF.
	!

  	double precision :: a
  	double precision :: b
  	double precision :: x

  	if ( x < a .or. b < x ) then
    	uniform_pdf = 0.0D+00
  	else
    	uniform_pdf = 1.0D+00 / ( b - a )
  	end if

	return
end Function

!-----------------------------------------------------
!-  RNORM
!-----------------------------------------------------
	double precision function normal_sample ( a, b)

	!*****************************************************************************80
	!
	!! NORMAL_SAMPLE samples the Normal PDF.
	!
	!  Discussion:
	!
	!    The Box-Muller method is used.
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    10 October 2004
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
	!    0.0D+00 < B.
	!
	!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
	!    generator.
	!
	!    Output, real ( kind = 8 ) X, a sample of the PDF.
	!
  	implicit none

  	double precision ::  a
  	double precision :: b
  	double precision :: x

	
  	x = normal_01_sample ()

  	normal_sample = a + b * x
  	
  	return
end FUNCTION

double precision function normal_01_sample ()
	!*****************************************************************************80
	!
	!NORMAL_01_SAMPLE samples the standard normal probability distribution.
	!
	!  Discussion:
	!
	!    The standard normal probability distribution function (PDF) has
	!    mean 0 and standard deviation 1.
	!
	!    The Box-Muller method is used, which is efficient, but
	!    generates two values at a time.
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    26 August 2013
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number	
	!    generator.
	!
	!    Output, real ( kind = 8 ) X, a sample of the standard normal PDF.
	!
  	implicit none
  	double precision:: r1
  	double precision:: r2
  	double precision, parameter :: r8_pi = 3.141592653589793D+00
  	!double precision:: x
  
  	call RANDOM_NUMBER(r1)
  	call RANDOM_NUMBER(r2)
  	!r1 = randn(0)
  	!r2 = randn(0)
  	  	
  	normal_01_sample = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  	return
end function

!------ Standard deviation calculation 
double precision function std(a,n)

	! Function to determine the standard deviation
    ! inputs are the vector of values and number of values included

    implicit none

    ! declare inputs
    integer, intent(in) :: n ! number of values in vector
    double precision, dimension(n), intent(in) :: a

    ! declare local variables
    double precision mean, sq_diff_sum, diff, variance
    integer i

    ! if no length has been returned then provide value which ensures crash (i.e.
    ! infinity)
    if (n == 0) then
        std=0.
        write(*,*) "no sample size has been provided to std function"
        return
    endif

    ! first calculate the mean
    mean=sum(a)/n
    
    ! ensure zero values
    diff=0. ; sq_diff_sum=0.

    ! calculate cumulative square difference
    !sq_diff_sum=(sum(a-mean))**2.
    do i=1,n
       diff=a(i)-mean
       sq_diff_sum=sq_diff_sum + diff**2.
    end do

    ! calculate the variance
    variance=sq_diff_sum/(n-1)
    
    ! return the standard deviation
    std=sqrt(variance)

    return

end function std

double precision function randn(option)

    ! from Numerical Receipes p271 Press et al., 1986 2nd Edition Chapter 7, Random
    ! Numbers
    ! function returns real random number between 0-1 based on an initial start
    ! point.
    ! The start point (default = -1) is reinitialised every time the model runs
    ! providing the same distribution each run
    ! To ensure random numbers each time use the sum of the current system time

    ! modified based on blooms C code to alter range of random numbers

    implicit none
    
    integer IA,IM,IQ,IR,NTAB,NDIV,option
    double precision AM,EPS,RNMX,const,r1,r2,pi
    parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-30,RNMX=1.-EPS)
    INTEGER j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/
    DATA iy /0/
    
    const=1.0
    pi=3.141592653589793

    if (option == 0) then
      if (idum < 0. .or. iy == 0) then
          idum=MAX(-idum,const)
          do j=(NTAB+8),1,-1
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum < 0.) idum=idum+IM
             if (j < NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum < 0.) idum=idum+IM
      j=1.+iy/NDIV
      iy=iv(j)
      iv(j)=idum

      ! output now
      randn=MIN(AM*iy, RNMX)
      return

    else

      if (idum < 0. .or. iy == 0) then
          idum=MAX(-idum,const)
          do j=(NTAB+8),1,-1
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum < 0.) idum=idum+IM
             if (j < NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum < 0.) idum=idum+IM
      j=1.+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      r1=MAX(MIN(DBLE(AM)*DBLE(iy), DBLE(RNMX)),DBLE(1e-30))

      if (idum < 0. .or. iy == 0) then
          idum=MAX(-idum,const)
          do j=(NTAB+8),1,-1
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum < 0.) idum=idum+IM
             if (j < NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum < 0.) idum=idum+IM
      j=1.+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      r2=MAX(MIN(DBLE(AM)*DBLE(iy), DBLE(RNMX)),DBLE(1e-30))

      ! output now
      randn=sqrt(-2.*log(r1)) * cos(2.*pi*r2)
      return

    endif

  end function
  
  double precision function par2nor(initial_par,min_par,max_par)

	!----------------------------------------------------------
    ! functions to normalised log parameter values and return them back to
    ! unlogged
    ! un-normalised value

    ! converting parameters on log scale between 0-1 for min/max values
    implicit none
    double precision initial_par, min_par, max_par

    if (max_par > 0. .and. min_par < 0.) then
        ! then normalise without logs and we cross zero
        par2nor=(initial_par-min_par)/(max_par-min_par)
    else
        par2nor=log(initial_par/min_par)/log(max_par/min_par)
    end if

    return

  end function par2nor
  !
  !---------------------and vise versa ------------------------------
  ! 
  double precision function nor2par(initial_par,min_par,max_par)

    ! converting values back from log scale 0-1 to 'real' numbers
    implicit none
    double precision initial_par, min_par, max_par

    if ( max_par > 0. .and. min_par < 0.) then
        ! then un-normalise without logs as we cross zero and logs wont work
        nor2par=min_par+(max_par-min_par)*initial_par
    else
        nor2par=min_par*((max_par/min_par)**initial_par)
    endif

    return

  end function nor2par
  
end module prob_functions_mod
