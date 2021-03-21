!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>  Kramers–Kronig transformation from the imaginary part of the dielectric function.
!>  See 6.1.3 in Fundamentals of Semiconductors By Cardona.

subroutine opt_kk_im(ehist, e_im, e_re, nhist)

! Written by Carlos Loia Reis. July 2020
! Modified, documentation, 20 September 2020. JLM
! Modified, tail of e_im, 17 December 2020. JLM
! copyright  Carlos Loia Reis/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)         ::  nhist                                  !<  grid size

  real(REAL64), intent(in)    ::  ehist(nhist)                           !<  frequencies
  real(REAL64), intent(in)    ::  e_im(nhist)                            !<  imaginary response

! output

  real(REAL64), intent(out)   :: e_re(nhist)                             !<  real responsse

! local allocatable arrays (for spline interpolation)

  real(REAL64), allocatable   ::  x(:)
  real(REAL64), allocatable   ::  y(:)
  real(REAL64), allocatable   ::  fout(:)

  real(REAL64), allocatable   ::  f(:)
  real(REAL64), allocatable   ::  fp(:)
  real(REAL64), allocatable   ::  fpp(:)
  real(REAL64), allocatable   ::  yp(:)
  real(REAL64), allocatable   ::  ypp(:)
  real(REAL64), allocatable   ::  w(:,:)

! local variables

  integer                     ::  n, km
  real(REAL64)                ::  pow

  real(REAL64)                ::  a1,b1,an,bn

  real(REAL64)                ::  ans1(1), ans2(1)

  real(REAL64)                ::  a
  integer                     ::  ierr, isx

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  real(REAL64), parameter     ::  UM = 1.0_REAL64
  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-8_REAL64

! counters

  integer :: i, j


! triplicates range

  n = 3*nhist

  allocate(x(n),y(n),fout(n))
  allocate(f(n), fp(n), fpp(n), yp(n), ypp(n), w(n,3))

  do i = 1,nhist
    x(i) = ehist(i)
  enddo
  do i = nhist+1,n
    x(i) = ehist(nhist) + (i-nhist)*( (ehist(nhist)-ehist(1)) / (nhist-1) )
  enddo

  do i = 1,nhist
    y(i) = e_im(i)
  enddo
  
! tries a power law based on last 10% of e_im

  km = max(1,nhist/10)

  if(e_im(nhist) < EPS .or. e_im(nhist-km) < EPS) then
  
    do i = nhist+1,n
      y(i) = ZERO
    enddo    

  else

    pow = log(e_im(nhist) / e_im(nhist-km)) / log(ehist(nhist-km) / ehist(nhist))
    if(pow < 2.0) pow = 2.0
    if(pow > 6.0) pow = 6.0

    do i = nhist+1,n
      y(i) = (x(nhist) / x(i))**pow * e_im(nhist)
    enddo

  endif

! Kramers-Kronig on extended range

  a1 = ZERO 
  b1 = ZERO 
  an = ZERO 
  bn = ZERO 

  isx = 0
  call splift (x,y,yp,ypp,n,w,ierr,isx,a1,b1,an,bn)

  do j = 2,n-1

    a = x(j)

    do i = 1,n

      if(i /= j) then

        f(i)  = x(i)*y(i) /( (-a*a+x(i)*x(i)) )

        fp(i) =       -2*x(i)*x(i)*y(i)/( (-a*a+x(i)*x(i))*(-a*a+x(i)*x(i)) )
        fp(i) = fp(i) +y(i)/( (-a*a+x(i)*x(i)))
        fp(i) = fp(i) +x(i)*yp(i)/( (-a*a+x(i)*x(i)))

        fpp(i) = 8*x(i)*x(i)*x(i)*y(i)/( (-a*a+x(i)*x(i))*(-a*a+x(i)*x(i))*(-a*a+x(i)*x(i)) )
        fpp(i) = fpp(i) -  6*x(i)*y(i)/( (-a*a+x(i)*x(i))*(-a*a+x(i)*x(i)) )
        fpp(i) = fpp(i) -  4*x(i)*x(i)*yp(i)/( (-a*a+x(i)*x(i))*(-a*a+x(i)*x(i)) )
        fpp(i) = fpp(i) +  2*yp(i)/( (-a*a+x(i)*x(i)) )
        fpp(i) = fpp(i) +  x(i)*ypp(i)/( (-a*a+x(i)*x(i)) )

      endif
    enddo

    call spliq(x,f,fp,fpp,n,x(1),x(j-1),1,ans1,ierr)
    call spliq(x,f,fp,fpp,n,x(j+1),x(n),1,ans2,ierr)

    fout(j) =  UM + (ans1(1)+ans2(1))*2/PI

  enddo

! extrapolation for 1 and n

  fout(1) = fout(2) - (fout(3)-fout(2))*(x(2)-x(1))/(x(3)-x(2))
  fout(n) = fout(n-1) + (fout(n-1)-fout(n-2))*(x(n)-x(n-1))/(x(n-1)-x(n-2))

  do i = 1,nhist
    e_re(i) = fout(i)
  enddo

  deallocate(f, fp, fpp, yp, ypp, w)
  deallocate(x, y, fout)

end subroutine opt_kk_im

