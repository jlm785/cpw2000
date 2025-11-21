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


!>  mgga exchange and LDA/mgga correlation.  Only potential (e.g. Tran-Blaha)
!>
!>  \author       Carlos Loia Reis
!>  \version      5.12
!>  \date         September 2015, 21 November 2025.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_vxc(author_x, author_c, rho, grho, lap_rho, tau,      &
                     epsx, epsc, vx, vc, ctb09 )

!

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author_x                        !<  type of exchange wanted.
  character(len = *), intent(in)     ::  author_c                        !<  type of correlation wanted.
  real(REAL64), intent(in)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64), intent(in)           ::  grho                            !<  gradient of charge density
  real(REAL64), intent(in)           ::  lap_rho                         !<  Laplacian of charge density
  real(REAL64), intent(in)           ::  tau                             !<  The kinetic energy density
  real(REAL64), intent(in)           ::  ctb09                           !<  TB constant

! output

  real(REAL64), intent(out)          ::  epsx                            !<  ZERO exchange energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  epsc                            !<  ZERO correlation energy density (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vx                              !<  Tran-Blaha exchange potential (hartree/bohr^3)
  real(REAL64), intent(out)          ::  vc                              !<  LDA correlation potential (hartree/bohr^3)

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-24_REAL64


  epsx = ZERO
  epsc = ZERO


  if(rho < EPS) then

    vx = ZERO
    vc = ZERO

  else

    if (author_x == 'tbl' .or. author_x == 'TBL' .or.                    &
        author_x == 'tb00' .or. author_x == 'TB09') then

      call xc_lda(author_c, rho, epsx, epsc, vx, vc )
      call xc_mgga_x_tb09(rho, grho, lap_rho, tau, ctb09, vx)

    else

      write(6,'("  STOPPED  in xc_mgga_vxc:   unknown exchange-correlation")')
      write(6,*) '     ',author_x,'  ', author_c

      stop

    endif

  endif

  return

end subroutine xc_mgga_vxc



!>  Calculates the Tran-Blaha exchange potential
!>  PhysRevLett, 102, 226401 (2009).
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.04
!>  \date         September 2015, 25 September 2022.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_x_tb09(rho, grho, lap_rho, tau, set_ctb09, vrho)

! Written September 2015.  CLR
! version 1.0 of xc, September 2015
! Modified, documentation, December 2019. JLM

  implicit none

  integer, parameter           :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           :: rho                         !<  charge density  (electrons per unit cell)
  real(REAL64), intent(in)           :: grho                        !<  gradient of charge density
  real(REAL64), intent(in)           :: lap_rho                     !<  Laplacian of charge density
  real(REAL64), intent(in)           :: tau                         !<  kinetic energy density
  real(REAL64), intent(in)           :: set_ctb09                   !<  TB constant

! output

  real(real64), intent(out)          :: vrho                        !<  Tran-Blaha exchange potential

! local variables

  real(REAL64)                 :: v_BR                              !  Becke and Roussel exchange
  real(REAL64)                 :: c_HEG                             !  Tran-Blaha correction

! parameters

  real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter   :: MIN_DENS = 5.0E-13_REAL64
  real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64


! Becke-Roussel exchange is with spin-density

  call xc_mgga_x_br89(rho/2, grho*grho/4, lap_rho/2, tau/2, v_BR)

  c_HEG  = (3*set_ctb09 - 2*UM)*sqrt((5*UM)/(12*UM)) / PI
  vrho  = set_ctb09*v_BR + c_HEG*sqrt(2*tau/rho)

  return

end subroutine xc_mgga_x_tb09




!>  Calculates the Becke-Roussel exchange potential
!>  Becke and Roussel PRA 39, 3761 (1989)
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.04
!>  \date         September 2015, 22 September 2022.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_x_br89(rho, sigma, lap_rho, tau, v_BR)

! Written September 2015.  CLR
! version 1.0 of xc, September 2015
! Modified, documentation, December 2019. JLM
! Simplified, 25 September 2022. JLM

  implicit none

  integer, parameter           :: REAL64 = selected_real_kind(12)

  real(REAL64), intent(in)                ::  rho                        !<  spin-density
  real(REAL64), intent(in)                ::  sigma                      !<  gradient squared of the spin density
  real(REAL64), intent(in)                ::  lap_rho                    !<  laplacian of the spin-density
  real(REAL64), intent(in)                ::  tau                        !<  TWICE the kinetic energy of the spin-density

  real(REAL64), intent(out)               ::  v_BR                       !<  Becke-Roussel exchange potential

  real(REAL64)                 :: br_Q
  real(REAL64)                 :: br_x

  real(REAL64)                :: rho1T, rho5T

  real(REAL64), parameter   :: BR_GAMMA = 0.8_REAL64                     !  adjusted parameter of the BR exchange

  real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter   :: UM = 1.0_REAL64
  real(REAL64), parameter   :: PITHRD = PI**(UM/3*UM)
  real(REAL64), parameter     :: MIN_TAU  = 1.0E-3_REAL64

  rho1T = rho**(UM/(3*UM))
  rho5T = rho*rho1T*rho1T

! beware that tau in Becke-Roussel paper is twice the kinetic energy density

  br_Q  = (lap_rho - 2*BR_GAMMA*2*tau + BR_GAMMA*sigma/(2*rho))/(6*UM)
  br_Q  = br_Q / rho5T

  call xc_mgga_num_br_find_x(br_Q,br_x)

  if ((br_x) > MIN_TAU) then
    v_BR =  exp(br_x/3) * (UM - exp(-br_x)*(UM + br_x/(2*UM))) / br_x
  else
    v_BR = UM/(2*UM)*(UM + (br_x/(3*UM))*(UM - br_x/(3*UM)*(UM - 11*br_x/(36*UM))))
  endif

  v_BR = -2*PITHRD*v_BR*rho1T

  return

end subroutine xc_mgga_x_br89




!>  Solves equation 21 of Becke and Roussel PRA 39, 3761 (1989)
!>  with Newton-Raphson.
!>
!>  \author       Carlos Loia Reis
!>  \version      5.04
!>  \date         September 2015, 22 September 2022.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_num_br_find_x(br_Q,br_x)
! libxc 1.2.0
! Modified, documentation, December 2019, 24 September 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  real(REAL64), intent(in)          :: br_Q                              !<  Q / rho_sigma^5/3 of eq. 21
  real(REAL64), intent(out)         :: br_x                              !<  solution x

  real(REAL64)                      :: res                               !   residual value
  real(REAL64)                      :: rhs                               !   right han sde of equation

  integer                           :: ierr                              !   1: success;  0: error

  real(REAL64), parameter   :: TOL = 5.0E-12_REAL64
  real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter   :: UM = 1.0_REAL64

  rhs = ((2*UM/(3*UM)) * PI**(2*UM/(3*UM))) / br_Q

! Newton-Raphson should always work

  call xc_mgga_num_br_newt_raph(br_x, rhs, TOL, res, ierr)

  if(ierr == 0) then

! Try bissection

     call xc_mgga_num_br_bisect(br_x, rhs, TOL, ierr)

     if(ierr == 0) then

!      This should never happen

       write(6,*)
       write(6,*) '   STOPPED in xc_becke_roussel_x: Convergence',       &
          ' not reached in Becke-Roussel functional'
       write(6,*)

       write(6,*) '  rhs:      ', rhs
       write(6,*) '  residual: ', res

       stop

     endif
  endif

  return

end subroutine xc_mgga_num_br_find_x




!>  Solves equation 21 of Becke and Roussel PRA 39, 3761 (1989)
!>  with Newton-Raphson.
!>
!>  \author       Carlos Loia Reis
!>  \version      5.034
!>  \date         September 2015, 22 September 2022.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_num_br_newt_raph(xc_br, a, tol, res, ierr)

! libxc 1.2.0
! Modified, documentation, December 2019. JLM
! Modified, starting point, subroutine, September 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  real(REAL64), intent(out)               :: xc_br                       !<  solution of equation

  real(REAL64), intent(in)                :: a                           !<  parameter "a"
  real(REAL64), intent(in)                :: tol                         !<  tolerance for successful solution
  real(REAL64), intent(out)               :: res                         !<  residual value
  integer, intent(out)                    :: ierr                        !<  1: success;  0: error

  integer                     :: icount
  real(REAL64)                :: x
  real(REAL64)                :: f

  real(REAL64)                ::  eee, fp

  real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  integer, parameter          :: max_iter = 50

  ierr = 0

! a is calculated as an inverse, so is never zero.
! the function diverges to +-infinity as a goes to +-0

  if(a == ZERO) then
    ierr = 1
    xc_br = ZERO
  else

! /* starting point */

    if(a < -0.75) then
      x = 1.5
    elseif(a < ZERO) then
      x = -2*a
    elseif(a < 0.15) then
      x = -1.5*log(a)
    elseif(a < 0.5) then
      x = 3
    else
      x =  2 + 1 / (2*a)
    endif

    icount = 0
    do icount=1, max_iter
       eee = exp(-2*x/(3*UM)) / a

       f  = x*eee - (x - 2*UM)
       fp = eee*(UM - 2*x/(3*UM)) - UM

       x  = x - f/fp

!      physical x is positive, should never occur
       x  = abs(x)

       res = abs(f)

       if(res <= tol) then

         ierr = 1
         xc_br = x

         exit

       endif

    enddo

  endif

  return

end subroutine xc_mgga_num_br_newt_raph



!>  Solves equation 21 of Becke and Roussel PRA 39, 3761 (1989)
!>  with Bissection.
!>
!>  \author       Carlos Loia Reis
!>  \version      5.04
!>  \date         September 2015, 22 September 2022.
!>  \copyright    GNU Public License v2

subroutine xc_mgga_num_br_bisect(xc_br, a, tol, ierr)

! libxc 1.2.0
! Modified, documentation, December 2019. JLM
! Modified, bug negative a, subroutine, September 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  real(REAL64), intent(out)               :: xc_br                       !<  solution of equation

  real(REAL64), intent(in)                :: a                           !<  parameter "a"
  real(REAL64), intent(in)                :: tol                         !<  tolerance for successful solution
  integer, intent(out)                    :: ierr                        !<  1: success;  0: error

  integer                     :: icount
  real(REAL64)                :: f, x, x1, x2

  real(REAL64)                :: eee

  real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  integer, parameter          :: max_iter = 500

  ierr = 0

! a is calculated as an inverse, so is never zero.
! the function diverges to +-infinity as a goes to +-0

  if(a == ZERO) then
    ierr = 1
    xc_br = ZERO
  else

!   /* starting interval */

    if(a > ZERO) then
      x1 = 2*UM + tol
      x2 = UM/a + 2*UM
    else
      x2 = 2*UM - tol
      x1 = ZERO
    endif

!   /* bisection */
    do icount=1,max_iter
      x   = (x1 + x2)/(2*UM)
!     xm2 = x - 2*UM
!     arg = 2*x/(3*UM)
      eee = exp(-2*x/(3*UM))
      f   = x*eee - a*(x - 2*UM)

      if(a > ZERO) then
        if(f > ZERO) x1 = x
        if(f < ZERO) x2 = x
      else
        if(f > ZERO) x2 = x
        if(f < ZERO) x1 = x
      endif

      if (abs(f) < tol) then

        ierr = 1
        xc_br = x

        exit

      endif
    enddo

  endif

  return

end subroutine xc_mgga_num_br_bisect

