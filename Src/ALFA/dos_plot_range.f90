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

!>  finds a reasonable scale for a density of states plot

  subroutine dos_plot_range(nx,egrid,np,emin,deltae,nhist,ezero,nxmax,mxdbnd)

! Written December 9, 2013.
! Modified, documentation, 19 September 2020. JLM
! Modified, egrid, 18 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)             ::  mxdbnd                             !<  size of number of bands
  integer, intent(in)             ::  nxmax(3)                           !<  maximum number of k-points in each direction in regular grid
  integer, intent(in)             ::  nx(3)                              !<  number of k-points in each direction in each regular grid
  real(REAL64), intent(in)  ::  egrid(nxmax(1),nxmax(2),nxmax(3),mxdbnd) !<  eigenvalues in Hartree in regular grid
  integer, intent(in)             ::  np                                 !<  approximate number of points for plot
  real(REAL64), intent(in)        ::  ezero                              !<  zero of energy

! output

  real(REAL64), intent(out)       ::  emin                               !<  minimum recomended energy for histogram 
  real(REAL64), intent(out)       ::  deltae                             !<  recomended energy step for histogram 
  integer, intent(out)            ::  nhist                              !<  recomended number of points for histogram

! local variables

  real(REAL64)  ::  emax, xmax, b

! counters

  integer   ::  i, j ,k, m

! constants

  real(REAL64), parameter :: HARTREE = 27.21138386_REAL64

! find emin and emax

  emin = egrid(1,1,1,1) - ezero
  emax = egrid(1,1,1,1) - ezero
  do j = 1,mxdbnd
    emax = max(emax,egrid(1,1,1,j) - ezero)
  enddo

  do k = 1,nx(3)
  do j = 1,nx(2)
  do i = 1,nx(1)

    xmax = egrid(i,j,k,1)
    do m = 1,mxdbnd
      xmax = max(xmax,egrid(i,j,k,m) - ezero)
      emin = min(emin,egrid(i,j,k,m) - ezero)
    enddo
    if(emax > xmax) emax = xmax

  enddo
  enddo
  enddo

! find scale for the plot

  deltae = HARTREE*(emax-emin) / np

  call plot_step(deltae,b)

  deltae = b/HARTREE

  emin = deltae*real(int(emin/deltae)-3)
  emax = emax + deltae
  nhist = nint((emax-emin)/deltae) + 1

  return
end subroutine dos_plot_range
