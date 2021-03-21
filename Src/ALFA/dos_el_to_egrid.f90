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

!>  unfolds the band information from a wedge
!>  into the regular grid

subroutine dos_el_to_egrid(el,nrk,nband,nx,kmap,egrid,nxmax,mxdbnd)

! Written November 19, 2013. jlm
! Modified May 26 2014. JLM
! Modified, documentation, 19 September 2020. JLM
! Modified egrid, 18 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)             ::  nrk                                !<  number of k-points
  integer, intent(in)             ::  mxdbnd                             !<  size of number of bands
  integer, intent(in)             ::  nxmax(3)                           !<  maximum number of k-points in each direction in regular grid
  integer, intent(in)             ::  nx(3)                              !<  number of k-points in each direction in regular grid
  real(REAL64), intent(in)        ::  el(mxdbnd,nrk)                     !<  eigenvalues in Hartree
  integer, intent(in)             ::  nband(nrk)                         !<  number of bands for each k-points
  integer, intent(in)             ::  kmap(3,nx(1),nx(2),nx(3))          !<  index of point in the irreducible zone corresponding to point in the regular mesh         

! output

  real(REAL64), intent(out) ::   egrid(nxmax(1),nxmax(2),nxmax(3),mxdbnd)   !<  eigenvalues in Hartree in regular grid

! local variables

  integer  ::  neig

! counters

  integer  ::  n, i, j, k

  neig = nband(1)
  do n = 1,nrk
    if(nband(n) < neig) neig = nband(n)
  enddo

  do n = 1,neig
    do k = 1,nx(3)
    do j = 1,nx(2)
    do i = 1,nx(1)
      egrid(i,j,k,n) = el(n,iabs(kmap(1,i,j,k)))
    enddo
    enddo
    enddo
  enddo

  return
end subroutine dos_el_to_egrid
