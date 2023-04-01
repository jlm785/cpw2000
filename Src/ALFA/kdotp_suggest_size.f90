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

!>  Suggests good sizes of kdotp matrices from symmetry and
!>  free electron bands.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         21 March 2023.
!>  \copyright    GNU Public License v2

subroutine kdotp_suggest_size(ztot, adot, ng, kgv, rkin, kdotpsize, nk,  &
    mxdgve)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  real(REAL64)                       ::  ztot                            !<  total charge density (electrons/cell)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  real(REAL64), intent(in)           ::  rkin(3)                         !<  k-vector (reciprocal lattice coordinates)

! output

  integer, intent(out)               ::  kdotpsize(10)                   !<  suggested sizes for kdotp problem
  integer, intent(out)               ::  nk                              !<  number of suggested valuex.

! allocatable arrays

  real(REAL64), allocatable          ::  efree(:)                        !  free electron eigenvalue
  real(REAL64), allocatable          ::  egap(:)                         !  free electron gaps
  integer, allocatable               ::  indx(:)                         !  gaps in increasing size

! local variables

  integer           ::  nfirst, nlast
  real(REAL64)      ::  bdot(3,3), vcell
  real(REAL64)      ::  r1, r2, r3

! constants

  real(REAL64), parameter     ::  EPS = 1.0E-6_REAL64

! counters

  integer    ::  n


! reasonable range

  nfirst = nint(ztot)
  nlast = nfirst*3

  allocate(efree(nlast-nfirst+1))
  allocate(egap(nlast-nfirst))
  allocate(indx(nlast-nfirst))

  call adot_to_bdot(adot,vcell,bdot)

  do n = nfirst,nlast
    r1 = rkin(1)+kgv(1,n)
    r2 = rkin(2)+kgv(2,n)
    r3 = rkin(3)+kgv(3,n)
    efree(n-nfirst+1) = r1*bdot(1,1)*r1 + r1*bdot(1,2)*r2 + r1*bdot(1,3)*r3 +    &
                        r2*bdot(2,1)*r1 + r2*bdot(2,2)*r2 + r2*bdot(2,3)*r3 +    &
                        r3*bdot(3,1)*r1 + r3*bdot(3,2)*r2 + r3*bdot(3,3)*r3
  enddo

  do n = nfirst,nlast-1
    egap(n-nfirst+1) = efree(n-nfirst+2) - efree(n-nfirst+1)
  enddo

  call sort(nlast-nfirst, egap, indx)

  nk = 0
  do n = nlast-1,nfirst,-1

    if(abs(egap(indx(n-nfirst+1))) > EPS) then
      nk = nk + 1
      kdotpsize(nk) = indx(n-nfirst+1) + nfirst-1
    endif

    if(nk == 10) exit

  enddo

  deallocate(efree)
  deallocate(egap)
  deallocate(indx)

  return

end subroutine kdotp_suggest_size
