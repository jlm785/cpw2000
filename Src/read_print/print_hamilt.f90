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

!>  Writes (part of) the hamiltonian
!>
!>  \author       José Luís Martins
!>  \version      5.12
!>  \date         Before November 2015. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine print_hamilt(icmplx, rkpt, ham, mtxd, nspin,                  &
    adot, isort, kgv,                                                    &
    mxddim, mxdgve, mxdham)

! Adapted frm print_eig, before November 2015
! Modified, documentation, 26 January 2020. JLM
! Added spin polarization option, indentation, 4 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdham                          !<  array dimension of hamiltonian

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  j-th component in lattice coordinates of the k-point
  complex(REAL64), intent(in)        ::  ham(mxdham,mxdham)              !<  hamiltonian matrix
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

! local varaibles

  real(REAL64)        ::  avec(3,3),bvec(3,3),rkcar(3)

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64

! counters

  integer    ::  i, j, k, jmax, nc, ii



  if(nspin /=1 .and. nspin /= 2) then
    write(6,*)
    write(6,*) '    STOPPED in print_hamilt, nspin = ', nspin
    write(6,*)

    STOP

  endif

! calculates coordinates

  call adot_to_avec_sym(adot,avec,bvec)

  write(6,*)
  write(6,*)

  do i = 1,mtxd

    if(nspin == 1) then
      ii = i
    else
      ii = (i+1)/2
    endif

    do k = 1,3
      rkcar(k) = ZERO
      do j=1,3
        rkcar(k) = rkcar(k) + bvec(k,j)*(rkpt(j) + kgv(j,isort(ii)))
      enddo
    enddo

    do j = 1,mtxd,8

      jmax = j+7
      if (jmax > mtxd) jmax = mtxd

      if(j == 1) then
        write(6,'(3f10.4,5x,8f12.5)') (rkcar(k),k=1,3), (real(ham(nc,i)),nc=j,jmax)
      else
        write(6,'(35x,8f12.5)') (real(ham(nc,i)),nc=j,jmax)
      endif
      if(icmplx /= 0) then
        write(6,'(35x,8f12.5)') (aimag(ham(nc,i)),nc=j,jmax)
      endif

    enddo

  enddo

  write(6,*)
  write(6,*)

  return

end subroutine print_hamilt
