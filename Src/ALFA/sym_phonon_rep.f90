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

!>  Given phonon modes, groups them in irreducible representations.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         5 December 2022.
!>  \copyright    GNU Public License v2

subroutine sym_phonon_rep(adot, qvec, nat_qe, freqsq, uvec, nrep, irep,  &
        irt_qe, rtau_qe, invs_qe,                                        &
        ntrans, mtrx, tnp)

! Adapted from other sym subroutines. 27 November 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nat_qe                          !<  total number of atoms

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  qvec(3)                         !<  wave-vector in lattice coordinates

  integer, intent(in)                ::  irt_qe(48,nat_qe)               !<  index of the rotated atom
  real(REAL64), intent(in)           ::  rtau_qe(3,48,nat_qe)            !<  rotated atom minus original, lattice coordinates.
  integer, intent(in)                ::  invs_qe(48)                     !<  index of the inverse operation

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the little group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the little group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the little group (unused).

! input and output

  real(REAL64), intent(inout)        ::  freqsq(3*nat_qe)                !<  eigenvalue associated with uvec. Modified order.
  complex(REAL64), intent(inout)     ::  uvec(3*nat_qe,3*nat_qe)         !<  proper modes in cartesian coordinates.  Modified order.

! output

  integer, intent(out)               ::  nrep                            !<  number of representations present
  integer, intent(out)               ::  irep(3*nat_qe)                  !<  representation of each mode

! allocatable arrays

  complex(REAL64), allocatable       ::  ctmp(:)
  integer, allocatable               ::  ideg(:)                         !  degeneracy of eigenvalue
  complex(REAL64), allocatable       ::  ruvec(:,:)

! local variables

  real(REAL64)         ::  tmp
  integer              ::  itmp, itdeg

  integer              ::  insert, ifirst, ilast

  real(REAL64)         ::  rot(3,3,48)                       !  rotational matrix in cartesian coordinates
  real(REAL64)         ::  tau(3,48)                         !  fractional translation vector in cartesian coordinates

  complex(REAL64)      ::  prod

  complex(REAL64), external   ::  zdotc

! constants

  real(REAL64), parameter  ::  EPS = 1.0E-4_REAL64                       !  tolerance for degeneracy

! counters

  integer    ::  i, j, n

! Eigenvalues are expected to be in increasing order.  Do an insertion sort just in case.
! It is fast in the expected case!

  allocate(ctmp(3*nat_qe))

  do j = 2,3*nat_qe
    if(freqsq(j) - freqsq(j-1) < -EPS*EPS) then
      tmp = freqsq(j)
      ctmp(:) = uvec(:,j)
      insert = 1
      do i = j-1,1,-1
        insert = i
        if(freqsq(i) - tmp < -EPS*EPS) then
          insert = i+1

          exit

        else
          freqsq(i+1) = freqsq(i)
          uvec(:,i+1) = uvec(:,i)
        endif
      enddo
      freqsq(insert) = tmp
      uvec(:,insert) = ctmp(:)
    endif
  enddo

! computes the eigenvalue degeneracy

  allocate(ideg(3*nat_qe))
  ideg(:) = 0

  ifirst = 1
  ilast = 1
  do j = 2,3*nat_qe
    if(abs(freqsq(j) - freqsq(j-1)) > EPS) then
      do i = ifirst,ilast
        ideg(i) = ilast-ifirst+1
      enddo
      ifirst = j
    endif
    ilast = j
  enddo
  do i = ifirst,ilast
    ideg(i) = ilast-ifirst+1
  enddo

! Do an insertion sort on degeneracy.  It is not optimal but it does not
! break the grouping by eigenvalue.

  do j = 2,3*nat_qe
    if(ideg(j) - ideg(j-1) < 0) then
      tmp = freqsq(j)
      ctmp(:) = uvec(:,j)
      itmp = ideg(j)
      insert = 1
      do i = j-1,1,-1
        insert = i
        if(ideg(i) - itmp <= 0) then
          insert = i+1

          exit

        else
          freqsq(i+1) = freqsq(i)
          uvec(:,i+1) = uvec(:,i)
          ideg(i+1) = ideg(i)
        endif
      enddo
      freqsq(insert) = tmp
      uvec(:,insert) = ctmp(:)
      ideg(insert) = itmp
    endif
  enddo

! starts grouping from representations

  call sym_cartesian_op(adot, rot, tau,                                  &
     ntrans, mtrx, tnp)

  allocate(ruvec(3*nat_qe,ntrans))

  irep(:) = 0
  nrep = 0

  do i = 1,3*nat_qe
    if(irep(i) == 0) then

!    new representation

      nrep = nrep + 1
      irep(i) = nrep

      call QE_rotate_mod(qvec, nat_qe, uvec(:,i), ruvec,                 &
            ntrans, rot, irt_qe, rtau_qe, invs_qe)

      do j = 1,3*nat_qe
        if(irep(j) == 0 .and. ideg(i) == ideg(j))then

!         could be the same representation

          do n = 1,ntrans
            prod = zdotc(3*nat_qe, ruvec(:,n), 1, uvec(:,j), 1)

            if(abs(prod) > eps) then
              irep(j) = nrep

              exit

            endif
          enddo

        endif
      enddo

    endif
  enddo

! Final sort by representation.  Again it may not be optimal but data is
! already somewhat ordered.

  do j = 2,3*nat_qe
    if(irep(j) - irep(j-1) < 0) then
      tmp = freqsq(j)
      itmp = irep(j)
      ctmp(:) = uvec(:,j)
      itdeg = ideg(j)
      insert = 1
      do i = j-1,1,-1
        insert = i
        if(irep(i) - itmp <= 0) then
          insert = i+1

          exit

        else
          freqsq(i+1) = freqsq(i)
          irep(i+1) = irep(i)
          uvec(:,i+1) = uvec(:,i)
          ideg(i+1) = ideg(i)
        endif
      enddo
      freqsq(insert) = tmp
      irep(insert) = itmp
      uvec(:,insert) = ctmp(:)
      ideg(insert) = itdeg
    endif
  enddo

  deallocate(ctmp)
  deallocate(ideg)

  deallocate(ruvec)

  return

end subroutine sym_phonon_rep

