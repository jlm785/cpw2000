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

!>  Given a q-vector and a crystal structure finds a basis
!>  for phonon modes according to the irreducible representations
!>  of the little group.  Based on Quantum Espresso.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         6 December 2022.
!>  \copyright    GNU Public License v2

subroutine sym_phonon_mode(qvec, vpol, nrep, irep, iseed,                &
       adot, ntype, natom, rat,                                          &
       ntrans, mtrx, tnp,                                                &
       ntrans_q, mtrx_q, tnp_q,                                          &
       mxdtyp, mxdatm, mxdnat)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdnat                          !<  array dimension of total number of atoms

  real(REAL64), intent(in)           ::  qvec(3)                         !<  wave-vector in lattice coordinates

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer,  intent(in)               ::  iseed                           !<  seed for the random numbers.

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! output

  complex(REAL64), intent(out)       ::  vpol(3,mxdatm,mxdtyp,mxdnat)    !<  proper modes in lattice coordinates.

  integer, intent(out)               ::  nrep                            !<  number of representations present
  integer, intent(out)               ::  irep(mxdnat)                    !<  representation of each mode

  integer, intent(out)               ::  ntrans_q                        !<  number of symmetry operations in the little group of q
  integer, intent(out)               ::  mtrx_q(3,3,48)                  !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the little group of q
  real(REAL64), intent(out)          ::  tnp_q(3,48)                     !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the little group of q

! local allocatable arrays

  integer, allocatable               ::  map_sr(:,:,:)                   !  map_sr(i,j,k) indicates to which atom the atom i, of type j is transformed by symmetry k
  integer, allocatable               ::  map_stau(:,:,:,:)               !  map_stau(:,i,j,k) indicates additional translation to bring rotated atom i to coincide with atom map_sr(i,j,k).  Lattice coordinates.

  real(REAL64), allocatable          ::  xau_qe(:,:)                     !  atomic positions
  integer, allocatable               ::  ityp_qe(:)                      !  type of atom

  integer, allocatable               ::  irt_qe(:,:)                     !  index of the rotated atom
  real(REAL64), allocatable          ::  rtau_qe(:,:,:)                  !  rotated atom minus original, lattice coordinates.
  integer                            ::  invs_qe(48)                     !  inverse of rotation

  complex(REAL64), allocatable       ::  wdyn_qe (:,:,:,:)               !  "dynamical" matrix
  complex(REAL64), allocatable       ::  phi_qe (:,:)                    !  "dynamical" matrix

  real(REAL64), allocatable          ::  freqsq(:)                       !  eigenvalue of "dynamical" matrix.
  complex(REAL64), allocatable       ::  uvec(:,:)                       !  eigenvector of "dynamical" matrix.

! local variables

  integer               ::  nat_qe                                       !  total number of atoms

  integer               ::  irotm_q                                      !  the rotation sending q -> -q
  logical               ::  lgam_q                                       !  if TRUE  q=0
  logical               ::  lminus_q                                     !  if TRUE there is an inversion symmetry
  integer               ::  t_rev_q(48)                                  !  1 if symmetry includes time reversal.  Here only 0 is used.

  real(REAL64)          ::  rot(3,3,48)                                  !  rotational matrix in cartesian coordinates
  real(REAL64)          ::  tau(3,48)                                    !  fractional translation vector in cartesian coordinates

  real(REAL64)          ::  avec(3,3), bvec(3,3)

  integer               ::  info, ipr, istatus

! constants

  real(REAL64), parameter  ::  TOL = 1.0E-7_REAL64                        !  tolerance

! counters

  integer             ::  n, j


! finds the little group of q

  call sym_little_group_q(qvec,                                          &
    ntrans_q, mtrx_q, tnp_q, lminus_q, irotm_q, lgam_q,                  &
    ntrans, mtrx, tnp)

! get cartesian components

  call sym_cartesian_op(adot, rot, tau,                                  &
      ntrans_q, mtrx_q, tnp_q)

  call sym_test(ipr, TOL, istatus,                                       &
       ntrans_q, mtrx_q, tnp_q,                                          &
       ntype, natom, rat, adot,                                          &
       mxdtyp,mxdatm)

  allocate(map_sr(mxdatm,mxdtyp,48))
  allocate(map_stau(3,mxdatm,mxdtyp,48))

  call sym_map_rat(adot, map_sr, map_stau,                               &
       ntrans_q, mtrx_q, tnp_q,                                          &
       ntype, natom, rat,                                                &
       mxdtyp, mxdatm)

  allocate(xau_qe(3,mxdatm*mxdtyp))
  allocate(ityp_qe(mxdatm*mxdtyp))
  allocate(irt_qe(48,mxdatm*mxdtyp))
  allocate(rtau_qe(3,48,mxdatm*mxdtyp))

  call sym_convto_qe(nat_qe, xau_qe, ityp_qe, irt_qe, rtau_qe, invs_qe,  &
       map_sr, map_stau,                                                 &
       ntrans_q, mtrx_q,                                                 &
       ntype, natom, rat,                                                &
       mxdtyp, mxdatm)

  allocate(wdyn_qe (3, 3, nat_qe, nat_qe))

  call QE_random_matrix_new (nat_qe, irt_qe, ntrans_q, lminus_q, irotm_q,      &
       wdyn_qe, lgam_q, iseed)


  t_rev_q = 0

  call QE_symdynph_gq_new( qvec, wdyn_qe, mtrx_q, invs_qe, rtau_qe, irt_qe, ntrans_q, &
       nat_qe, irotm_q, lminus_q, t_rev_q)

  call adot_to_avec_sym(adot,avec,bvec)

  do n = 1, nat_qe
  do j = 1, nat_qe
      call QE_trntnsc( wdyn_qe(1,1,n,j), avec, bvec, 1 )
  enddo
  enddo

  allocate(phi_qe (3*nat_qe, 3*nat_qe))

  call QE_compact_dyn(nat_qe, phi_qe, wdyn_qe)

  allocate(freqsq(3*nat_qe))
  allocate(uvec(3*nat_qe,3*nat_qe))

  call diag_c16(3*nat_qe, phi_qe, freqsq, uvec, 3*nat_qe, info)

! group by representation

  call sym_phonon_rep(adot, qvec, nat_qe, freqsq, uvec, nrep, irep,      &
                irt_qe, rtau_qe, invs_qe,                                &
                ntrans_q, mtrx_q, tnp_q)

  call sym_convfrom_qe(adot, nat_qe, uvec, vpol,                         &
            ntype, natom,                                                &
            mxdtyp, mxdatm)

  deallocate(map_sr)
  deallocate(map_stau)

  deallocate(xau_qe)
  deallocate(ityp_qe)
  deallocate(irt_qe)
  deallocate(rtau_qe)

  deallocate(wdyn_qe)
  deallocate(phi_qe)

  deallocate(freqsq)
  deallocate(uvec)

  return

end subroutine sym_phonon_mode

