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

!>  Converts data for crystal structure from cpw2000 conventions
!>  to quantum espresso conventions.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         27 November 2022.
!>  \copyright    GNU Public License v2

subroutine sym_convto_qe(nat_qe, xau_qe, ityp_qe, irt_qe, rtau_qe, invs_qe,       &
      map_sr, map_stau,                                                  &
      ntrans, mtrx,                                                      &
      ntype, natom, rat,                                                 &
      mxdtyp, mxdatm)

! Adapted from other sym subroutines. 27 November 2022. JLM


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group


  integer, intent(in)                ::  map_sr(mxdatm,mxdtyp,48)        !<  map_sr(i,j,k) indicates to which atom the atom i, of type j is transformed by symmetry k
  integer, intent(in)                ::  map_stau(3,mxdatm,mxdtyp,48)    !<  map_stau(:,i,j,k) indicates additional translation to bring rotated atom i to coincide with atom map_sr(i,j,k).  Lattice coordinates.

! output

  integer, intent(out)               ::  nat_qe                          !<  total number of atoms
  real(REAL64), intent(out)          ::  xau_qe(3,mxdatm*mxdtyp)         !<  atomic positions lattice coordinates
  integer, intent(out)               ::  ityp_qe(mxdatm*mxdtyp)          !<  type of atom

  integer, intent(out)               ::  irt_qe(48,mxdatm*mxdtyp)        !<  index of the rotated atom
  real(REAL64), intent(out)          ::  rtau_qe(3,48,mxdatm*mxdtyp)     !<  rotated atom minus original, lattice coordinates.
  integer, intent(out)               ::  invs_qe(48)                     !<  index of the inverse operation

! local variables

  integer, allocatable     ::  imap(:,:)
  integer                  ::  iprod(3,3), isum

! counters

  integer    ::  i, j, n, nt, k, nq


! initialization

  allocate(imap(mxdatm,mxdtyp))

  ityp_qe(:) = 0
  irt_qe(:,:) = 0

! atom positions and type

  nat_qe = 0
  do nt = 1,ntype
    do i = 1,natom(nt)
      nat_qe = nat_qe + 1
      ityp_qe(nat_qe) = nt
      do k = 1,3
        xau_qe(k,nat_qe) = rat(k,i,nt)
        imap(i,nt) = nat_qe
      enddo
    enddo
  enddo

! rotated atoms

  nq = 0
  do nt = 1,ntype
    do i = 1,natom(nt)
      nq = nq + 1
      do n = 1,ntrans
        irt_qe(n,nq) = imap(map_sr(i,nt,n),nt)
        do k = 1,3
          rtau_qe(k,n,nq) = rat(k,map_sr(i,nt,n),nt) - rat(k,i,nt) - map_stau(k,i,nt,n)
        enddo
      enddo
    enddo
  enddo

  deallocate(imap)

! identify the inverse operation

  do n = 1,ntrans
    invs_qe(n) = 0

    do nq = 1,ntrans

    do i = 1,3
      do j = 1,3
        iprod(j,i) = 0
        do k = 1,3
          iprod(j,i) = iprod(j,i) + mtrx(j,k,nq)*mtrx(k,i,n)
        enddo
      enddo
      enddo
      isum = 0
      do i = 1,3
      do j = 1,3
        if(i == j) then
          isum = isum + iabs(iprod(i,i)-1)
        else
          isum = isum + iabs(iprod(i,j))
        endif
      enddo
      enddo
      if(isum == 0) then
        invs_qe(n) = nq

        exit

      endif

    enddo

  enddo

  return

end subroutine sym_convto_qe

