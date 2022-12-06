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

!>  Creates the mapping of to an atom is transformed into
!>  by a symmetry operation, and the resspective coordinates.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         27 November 2022.
!>  \copyright    GNU Public License v2

subroutine sym_map_rat(adot, map_sr, map_stau,                           &
      ntrans, mtrx, tnp,                                                 &
      ntype, natom, rat,                                                 &
      mxdtyp, mxdatm)

! Adapted from other sym subroutines. 27 November 2022. JLM


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space


  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! output

  integer, intent(out)               ::  map_sr(mxdatm,mxdtyp,48)        !<  map_sr(i,j,k) indicates to which atom the atom i, of type j is transformed by symmetry k
  integer, intent(out)               ::  map_stau(3,mxdatm,mxdtyp,48)    !<  map_stau(:,i,j,k) indicates additional translation to bring rotated atom i to coincide with atom map_sr(i,j,k).  Lattice coordinates

! local variables

  real(REAL64)         ::  vcell, bdot(3,3)
  real(REAL64)         ::  rt(3)
  real(REAL64)         ::  distsq, dif(3)

  integer              ::  irotdir(3,3), idet

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-8_REAL64                        !  tolerance

! counters

  integer    ::  i, j, n, nt, j1, j2, k

  call adot_to_bdot(adot,vcell,bdot)

  do i = 1,3
  do j = 1,3
    bdot(i,j) = bdot(i,j)/ (2*PI*2*PI)
  enddo
  enddo

  do n = 1,ntrans

!   computes the minors of mtrx

    irotdir(1,1) = mtrx(2,2,n)*mtrx(3,3,n) - mtrx(3,2,n)*mtrx(2,3,n)
    irotdir(2,1) = mtrx(3,2,n)*mtrx(1,3,n) - mtrx(1,2,n)*mtrx(3,3,n)
    irotdir(3,1) = mtrx(1,2,n)*mtrx(2,3,n) - mtrx(2,2,n)*mtrx(1,3,n)
    irotdir(1,2) = mtrx(2,3,n)*mtrx(3,1,n) - mtrx(3,3,n)*mtrx(2,1,n)
    irotdir(2,2) = mtrx(3,3,n)*mtrx(1,1,n) - mtrx(1,3,n)*mtrx(3,1,n)
    irotdir(3,2) = mtrx(1,3,n)*mtrx(2,1,n) - mtrx(2,3,n)*mtrx(1,1,n)
    irotdir(1,3) = mtrx(2,1,n)*mtrx(3,2,n) - mtrx(3,1,n)*mtrx(2,2,n)
    irotdir(2,3) = mtrx(3,1,n)*mtrx(1,2,n) - mtrx(1,1,n)*mtrx(3,2,n)
    irotdir(3,3) = mtrx(1,1,n)*mtrx(2,2,n) - mtrx(2,1,n)*mtrx(1,2,n)

!   determinant

    idet = irotdir(1,1)*mtrx(1,1,n) + irotdir(2,1)*mtrx(2,1,n) +         &
           irotdir(3,1)*mtrx(3,1,n)

!   paranoid check

    if(iabs(idet) /= 1) then
      write(6,*) '    STOPPED in sym_map_rat'
      write(6,*) '    determinant = ',idet

      stop

    endif

!   inverse transpose

    do j = 1,3
    do i = 1,3
      irotdir(i,j) = irotdir(i,j) / idet
    enddo
    enddo

!   loop over types of atoms and atoms

    do nt = 1,ntype
      do j1 = 1,natom(nt)
        map_sr(j1,nt,n) = 0
        do k = 1,3
          map_stau(k,j1,nt,n) = 0
        enddo
        do k = 1,3
          rt(k) = irotdir(k,1)*rat(1,j1,nt) + irotdir(k,2)*rat(2,j1,nt) +  &
                  irotdir(k,3)*rat(3,j1,nt) + tnp(k,n) / (2*PI)
        enddo

!       identifies the rotated atom

        do j2 = 1,natom(nt)
          distsq = zero
          do k = 1,3
            dif(k) = rt(k) - rat(k,j2,nt)
            dif(k) = dif(k) - UM*nint(dif(k))
            distsq = distsq + dif(k)*dif(k)
          enddo
          if(distsq < EPS) then
            map_sr(j1,nt,n) = j2
            do k = 1,3
              map_stau(k,j1,nt,n) = nint(rat(k,j2,nt) - rt(k))
            enddo

            exit

          endif
        enddo

      enddo
    enddo

  enddo

  return

end subroutine sym_map_rat

