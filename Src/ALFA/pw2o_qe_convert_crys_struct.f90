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

!>  Converts the crystal structure in the cow2000 format to the
!>  format used in the pwscf.in file of QuantumEspresso
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         24 February 2025.
!>  \copyright    GNU Public License v2

subroutine pw2o_qe_convert_crys_struct(nat, iqe, rat_qe, avec_qe,        &
    aa, bb, cc, cosbc, cosac, cosab, ibravais, avec,                     &
    adot, ntype, natom, rat, alatt,                                      &
    mxdtyp, mxdatm, mxdnat)

! Extracted from the write_pwscf_in subroutine. 24 February 2025. JLM

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdnat                          !<  array dimension of rat_qe

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  real(REAL64), intent(in)           ::  alatt                           !<  lattice constant

! output:

  integer, intent(out)               ::  nat                             !<  total number of atoms
  integer, intent(out)               ::  iqe                             !<  type of Bravais lattice for QE
  real(REAL64), intent(out)          ::  rat_qe(3,mxdnat)                !<  atom positions for QE
  real(REAL64), intent(out)          ::  avec_qe(3,3)                    !<  lattice vectors for QE
  real(REAL64), intent(out)          ::  aa, bb, cc                      !<  lattice constants for QE
  real(REAL64), intent(out)          ::  cosbc, cosac, cosab             !<  cosines of lattice angles for QE
  integer, intent(out)               ::  ibravais                        !<  Bravais lattice ic cpw2000
  real(REAL64), intent(out)          ::  avec(3,3)                       !<  lattice vectors (with conventional orientation)

! allocatable array

  real(REAL64), allocatable          ::  rat_car(:,:,:)

! metric variables

  real(REAL64)                       ::  adotnig(3,3), adotsym(3,3)
  integer                            ::  mtotal(3,3)
  real(REAL64)                       ::  bvec(3,3)
  real(REAL64)                       ::  aconv(3,3), avecnig(3,3)

  real(REAL64)                       ::  bvec_qe(3,3)
  real(REAL64)                       ::  vcellqe

! local:

  integer               ::  ipr

! parameter

  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter ::  TOL = 1.0E-6_REAL64

! counters

  integer       ::  nt, i, j, k



! lattice vectors

  call adot_to_avec_sym(adot,avec,bvec)

  ipr = 0
  call metric_ident(adot, adotnig, adotsym, ibravais, mtotal,            &
                      avec, aconv, avecnig, TOL, ipr)

! Finds total number of atoms, converts to cartesian

  nat = 0
  do nt = 1,ntype
    nat = nat + natom(nt)
  enddo

  if(mxdnat < nat) then
    write(6,*)
    write(6,*) '   STOPPED in qe_convert_crystal_structure'
    write(6,*) '   Incorrect total number of atoms', mxdnat, nat
    write(6,*)

    stop

  endif

  allocate(rat_car(3,mxdatm,mxdtyp))

  do nt = 1,ntype
  do j = 1,natom(nt)
  do k = 1,3
    rat_car(k,j,nt) = avec(k,1)*rat(1,j,nt) + avec(k,2)*rat(2,j,nt)      &
                                            + avec(k,3)*rat(3,j,nt)
  enddo
  enddo
  enddo


! converts convention of bravais lattice

  iqe = 0
  aa = alatt
  if(ibravais == 1) then

    iqe = 1
    aa = aconv(1,1)

    avec_qe(1,1) = aa
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = ZERO

    avec_qe(1,2) = ZERO
    avec_qe(2,2) = aa
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = ZERO
    avec_qe(2,3) = ZERO
    avec_qe(3,3) = aa

  elseif(ibravais == 2) then

    iqe = 2
    aa = aconv(1,1)

    avec_qe(1,1) =-aa / 2
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = aa / 2

    avec_qe(1,2) = ZERO
    avec_qe(2,2) = aa / 2
    avec_qe(3,2) = aa / 2

    avec_qe(1,3) =-aa / 2
    avec_qe(2,3) = aa / 2
    avec_qe(3,3) = ZERO

  elseif(ibravais == 3) then

    iqe = 3
    aa = aconv(1,1)

    avec_qe(1,1) = aa / 2
    avec_qe(2,1) = aa / 2
    avec_qe(3,1) = aa / 2

    avec_qe(1,2) =-aa / 2
    avec_qe(2,2) = aa / 2
    avec_qe(3,2) = aa / 2

    avec_qe(1,3) =-aa / 2
    avec_qe(2,3) =-aa / 2
    avec_qe(3,3) = aa / 2

  elseif(ibravais == 4) then

    iqe = 6
    aa = aconv(1,1)
    cc = aconv(3,3)

    avec_qe(1,1) = aa
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = ZERO

    avec_qe(1,2) = ZERO
    avec_qe(2,2) = aa
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = ZERO
    avec_qe(2,3) = ZERO
    avec_qe(3,3) = cc

  elseif(ibravais == 5) then

    iqe = 7
    aa = aconv(1,1)
    cc = aconv(3,3)

    avec_qe(1,1) = aa / 2
    avec_qe(2,1) =-aa / 2
    avec_qe(3,1) = cc / 2

    avec_qe(1,2) = aa / 2
    avec_qe(2,2) = aa / 2
    avec_qe(3,2) = cc / 2

    avec_qe(1,3) =-aa / 2
    avec_qe(2,3) =-aa / 2
    avec_qe(3,3) = cc / 2

  elseif(ibravais == 6) then

    iqe = 8
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avec_qe(1,1) = aa
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = ZERO

    avec_qe(1,2) = ZERO
    avec_qe(2,2) = bb
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = ZERO
    avec_qe(2,3) = ZERO
    avec_qe(3,3) = cc

  elseif(ibravais == 7) then

    iqe = 9
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avec_qe(1,1) = aa / 2
    avec_qe(2,1) = bb / 2
    avec_qe(3,1) = ZERO

    avec_qe(1,2) =-aa / 2
    avec_qe(2,2) = bb / 2
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = ZERO
    avec_qe(2,3) = ZERO
    avec_qe(3,3) = cc

  elseif(ibravais == 8) then

    iqe = 10
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avec_qe(1,1) = aa / 2
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = cc / 2

    avec_qe(1,2) = aa / 2
    avec_qe(2,2) = bb / 2
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = ZERO
    avec_qe(2,3) = bb / 2
    avec_qe(3,3) = cc / 2

  elseif(ibravais == 9) then

    iqe = 11
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avec_qe(1,1) = aa / 2
    avec_qe(2,1) = bb / 2
    avec_qe(3,1) = cc / 2

    avec_qe(1,2) =-aa / 2
    avec_qe(2,2) = bb / 2
    avec_qe(3,2) = cc / 2

    avec_qe(1,3) =-aa / 2
    avec_qe(2,3) =-bb / 2
    avec_qe(3,3) = cc / 2

  elseif(ibravais == 10) then

    iqe = 4
    aa = aconv(1,1)
    cc = aconv(3,3)

    avec_qe(1,1) = aa
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = ZERO

    avec_qe(1,2) =-aa / 2
    avec_qe(2,2) = aa*sqrt(3*UM) / 2
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = ZERO
    avec_qe(2,3) = ZERO
    avec_qe(3,3) = cc

  elseif(ibravais == 11) then

    iqe = 5
    aa = sqrt((3*aconv(1,1)*aconv(1,1)+aconv(3,3)*aconv(3,3)) / (9*UM))
    cosbc = (2*aconv(3,3)*aconv(3,3)-3*aconv(1,1)*aconv(1,1)) / (18*aa*aa)

    avec_qe(1,1) = aa*sqrt(UM - cosbc)/sqrt(2*UM)
    avec_qe(2,1) =-aa*sqrt(UM - cosbc)/sqrt(6*UM)
    avec_qe(3,1) = aa*sqrt(UM + 2*cosbc)/sqrt(3*UM)

    avec_qe(1,2) = ZERO
    avec_qe(2,2) = 2*aa*sqrt(UM - cosbc)/sqrt(6*UM)
    avec_qe(3,2) = aa*sqrt(UM + 2*cosbc)/sqrt(3*UM)

    avec_qe(1,3) =-aa*sqrt(UM - cosbc)/sqrt(2*UM)
    avec_qe(2,3) =-aa*sqrt(UM - cosbc)/sqrt(6*UM)
    avec_qe(3,3) = aa*sqrt(UM + 2*cosbc)/sqrt(3*UM)

  elseif(ibravais == 12) then

    iqe = -12
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = sqrt(aconv(1,3)*aconv(1,3)+aconv(3,3)*aconv(3,3))
    cosac = aconv(1,3) / cc

    avec_qe(1,1) = aa
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = ZERO

    avec_qe(1,2) = ZERO
    avec_qe(2,2) = bb
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = aconv(1,3)
    avec_qe(2,3) = ZERO
    avec_qe(3,3) = aconv(3,3)

  elseif(ibravais == 13) then

    iqe = -13
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = sqrt(aconv(1,3)*aconv(1,3)+aconv(3,3)*aconv(3,3))
    cosac = aconv(1,3) / cc

    avec_qe(1,1) = aa / 2
    avec_qe(2,1) = bb / 2
    avec_qe(3,1) = ZERO

    avec_qe(1,2) =-aa / 2
    avec_qe(2,2) = bb / 2
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = aconv(1,3)
    avec_qe(2,3) = ZERO
    avec_qe(3,3) = aconv(3,3)

  elseif(ibravais == 14) then

    iqe = 14
    aa = aconv(1,1)
    bb = sqrt(aconv(1,2)*aconv(1,2)+aconv(2,2)*aconv(2,2))
    cc = sqrt(aconv(1,3)*aconv(1,3)+aconv(2,3)*aconv(2,3)+aconv(3,3)*aconv(3,3))
    cosbc = aconv(2,3) / cc
    cosac = aconv(1,3) / cc
    cosab = aconv(1,2) / bb

    avec_qe(1,1) = aconv(1,1)
    avec_qe(2,1) = ZERO
    avec_qe(3,1) = ZERO

    avec_qe(1,2) = aconv(1,2)
    avec_qe(2,2) = aconv(2,2)
    avec_qe(3,2) = ZERO

    avec_qe(1,3) = aconv(1,3)
    avec_qe(2,3) = aconv(2,3)
    avec_qe(3,3) = aconv(3,3)

  endif

! converts cartesian atomic positions to QE lattice coordinates

  bvec_qe(1,1) = avec_qe(2,2)*avec_qe(3,3) - avec_qe(3,2)*avec_qe(2,3)
  bvec_qe(2,1) = avec_qe(3,2)*avec_qe(1,3) - avec_qe(1,2)*avec_qe(3,3)
  bvec_qe(3,1) = avec_qe(1,2)*avec_qe(2,3) - avec_qe(2,2)*avec_qe(1,3)
  bvec_qe(1,2) = avec_qe(2,3)*avec_qe(3,1) - avec_qe(3,3)*avec_qe(2,1)
  bvec_qe(2,2) = avec_qe(3,3)*avec_qe(1,1) - avec_qe(1,3)*avec_qe(3,1)
  bvec_qe(3,2) = avec_qe(1,3)*avec_qe(2,1) - avec_qe(2,3)*avec_qe(1,1)
  bvec_qe(1,3) = avec_qe(2,1)*avec_qe(3,2) - avec_qe(3,1)*avec_qe(2,2)
  bvec_qe(2,3) = avec_qe(3,1)*avec_qe(1,2) - avec_qe(1,1)*avec_qe(3,2)
  bvec_qe(3,3) = avec_qe(1,1)*avec_qe(2,2) - avec_qe(2,1)*avec_qe(1,2)

  vcellqe = bvec_qe(1,1)*avec_qe(1,1) + bvec_qe(2,1)*avec_qe(2,1) +          &
            bvec_qe(3,1)*avec_qe(3,1)

  do j = 1,3
  do k = 1,3
    bvec_qe(k,j) = bvec_qe(k,j) / vcellqe
  enddo
  enddo

  i = 0
  do nt=1,ntype
  do j = 1,natom(nt)
    i  = i+1
    do k = 1,3
      rat_qe(k,i) = bvec_qe(1,k)*rat_car(1,j,nt) + bvec_qe(2,k)*rat_car(2,j,nt)      &
                                                + bvec_qe(3,k)*rat_car(3,j,nt)
    enddo
  enddo
  enddo

  deallocate(rat_car)

  return

end subroutine pw2o_qe_convert_crys_struct
