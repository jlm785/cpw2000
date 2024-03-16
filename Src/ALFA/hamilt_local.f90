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

!>  Adds the contribution of the local potential to the hamiltonian.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         15 March 2024,
!>  \copyright    GNU Public License v2

subroutine hamilt_local(mtxd, isort, qmod, hamk,                         &
    ng, kgv, phase, conj, inds, kmax, indv,                              &
    veff,                                                                &
    mxdgve, mxdnst, mxdcub, mxddim, mxdsml)

! Extracted from hamilt_kb.  15 March 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer , intent(in)               ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdsml                          !<  array dimension of hamiltonian


  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(in)           ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  real part of the phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of kgv(i,n)
  integer, intent(in)                ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax

  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  real part of the ionic potential (hartree) for the prototype g-vector in star j

! input and output

  complex(REAL64), intent(inout)     ::  hamk(mxdsml,mxdsml)             !<  hamiltonian

! local varaibles

  real(REAL64)      ::  qi,qj
  complex(REAL64)   ::  hloc

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64

! counters

  integer   ::  i, j
  integer   ::  im, iadd, jadd, kadd
  integer   ::  imx, jmx, kmx
  integer   ::  kxd, kyd, kzd


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  do i = 1,mtxd
    hamk(i,i) = hamk(i,i) + veff(1)
  enddo

! set up non diagonal part of matrix

  imx = 2*kmax(1) + 1
  jmx = 2*kmax(2) + 1
  kmx = 2*kmax(3) + 1

! start loop over columns for potential

  do i = 1,mtxd
    iadd = isort(i)
    qi = qmod(i)

!   start loop over rows (off diagonal elements)

    if (i > 1) then
      im = i-1
      do j = 1,im
        jadd = isort(j)
        qj = qmod(j)
        kxd = kgv(1,iadd)-kgv(1,jadd)+kmax(1)+1
        if (kxd >= 1 .and. kxd <= imx) then
        kyd = kgv(2,iadd)-kgv(2,jadd)+kmax(2)+1
        if (kyd >= 1 .and. kyd <= jmx) then
        kzd = kgv(3,iadd)-kgv(3,jadd)+kmax(3)+1
        if (kzd >= 1 .and. kzd <= kmx) then

!         calculate address for arrays vloc(r)(i) and s(r)(i)

          kadd = ((kxd-1)*jmx+kyd-1)*kmx+kzd
          kadd = indv(kadd)
          if (kadd /= 0) then

!           local part of potential

            hloc = veff(inds(kadd))*conjg(phase(kadd))
            if(conj(kadd) < ZERO) hloc = conjg(hloc)
            hamk(i,j) = hamk(i,j) + hloc
            hamk(j,i) = hamk(j,i) + conjg(hloc)

          endif
        endif
        endif
        endif

      enddo
    endif

  enddo

  return

end subroutine hamilt_local
