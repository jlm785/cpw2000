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

subroutine sym_convfrom_qe(adot, nat_qe, uvec, vpol,                     &
      ntype, natom,                                                      &
      mxdtyp, mxdatm)

! Adapted from other sym subroutines. 27 November 2022. JLM
! Adapts some code from Quantum Espresso set_iss_new.f90


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  integer, intent(in)                ::  nat_qe                          !<  total number of atoms

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i

! input and output

  complex(REAL64), intent(inout)     ::  uvec(3*nat_qe,3*nat_qe)         !<  proper modes in cartesian coordinates.  Modified just by a phase.

! output

  complex(REAL64), intent(out)       ::  vpol(3,mxdatm,mxdtyp,3*nat_qe)  !<  proper modes in reciprocal coordiantes

! local variables

  real(REAL64)        ::  avec(3,3), bvec(3,3)
  real(REAL64)        ::  xmax
  complex(REAL64)     ::  fase
  integer             ::  jmax

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-9_REAL64

! counters

  integer    ::  i, j, nt, k, nq, imode


  call adot_to_avec_sym(adot,avec,bvec)

! We adjust the phase of each mode in such a way that the first
! non zero element is real (from QE)

  do imode = 1, 3*nat_qe

    xmax = ZERO
    jmax = 0
    do j = 1, 3*nat_qe
      if(abs(uvec(j, imode)) > xmax) then
        jmax = j
        xmax = abs(uvec(j, imode))
      endif
    enddo
    if(xmax  > EPS) then
      fase = uvec(jmax, imode) / xmax
      do j = 1, 3*nat_qe
        uvec(j, imode) = uvec(j, imode) * CONJG(fase)
       enddo
    endif

  enddo

! converts

  nq = 0
  do nt = 1,ntype
  do j = 1,natom(nt)
    nq = nq + 1
    do imode = 1,3*nat_qe
      do i = 1,3
        vpol(i,j,nt,imode) = C_ZERO
        do k = 1,3
          vpol(i,j,nt,imode) = vpol(i,j,nt,imode) + uvec(3*(nq-1)+k,imode)*bvec(k,i)
        enddo
        vpol(i,j,nt,imode) = vpol(i,j,nt,imode) / (2*PI)
      enddo
    enddo
  enddo
  enddo

! paranoid test

  if(nq /= nat_qe) then
    write(6,*)
    write(6,*) '   stopped in sym_convfrom_qe'
    write(6,*) '   nq, nat_qe = ', nq, nat_qe

    STOP

  endif

  return

end subroutine sym_convfrom_qe

