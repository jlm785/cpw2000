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

!>  Adds the contribution of a non-local Kleinman-Bylander
!>  pseudopotential to the hamiltonian matrix
!>  Alternative method using proj_nl_kb
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         15 March 2024,
!>  \copyright    GNU Public License v2

subroutine hamilt_kb_non_local(rkpt, mtxd, isort, hamk,                  &
    ng, kgv,                                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdgve, mxdlqp, mxddim, mxdsml)

! Written 15 March 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdsml                          !<  array dimension of hamiltonian


  real(REAL64), intent(in)           ::  rkpt(3)                         !<  j-th component in lattice coordinates of the k-point
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms

  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

! input and output

  complex(REAL64), intent(inout)     ::  hamk(mxdsml,mxdsml)             !<  hamiltonian

! local allocatable arrays

  real(REAL64), allocatable           ::  xnlkb(:)
  complex(REAL64), allocatable        ::  anlga(:,:)                     !  KB projectors
  complex(REAL64), allocatable        ::  xanl(:,:)                      !  xnl*anl

! local variables

  integer           ::  nanl, nanlso, nanlspin
  integer           ::  mxdanl

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer   ::  i, j


  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
        mxdtyp)

  if(nanl > 0) then

    mxdanl = nanl

    allocate(xnlkb(mxdanl))
    allocate(anlga(mxddim,mxdanl))
    allocate(xanl(mxddim,mxdanl))

    call proj_nl_kb_c16(rkpt, mtxd, isort, nanl,                         &
        ng, kgv,                                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        anlga, xnlkb,                                                    &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

    do i = 1,nanl
      do j = 1,mtxd
        xanl(j,i) = anlga(j,i)*xnlkb(i)
      enddo
    enddo

    call zgemm('n', 'c', mtxd, mtxd, nanl, C_UM, xanl, mxddim,           &
                   anlga, mxddim, C_UM, hamk, mxdsml)

    deallocate(xnlkb)
    deallocate(anlga)
    deallocate(xanl)

  endif

  return

end subroutine hamilt_kb_non_local
