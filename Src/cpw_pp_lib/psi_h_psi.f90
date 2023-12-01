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

!>  Calculates <psi|H|psi> = <psi|v_loc|psi> + <psi|v_nl|psi> + <psi|0.5 nabla^2|psi>
!>
!>  \author       José Luís Martins
!>  \version      5.09
!>  \date         7 February 2014. 30 November 2023.
!>  \copyright    GNU Public License v2

subroutine psi_h_psi(rkpt, neig, psi, mtxd, isort, ekpg,                &
    hloc, hkin, hnl,                                                     &
    ng, kgv,                                                            &
    vscr, kmscr, nqnl, delqnl, vkb, nkb,                                    &
    ntype, natom, rat, adot,                                              &
    mxdtyp, mxdatm, mxdgve, mxddim, mxdbnd, mxdlqp, mxdscr)


! Written 7 February 2014. JLM
! modified, vkb dimensions, March 31, 2014. jlm
! Modified 8 November 2015. Compatibility new libpw. JLM
! Modified, documentation, 21 February 2020. JLM
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! Modified, nanlspin, 30 November 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (guess on input)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

! output

  complex(REAL64), intent(out)       ::  hloc(mxdbnd,mxdbnd)             !<  local component of reduced hamiltonian
  complex(REAL64), intent(out)       ::  hkin(mxdbnd,mxdbnd)             !<  kinetic component of reduced hamiltonian
  complex(REAL64), intent(out)       ::  hnl(mxdbnd,mxdbnd)              !<  non-local pseudopotential component of reduced hamiltonian

! allocatable arrays

  complex(REAL64), allocatable       ::  anlga(:,:)                      !  KB projectors without spin-orbit
  real(REAL64), allocatable          ::  xnlkb(:)                        !  KB normalization without spin-orbit

! local variables

  integer             ::  nanl, nanlso, nanlspin                         !  number of KB projectors
  integer             ::  mxdanl                                         !  array dimension of number of projectors


  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
       mxdtyp)

  mxdanl = nanl

  allocate(anlga(mxddim,mxdanl))
  allocate(xnlkb(mxdanl))

  call proj_nl_kb_c16(rkpt, mtxd, isort, nanl,                           &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb,                                                      &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

  call psi_loc_psi(mtxd, neig, psi, hloc, ng, kgv, isort,                &
      vscr, kmscr,                                                       &
      mxddim, mxdbnd, mxdgve, mxdscr)

  call psi_kin_psi(mtxd, neig, psi, hkin, ekpg,                          &
      mxddim, mxdbnd)

  call psi_vnl_psi(mtxd, neig, psi, hnl, anlga, xnlkb, nanl,             &
      mxddim, mxdbnd, mxdanl)

  deallocate(anlga)
  deallocate(xnlkb)

  return

end subroutine psi_h_psi
