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

!>  calculates the spin-orbit perturbation hamiltonian
!>  for a group of eigen-vectors
!>
!>  \author       José Luís Martins
!>  \version      5.09
!>  \date         4 February 2014. 30 November 2023.
!>  \copyright    GNU Public License v2

subroutine spin_orbit_perturb(rkpt, mtxd, isort,                         &
    neig, psi, ei, ei_so, psi_so, lpsiso,                                &
    ng, kgv,                                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)


! Written 4 February 2014.  JLM
! modified, dimensions vkb, March 31, 2014. jlm
! modified, less memory, faster, 13 January 2020. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction
  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  neig                            !<  number of eigenvectors required (maybe modified on output)
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (guess on input)
  real(REAL64), intent(in)           ::  ei(mxddim)                      !<  eigenvalue no. i. (hartree)

  logical, intent(in)                ::  lpsiso                          !<  If true calculates the spin-orbit perturbed wave-functions

! output

  real(REAL64), intent(out)          ::  ei_so(2*mxdbnd)                 !<  spin-orbit eigenvalues
  complex(REAL64), intent(out)       ::  psi_so(2*mxddim,2*mxdbnd)       !<  vectors in the spin-orbit form

! main local variables

  integer       ::  mxdanl                                               !  array dimension of number of projectors
  integer       ::  mxdaso                                               !  array dimension of number of projectors with spin-orbit
  integer       ::  nanl                                                 !  half of number of projectors without spin
  integer       ::  nanlso, nanlspin                                     !  number of projectors witn spin

  integer            ::  info

! allocatable local arrays

  complex(REAL64), allocatable       ::  anlspin(:,:)               !  KB projectors without spin-orbit
  real(REAL64), allocatable          ::  xnlkbspin(:)               !  KB normalization without spin-orbit

  complex(REAL64), allocatable       ::  anlsop(:,:)                !  KB projectors with spin-orbit
  complex(REAL64), allocatable       ::  anlsom(:,:)                !  KB projectors with spin-orbit
  real(REAL64), allocatable          ::  xnlkbso(:)                 !  KB normalization with spin-orbit

  complex(REAL64), allocatable       ::  h_so(:,:)                  !  spin-orbit hamiltonian
  complex(REAL64), allocatable       ::  hnl(:,:)                   !  non-local component of hamiltonian

  complex(REAL64), allocatable       ::  psi_in(:,:)                !  input vectors in the spin-orbit form
  complex(REAL64), allocatable       ::  vec_so(:,:)                !  spin-orbit vectors in the psi_so basis set

  complex(REAL64), allocatable       ::  hnlhalf(:,:)               !  non-local component of hamiltonian

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer         ::  i, j

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
       mxdtyp)

  mxdanl = nanl
  mxdaso = nanlso

  allocate(anlspin(mxddim,mxdanl))
  allocate(xnlkbspin(mxdanl))
  allocate(xnlkbso(mxdaso))

  allocate(anlsop(mxddim,mxdaso))
  allocate(anlsom(mxddim,mxdaso))

  call proj_nl_kb_so_pert_c16(rkpt, mtxd, isort, nanl, nanlso,           &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlspin, anlsop, anlsom, xnlkbspin, xnlkbso,                       &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdaso, mxdgve)

  allocate(h_so(2*neig,2*neig))
  allocate(hnl(2*neig,2*neig))

  allocate(hnlhalf(neig,neig))

  do i=1,2*neig
  do j=1,2*neig
    h_so(j,i) = C_ZERO
  enddo
  enddo

  do i=1,neig
    h_so(i,i) = cmplx(ei(i),ZERO,REAL64)
    h_so(i+neig,i+neig) = cmplx(ei(i),ZERO,REAL64)
  enddo

  call psi_vnl_sp_psi(mtxd, neig, nanlso, psi, hnl,                      &
      anlsop, anlsom, xnlkbso,                                           &
      mxddim, neig, mxdaso)

  do i=1,2*neig
  do j=1,2*neig
    h_so(j,i) = h_so(j,i)  + hnl(j,i)
  enddo
  enddo

  call psi_vnl_psi(mtxd, neig, psi, hnlhalf, anlspin, xnlkbspin,         &
      nanl, mxddim, neig, mxdanl)

  do i=1,neig
  do j=1,neig
    h_so(j,i) = h_so(j,i) - hnlhalf(j,i)
  enddo
  enddo
  do i=1,neig
  do j=1,neig
    h_so(j+neig,i+neig) = h_so(j+neig,i+neig) - hnlhalf(j,i)
  enddo
  enddo

  allocate(vec_so(2*neig,2*neig))

  call diag_c16(2*neig, h_so, ei_so, vec_so, 2*neig, info)


  if(info /= 0) stop


  if(lpsiso) then

    allocate(psi_in(2*mxddim,2*neig))

    do i=1,neig
    do j=1,mtxd
      psi_in(2*j-1,i     ) = psi(j,i)
      psi_in(2*j  ,i     ) = C_ZERO
      psi_in(2*j-1,i+neig) = C_ZERO
      psi_in(2*j  ,i+neig) = psi(j,i)
    enddo
    enddo

    call zgemm('n', 'n', 2*mtxd, 2*neig, 2*neig, C_UM, psi_in,           &
           2*mxddim, vec_so, 2*neig, C_ZERO, psi_so, 2*mxddim)

    deallocate(psi_in)

  endif

  deallocate(anlspin)
  deallocate(xnlkbspin)
  deallocate(xnlkbso)

  deallocate(anlsop)
  deallocate(anlsom)

  deallocate(h_so)
  deallocate(vec_so)
  deallocate(hnl)

  deallocate(hnlhalf)

  return

end subroutine spin_orbit_perturb
