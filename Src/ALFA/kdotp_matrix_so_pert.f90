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

!>  Calculates the matrices for the k.p method or oscillator strength
!>  with spin-orbit. Perturbative version
!>
!>  \author       José Luís Martins
!>  \version      5.09
!>  \date         February 29, 2020. 30 November 2023.
!>  \copyright    GNU Public License v2

subroutine kdotp_matrix_so_pert(mtxd, neig, psi, ei, rkpt, isort, nder,  &
      hso0, dhso0drk, d2hso0drk2,                                        &
      ng, kgv,                                                           &
      ntype, natom, rat, adot,                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

! Written February 29, 2020, from previous code,
! spin_orbir_perturb and kdotp_matrix_so.  JLM
! Modified, nanlspin, indentation, 30 November 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalues (Hartree) for rk0

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction

  integer, intent(in)                ::  nder                            !<  order of derivatives to be calculated. kdotp nder =2, oscillator nder = 1.

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

! output

  complex(REAL64), intent(out)       ::  hso0(2*mxdbnd,2*mxdbnd)         !<  <Psi|H|Psi> without spin-orbit
  complex(REAL64), intent(out)       ::  dhso0drk(2*mxdbnd,2*mxdbnd,3)   !<  d <Psi|H|Psi> d k
  complex(REAL64), intent(out)       ::  d2hso0drk2(2*mxdbnd,2*mxdbnd,3,3)  !<  d^2 <Psi|H|Psi> d k^2

! local allocatable arrays

  complex(REAL64), allocatable       ::  pmat(:,:,:)                     ! <Psi|p_j|Psi>

  complex(REAL64), allocatable       ::  anlnoso(:,:)                    !  KB projectors without spin-orbit but spin representation
  real(REAL64), allocatable          ::  xnlkbnoso(:)                    !  KB normalization without spin-orbit but spin representation

  complex(REAL64), allocatable       ::  anlsop(:,:)                     !  KB projectors with spin-orbit
  complex(REAL64), allocatable       ::  anlsom(:,:)                     !  KB projectors with spin-orbit
  real(REAL64), allocatable          ::  xnlkbso(:)                      !  KB normalization without spin-orbit

  complex(REAL64), allocatable       ::  danlnosodrk(:,:,:)              !  d anlnoso / d rkpt
  complex(REAL64), allocatable       ::  d2anlnosodrk2(:,:,:,:)          !  d^2 anlnoso / d rkpt^2

  complex(REAL64), allocatable       ::  danlsopdrk(:,:,:)               !  d anlsop / d rkpt
  complex(REAL64), allocatable       ::  d2anlsopdrk2(:,:,:,:)           !  d^2 anlsop / d rkpt^2
  complex(REAL64), allocatable       ::  danlsomdrk(:,:,:)               !  d anlsom / d rkpt
  complex(REAL64), allocatable       ::  d2anlsomdrk2(:,:,:,:)           !  d^2 anlsom / d rkpt^2

  complex(REAL64), allocatable       ::  vnlso0(:,:)                     !  <Psi|V_NL|Psi>
  complex(REAL64), allocatable       ::  dvnlso0drk(:,:,:)               !  d <Psi|V_NL|Psi> d k
  complex(REAL64), allocatable       ::  d2vnlso0drk2(:,:,:,:)           !  d^2 <Psi|V_NL|Psi> d k^2

  complex(REAL64), allocatable       ::  vnlhalf(:,:)                    !  non-local component of hamiltonian

! local variables

  integer           ::  mxdanl         !  array dimension of number of projectors
  integer           ::  mxdaso         !  array dimension of number of projectors with spin-orbit
  integer           ::  nanl, nanlso, nanlnoso   !  number of KB projectors
  real(REAL64)      ::  vcell, bdot(3,3)

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer    ::  i, j, m, n

  if(nder < 0 .or. nder > 2) then
    write(6,*)
    write(6,*) '  Stopped in kdotp_matrix:  nder = ',nder
    write(6,*)

    stop

  elseif(nder == 0) then
    write(6,*)
    write(6,*) '  WARNING in kdotp_matrix. nder = 0',                    &
                   ' does not make sense'
    write(6,*)
  endif



  call adot_to_bdot(adot, vcell, bdot)

  allocate(pmat(3,mxdbnd,mxdbnd))

  call psi_p_psi(mtxd, neig, psi, pmat, rkpt, isort, ng, kgv, 1,         &
      mxddim, mxdbnd, mxdgve)

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlnoso,        &
      mxdtyp)

  mxdanl = nanl
  mxdaso = nanlso

  allocate(anlnoso(mxddim,mxdanl))
  allocate(xnlkbnoso(mxdanl))
  allocate(xnlkbso(mxdaso))

  allocate(danlnosodrk(mxddim,mxdanl,3))
  allocate(d2anlnosodrk2(mxddim,mxdanl,3,3))

  allocate(anlsop(mxddim,mxdaso))
  allocate(anlsom(mxddim,mxdaso))
  allocate(danlsopdrk(mxddim,mxdaso,3))
  allocate(danlsomdrk(mxddim,mxdaso,3))
  allocate(d2anlsopdrk2(mxddim,mxdaso,3,3))
  allocate(d2anlsomdrk2(mxddim,mxdaso,3,3))

  call proj_nl_kb_so_pert_der_c16(rkpt, mtxd, isort,                     &
      nanl, nanlso, nder,                                                &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlnoso, anlsop, anlsom, xnlkbnoso, xnlkbso,                       &
      danlnosodrk, danlsopdrk, danlsomdrk,                               &
      d2anlnosodrk2, d2anlsopdrk2, d2anlsomdrk2,                         &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdaso, mxdgve)

  allocate(vnlso0(2*mxdbnd,2*mxdbnd))
  allocate(dvnlso0drk(2*mxdbnd,2*mxdbnd,3))
  allocate(d2vnlso0drk2(2*mxdbnd,2*mxdbnd,3,3))

  allocate(vnlhalf(neig,neig))

! initialize

  do i=1,2*neig
  do j=1,2*neig
    hso0(j,i) = C_ZERO
  enddo
  enddo

! eigenvalues without spin-orbit

  do i=1,neig
    hso0(i,i) = cmplx(ei(i),ZERO,REAL64)
    hso0(i+neig,i+neig) = cmplx(ei(i),ZERO,REAL64)
  enddo

! add non-local spin-orbit pseudopotential

  call psi_vnl_sp_psi_der(mtxd, neig, nanlso, psi, 2,                    &
      vnlso0, dvnlso0drk, d2vnlso0drk2,                                  &
      anlsop, anlsom, xnlkbso, danlsopdrk, danlsomdrk,                   &
      d2anlsopdrk2, d2anlsomdrk2,                                        &
      mxddim, mxdbnd, mxdaso)

  do i=1,2*neig
  do j=1,2*neig
    hso0(j,i) = hso0(j,i)  + vnlso0(j,i)
  enddo
  enddo

! remove non-local pseudopotential without spin-orbit

  call psi_vnl_psi(mtxd, neig, psi, vnlhalf, anlnoso, xnlkbnoso,         &
      nanl, mxddim, neig, mxdanl)

  do i=1,neig
  do j=1,neig
    hso0(j,i) = hso0(j,i) - vnlhalf(j,i)
  enddo
  enddo
  do i=1,neig
  do j=1,neig
    hso0(j+neig,i+neig) = hso0(j+neig,i+neig) - vnlhalf(j,i)
  enddo
  enddo


  do i=1,2*neig
  do j=1,2*neig
    do m=1,3
    dhso0drk(j,i,m) = dvnlso0drk(j,i,m)
    enddo
  enddo
  enddo

  do i=1,neig
  do j=1,neig
    do m=1,3
    do n=1,3
      dhso0drk(j     ,i     ,m)  = dhso0drk(j     ,i     ,m) +           &
                                          bdot(m,n)*pmat(n,j,i)
      dhso0drk(j+neig,i+neig,m)  = dhso0drk(j+neig,i+neig,m) +           &
                                          bdot(m,n)*pmat(n,j,i)
    enddo
    enddo
  enddo
  enddo

  if(nder == 2) then

    do i=1,2*neig
    do j=1,2*neig
      do n=1,3
      do m=1,3
        d2hso0drk2(j,i,m,n) = d2vnlso0drk2(j,i,m,n)
      enddo
      enddo
    enddo
    enddo

    do j=1,2*neig
      do n=1,3
      do m=1,3
        d2hso0drk2(j,j,m,n) = d2hso0drk2(j,j,m,n) + bdot(m,n)
      enddo
      enddo
    enddo

  endif

  deallocate(pmat)

  deallocate(anlnoso)
  deallocate(xnlkbnoso)

  deallocate(anlsop)
  deallocate(anlsom)
  deallocate(xnlkbso)

  deallocate(danlnosodrk)
  deallocate(d2anlnosodrk2)

  deallocate(danlsopdrk)
  deallocate(danlsomdrk)
  deallocate(d2anlsopdrk2)
  deallocate(d2anlsomdrk2)

  deallocate(vnlso0)
  deallocate(dvnlso0drk)
  deallocate(d2vnlso0drk2)

  deallocate(vnlhalf)

  return

end subroutine kdotp_matrix_so_pert

