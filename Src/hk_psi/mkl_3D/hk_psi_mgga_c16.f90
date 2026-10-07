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

!>  Calculates the product of the generalized Kohn-Sham hamiltonian
!>  of a meta-GGA times neig wavevectors.
!>
!>    H = -1/2 nabla^2 + V_NL + vscr(r) - 1/2 nabla . vtau(r) nabla
!>
!>  where vtau = d (rho eps_xc) / d tau.  The local potential and the
!>  vtau term are dealt with fast Fourier transforms.  For the vtau term
!>
!>    < k+G | -1/2 nabla . vtau nabla | k+G' > = 1/2 (k+G).(k+G') vtau(G-G')
!>
!>  is calculated with the reciprocal lattice components q_j of k+G
!>  as 1/2 sum_ij q_i FFT[ vtau sum_j bdot_ij FFT^-1[ q_j psi ] ].
!>
!>  This is the mkl 3D FFT version. There are generic and special versions for
!>  combinations of libraries and CPU/GPU
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         7 October 2026.
!>  \copyright    GNU Public License v2

subroutine hk_psi_mgga_c16(mtxd, neig, psi, hpsi, lnewanl,               &
    ng, kgv, rkpt, adot,                                                 &
    ekpg, isort, vscr, vtau, kmscr,                                      &
    anlga, xnlkb, nanl,                                                  &
    mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

! Written 7 October 2026 from the generic hk_psi_mgga_c16 and the mkl hk_psi_c16. JLM+claude

  use mkl_dfti

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  real(REAL64), intent(in)           ::  vtau(mxdscr)                    !<  d (rho eps_xc) / d tau in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanl                            !<  number of projectors
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  Kleinman-Bylander normalization
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector

! input and output

  logical, intent(inout)             ::  lnewanl                         !<  indicates that anlga has been recalculated (not used in default implementation)

! output

  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  |hpsi> =  H |psi>

! local allocatable arrays

  integer,allocatable          ::  ipoint(:)
  real(REAL64),allocatable     ::  qkpg(:,:)                             !  reciprocal lattice components of k+G
  complex(REAL64),allocatable  ::  chd(:)
  complex(REAL64),allocatable  ::  phi(:,:)                              !  FFT^-1 of q_j psi

! local variables

  integer                ::  mxdfft                                      !  array dimension for fft transform
  integer                ::  mxdwrk                                      !  array dimension for fft transform workspace

  integer         ::  kd1, kd2, kd3
  integer         ::  n1, n2, n3, id
  integer         ::  ntot,it
  integer         ::  nsfft(3)

  integer         ::  status
  integer         ::  strides_in(4)
  real(REAL64)    ::  fac

  real(REAL64)    ::  vcell, bdot(3,3)
  complex(REAL64) ::  ph1, ph2, ph3

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer   ::   i, j, n, k1, k2, k3

  type(dfti_descriptor), pointer :: plan


! find n for fast fourier transform
! ni is the number of points used in direction i.

  if(ng > mxdgve) stop

  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  if(mxdwrk > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_mgga_c16.  mxdwrk = ",i8,            &
           &  " is greater than mxdscr = ",2i8)') mxdwrk, mxdscr, ng

    stop

  endif

  n1 = kmscr(4)
  n2 = kmscr(5)
  n3 = kmscr(6)
  id = kmscr(7)
  ntot = id * n2 * n3

  if(ntot > mxdfft) then
    mxdfft = ntot
    write(6,*)
    write(6,*)  '   WARNING   in hk_psi_mgga_c16:   mxdwrk may be',      &
      '  incorrectly calculated ', lnewanl
    write(6,*)
  endif

  allocate(chd(mxdfft))
  allocate(phi(mxdfft,3))

  call adot_to_bdot(adot, vcell, bdot)

! fills array ipoint and the components of k+G

  allocate(ipoint(mtxd))
  allocate(qkpg(3,mtxd))

  kd1 = 0
  kd2 = 0
  kd3 = 0
  do i=1,mtxd
    it = isort(i)
    do j = 1,3
      qkpg(j,i) = rkpt(j) + kgv(j,it)
    enddo
    k1 = kgv(1,it)
    if (iabs(k1) > kd1) kd1 = iabs(k1)
    if (k1 < 0) k1 = n1 + k1
    k2 = kgv(2,it)
    if (iabs(k2) > kd2) kd2 = iabs(k2)
    if (k2 < 0) k2 = n2 + k2
    k3 = kgv(3,it)
    if (iabs(k3) > kd3) kd3 = iabs(k3)
    if (k3 < 0) k3 = n3 + k3
    ipoint(i) = (k3*n2 + k2)*id + k1 + 1
  enddo

! no wraparound for eigenvectors

  if(kd1 > (n1-1)/2 .or. kd2 > (n2-1)/2 .or. kd3 > (n3-1)/2) then

    write(6,'("   STOPPED in hk_psi_mgga_c16:   dimension of k ",        &
         &  "index= ",3i5," exceeds ",3i5)') kd1,kd2,kd3,                &
                     (n1-1)/2,(n2-1)/2,(n3-1)/2
    write(6,*) ' You are probably using the dual approximation'
    write(6,*) ' with a k-point far away from the Brillouin zone'

    stop

  endif

! prepare FFT

  strides_in(1) = 0
  strides_in(2) = 1
  strides_in(3) = id
  strides_in(4) = id*n2

  fac = UM / (n1*n2*n3)

  Status = DftiCreateDescriptor(plan, DFTI_DOUBLE, DFTI_COMPLEX, 3, nsfft)
  Status = DftiSetValue(plan, DFTI_INPUT_STRIDES, strides_in)
  Status = DftiSetValue(plan, DFTI_FORWARD_SCALE, fac)

  Status = DftiCommitDescriptor(plan)

! local potential, vtau term and kinetic energy

  do n=1,neig

!   local potential

!$omp parallel do default(shared) private(i)
    do i=1,ntot
       chd(i) = C_ZERO
    enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
    do i=1,mtxd
      chd(ipoint(i)) = psi(i,n)
    enddo
!$omp end parallel do

    Status = DftiComputeBackward(plan, chd)

!$omp parallel do default(shared) private(i)
    do i=1,ntot
      chd(i) = vscr(i)*chd(i)
    enddo
!$omp end parallel do

    Status = DftiComputeForward(plan, chd)

!$omp parallel do default(shared) private(i)
    do i=1,mtxd
      hpsi(i,n) = chd(ipoint(i)) + ekpg(i)*psi(i,n)
    enddo
!$omp end parallel do

!   vtau term.  phi_j = FFT^-1[ q_j psi ]

    do j = 1,3

!$omp parallel do default(shared) private(i)
      do i=1,ntot
        phi(i,j) = C_ZERO
      enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
      do i=1,mtxd
        phi(ipoint(i),j) = qkpg(j,i)*psi(i,n)
      enddo
!$omp end parallel do

      Status = DftiComputeBackward(plan, phi(:,j))

    enddo

!   chi_i = vtau sum_j bdot_ij phi_j, stored in phi

!$omp parallel do default(shared) private(i, ph1, ph2, ph3)
    do i=1,ntot
      ph1 = phi(i,1)
      ph2 = phi(i,2)
      ph3 = phi(i,3)
      phi(i,1) = vtau(i)*(bdot(1,1)*ph1 + bdot(1,2)*ph2 + bdot(1,3)*ph3)
      phi(i,2) = vtau(i)*(bdot(2,1)*ph1 + bdot(2,2)*ph2 + bdot(2,3)*ph3)
      phi(i,3) = vtau(i)*(bdot(3,1)*ph1 + bdot(3,2)*ph2 + bdot(3,3)*ph3)
    enddo
!$omp end parallel do

    do j = 1,3

      Status = DftiComputeForward(plan, phi(:,j))

!$omp parallel do default(shared) private(i)
      do i=1,mtxd
        hpsi(i,n) = hpsi(i,n) + (UM/2)*qkpg(j,i)*phi(ipoint(i),j)
      enddo
!$omp end parallel do

    enddo

!   end of loop over eigenvectors

  enddo

  Status = DftiFreeDescriptor(plan)

! non local potential

  call hk_psi_nl_c16(mtxd, neig, psi, hpsi, anlga, xnlkb, nanl, .TRUE.,  &
      mxddim,mxdbnd,mxdanl)


  deallocate(chd)
  deallocate(phi)
  deallocate(ipoint)
  deallocate(qkpg)

  return

end subroutine hk_psi_mgga_c16
