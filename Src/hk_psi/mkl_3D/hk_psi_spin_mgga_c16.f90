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
!>  of a meta-GGA times neig spin-wavevectors.
!>
!>    H = -1/2 nabla^2 + V_NL + V(r) - 1/2 nabla . T(r) nabla
!>
!>  where V and T are 2x2 matrices in spin space with components
!>  (1, sigma_z, sigma_x, sigma_y) given by vscr_sp and vtau_sp,
!>  vtau = d (rho eps_xc) / d tau.  For nsp = 2, vtau_sp(:,1) =
!>  (vtau_up + vtau_dn)/2 and vtau_sp(:,2) = (vtau_up - vtau_dn)/2.
!>  The vtau term is calculated with the reciprocal lattice components q_j
!>  of k+G as 1/2 sum_ij q_i FFT[ T sum_j bdot_ij FFT^-1[ q_j psi ] ].
!>
!>  This is the mkl 3D FFT spin version. There are special versions for other
!>  combinations of libraries and CPU/GPU
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         7 October 2026.
!>  \copyright    GNU Public License v2

subroutine hk_psi_spin_mgga_c16(mtxd, neig, psi_sp, hpsi_sp, lnewanl,    &
    ng, kgv, rkpt, adot,                                                 &
    ekpg, isort, vscr_sp, vtau_sp, kmscr, nsp,                           &
    anlsp, xnlkbsp, nanlsp,                                              &
    mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

! Written 7 October 2026 from the generic hk_psi_spin_mgga_c16 and the mkl hk_psi_spin_c16. JLM+claude

  use mkl_dfti

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdasp                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr_sp
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size, does not include spin)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  vscr_sp(mxdscr,mxdnsp)          !<  screened potential in the fft real space mesh (1, sigma_z, sigma_x, sigma_y)
  real(REAL64), intent(in)           ::  vtau_sp(mxdscr,mxdnsp)          !<  d (rho eps_xc) / d tau in the fft real space mesh (1, sigma_z, sigma_x, sigma_y)
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size
  integer, intent(in)                ::  nsp                             !<  number of spin components of the potentials (1,2,4)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanlsp                          !<  number of projectors
  complex(REAL64), intent(in)        ::  anlsp(2*mxddim,mxdasp)          !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkbsp(mxdasp)                 !<  Kleinman-Bylander normalization
  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim,mxdbnd)         !<  wavevector

! input and output

  logical, intent(inout)             ::  lnewanl                         !<  indicates that anlsp has been recalculated (not used in default implementation)

! output

  complex(REAL64), intent(out)       ::  hpsi_sp(2*mxddim,mxdbnd)        !<  |hpsi_sp> =  H |psi_sp>

! local allocatable arrays

  integer,allocatable          ::  ipoint(:)
  real(REAL64),allocatable     ::  qkpg(:,:)                             !  reciprocal lattice components of k+G
  complex(REAL64),allocatable  ::  chd_m(:)                              !  s = -1/2
  complex(REAL64),allocatable  ::  chd_p(:)                              !  s =  1/2
  complex(REAL64),allocatable  ::  phi_m(:,:)                            !  FFT^-1 of q_j psi, s = -1/2
  complex(REAL64),allocatable  ::  phi_p(:,:)                            !  FFT^-1 of q_j psi, s =  1/2

! local variables

  integer                ::  mxdfft                                      !  array dimension for fft transform
  integer                ::  mxdwrk                                      !  array dimension for fft transform workspace

  complex(REAL64)        ::  tmp_m, tmp_p
  complex(REAL64)        ::  gr_m(3), gr_p(3)                            !  sum_j bdot_ij phi_j

  integer         ::  kd1, kd2, kd3
  integer         ::  n1, n2, n3, id
  integer         ::  ntot,it
  integer         ::  nsfft(3)

  integer         ::  status
  integer         ::  strides_in(4)
  real(REAL64)    ::  fac

  real(REAL64)    ::  vcell, bdot(3,3)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)

! counters

  integer   ::   i, j, n, k1, k2, k3

  type(dfti_descriptor), pointer :: plan


! paranoid checks

  if(ng > mxdgve) stop

  if(nsp > mxdnsp) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_spin_mgga_c16.  nsp = ",i8,          &
           &  " is greater than mxdnsp = ",i8)') nsp, mxdnsp

    stop

  endif

  if(nsp /=1 .and. nsp /=2 .and. nsp /= 4) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_spin_mgga_c16.  nsp = ",i8,          &
           &  " is an incorrect value, should be 1,2,4")') nsp

    stop

  endif

! find n for fast fourier transform
! ni is the number of points used in direction i.

  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  if(mxdwrk > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_spin_mgga_c16.  mxdwrk = ",i8,       &
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
    write(6,*)  '   WARNING   in hk_psi_spin_mgga_c16:   mxdwrk may be', &
      '  incorrectly calculated ', lnewanl
    write(6,*)
  endif

  allocate(chd_m(mxdfft), chd_p(mxdfft))
  allocate(phi_m(mxdfft,3), phi_p(mxdfft,3))

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

    write(6,'("   STOPPED in hk_psi_spin_mgga_c16:   dimension of k ",   &
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

  do n = 1,neig

!   local potential

!$omp parallel do default(shared) private(i)
    do i = 1,ntot
       chd_p(i) = C_ZERO
       chd_m(i) = C_ZERO
    enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
    do i = 1,mtxd
      chd_p(ipoint(i)) = psi_sp(2*i-1,n)
      chd_m(ipoint(i)) = psi_sp(2*i  ,n)
    enddo
!$omp end parallel do

    Status = DftiComputeBackward(plan, chd_p)
    Status = DftiComputeBackward(plan, chd_m)

    if(nsp == 1) then

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd_p(i) = vscr_sp(i,1)*chd_p(i)
        chd_m(i) = vscr_sp(i,1)*chd_m(i)
      enddo
!$omp end parallel do

    elseif(nsp == 2) then

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd_p(i) = (vscr_sp(i,1)+vscr_sp(i,2))*chd_p(i)
        chd_m(i) = (vscr_sp(i,1)-vscr_sp(i,2))*chd_m(i)
      enddo
!$omp end parallel do

    else

!$omp parallel do default(shared) private(i, tmp_m, tmp_p)
      do i = 1,ntot
        tmp_p = ( vscr_sp(i,1) + vscr_sp(i,2) )*chd_p(i) +               &
                ( vscr_sp(i,3) - C_I*vscr_sp(i,4) )*chd_m(i)
        tmp_m = ( vscr_sp(i,1) - vscr_sp(i,2) )*chd_m(i) +               &
                ( vscr_sp(i,3) + C_I*vscr_sp(i,4) )*chd_p(i)
        chd_p(i) = tmp_p
        chd_m(i) = tmp_m
      enddo
!$omp end parallel do

    endif

    Status = DftiComputeForward(plan, chd_p)
    Status = DftiComputeForward(plan, chd_m)

!$omp parallel do default(shared) private(i)
    do i = 1,mtxd
      hpsi_sp(2*i-1,n) = chd_p(ipoint(i)) + ekpg(i)*psi_sp(2*i-1,n)
      hpsi_sp(2*i  ,n) = chd_m(ipoint(i)) + ekpg(i)*psi_sp(2*i  ,n)
    enddo
!$omp end parallel do

!   vtau term.  phi_j = FFT^-1[ q_j psi ]

    do j = 1,3

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        phi_p(i,j) = C_ZERO
        phi_m(i,j) = C_ZERO
      enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        phi_p(ipoint(i),j) = qkpg(j,i)*psi_sp(2*i-1,n)
        phi_m(ipoint(i),j) = qkpg(j,i)*psi_sp(2*i  ,n)
      enddo
!$omp end parallel do

      Status = DftiComputeBackward(plan, phi_p(:,j))
      Status = DftiComputeBackward(plan, phi_m(:,j))

    enddo

!   chi_i = T sum_j bdot_ij phi_j, stored in phi

    if(nsp == 1) then

!$omp parallel do default(shared) private(i, j, gr_p, gr_m)
      do i = 1,ntot
        do j = 1,3
          gr_p(j) = bdot(j,1)*phi_p(i,1) + bdot(j,2)*phi_p(i,2) + bdot(j,3)*phi_p(i,3)
          gr_m(j) = bdot(j,1)*phi_m(i,1) + bdot(j,2)*phi_m(i,2) + bdot(j,3)*phi_m(i,3)
        enddo
        do j = 1,3
          phi_p(i,j) = vtau_sp(i,1)*gr_p(j)
          phi_m(i,j) = vtau_sp(i,1)*gr_m(j)
        enddo
      enddo
!$omp end parallel do

    elseif(nsp == 2) then

!$omp parallel do default(shared) private(i, j, gr_p, gr_m)
      do i = 1,ntot
        do j = 1,3
          gr_p(j) = bdot(j,1)*phi_p(i,1) + bdot(j,2)*phi_p(i,2) + bdot(j,3)*phi_p(i,3)
          gr_m(j) = bdot(j,1)*phi_m(i,1) + bdot(j,2)*phi_m(i,2) + bdot(j,3)*phi_m(i,3)
        enddo
        do j = 1,3
          phi_p(i,j) = (vtau_sp(i,1)+vtau_sp(i,2))*gr_p(j)
          phi_m(i,j) = (vtau_sp(i,1)-vtau_sp(i,2))*gr_m(j)
        enddo
      enddo
!$omp end parallel do

    else

!$omp parallel do default(shared) private(i, j, gr_p, gr_m)
      do i = 1,ntot
        do j = 1,3
          gr_p(j) = bdot(j,1)*phi_p(i,1) + bdot(j,2)*phi_p(i,2) + bdot(j,3)*phi_p(i,3)
          gr_m(j) = bdot(j,1)*phi_m(i,1) + bdot(j,2)*phi_m(i,2) + bdot(j,3)*phi_m(i,3)
        enddo
        do j = 1,3
          phi_p(i,j) = ( vtau_sp(i,1) + vtau_sp(i,2) )*gr_p(j) +         &
                       ( vtau_sp(i,3) - C_I*vtau_sp(i,4) )*gr_m(j)
          phi_m(i,j) = ( vtau_sp(i,1) - vtau_sp(i,2) )*gr_m(j) +         &
                       ( vtau_sp(i,3) + C_I*vtau_sp(i,4) )*gr_p(j)
        enddo
      enddo
!$omp end parallel do

    endif

    do j = 1,3

      Status = DftiComputeForward(plan, phi_p(:,j))
      Status = DftiComputeForward(plan, phi_m(:,j))

!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        hpsi_sp(2*i-1,n) = hpsi_sp(2*i-1,n) + (UM/2)*qkpg(j,i)*phi_p(ipoint(i),j)
        hpsi_sp(2*i  ,n) = hpsi_sp(2*i  ,n) + (UM/2)*qkpg(j,i)*phi_m(ipoint(i),j)
      enddo
!$omp end parallel do

    enddo

!   end of loop over eigenvectors

  enddo

  Status = DftiFreeDescriptor(plan)

! non local potential

  call hk_psi_nl_c16(2*mtxd, neig, psi_sp, hpsi_sp, anlsp, xnlkbsp, nanlsp, .TRUE.,  &
      2*mxddim, mxdbnd, mxdasp)


  deallocate(chd_m, chd_p)
  deallocate(phi_m, phi_p)
  deallocate(ipoint)
  deallocate(qkpg)

  return

end subroutine hk_psi_spin_mgga_c16
