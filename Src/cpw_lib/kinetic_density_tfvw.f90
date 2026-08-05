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

!>  Calculates the kinetic energy density in the
!>  Thomas-Fermi-von Weizsaker approximation
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         17 May 2026.
!>  \copyright    GNU Public License v2

subroutine kinetic_density_tfvw(ipr, adot, den, tau_tfvw,                &
    ng, kgv, phase, conj, ns, inds, kmax, mstar,                         &
    mxdgve, mxdnst)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  ipr                             !<  contrlos printing
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  elecron density in prototype G-vector

! output

  complex(REAL64), intent(out)       ::  tau_tfvw(mxdnst)                !<  "kinetic energy density" in prototype G-vector

! local allocatable arrays

  real(REAL64), allocatable          ::  rhomsh(:,:,:)                   !  density on regular mesh in real space
  real(REAL64), allocatable          ::  taumsh(:,:,:)                   !  "kinetic energy density" on regular mesh in real space

! local variables

  integer                            ::  mxdfft                          !  array dimension for rhomsh
  integer                            ::  id,n1,n2,n3                     !  packing of rhomsh(id,n2,n3), id >= n1
  integer                            ::  nsfft(3), ntot
  integer                            ::  mxdwrk
  real(REAL64)                       ::  rho
  real(REAL64)                       ::  tauunif, tausingle
  real(REAL64)                       ::  d_tauunif_dr, d_tausingle_dr, d_tausingle_dgr
  real(REAL64)                       ::  grho

  integer, parameter  ::  mxdnn = 3                                      !  Lagrange interpolation uses 2*mxdnn+1 points
  real(REAL64)        ::  dgdm(-mxdnn:mxdnn), drhocon(3)
  integer             ::  if1, if2

  real(REAL64)        ::  vcell, bdot(3,3)


! parameters

  real(REAL64), parameter     :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer      ::  i1, i2, i3
  integer      ::  nn, in, jn
  integer      ::  i, J


! initial stuff

  call adot_to_bdot(adot,vcell,bdot)

! correct for the 2*PI

  do j = 1,3
  do i = 1,3
    bdot(i,j) = bdot(i,j) / (2*PI*2*PI)
  enddo
  enddo

  call size_fft(kmax, nsfft, mxdfft, mxdwrk)

  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)
  id = nsfft(1) + 1
  ntot = id * n2 * n3

  allocate(rhomsh(id,n2,n3))
  allocate(taumsh(id,n2,n3))

  rhomsh(:,:,:) = ZERO
  taumsh(:,:,:) = ZERO

! Weights for lagrange interpolation formula

  nn = min(n1,n2,n3) / 2
  nn = min(nn,mxdnn)
  do in = -nn,nn
    if1 = 1
    if2 = 1
    do jn = -nn,nn
      if (jn /= in .and. jn /= 0) if1 = if1 * (0  - jn)
      if (jn /= in)               if2 = if2 * (in - jn)
    enddo
    dgdm(in) = (if1*UM) / (if2*UM)
  enddo
  dgdm(0) = ZERO

! sets on real mesh

  call gvec_mesh_set(ipr, 'density-tau', adot, den,                      &
    rhomsh, id,n1,n2,n3, .TRUE.,                                         &
    ng, kgv, phase, conj, inds, kmax,                                    &
    mxdgve, mxdnst, ntot)

! calculates tau_tfvw

  do i3 = 1,n3
  do i2 = 1,n2
  do i1 = 1,n1

    call xc_cell_deriv(rhomsh, i1,i2,i3, id,n2, n1,n2,n3,                &
          nn, dgdm, bdot, rho, grho, drhocon,                            &
          mxdnn)

    call xc_tau(rho, grho, tauunif, tausingle,                           &
                  d_tauunif_dr, d_tausingle_dr, d_tausingle_dgr)
    taumsh(i1,i2,i3) = tauunif  + tausingle

  enddo
  enddo
  enddo

! back to stars

  tau_tfvw = C_ZERO

  call gvec_mesh_unset(ipr, 'tau-tfvw', adot, tau_tfvw,                  &
    taumsh, id,n1,n2,n3, .TRUE.,                                         &
    ng, kgv, phase, conj, ns, inds, kmax, mstar,                         &
    mxdgve, mxdnst, ntot)

  return

end subroutine kinetic_density_tfvw
