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
!>  Thomas-Fermi-von Weizsaker approximation, and its
!>  derivative with respect to the reciprocal space metric
!>  on the real space mesh of the potential.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         17 May 2026. 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine tau_tfvw_stress(ipr, adot, den, kmscr,                        &
    tau_tfvw, dtau_tfvw_dbdot,                                           &
    ng, kgv, phase, conj, ns, inds, kmax, mstar,                         &
    mxdgve, mxdnst, mxdscr)

! Written 17 May 2026. JLM
! Derivative with respect to bdot on the mesh of the potential.
! Removed debugging writes. 1 October 2026. JLM+claude
! renamed kinetic_density_tfvw to tau_tfvw_stress. 3 October 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdscr                          !<  array dimension of dtau_tfvw_dbdot (vscr)

  integer, intent(in)                ::  ipr                             !<  controls printing
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  electron density in prototype G-vector

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

! output

  complex(REAL64), intent(out)       ::  tau_tfvw(mxdnst)                !<  "kinetic energy density" in prototype G-vector
  real(REAL64), intent(out)          ::  dtau_tfvw_dbdot(3,3,mxdscr)     !<  (1/V) d (V tau_tfvw) / d bdot at fixed density in lattice coordinates, on the real space mesh of the potential (hartree/bohr^3).

! local allocatable arrays

  real(REAL64), allocatable          ::  rhomsh(:)                       !  density on regular mesh in real space, packed (id,n2,n3)
  real(REAL64), allocatable          ::  taumsh(:)                       !  "kinetic energy density" on regular mesh in real space
  complex(REAL64), allocatable       ::  taug(:)                         !  "kinetic energy density" for each G-vector

! local variables

  integer                            ::  id, n1, n2, n3                  !  packing of rhomsh(id,n2,n3), id >= n1
  integer                            ::  ntot
  real(REAL64)                       ::  rho
  real(REAL64)                       ::  tauunif, tausingle
  real(REAL64)                       ::  d_tauunif_dr, d_tausingle_dr, d_tausingle_dgr
  real(REAL64)                       ::  grho
  real(REAL64)                       ::  dcov(3)                         !  grho * adot * drhocon = d rho / d x_i (x lattice coordinates)
  real(REAL64)                       ::  fac

  integer, parameter  ::  mxdnn = 3                                      !  Lagrange interpolation uses 2*mxdnn+1 points
  real(REAL64)        ::  dgdm(-mxdnn:mxdnn), drhocon(3)
  integer             ::  if1, if2

  real(REAL64)        ::  vcell, bdot(3,3)                               !  cell volume, metric in reciprocal space
  real(REAL64)        ::  adotm1(3,3)                                    !  inverse of adot, bdot / (2 pi)^2

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer      ::  i1, i2, i3
  integer      ::  nn, in, jn
  integer      ::  i, j
  integer      ::  ind


  n1 = kmscr(4)
  n2 = kmscr(5)
  n3 = kmscr(6)
  id = kmscr(7)
  ntot = id * n2 * n3

  if(ntot > mxdscr) then
    write(6,*)
    write(6,*) '   STOPPED in tau_tfvw_stress:  mesh ', ntot,            &
               ' larger than mxdscr = ', mxdscr

    stop

  endif

  dtau_tfvw_dbdot(:,:,:) = ZERO

  call adot_to_bdot(adot, vcell, bdot)

! xc_cell_deriv uses the inverse of adot

  do j = 1,3
  do i = 1,3
    adotm1(i,j) = bdot(i,j) / (4*PI*PI)
  enddo
  enddo

  allocate(rhomsh(mxdscr))
  allocate(taumsh(mxdscr))

  rhomsh(:) = ZERO
  taumsh(:) = ZERO

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
    mxdgve, mxdnst, mxdscr)

! calculates tau_tfvw and its derivative with respect to bdot (with the 2 pi factors)
! for a variation of the metric at fixed density in lattice coordinates
! (fixed number of electrons in each volume element in lattice coordinates).
!
! As in tau_by_fft_stress, the derivative is of the quantity per cell,
!   dtau_dbdot = (1/V) d (V tau) / d bdot,
! because xc_cell adds separately the 1/V dependence of tau
! (the taumsh*adotm1 term).  With V = (2 pi)^3 / sqrt(det bdot),
!   d V / d bdot_ij = - (V/2) adot_ij / (2 pi)^2,
!   d rho / d bdot_ij = (rho/2) adot_ij / (2 pi)^2.
!
! von Weizsaker:  V tausingle = V grho^2 / (8 rho) depends on the metric only through
! grho^2 = sum_ij dcov_i bdot_ij dcov_j / (2 pi)^2  (dcov = d rho / d x_i, x lattice coordinates),
!   d tausingle / d bdot_ij = d_tausingle_dgr * dcov_i dcov_j / (2 grho (2 pi)^2),
! with dcov = grho * adot * drhocon.
!
! Thomas-Fermi:  V tauunif is proportional to V rho^(5/3), that is to V^(-2/3), so
!   (1/V) d (V tauunif) / d bdot_ij = (-2/3) (tauunif/V) d V / d bdot_ij
!                                   = (tauunif/3) adot_ij / (2 pi)^2.
! Equivalently, d tauunif / d rho * d rho / d bdot_ij = (5/6) tauunif adot_ij / (2 pi)^2
! minus the (1/2) tauunif adot_ij / (2 pi)^2 of the 1/V dependence.
!
! V tau is homogeneous of degree one in bdot, so
!   sum_ij bdot_ij dtau_tfvw_dbdot_ij = tauunif + tausingle.

  do i3 = 1,n3
  do i2 = 1,n2
  do i1 = 1,n1

    call xc_cell_deriv(rhomsh, i1,i2,i3, id,n2, n1,n2,n3,                &
          nn, dgdm, adotm1, rho, grho, drhocon,                          &
          mxdnn)

    call xc_tau(rho, grho, tauunif, tausingle,                           &
                  d_tauunif_dr, d_tausingle_dr, d_tausingle_dgr)
    ind = i1 + (i2-1)*id + (i3-1)*id*n2

    taumsh(ind) = tauunif  + tausingle

    do i = 1,3
      dcov(i) = grho * (adot(i,1)*drhocon(1) + adot(i,2)*drhocon(2) +    &
                        adot(i,3)*drhocon(3))
    enddo

    if(grho > ZERO) then
      fac = d_tausingle_dgr / (2*grho*2*PI*2*PI)
    else
      fac = ZERO
    endif

    do i = 1,3
    do j = 1,3
      dtau_tfvw_dbdot(j,i,ind) = fac * dcov(j)*dcov(i) +                 &
                        tauunif * adot(j,i) / (3*2*PI*2*PI)
    enddo
    enddo

  enddo
  enddo
  enddo

! back to stars for the kinetic energy density.
! gvec_mesh_unset only accepts the mesh of size_fft(kmax),
! so it is done in two steps for the (smaller) mesh of the potential.

  allocate(taug(mxdgve))

  call gvec_mesh_unset_nostar(ipr, 'tau-tfvw', adot, taug,               &
    taumsh, id,n1,n2,n3, .TRUE.,                                         &
    ng, kgv, kmax,                                                       &
    mxdgve, mxdscr)

  tau_tfvw(:) = C_ZERO

  call gvec_star_of_g_fold(tau_tfvw, taug, .FALSE.,                      &
    ng, phase, conj, ns, inds, mstar,                                    &
    mxdgve, mxdnst)

  deallocate(taug)
  deallocate(rhomsh)
  deallocate(taumsh)

  return

end subroutine tau_tfvw_stress
