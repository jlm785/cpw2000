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

!>  Loop over the k-points to calculate the charge density
!>  and if required the "kinetic energy density".
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         12 May 2026. 1 October 2026.
!>  \copyright    GNU Public License v2


subroutine cpw_scf_loop_rho(ektot, ekl, lxctau, lxcmgga, tau,            &
      kmscr, dtau_dbdot, mxdscr,                                         &
      dims_, crys_, spaceg_, recip_, kpoint_, hamallk_, psiallk_,        &
      chdens_, filename_)

! extracted from cpw_scf (it was too long), 12 may 2026. JLM
! renamed occp to occ_x_wgk.  28 September 2026. JLM+claude
! tau_by_fft_stress. 1 October 2026. JLM+claude


  use cpw_variables

  implicit none

  type(dims_t)                       ::  dims_                           !<  array dimensions

  type(crys_t)                       ::  crys_                           !<  crystal structure

  type(spaceg_t)                     ::  spaceg_                         !<  space group information

  type(recip_t)                      ::  recip_                          !<  reciprocal space information

  type(kpoint_t)                     ::  kpoint_                         !<  k-point data

  type(hamallk_t)                    ::  hamallk_                        !<  hamiltonian size and indexation for all k-points

  type(psiallk_t)                    ::  psiallk_                        !<  psi for all k-points

  type(chdens_t)                     ::  chdens_                         !<  charge densities

  type(filename_t)                   ::  filename_                       !<  filenames

! input

  real(REAL64), intent(in)           ::  ekl(dims_%mxdnrk*dims_%mxdbnd)  !<  kinetic energy of wave-function j, for all the k-points
  logical, intent(in)                ::  lxctau                          !<  calculate tau
  logical, intent(in)                ::  lxcmgga                         !<  meta-GGA, also calculates dtau_dbdot (used in the stress)
  integer, intent(in)                ::  mxdscr                          !<  array dimension of dtau_dbdot (vscr)
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

! output

  real(REAL64), intent(out)          ::  ektot                           !<  electron kinetic energy (Hartree)
  complex(REAL64), intent(out)       ::  tau(dims_%mxdnst)               !<  "kinetic energy density"
  real(REAL64), intent(out)          ::  dtau_dbdot(3,3,mxdscr)          !<  d tau / d bdot on the real space mesh of the potential (hartree/bohr^3), only for meta-GGA

! allocatable arrays

  real(REAL64), allocatable          ::  occ_x_wgk(:)                    !  occupation*k-weight*spin deg. of eigenvector j
  complex(REAL64), allocatable       ::  denk(:)

  complex(REAL64), allocatable       ::  tauk(:)                         !  contribution of k-point to tau
  real(REAL64), allocatable          ::  dtau_k_dbdot(:,:,:)             !  contribution of k-point to dtau_dbdot

! local variables

  integer                ::  neig, mtxd
  real(REAL64)           ::  rkpt(3)
  integer                ::  irkpsi                                      !  used for saving to disk

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64

! counters

  integer       ::  irk, iel, i, j

  allocate(denk(dims_%mxdnst))
  allocate(occ_x_wgk(dims_%mxdbnd))

  if(lxctau) then
    allocate(tauk(dims_%mxdnst))
    if(lxcmgga) allocate(dtau_k_dbdot(3,3,mxdscr))
  endif

  if(lxctau .and. lxcmgga) dtau_dbdot(:,:,:) = ZERO

  ektot = zero
  iel = 0
  do irk = 1,kpoint_%nrk

!   reads psi if saved to disk

    if(filename_%itape_save_psi > 9) then
      read(filename_%itape_save_psi,rec = irk) psiallk_%psi_allk
      irkpsi = 1
    else
      irkpsi = irk
    endif

    rkpt(1) = kpoint_%rk(1,irk)
    rkpt(2) = kpoint_%rk(2,irk)
    rkpt(3) = kpoint_%rk(3,irk)

    neig = kpoint_%nband(irk)
    mtxd = hamallk_%mtxd_allk(irk)

!   adds to sum of occupied eigenvalues and kinetic energy

    do j = 1,neig
      iel = iel + 1
      occ_x_wgk(j) = 2*kpoint_%wgk(irk)*psiallk_%occ_allk(iel)
      ektot = ektot + occ_x_wgk(j)*ekl(iel)
    enddo

!   adds to total charge density


    call charge_by_fft(mtxd, neig, occ_x_wgk,                            &
        hamallk_%isort_allk(:,irk), psiallk_%psi_allk(:,:,irkpsi), denk, &
        recip_%ng, recip_%kgv, recip_%phase , recip_%conj, recip_%ns,    &
        recip_%inds, recip_%kmax, recip_%mstar,                          &
        dims_%mxddim, dims_%mxdbnd, dims_%mxdgve, dims_%mxdnst)

    do i=1,recip_%ns
      chdens_%den(i) = chdens_%den(i) + denk(i)
    enddo

!   adds to the kinetic energy density.  For meta-GGA also d tau / d bdot,
!   that requires the sum over the full Brillouin zone mesh.

    if(lxctau) then

      if(lxcmgga) then

        call tau_by_fft_stress(tauk, .TRUE., kmscr, dtau_k_dbdot,        &
            kpoint_%nx, kpoint_%ny, kpoint_%nz,                          &
            kpoint_%sx, kpoint_%sy, kpoint_%sz,                          &
            mtxd, neig, occ_x_wgk,                                       &
            hamallk_%isort_allk(:,irk), psiallk_%psi_allk(:,:,irkpsi),   &
            rkpt, irk, crys_%adot,                                       &
            spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                   &
            kpoint_%nrk, kpoint_%rk, kpoint_%wgk, kpoint_%kmap,          &
            recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns, &
            recip_%inds, recip_%kmax, recip_%mstar,                      &
            dims_%mxddim, dims_%mxdbnd, dims_%mxdnrk,                    &
            dims_%mxdgve, dims_%mxdnst, mxdscr)

        dtau_dbdot(:,:,:) = dtau_dbdot(:,:,:) + dtau_k_dbdot(:,:,:)

      else

        call tau_by_fft(tauk, mtxd, neig, occ_x_wgk,                     &
            hamallk_%isort_allk(:,irk), psiallk_%psi_allk(:,:,irkpsi),   &
            rkpt, crys_%adot,                                            &
            recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns, &
            recip_%inds, recip_%kmax, recip_%mstar,                      &
            dims_%mxddim, dims_%mxdbnd, dims_%mxdgve, dims_%mxdnst)

      endif

      do i = 1,recip_%ns
        tau(i) = tau(i) + tauk(i)
      enddo

    endif

  enddo

  deallocate(denk)
  deallocate(occ_x_wgk)

  if(lxctau) then
    deallocate(tauk)
    if(lxcmgga) deallocate(dtau_k_dbdot)
  endif

  return

end subroutine cpw_scf_loop_rho
