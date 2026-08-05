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

!>  Loop over the k-points to calculate the wavefunctions
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         12 May 2026.
!>  \copyright    GNU Public License v2

subroutine cpw_scf_loop_psi(iprglob, iter, minifail,                     &
      flgaopw,  iguess,  lkpg,                                           &
      kmscr, vscr, ekl,                                                  &
      dims_, crys_, flags_, pwexp_, recip_, acc_, strfac_,               &
      vcomp_, pseudo_, atorb_, kpoint_, hamallk_, psiallk_, filename_,   &
      mxdscr)

! extracted from cpw_scf (it was too long), 12 may 2026. JLM


  use cpw_variables

  implicit none


  type(dims_t)                       ::  dims_                           !<  array dimensions

  type(crys_t)                       ::  crys_                           !<  crystal structure

  type(flags_t)                      ::  flags_                          !<  computational flags

  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices

  type(recip_t)                      ::  recip_                          !<  reciprocal space information

  type(acc_t)                        ::  acc_                            !<  accuracy parameters

  type(strfac_t)                     ::  strfac_                         !<  structure factors

  type(vcomp_t)                      ::  vcomp_                          !<  Componemts of local potential

  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)

  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space

  type(kpoint_t)                     ::  kpoint_                         !<  k-point data

  type(hamallk_t)                    ::  hamallk_                        !<  hamiltonian size and indexation for all k-points

  type(psiallk_t)                    ::  psiallk_                        !<  psi for all k-points

  type(filename_t)                   ::  filename_                       !<  filenames

! input

  integer, intent(in)                ::  mxdscr                          !<  array dimension for vscr

  integer, intent(in)                ::  iter                            !<  iteration number

  character(len=2), intent(in)       ::  flgaopw                         !<  type of calculation 'AO' or 'PW
  integer, intent(in)                ::  iprglob                         !<  Global printing level

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the FFT real space mesh

  logical, intent(in)                ::  lkpg                            !<  If true use the previous G-vectors (same mtxd and isort)

! input and output

  integer, intent(inout)             ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential FFT mesh and fft mesh size

  integer, intent(inout)             ::  iguess                          !<  tells if guess eigenvectors are available

! output

  real(REAL64), intent(out)          ::  ekl(dims_%mxdnrk*dims_%mxdbnd)  !<  kinetic energy of wave-function j, for all the k-points
  integer, intent(out)               ::  minifail                        !<  minimum value of non-zero ifail, or 100 if all ifail is 0


! allocatable arrays

  complex(REAL64), allocatable       ::  hpsi(:,:)

  real(REAL64), allocatable          ::  hdiag(:)                        !  Hamiltonian diagonal for k+G-vector i
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+G-vector i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (Hartree) of k+G-vector  i

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue j
  real(REAL64), allocatable          ::  ekn(:)                          !  kinetic energy of wave-function j

! local variables

  integer                ::  neig, mtxd
  real(REAL64)           ::  rkpt(3)
  integer                ::  irkpsi                                      !  used for saving to disk

  integer                ::  ipr, nrka
  character(len=4)       ::  diag_type                                   !  selects diagonalization, 'pw  ','ao  ','aojc'
  integer                ::  nocc

  integer                ::  ifail                                       !  if ifail=0 the subroutine was successfull. Otherwise ifail indicates the number of correct digits.

  integer                ::  ickin
  character(len=5)       ::  labelk

  real(REAL64)           ::  veffr1

! counters

  integer       ::  irk, iel, j


  allocate(hdiag(dims_%mxddim))
  allocate(qmod(dims_%mxddim))
  allocate(ekpg(dims_%mxddim))

  allocate(hpsi(dims_%mxddim,dims_%mxdbnd))
  allocate(ei(dims_%mxdbnd))

  allocate(ekn(dims_%mxdbnd))

  minifail = 100
  iel = 0

  do irk = 1,kpoint_%nrk

!   loop over k-points

    rkpt(1) = kpoint_%rk(1,irk)
    rkpt(2) = kpoint_%rk(2,irk)
    rkpt(3) = kpoint_%rk(3,irk)
!
    neig = kpoint_%nband(irk)

!   calculates hamiltonian and diagonalizes
!   ***************************************

!   past first iteration one has a guess of the eigenvalues

    if(filename_%itape_save_psi > 9) then
      irkpsi = 1
    else
      irkpsi = irk
    endif

    if(iter > 1) then
      iguess = 1
      if(filename_%itape_save_psi > 9) then
        read(filename_%itape_save_psi,rec = irk) psiallk_%psi_allk
      endif
    endif


    mtxd = hamallk_%mtxd_allk(irk)

    if(flgaopw == 'PW') then

      ipr = 0
      if(iprglob > 2) ipr = 1

      diag_type = 'pw  '
      nocc = neig

      call h_kb_dia_all(diag_type, pwexp_%emax, rkpt, neig, nocc,        &
          flags_%flgpsd, ipr, ifail, acc_%icdiagmax,                     &
          iguess, acc_%epspsi,                                           &
          recip_%ng, recip_%kgv, recip_%phase, recip_%conj,              &
          recip_%ns, recip_%inds, recip_%kmax,                           &
          recip_%indv, recip_%ek,                                        &
          strfac_%sfact, vcomp_%veff, strfac_%icmplx,                    &
          pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,            &
          crys_%ntype,crys_%natom,crys_%rat,crys_%adot,                  &
          mtxd, hdiag, hamallk_%isort_allk(:,irk),                       &
          qmod, ekpg, lkpg,                                              &
          psiallk_%psi_allk(:,:,irkpsi), hpsi, ei,                       &
          vscr, kmscr,                                                   &
          atorb_%latorb, atorb_%norbat, atorb_%nqwf,                     &
          atorb_%delqwf, atorb_%wvfao, atorb_%lorb,                      &
          dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,        &
          dims_%mxdcub, dims_%mxdlqp, dims_%mxddim, dims_%mxdbnd,        &
          mxdscr, dims_%mxdlao)

      if(ifail /= 0) minifail = min(ifail,minifail)

      if(ifail < -3) then
        write(6,*)
        write(6,'("   Stopped in cpw_scf_loop_psi:  cycle is diverging", &
           &    " negative number of accuracy digits",i5)') ifail

        stop

      endif

    elseif(flgaopw == 'AO') then

      veffr1 = real(vcomp_%veff(1),REAL64)
      nocc = neig

      if(flags_%flgscf == 'AOJCPW') diag_type = 'aojc'
      if(flags_%flgscf == 'AOJC  ') diag_type = 'aojc'
      if(flags_%flgscf == 'AO    ') diag_type = 'ao  '

      call h_kb_dia_all(diag_type, pwexp_%emax, rkpt, neig, nocc,        &
          flags_%flgpsd, ipr, ifail, acc_%icdiagmax,                     &
          iguess, acc_%epspsi,                                           &
          recip_%ng, recip_%kgv, recip_%phase, recip_%conj,              &
          recip_%ns, recip_%inds, recip_%kmax,                           &
          recip_%indv, recip_%ek,                                        &
          strfac_%sfact, vcomp_%veff, strfac_%icmplx,                    &
          pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,            &
          crys_%ntype,crys_%natom,crys_%rat,crys_%adot,                  &
          mtxd, hdiag, hamallk_%isort_allk(:,irk),                       &
          qmod, ekpg, lkpg,                                              &
          psiallk_%psi_allk(:,:,irkpsi), hpsi, ei,                       &
          vscr, kmscr,                                                   &
          atorb_%latorb, atorb_%norbat, atorb_%nqwf,                     &
          atorb_%delqwf, atorb_%wvfao, atorb_%lorb,                      &
          dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,        &
          dims_%mxdcub, dims_%mxdlqp, dims_%mxddim, dims_%mxdbnd,        &
          mxdscr, dims_%mxdlao)

    else

      write(6,*)
      write(6,'("     STOPPED in scf_kb_loop_psi_c16:   unknown type",   &
         &      " of basis   ",a2)') flgaopw

      stop

    endif


    hamallk_%mtxd_allk(irk) = mtxd
    kpoint_%nband(irk) = neig


!   end of diagonalization


!   calculates the kinetic energy

    call kinetic_energy(neig, mtxd, ekpg, psiallk_%psi_allk(:,:,irkpsi), ekn,  &
        dims_%mxddim,dims_%mxdbnd)

!   prints the eigensolutions

    ipr = 0
    if(iprglob == 3) ipr = 1
    if(iprglob == 4) ipr = 2

    nrka = -1

    ickin = 1

    call print_eig(ipr, irk, labelk, nrka, rkpt,                         &
        mtxd, ickin, neig, psiallk_%psi_allk(:,:,irkpsi),                &
        crys_%adot, ei, ekn, hamallk_%isort_allk(:,irk), recip_%kgv,     &
        dims_%mxddim, dims_%mxdbnd, dims_%mxdgve)

!   stores eigenvalues and kinetic energies

    do j=1,neig
      iel = iel + 1
      psiallk_%eig_allk(iel) = ei(j)
      ekl(iel) = ekn(j)
    enddo

!   stores the wavefunctions

    if(filename_%itape_save_psi > 9) then
      write(filename_%itape_save_psi,rec = irk) psiallk_%psi_allk
    endif

!   end of loop over k-points

  enddo

  deallocate(hpsi)

  deallocate(ei)

  deallocate(hdiag)
  deallocate(qmod)
  deallocate(ekpg)

  return

end subroutine cpw_scf_loop_psi

