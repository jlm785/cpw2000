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

!>  This subroutine suggests a reference energy for a band structure plot
!>  based on the information from the k-point path.
!>  Another path may give a different value...
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         September 4, 201, 13 August 2025.
!>  \copyright    GNU Public License v2

subroutine out_band_eref(neig, nrk, rk, ztot, efermi, ispin, ivc,        &
             e_of_k, eref,nocc)

! Written September 4, 2014. JLM
! Modified, writes information, 7 November 2018. JLM
! Modified, documentation, August 1, 2019. JLM
! copyright  Jose Luis Martins/INESC-MN
! Modified, efermi, 29 November 2021. JLM
! Modified, extra information about local of CBM and VBM, 13 August 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  neig                            !<  number of eigenvectors (without spin)
  integer, intent(in)                ::  nrk                             !<  number of k-points in path

  real(REAL64), intent(in)           ::  rk(3,nrk)                            !<  k-point to be plotted in lattice coordinates
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree
  integer, intent(in)                ::  ispin                           !<  spin degeneracy (must be 1 or 2)
  integer, intent(in)                ::  ivc                             !<  1=CBM, 2=VBM, 3=midgap=(CBM+VBM)/2

  real(REAL64), intent(in)           ::  e_of_k(2*neig/ispin,nrk)        !<  band energies of k-point in plot

! output

  real(REAL64), intent(out)          ::  eref                            !<  reference energy for plot
  integer, intent(out)               ::  nocc                            !<  number of occupied bands (in semiconductors)

! allocatable arrays

  integer, allocatable               ::  indx(:)                         !  sorting index


! local variables

  real(REAL64)          ::  evbm           !  energy of valence band maximum
  real(REAL64)          ::  ecbm           !  energy of conduction band minimum
  logical               ::  lmetal         !  metallic case
  character(len=16)     ::  cso            !  for spin-orbit
  real(REAL64)          ::  rk_vbm(3)      !  k-point of valence band maximum
  real(REAL64)          ::  rk_cbm(3)      !  k-point of conduction band minimum
  real(REAL64)          ::  gam_vb, gam_cb !  eigenvalues at gamma


! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64
  real(REAL64), parameter  ::  EV = 27.21138505_REAL64
  real(REAL64), parameter  ::  EPS = 1.0E-6_REAL64

! counters

  integer    ::  j, n, k, irk


  if(ispin /=1 .and. ispin /= 2) then

!   inconsistent input, returns zero

    eref = ZERO
    nocc = 0

  else

    if(ispin == 2) then
      cso = ' no spin-orbit  '
    else
      cso = ' with spin-orbit'
    endif

    lmetal = .FALSE.
    if(ispin == 2 .and. mod(nint(ztot),2) == 1) lmetal = .TRUE.

    if(.not. lmetal) then

!     could be metal or insulator

      n = min(nint(ztot/ispin + 0.01),2*neig/ispin-1)

      if(n == 0) then

!       Exceptional case ztot = 0 or just one band...

        evbm = e_of_k(1,1)

        if(nrk > 1) then

          do irk = 2,nrk

            evbm = max(e_of_k(1,irk),evbm)

          enddo

        endif

        ecbm = evbm

      else

        allocate(indx(2*neig/ispin))

        call sort(2*neig/ispin, e_of_k(:,1), indx)

        evbm = e_of_k(indx(n),1)
        ecbm = e_of_k(indx(n+1),1)

        rk_vbm(:) = rk(:,1)
        rk_cbm(:) = rk(:,1)

        if(nrk > 1) then

          do irk = 2,nrk

            call sort(2*neig/ispin, e_of_k(1,irk), indx)

!            evbm = max(e_of_k(indx(n),irk),evbm)
!            ecbm = min(e_of_k(indx(n+1),irk),ecbm)

            if(e_of_k(indx(n),irk) > evbm) then
              evbm = e_of_k(indx(n),irk)
              rk_vbm(:) = rk(:,irk)
            endif

            if(e_of_k(indx(n+1),irk) < ecbm) then
              ecbm = e_of_k(indx(n+1),irk)
              rk_cbm(:) = rk(:,irk)
            endif

!           gets information for Gamma point

            if(abs(rk(1,irk)) < EPS .and. abs(rk(2,irk)) < EPS .and.     &
               abs(rk(3,irk)) < EPS) then
              gam_vb = e_of_k(indx(n),irk)
              gam_cb = e_of_k(indx(n+1),irk)
            endif

          enddo

        endif

        deallocate(indx)

      endif

      if(ecbm > evbm) then

        if(ivc == 1) then
          eref = evbm
        elseif(ivc == 2) then
          eref = ecbm
        else
          eref = (evbm+ecbm)/2
        endif
        nocc = n

        write(6,*)
        write(6,'(i8,"   Occupied bands in apparent ",                   &
           &     "semiconductor/insulator",a16)') nocc,cso
        write(6,'(f12.6,"   shift applied to bands (eV) ",a16)')         &
              -eref*EV,cso
        write(6,'(2f12.6,"   valence band maximum and conduction",       &
           &     " band minimum (eV) ",a16)') evbm*EV,ecbm*EV,cso
        write(6,'(3f8.3,8x,3f8.3,"   at these k-points")') rk_vbm(:), rk_cbm(:)
        write(6,*)
        write(6,'(2f12.6,"   energies at gamma for bands ",2i5 )')       &
            gam_vb*EV, gam_cb*EV, n, n+1
        write(6,*)

      else

        lmetal = .TRUE.

        write(6,*)
        write(6,'(i8,"   Average occupied bands in a metal ",            &
           &     "or semimetal",a16)') n,cso
        write(6,'(f12.6,"   Maximum energy in band ",i5)') evbm*EV, n
        write(6,'(f12.6,"   Minimum energy in band ",i5)') ecbm*EV, n+1
        write(6,*)

      endif
    endif

    if(lmetal) then

      n = min(nint((nrk*ztot)/ispin + 0.01),2*neig*nrk/ispin)

      allocate(indx(2*neig*nrk/ispin))

      call sort(2*neig*nrk/ispin, e_of_k, indx)

      k = (indx(n)-1)/(2*neig/ispin) + 1
      j = mod((indx(n)-1),(2*neig/ispin)) + 1
      eref = e_of_k(j,k)
      nocc = 0

!     uses efermi if not too different from the current estimate

      if(abs(eref - efermi) < 2/EV) then
        eref = efermi
        write(6,*)
        write(6,'(f12.6,"   E_F, from self-consistent calculation", a16)') eref*EV, cso
        write(6,'(f12.6,"   shift applied to bands (eV) ",a16)') -eref*EV, cso
        write(6,*)
      else
        write(6,*)
        write(6,'(f12.6,"   E_F, estimate of Fermi energy (eV) ", a16)') eref*EV, cso
        write(6,'(f12.6,"   shift applied to bands (eV) ",a16)') -eref*EV, cso
        write(6,*)
      endif

      deallocate(indx)

    endif

  endif

  return

  end subroutine out_band_eref
