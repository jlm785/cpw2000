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

!>  Computes the exchange correlation potential and energy (Hartree)
!>  for the Tran-Blaha meta-GGA
!>  from the charge density (in electrons per unit cell), density laplacian
!>  and "kinetic energy" density.
!>  Adapted from older LDA code.
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.13
!>  \date         23 February 1999, 1 October 2026.
!>  \copyright    GNU Public License v2

subroutine xc_cell(author, adot, tblaha, lkincalc, id1,id2, n1,n2,n3,    &
        rhomsh, taumsh, dtau_dbdot, rholapmsh,                           &
        exc, vxc, rhovxc, strxc)

! Written 23 February 1999. jlm
! Modified for MMGA Tran-Blaha. CLR
! Modified for f90. 18 December 2016.  JLM
! Modified, documentation, December 2019. JLM
! sigma, twotau, 29 September 2022. JLM
! indentation, avois noise with Tran-Blaha and slabs, 6 October 2025. JLM
! Changed name of old xc_mgga to xc_mgga_vxc in preparaation for new functionals. 21 November 2025. JLM
! Kinetic energy density not twice (twotau -> tau). 25 November 2025. JLM
! Thomas-Fermi-von Weizsaker for tau. 12 May 2026. JLM
! replace dble, 18 August 2026. JLM
! d_taumsh_dgij(3,3,mesh) as in tau_by_fft_stress. 29 September 2026. JLM+claude
! stress contribution of the tau correction taumsh and d_taumsh_dgij. 29 September 2026. JLM+claude
! renamed d_taumsh_dgij to dtau_dbdot. 1 October 2026. JLM+claude

! WARNING choice of correlation for Tran-Blaha is hard coded as Perdew-Zunger
! WARNING correction for slab for Tran-Blaha are hard coded.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len = *), intent(in)     ::  author                          !<  type of xc wanted (ca=pz , pw92 , vwn, wi, pbe)
  real(REAL64), intent(in)           ::  tblaha                          !<  Tran-Blaha constant, if negative calculates it...
  logical, intent(in)                ::  lkincalc                        !<  Indicates that the kinetic energy density has been calculated.

  integer, intent(in)                ::  id1, id2                        !<  first and second dimension of the fft array
  integer, intent(in)                ::  n1, n2, n3                      !<  fft dimensions in directions 1,2,3
  real(REAL64), intent(in)           ::  rhomsh(id1,id2,n3)              !<  charge density (1/bohr^3)
  real(REAL64), intent(in)           ::  taumsh(id1,id2,n3)              !<  kinetic energy density (Hartree/bohr^3) [total for lxcmggavxc, correction for lxcmgga]
  real(REAL64), intent(in)           ::  dtau_dbdot(3,3,id1,id2,n3)   !<  (1/V) d (V taumsh) / d bdot on the mesh (Hartree/bohr^3), bdot with the 2 pi factors [only correction for lxcmgga]
  real(REAL64), intent(in)           ::  rholapmsh(id1,id2,n3)           !<  Laplacian of charge density (1/bohr^5)
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space (covariant components)

! output

  real(REAL64), intent(out)          ::  exc                             !<  the total exchange correlation energy given by the integral of the density times epsilon xc (Hartree).
  real(REAL64), intent(out)          ::  vxc(id1,id2,n3)                 !<  exchange-correlation potential vxc (Hartree).
  real(REAL64), intent(out)          ::  rhovxc                          !<  integral of the density times vxc (Hartree)
  real(REAL64), intent(out)          ::  strxc(3,3)                      !<  contribution of xc to the stress tensor (contravariant,a.u.)

! local variables

  real(REAL64)        ::  bdot(3,3), vcell                               !  metric in reciprocal space, cell volume
  real(REAL64)        ::  adotm1(3,3)                                    !  inverse of adot, bdot / (2 pi)^2
  real(REAL64)        ::  rho, epsx, epsc, vx, vc
  real(REAL64)        ::  strgga(3,3)                                    !  contribution to stress
  real(REAL64)        ::  grho,dexdr,decdr,dexdgr,decdgr
  real(REAL64)        ::  dexdtau, decdtau

  integer, parameter  ::  mxdnn = 3                                      !  Lagrange interpolation uses 2*mxdnn+1 points
  real(REAL64)        ::  dgdm(-mxdnn:mxdnn), drhocon(3)

  real(REAL64)        ::  coef

  real(REAL64)        ::  tau, rholap
  real(REAL64)        ::  tb09_integral, tb09_const_c

  real(REAL64)        ::  rhomax                                         !  maximum value of density
  real(REAL64)        ::  epsx_lda, epsc_lda, vx_lda, vc_lda             !  LDA XC
  real(REAL64)        ::  cprop                                          !  mix of TBL and LDA

  logical             ::  lxclda, lxcgga, lxcmgga, lxcmggavxc            !  family of xc functionals
  logical             ::  lxcgrad, lxclap, lxctau, lxctb09, lxccalc      !  properties of xc functionals

  real(REAL64)        ::  tauunif                         !  Thomas-Fermi kinetic energy density (hartree/bohr^3)
  real(REAL64)        ::  tausingle                       !  von Weizsaker correction to tauunif (hartree/bohr^3)

  real(REAL64)        ::  d_tauunif_dr                    !  d tauunif / d rho
  real(REAL64)        ::  d_tausingle_dr                  !  d tausingle / d rho
  real(REAL64)        ::  d_tausingle_dgr                 !  d tausingle / d grho
  real(REAL64)        ::  tautfvw                         !  tauunif + tausingle
  real(REAL64)        ::  d_exc_dgr                       !  d E_xc / d grad_rho
  real(REAL64)        ::  d_exc_dtau                      !  d E_xc / d tau

  real(REAL64)        ::  d_tau_dgr                       !  d tau / d grad_rho
  real(REAL64)        ::  d_tau_dr                        !  d tau / d rho
  real(REAL64)        ::  bdprod(3,3)                     !  dtau_dbdot * adotm1 at a mesh point

! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-18_REAL64
  real(REAL64), parameter  :: RHOEPS = 1.E-3_REAL64

! TB parameters

  real(REAL64), parameter  :: TB09_ALPHA = -0.012_REAL64
  real(REAL64), parameter  :: TB09_BETA = 1.023_REAL64

! counters

  integer       ::  i1, i2, i3
  integer       ::  i, j
  integer       ::  if1,if2
  integer       ::  in,jn,ip,nn


! initial stuff

  exc = ZERO
  rhovxc = ZERO
  do i3 = 1,n3
  do i2 = 1,n2
  do i1 = 1,n1
      vxc(i1,i2,i3) = ZERO
  enddo
  enddo
  enddo
  do i = 1,3
  do j = 1,3
    strxc(i,j) = ZERO
  enddo
  enddo


  if(n1 <= 0 .or. n2 <= 0 .or. n3 <= 0) return


  call xc_author_family(author, lxclda, lxcgga, lxcmgga, lxcmggavxc)

  call xc_author_info(author, lxcgrad, lxclap, lxctau, lxctb09, lxccalc)

  rhomax = rhomsh(1,1,1)
  do i3 = 1,n3
  do i2 = 1,n2
  do i1 = 1,n1
    if(rhomsh(i1,i2,i3) > rhomax) rhomax = rhomsh(i1,i2,i3)
  enddo
  enddo
  enddo

  if(n1 > id1 .or. n2 > id2) then
    write(6,'("   STOPPED in xc_cell:   wrong dimensions ",5i5)')  n1,n2,n3, id1,id2

    stop

  endif

  call adot_to_bdot(adot,vcell,bdot)

  do j = 1,3
  do i = 1,3
    adotm1(i,j) = bdot(i,j) / (4*PI*PI)
  enddo
  enddo

  do i = 1,3
  do j = 1,3
    strgga(i,j) = ZERO
  enddo
  enddo

  if(lxcgrad) then

    if(n1 < 3 .or. n2 < 3 .or. n3 < 3) then
      write(6,'("   STOPPED in xc_cell:   dimensions too small ",3i5)') n1,n2,n3

      stop

    endif

!   Weights for lagrange interpolation formula

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

  endif

! calculates the Tran-Blaha constant

  if(lxctb09) then

    if(tblaha > ZERO) then

      tb09_const_c = tblaha

    else

      tb09_integral = ZERO

      do i3=1,n3
      do i2=1,n2
      do i1=1,n1

        call xc_cell_deriv(rhomsh, i1,i2,i3, id1,id2, n1,n2,n3,          &
            nn, dgdm, adotm1, rho, grho, drhocon,                        &
            mxdnn)

        tb09_integral = tb09_integral + (grho/rho)

      enddo
      enddo
      enddo

      tb09_integral = tb09_integral / (n1*n2*n3)
      tb09_const_c = TB09_ALPHA + TB09_BETA*sqrt(tb09_integral)

      write(6,*)
      write(6,*) '   tb09_integral is: ',  tb09_integral
      write(6,*) '   tb09_const_c is:  ',  tb09_const_c
      write(6,*)

    endif

  endif

! LDA exchange and correlation

  if(lxclda) then

    do i3 = 1,n3
    do i2 = 1,n2
    do i1 = 1,n1
      rho = rhomsh(i1,i2,i3)

      call xc_lda( author, rho, epsx, epsc, vx, vc )

      exc = exc + rho * (epsx + epsc)
      vxc(i1,i2,i3) = vx + vc
      rhovxc = rhovxc + rho*(vx + vc)

    enddo
    enddo
    enddo

  elseif(lxcgga) then

!   gga exchange and correlation

    do i3 = 1,n3
    do i2 = 1,n2
    do i1 = 1,n1

      call xc_cell_deriv(rhomsh, i1,i2,i3, id1,id2, n1,n2,n3,            &
          nn, dgdm, adotm1, rho, grho, drhocon,                          &
          mxdnn)

      call xc_gga( author, rho, grho,                                    &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr )

      exc = exc + rho * (epsx + epsc)
      vxc(i1,i2,i3) = vxc(i1,i2,i3) + dexdr + decdr
      coef = (dexdgr + decdgr) * grho

      do i=1,3
      do j=1,3
        strgga(j,i) = strgga(j,i) + coef*drhocon(i)*drhocon(j)
      enddo
      enddo

      do in = -nn,nn
        ip = i1 + in
        ip = mod(ip+n1-1,n1) + 1
        vxc(ip,i2,i3) = vxc(ip,i2,i3) + n1*(dexdgr + decdgr)*drhocon(1)*dgdm(in)
      enddo
      do in = -nn,nn
        ip = i2 + in
        ip = mod(ip+n2-1,n2) + 1
        vxc(i1,ip,i3) = vxc(i1,ip,i3) + n2*(dexdgr + decdgr)*drhocon(2)*dgdm(in)
      enddo
      do in = -nn,nn
        ip = i3 + in
        ip = mod(ip+n3-1,n3) + 1
        vxc(i1,i2,ip) = vxc(i1,i2,ip) + n3*(dexdgr + decdgr)*drhocon(3)*dgdm(in)
      enddo
    enddo
    enddo
    enddo

    do i3=1,n3
    do i2=1,n2
    do i1=1,n1
      rhovxc = rhovxc + rhomsh(i1,i2,i3)*vxc(i1,i2,i3)
    enddo
    enddo
    enddo

  elseif(lxcmgga) then

    if(lkincalc) then

!     meta-gga exchange and correlation with both potential and energy

      do i3 = 1,n3
      do i2 = 1,n2
      do i1 = 1,n1

        call xc_cell_deriv(rhomsh, i1,i2,i3, id1,id2, n1,n2,n3,          &
            nn, dgdm, adotm1, rho, grho, drhocon,                        &
            mxdnn)

!       uses an approximate expression for tau and its derivatives

        call xc_tau(rho, grho, tauunif, tausingle,                       &
                  d_tauunif_dr, d_tausingle_dr, d_tausingle_dgr)

        tautfvw = tauunif + tausingle
        tau = taumsh(i1,i2,i3) + tautfvw
        rholap = ZERO

        call xc_mgga( author, rho, grho, tau,                            &
                     epsx, epsc, dexdr, decdr, dexdgr, decdgr,           &
                     dexdtau, decdtau  )

        exc = exc + rho * (epsx + epsc)

        d_exc_dtau = dexdtau + decdtau

        d_tau_dr = d_tauunif_dr + d_tausingle_dr
        d_tau_dgr = d_tausingle_dgr

        d_exc_dgr = dexdgr + decdgr

        vxc(i1,i2,i3) = vxc(i1,i2,i3) + dexdr + decdr + d_exc_dtau*d_tau_dr

        do i = 1,3
        do j = 1,3
          bdprod(j,i) = dtau_dbdot(j,1,i1,i2,i3)*adotm1(1,i) +        &
                        dtau_dbdot(j,2,i1,i2,i3)*adotm1(2,i) +        &
                        dtau_dbdot(j,3,i1,i2,i3)*adotm1(3,i)
        enddo
        enddo

        do i=1,3
        do j=1,3
          strgga(j,i) = strgga(j,i) + drhocon(j)*drhocon(i) * d_exc_dgr * grho
          strgga(j,i) = strgga(j,i) + drhocon(j)*drhocon(i) * d_exc_dtau * d_tau_dgr * grho

!         contribution of the correction taumsh to tau.  It scales as 1/volume
!         and V taumsh depends on the metric through dtau_dbdot = (1/V) d (V taumsh) / d bdot
!         (bdot with the 2 pi factors), so that -2 d E / d adot has
!         d_exc_dtau * ( taumsh adotm1 + 8 pi^2 adotm1 dtau_dbdot adotm1 ).

          strgga(j,i) = strgga(j,i) + d_exc_dtau * (                     &
                   taumsh(i1,i2,i3) * adotm1(j,i) +                      &
                   8*PI*PI * ( adotm1(j,1)*bdprod(1,i) +                 &
                               adotm1(j,2)*bdprod(2,i) +                 &
                               adotm1(j,3)*bdprod(3,i) ) )
        enddo
        enddo


        d_exc_dgr = d_exc_dgr + d_exc_dtau*d_tau_dgr

        do in = -nn,nn
          ip = i1 + in
          ip = mod(ip+n1-1,n1) + 1
          vxc(ip,i2,i3) = vxc(ip,i2,i3) + n1*d_exc_dgr*drhocon(1)*dgdm(in)
        enddo
        do in = -nn,nn
          ip = i2 + in
          ip = mod(ip+n2-1,n2) + 1
          vxc(i1,ip,i3) = vxc(i1,ip,i3) + n2*d_exc_dgr*drhocon(2)*dgdm(in)
        enddo
        do in = -nn,nn
          ip = i3 + in
          ip = mod(ip+n3-1,n3) + 1
          vxc(i1,i2,ip) = vxc(i1,i2,ip) + n3*d_exc_dgr*drhocon(3)*dgdm(in)
        enddo

      enddo
      enddo
      enddo

      do i3=1,n3
      do i2=1,n2
      do i1=1,n1
        rhovxc = rhovxc + rhomsh(i1,i2,i3)*vxc(i1,i2,i3)
      enddo
      enddo
      enddo

    else

!     kinetic energy has not been calculated, default to LDA

      do i3 = 1,n3
      do i2 = 1,n2
      do i1 = 1,n1
        rho = rhomsh(i1,i2,i3)

        call xc_lda( 'PW92', rho, epsx, epsc, vx, vc )

        exc = exc + rho * (epsx + epsc)
        vxc(i1,i2,i3) = vx + vc
        rhovxc = rhovxc + rho*(vx + vc)

      enddo
      enddo
      enddo

    endif


  elseif(lxcmggavxc) then

!   meta-gga exchange and correlation with only potential

    if(lkincalc) then

      do i3=1,n3
      do i2=1,n2
      do i1=1,n1

        call xc_cell_deriv(rhomsh, i1,i2,i3, id1,id2, n1,n2,n3,          &
            nn, dgdm, adotm1, rho, grho, drhocon,                        &
            mxdnn)

        tau = taumsh(i1,i2,i3)
        rholap = rholapmsh(i1,i2,i3)

!       avoids unphysical values due to noise at low densities
!       useful for slabs

        if(rho > RHOEPS*rhomax) then

          call xc_mgga_vxc('TB09','pz', rho, grho, rholap, tau,          &
                             epsx, epsc, vx, vc, tb09_const_c )

          vxc(i1,i2,i3) = vx + vc

        else

!         low density unstable region

          if(tau/rho > 1.0) tau = rho
          if(grho/rho > 2.5) grho = 2.5*rho

          call xc_mgga_vxc('TB09','pz', rho, grho, rholap, tau,          &
                             epsx, epsc, vx, vc, tb09_const_c )

          call xc_lda('pz', rho, epsx_lda, epsc_lda, vx_lda, vc_lda )

          cprop = (UM + COS(PI*rho/(RHOEPS*rhomax)))/2
          vxc(i1,i2,i3) = (UM-cprop)*(vx + vc) + cprop*(vx_lda+vc_lda)

          if(vxc(i1,i2,i3) >  20.0) vxc(i1,i2,i3) = 20.0
          if(vxc(i1,i2,i3) < -20.0) vxc(i1,i2,i3) =-20.0

        endif

      enddo
      enddo
      enddo

    else

!     kinetic energy has not been calculated, default to LDA

      do i3 = 1,n3
      do i2 = 1,n2
      do i1 = 1,n1
        rho = rhomsh(i1,i2,i3)

        call xc_lda( 'CA', rho, epsx, epsc, vx, vc )

        exc = exc + rho * (epsx + epsc)
        vxc(i1,i2,i3) = vx + vc
        rhovxc = rhovxc + rho*(vx + vc)

      enddo
      enddo
      enddo

    endif

  else

    write(6,'("    STOPPED in xc_cell:   unknown correlation")')
    write(6,*) '     ',author

    stop

  endif

! scales by volume factor

  coef = UM / (n1*n2*n3)
  exc  = coef* exc  * vcell
  rhovxc  = coef* rhovxc * vcell
  do i=1,3
  do j=1,3
    strgga(i,j) = coef * strgga(i,j) * vcell
  enddo
  enddo

! calculates stress

  do i=1,3
  do j=1,3
    strxc(i,j) = (rhovxc - exc) * adotm1(j,i) + strgga(i,j)
  enddo
  enddo

  return

end subroutine xc_cell

subroutine xc_cell_deriv(rhomsh, i1,i2,i3, id1,id2, n1,n2,n3,            &
       nn, dgdm, adotm1, rho, grho, drhocon,                             &
       mxdnn)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdnn                           !<  Lagrange interpolation uses at most 2*mxdnn+1 points

  integer, intent(in)                ::  id1, id2                        !<  first and second dimensions of the fft array
  integer, intent(in)                ::  n1, n2, n3                      !<  fft dimensions in directions 1,2,3

  integer, intent(in)                ::  i1, i2, i3                      !<  target point in the array
  real(REAL64), intent(in)           ::  rhomsh(id1,id2,n3)              !<  charge density (1/bohr^3)
  real(REAL64), intent(in)           ::  adotm1(3,3)                     !<  inverse of the metric adot, bdot / (2 pi)^2

  integer, intent(in)                ::  nn                              !<  Lagrange interpolation used 2*mxdnn+1 points
  real(REAL64), intent(in)           ::  dgdm(-mxdnn:mxdnn)              !<  Lagrange interpolation coefficients

! output

  real(REAL64), intent(out)          ::  rho                             !<  charge density
  real(REAL64), intent(out)          ::  grho                            !<  absolute value of gradient of rho
  real(REAL64), intent(out)          ::  drhocon(3)                      !<  components of gradient of rho (contravariant)

! allocatable arrays

  real(REAL64), allocatable          ::  drhodm(:)                       !  to be thread safe

! local variables

  integer           ::  ip

! counters

  integer           ::  in

! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: EPS = 1.0E-18_REAL64


  allocate(drhodm(3))

  rho = rhomsh(i1,i2,i3)

! calculates gradient of rho

  drhodm(1) = ZERO
  do in = -nn,nn
    ip = i1 + in
    ip = mod(ip+n1-1,n1) + 1
    drhodm(1) = drhodm(1) + dgdm(in)*rhomsh(ip,i2,i3)
  enddo
  drhodm(1) = n1*drhodm(1)

  drhodm(2) = ZERO
  do in = -nn,nn
    ip = i2 + in
    ip = mod(ip+n2-1,n2) + 1
    drhodm(2) = drhodm(2) + dgdm(in)*rhomsh(i1,ip,i3)
  enddo
  drhodm(2) = n2*drhodm(2)

  drhodm(3) = ZERO
  do in = -nn,nn
    ip = i3 + in
    ip = mod(ip+n3-1,n3) + 1
    drhodm(3) = drhodm(3) + dgdm(in)*rhomsh(i1,i2,ip)
  enddo
  drhodm(3) = n3*drhodm(3)

  drhocon(1) = adotm1(1,1)*drhodm(1) + adotm1(1,2)*drhodm(2) + adotm1(1,3)*drhodm(3)
  drhocon(2) = adotm1(2,1)*drhodm(1) + adotm1(2,2)*drhodm(2) + adotm1(2,3)*drhodm(3)
  drhocon(3) = adotm1(3,1)*drhodm(1) + adotm1(3,2)*drhodm(2) + adotm1(3,3)*drhodm(3)
  grho = drhodm(1)*drhocon(1) + drhodm(2)*drhocon(2) + drhodm(3)*drhocon(3)

  if(grho < EPS) then
    grho = ZERO
    drhocon(1) = ZERO
    drhocon(2) = ZERO
    drhocon(3) = ZERO
  else
    grho = sqrt(grho)
    drhocon(1) = drhocon(1) / grho
    drhocon(2) = drhocon(2) / grho
    drhocon(3) = drhocon(3) / grho
  endif

  deallocate(drhodm)

  return

end subroutine xc_cell_deriv
