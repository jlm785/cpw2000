!>  Compares the exchange-correlation subroutines of cpw2000 (Src/ALFA)
!>  with libxc on a grid of reduced variables:
!>  rho, s = |grad rho| / (2 kF rho), alpha = (tau - tau_W) / tau_unif
!>  and q = lap rho / (4 (3 pi^2)^(2/3) rho^(5/3)).
!>  For each functional and quantity prints the largest error
!>  |cpw - libxc| / max(|libxc|, ref), where ref is the natural scale
!>  of the quantity built from the LDA exchange at that density,
!>  and the point where it occurs, for alpha > 0.  The last two columns
!>  are the largest errors on the alpha = 0 boundary (s > 0), where the
!>  derivatives may be one sided, and for alpha < 0 (tau < tau_W, which
!>  only appears through numerical noise), where cpw2000 and libxc follow
!>  different conventions.  Points with tau <= 0 are skipped.
!>  The meta-GGA with energy are repeated with the libxc enforcement of
!>  the Fermi hole curvature (sigma <= 8 rho tau) switched off.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         10 November 2026.
!>  \copyright    GNU Public License v2

program test_xc_grid

! Written 10 November 2026. claude

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! grid

  integer, parameter          ::  NRHO = 37                              !  rho = 10^(-6 + (i-1)/4)
  integer, parameter          ::  NS = 9
  integer, parameter          ::  NALPHA = 15
  integer, parameter          ::  NQ = 4
  integer, parameter          ::  NPT = NRHO*NS*NALPHA*NQ

  real(REAL64), parameter     ::  SGRID(NS) = [ 0.0_REAL64, 0.05_REAL64, 0.2_REAL64, 0.5_REAL64,  &
                        1.0_REAL64, 1.5_REAL64, 2.0_REAL64, 3.0_REAL64, 5.0_REAL64 ]
  real(REAL64), parameter     ::  AGRID(NALPHA) = [ -0.5_REAL64, -0.1_REAL64, -0.01_REAL64,        &
                        0.0_REAL64, 0.001_REAL64, 0.01_REAL64,                                &
                        0.1_REAL64, 0.5_REAL64, 0.9_REAL64, 1.0_REAL64, 1.1_REAL64,           &
                        2.0_REAL64, 2.6_REAL64, 5.0_REAL64, 10.0_REAL64 ]
  real(REAL64), parameter     ::  QGRID(NQ) = [ -0.5_REAL64, 0.0_REAL64, 0.5_REAL64, 2.0_REAL64 ]

  real(REAL64)                ::  rho(NPT), grho(NPT), sigma(NPT)
  real(REAL64)                ::  tau(NPT), lap(NPT)
  real(REAL64)                ::  sval(NPT), aval(NPT), qval(NPT)

! reference scales

  real(REAL64)                ::  eref(NPT), gref(NPT), tref(NPT)

! cpw2000 results

  real(REAL64)                ::  ex(NPT), ec(NPT), vx(NPT), vc(NPT)
  real(REAL64)                ::  dexdgr(NPT), decdgr(NPT)
  real(REAL64)                ::  dexdtau(NPT), decdtau(NPT)

! libxc results

  real(REAL64)                ::  exl(NPT), ecl(NPT), vxl(NPT), vcl(NPT)
  real(REAL64)                ::  vsx(NPT), vsc(NPT), vlx(NPT), vlc(NPT)
  real(REAL64)                ::  vtx(NPT), vtc(NPT)

  real(REAL64)                ::  ctb09
  logical                     ::  lfhc                                   !  libxc enforces Fermi hole curvature
  integer                     ::  ifhc
  real(REAL64)                ::  kf, tauw, tauu, e1, e2, v2
  integer                     ::  i, ir, is, ia, iq
  integer(8)                  ::  np

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

  np = NPT

! builds the grid

  i = 0
  do ir = 1,NRHO
  do is = 1,NS
  do ia = 1,NALPHA
  do iq = 1,NQ
    i = i + 1
    rho(i) = 10.0_REAL64**(-6 + (ir-1)*0.25_REAL64)
    kf = (3*PI*PI*rho(i))**(UM/3)
    grho(i) = 2*kf*rho(i)*SGRID(is)
    sigma(i) = grho(i)*grho(i)
    tauw = sigma(i) / (8*rho(i))
    tauu = (3*UM/10)*kf*kf*rho(i)
    tau(i) = tauw + AGRID(ia)*tauu
    lap(i) = 4*(3*PI*PI)**(2*UM/3)*rho(i)**(5*UM/3)*QGRID(iq)
    sval(i) = SGRID(is)
    aval(i) = AGRID(ia)
    qval(i) = QGRID(iq)
    eref(i) = (3*UM/4)*kf/PI
    gref(i) = eref(i) / (2*kf)
    tref(i) = eref(i)*rho(i) / tauu
  enddo
  enddo
  enddo
  enddo

  write(6,*)
  write(6,'("  cpw2000 vs libxc on ",i6," points")') NPT
  write(6,'("  rho from 1e-6 to 1e3, s from 0 to 5, alpha from 0 to 10, q from -0.5 to 2")')
  write(6,*)
  write(6,'("  functional   quantity    max err(a>0)    rho         s       alpha      q", &
          & "           cpw                libxc         err(a=0)    err(a<0)")')
  write(6,*)

! LDA

  do i = 1,NPT
    call xc_lda('CA', rho(i), ex(i), ec(i), vx(i), vc(i))
  enddo
  call lib_lda(1, NPT, rho, exl, vxl)
  call lib_lda(9, NPT, rho, ecl, vcl)
  call report('LDA_X', 'eps_x', ex, exl, eref)
  call report('LDA_X', 'v_x', vx, vxl, eref)
  call report('LDA CA/PZ', 'eps_c', ec, ecl, eref)
  call report('LDA CA/PZ', 'v_c', vc, vcl, eref)

  do i = 1,NPT
    call xc_lda('VWN', rho(i), ex(i), ec(i), vx(i), vc(i))
  enddo
  call lib_lda(7, NPT, rho, ecl, vcl)
  call report('LDA VWN', 'eps_c', ec, ecl, eref)
  call report('LDA VWN', 'v_c', vc, vcl, eref)

  do i = 1,NPT
    call xc_lda('WI', rho(i), ex(i), ec(i), vx(i), vc(i))
  enddo
  call lib_lda(2, NPT, rho, ecl, vcl)
  call report('LDA WI', 'eps_c', ec, ecl, eref)
  call report('LDA WI', 'v_c', vc, vcl, eref)

  do i = 1,NPT
    call xc_lda('PW92', rho(i), ex(i), ec(i), vx(i), vc(i))
  enddo
  call lib_lda(13, NPT, rho, ecl, vcl)
  call report('LDA PW92', 'eps_c', ec, ecl, eref)
  call report('LDA PW92', 'v_c', vc, vcl, eref)
  write(6,*)

! GGA

  do i = 1,NPT
    call xc_gga('PBE', rho(i), grho(i), ex(i), ec(i),                   &
                vx(i), vc(i), dexdgr(i), decdgr(i))
  enddo

  call lib_gga(320, 0, NPT, rho, sigma, exl, vxl, vsx)
  vsx = 2*vsx*grho
  call report('PBE_X(320)', 'eps_x', ex, exl, eref)
  call report('PBE_X(320)', 'v_x', vx, vxl, eref)
  call report('PBE_X(320)', 'dex/dgrho', dexdgr, vsx, gref)

  call lib_gga(130, 0, NPT, rho, sigma, ecl, vcl, vsc)
  vsc = 2*vsc*grho
  call report('PBE_C libxc', 'eps_c', ec, ecl, eref)
  call report('PBE_C libxc', 'v_c', vc, vcl, eref)
  call report('PBE_C libxc', 'dec/dgrho', decdgr, vsc, gref)

  call lib_gga(130, 1, NPT, rho, sigma, ecl, vcl, vsc)
  vsc = 2*vsc*grho
  call report('PBE_C paper', 'eps_c', ec, ecl, eref)
  call report('PBE_C paper', 'v_c', vc, vcl, eref)
  call report('PBE_C paper', 'dec/dgrho', decdgr, vsc, gref)
  write(6,*)

! meta-GGA with energy, first with the libxc default enforcement
! of the Fermi hole curvature, then without

  do ifhc = 1,2

  lfhc = (ifhc == 1)
  if(lfhc) then
    write(6,'("  libxc enforces the Fermi hole curvature (default)")')
  else
    write(6,'("  libxc does not enforce the Fermi hole curvature")')
  endif

  do i = 1,NPT
    call xc_mgga('LAK', rho(i), grho(i), tau(i),                         &
                 ex(i), ec(i), vx(i), vc(i), dexdgr(i), decdgr(i),       &
                 dexdtau(i), decdtau(i))
  enddo
  call lib_mgga(342, NPT, rho, sigma, lap, tau, exl, vxl, vsx, vlx, vtx)
  call lib_mgga(345, NPT, rho, sigma, lap, tau, ecl, vcl, vsc, vlc, vtc)
  call report_mgga('LAK')

  do i = 1,NPT
    call xc_mgga('TASK', rho(i), grho(i), tau(i),                        &
                 ex(i), ec(i), vx(i), vc(i), dexdgr(i), decdgr(i),       &
                 dexdtau(i), decdtau(i))
  enddo
  call lib_mgga(707, NPT, rho, sigma, lap, tau, exl, vxl, vsx, vlx, vtx)
  call lib_lda(13, NPT, rho, ecl, vcl)
  vsc = ZERO
  vtc = ZERO
  call report_mgga('TASK+PW92')

  do i = 1,NPT
    call xc_mgga('R2SC', rho(i), grho(i), tau(i),                        &
                 ex(i), ec(i), vx(i), vc(i), dexdgr(i), decdgr(i),       &
                 dexdtau(i), decdtau(i))
  enddo
  call lib_mgga(497, NPT, rho, sigma, lap, tau, exl, vxl, vsx, vlx, vtx)
  call lib_mgga(498, NPT, rho, sigma, lap, tau, ecl, vcl, vsc, vlc, vtc)
  call report_mgga('r2SCAN')
  write(6,*)

  enddo

  lfhc = .TRUE.

! potential only meta-GGA

  ctb09 = 1.09_REAL64
  do i = 1,NPT
    call xc_mgga_vxc('TB09', 'PW92', rho(i), grho(i), lap(i), tau(i),    &
                     e1, e2, vx(i), v2, ctb09)
  enddo
  call lib_tb09(ctb09, NPT, rho, sigma, lap, tau, vxl)
  call report('TB09 c=1.09', 'v_x', vx, vxl, eref)

! pure Becke-Roussel hole potential, v_BR = 3/2 v_TB09(c = 2/3)

  ctb09 = 2*UM/3
  do i = 1,NPT
    call xc_mgga_vxc('BR', 'PW92', rho(i), grho(i), lap(i), tau(i),      &
                     e1, e2, vx(i), v2, ctb09)
  enddo
  call lib_tb09(ctb09, NPT, rho, sigma, lap, tau, vxl)
  vxl = 3*vxl/2
  call report('BR', 'v_x', vx, vxl, eref)

! the old comparison with the BR89 energy functional

  call lib_mgga(586, NPT, rho, sigma, lap, tau, exl, vxl, vsx, vlx, vtx)
  call report('BR vs 586', 'v_x', vx, vxl, eref)
  write(6,*)

  stop

  contains

  subroutine report(label, qname, a, b, ref)

    character(len=*), intent(in)       ::  label, qname
    real(REAL64), intent(in)           ::  a(NPT), b(NPT), ref(NPT)

    real(REAL64)                       ::  err, errmax, errbnd, errneg
    integer                            ::  j, jmax, nnan

    errmax = -UM
    errbnd = ZERO
    errneg = ZERO
    jmax = 1
    nnan = 0
    do j = 1,NPT
      if(tau(j) <= ZERO) then
        cycle
      elseif(a(j) /= a(j) .or. b(j) /= b(j)) then
        nnan = nnan + 1
      else
        err = abs(a(j) - b(j)) / max(abs(b(j)), ref(j))
        if(aval(j) > ZERO) then
          if(err > errmax) then
            errmax = err
            jmax = j
          endif
        elseif(aval(j) < ZERO) then
          errneg = max(err,errneg)
        elseif(sval(j) > ZERO) then
          errbnd = max(err,errbnd)
        endif
      endif
    enddo

    write(6,'(2x,a12,1x,a10,es12.2,es12.2,3f8.3,2es19.10,2es12.2)') label, qname, errmax,    &
           rho(jmax), sval(jmax), aval(jmax), qval(jmax), a(jmax), b(jmax), errbnd, errneg
    if(nnan > 0) write(6,'(2x,"      WARNING:  ",i6," NaN values")') nnan

  end subroutine report


  subroutine report_mgga(label)

    character(len=*), intent(in)       ::  label

    call report(label, 'eps_x', ex, exl, eref)
    call report(label, 'eps_c', ec, ecl, eref)
    call report(label, 'v_x', vx, vxl, eref)
    call report(label, 'v_c', vc, vcl, eref)
    vsx = 2*vsx*grho
    vsc = 2*vsc*grho
    call report(label, 'dex/dgrho', dexdgr, vsx, gref)
    call report(label, 'dec/dgrho', decdgr, vsc, gref)
    call report(label, 'dex/dtau', dexdtau, vtx, tref)
    call report(label, 'dec/dtau', decdtau, vtc, tref)

  end subroutine report_mgga


  subroutine lib_lda(id, n, r, e, v)

    integer, intent(in)                ::  id, n
    real(REAL64), intent(in)           ::  r(n)
    real(REAL64), intent(out)          ::  e(n), v(n)

    type(xc_f03_func_t)                ::  func
    integer(8)                         ::  nn

    nn = n
    call xc_f03_func_init(func, id, XC_UNPOLARIZED)
    call xc_f03_lda_exc_vxc(func, nn, r, e, v)
    call xc_f03_func_end(func)

  end subroutine lib_lda


  subroutine lib_gga(id, iparam, n, r, s, e, v, vs)

    integer, intent(in)                ::  id, iparam, n
    real(REAL64), intent(in)           ::  r(n), s(n)
    real(REAL64), intent(out)          ::  e(n), v(n), vs(n)

    type(xc_f03_func_t)                ::  func
    integer(8)                         ::  nn

    nn = n
    call xc_f03_func_init(func, id, XC_UNPOLARIZED)
    if(iparam == 1) then
!     beta and gamma of the original PBE paper, as in xc_gga_c_pbe
      call xc_f03_func_set_ext_params_name(func, '_beta', 0.066725_REAL64)
      call xc_f03_func_set_ext_params_name(func, '_gamma', 0.031091_REAL64)
    endif
    call xc_f03_gga_exc_vxc(func, nn, r, s, e, v, vs)
    call xc_f03_func_end(func)

  end subroutine lib_gga


  subroutine lib_mgga(id, n, r, s, l, t, e, v, vs, vl, vt)

    integer, intent(in)                ::  id, n
    real(REAL64), intent(in)           ::  r(n), s(n), l(n), t(n)
    real(REAL64), intent(out)          ::  e(n), v(n), vs(n), vl(n), vt(n)

    type(xc_f03_func_t)                ::  func
    integer(8)                         ::  nn

    nn = n
    call xc_f03_func_init(func, id, XC_UNPOLARIZED)
    if(.not. lfhc) call xc_f03_func_set_fhc_enforcement(func, .FALSE.)
    call xc_f03_mgga_exc_vxc(func, nn, r, s, l, t, e, v, vs, vl, vt)
    call xc_f03_func_end(func)

  end subroutine lib_mgga


  subroutine lib_tb09(c, n, r, s, l, t, v)

    real(REAL64), intent(in)           ::  c
    integer, intent(in)                ::  n
    real(REAL64), intent(in)           ::  r(n), s(n), l(n), t(n)
    real(REAL64), intent(out)          ::  v(n)

    type(xc_f03_func_t)                ::  func
    real(REAL64)                       ::  vs(n), vl(n), vt(n)
    integer(8)                         ::  nn

    nn = n
    call xc_f03_func_init(func, 208, XC_UNPOLARIZED)
    call xc_f03_func_set_ext_params_name(func, 'c', c)
    call xc_f03_func_set_ext_params_name(func, '_alpha', ZERO)
    call xc_f03_mgga_vxc(func, nn, r, s, l, t, v, vs, vl, vt)
    call xc_f03_func_end(func)

  end subroutine lib_tb09

end program test_xc_grid
