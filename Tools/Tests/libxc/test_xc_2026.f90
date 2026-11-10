program test_xc

  use xc_f03_lib_m

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  real(REAL64)           ::  rho                             !<  electron density (1/bohr^3)
  real(REAL64)           ::  grho                            !<  gradient of the electron density (1/bohr^4)
  real(REAL64)           ::  lap_rho                         !<  laplacian of rho
  real(REAL64)           :: tau                              !<  kinetic energy density

! output

  real(REAL64)          ::  epsx, epsx_old, epsx_lib                     !<  exchange energy density (hartree/bohr^3)
  real(REAL64)          ::  epsc, epsc_old, epsc_lib                     !<  correlation energy density (hartree/bohr^3)
  real(REAL64)          ::  vx, vx_old, vx_lib                           !<  exchange potential (hartree/bohr^3)
  real(REAL64)          ::  vc, vc_old, vc_lib                           !<  correlation potential (hartree/bohr^3)

  real(REAL64)          ::  dexdr, dexdr_old, dexdr_lib
  real(REAL64)          ::  decdr, decdr_old, decdr_lib
  real(REAL64)          ::  dexdgr, dexdgr_old, dexdgr_lib
  real(REAL64)          ::  decdgr, decdgr_old, decdgr_lib
  real(REAL64)          ::  dexdtau, decdtau
  real(REAL64)          ::  dexdtau_lib, decdtau_lib

! libxc interface stuff

  TYPE(xc_f03_func_t)        ::  xc_func
!  TYPE(xc_f03_func_info_t)   ::  xc_info
  integer(8), parameter      ::  np = 1
  integer                    ::  func_id
  real(REAL64)               ::  rhotmp(1)
  real(REAL64)               ::  sigtmp(1)
  real(REAL64)               ::  etmp(1)
  real(REAL64)               ::  vtmp(1)
  real(REAL64)               ::  dedsigtmp(1)
  real(REAL64)               ::  tautmp(1)
  real(REAL64)               ::  laptmp(1)
  real(REAL64)               ::  dedgtmp(1)
  real(REAL64)               ::  vtautmp(1)
  real(REAL64)               ::  vlaptmp(1)
  real(REAL64)               ::  set_ctb09(2)

!   rho = 5.30_REAL64
!   grho = 1.7_REAL64
!   lap_rho = 0.4_REAL64
!   tau = 0.47_REAL64

  rho = 0.2715331182_REAL64
  grho = 0.2048515062_REAL64
  lap_rho = 0.0_REAL64
  tau = 0.4239425252_REAL64

  rhotmp(1) = rho
  sigtmp(1) = grho*grho
  tautmp(1) = tau

  write(6,*)
  write(6,'(10x,"E_x",13x,"E_c",6x,"d E_x / d rho",5x,"d E_c / d rho",   &
        &       3x,"d E_x / d grho",3x,"d E_c / d grho")')


  write(6,*)
  write(6,*)  '  testing CA'
  write(6,*)

  call xc_lda('CA', rho, epsx, epsc, vx, vc)
!  call xc_lda_old('CA', rho, epsx_old, epsc_old, vx_old, vc_old)

  func_id = 1

   call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
     call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
     epsx_lib = etmp(1)
     vx_lib = vtmp(1)
   call xc_f03_func_end(xc_func)

  func_id = 9

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
    epsc_lib = etmp(1)
    vc_lib = vtmp(1)
  call xc_f03_func_end(xc_func)

  write(6,'(4f16.10)') epsx, epsc, vx, vc
!  write(6,'(4f16.10,"    OLD")') epsx_old-epsx, epsc_old-epsc,           &
!          vx_old-vx, vc_old-vc
  write(6,'(4f16.10,"    LIBXC")') epsx_lib-epsx, epsc_lib-epsc,         &
          vx_lib-vx, vc_lib-vc


  write(6,*)
  write(6,*)  '  testing VWN'
  write(6,*)

  call xc_lda('VWN', rho, epsx, epsc, vx, vc)
!  call xc_lda_old('VWN', rho, epsx_old, epsc_old, vx_old, vc_old)

  func_id = 1

   call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
     call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
     epsx_lib = etmp(1)
     vx_lib = vtmp(1)
   call xc_f03_func_end(xc_func)

  func_id = 7

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
    epsc_lib = etmp(1)
    vc_lib = vtmp(1)
  call xc_f03_func_end(xc_func)

  write(6,'(4f16.10)') epsx, epsc, vx, vc
!  write(6,'(4f16.10,"    OLD")') epsx_old-epsx, epsc_old-epsc,           &
!          vx_old-vx, vc_old-vc
  write(6,'(4f16.10,"    LIBXC")') epsx_lib-epsx, epsc_lib-epsc,         &
          vx_lib-vx, vc_lib-vc


  write(6,*)
  write(6,*)  '  testing WI'
  write(6,*)

  call xc_lda('WI', rho, epsx, epsc, vx, vc)
!  call xc_lda_old('WI', rho, epsx_old, epsc_old, vx_old, vc_old)

  func_id = 1

   call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
     call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
     epsx_lib = etmp(1)
     vx_lib = vtmp(1)
   call xc_f03_func_end(xc_func)

  func_id = 2

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
    epsc_lib = etmp(1)
    vc_lib = vtmp(1)
  call xc_f03_func_end(xc_func)

  write(6,'(4f16.10)') epsx, epsc, vx, vc
!  write(6,'(4f16.10,"    OLD")') epsx_old-epsx, epsc_old-epsc,           &
!          vx_old-vx, vc_old-vc
  write(6,'(4f16.10,"    LIBXC")') epsx_lib-epsx, epsc_lib-epsc,         &
          vx_lib-vx, vc_lib-vc


  write(6,*)
  write(6,*)  '  testing PW92'
  write(6,*)

  call xc_lda('PW92', rho, epsx, epsc, vx, vc)
!  call xc_lda_old('PW92', rho, epsx_old, epsc_old, vx_old, vc_old)

  func_id = 1

   call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
     call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
     epsx_lib = etmp(1)
     vx_lib = vtmp(1)
   call xc_f03_func_end(xc_func)

  func_id = 13

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_lda_exc_vxc(xc_func, np, rhotmp, etmp, vtmp)
    epsc_lib = etmp(1)
    vc_lib = vtmp(1)
  call xc_f03_func_end(xc_func)

  write(6,'(4f16.10)') epsx, epsc, vx, vc
!  write(6,'(4f16.10,"    OLD")') epsx_old-epsx, epsc_old-epsc,           &
!          vx_old-vx, vc_old-vc
  write(6,'(4f16.10,"    LIBXC")') epsx_lib-epsx, epsc_lib-epsc,         &
          vx_lib-vx, vc_lib-vc



  write(6,*)
  write(6,*)  '  testing PBE'
  write(6,*)

  call xc_gga('PBE', rho, grho, epsx, epsc,                              &
              dexdr, decdr, dexdgr, decdgr )


!  call xc_gga_old('PBE', rho, grho, epsx_old, epsc_old,                  &
!                  dexdr_old, decdr_old, dexdgr_old, decdgr_old )

  func_id = 130

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_gga_exc_vxc(xc_func, np, rhotmp, sigtmp, etmp, vtmp, dedsigtmp)
    epsc_lib = etmp(1)
    decdr_lib = vtmp(1)
    decdgr_lib = 2*dedsigtmp(1)*grho
  call xc_f03_func_end(xc_func)

  func_id = 320

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_gga_exc_vxc(xc_func, np, rhotmp, sigtmp, etmp, vtmp, dedsigtmp)
    epsx_lib = etmp(1)
    dexdr_lib = vtmp(1)
    dexdgr_lib = 2*dedsigtmp(1)*grho
  call xc_f03_func_end(xc_func)

  write(6,'(6f16.10)') epsx, epsc, dexdr, decdr, dexdgr, decdgr
!  write(6,'(6f16.10,"    OLD")') epsx_old-epsx, epsc_old-epsc,           &
!          dexdr_old-dexdr, decdr_old-decdr,                              &
!          dexdgr_old-dexdgr, decdgr_old-decdgr
  write(6,'(6f16.10,"    LIBXC")') epsx_lib-epsx, epsc_lib-epsc,         &
          dexdr_lib-dexdr, decdr_lib-decdr,                              &
          dexdgr_lib-dexdgr, decdgr_lib-decdgr



  write(6,*)
  write(6,*)  '  testing TB09'
  write(6,*)

  laptmp(1) = lap_rho
  tautmp(1) = tau
  set_ctb09(1) = 1.09_REAL64
  set_ctb09(2) = 0.0_REAL64

! call xc_mgga_old('PW92',rho, grho, lap_rho, 2*tau,                    &
!            epsx_old, epsc_old, vx_old, vc_old,  1.09)
!  call xc_tran_blaha_old(rho, grho, lap_rho, 2*tau, 1.09, vx_old)
  call xc_mgga_vxc('TB09', 'PW92', rho, grho, lap_rho, tau,                    &
                     epsx, epsc, vx, vc, 1.09_REAL64 )


  func_id = 208                             ! tb09

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    if(func_id == 208) call xc_f03_func_set_ext_params(xc_func, set_ctb09)
    call xc_f03_mgga_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp, vtmp, dedgtmp, vlaptmp, vtautmp )
    vx_lib = vtmp(1)
  call xc_f03_func_end(xc_func)

  write(6,'(f16.10)') vx
!  write(6,'(f16.10,"    OLD")') vx_old
  write(6,'(f16.10,"    LIBXC")') vx_lib


  write(6,*)
  write(6,*)  '  testing LAK'
  write(6,*)

!  call xc_mgga_x_lak(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau )
  call xc_mgga( 'LAK', rho, grho, tau,                                   &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr,             &
                   dexdtau, decdtau )

  func_id = 342                             ! LAK exchange

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_mgga_exc_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp,     &
    etmp, vtmp, dedgtmp, vlaptmp, vtautmp)
    epsx_lib = etmp(1)
    dexdr_lib = vtmp(1)
    dexdgr_lib = 2*dedgtmp(1)*grho
    dexdtau_lib = vtautmp(1)
  call xc_f03_func_end(xc_func)

!  call xc_mgga_c_lak(rho, grho, tau, epsc, decdr, decdgr, decdtau )

  func_id = 345                             ! LAK correlation

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_mgga_exc_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp,     &
    etmp, vtmp, dedgtmp, vlaptmp, vtautmp)
    epsc_lib = etmp(1)
    decdr_lib = vtmp(1)
    decdgr_lib = 2*dedgtmp(1)*grho
    decdtau_lib = vtautmp(1)
  call xc_f03_func_end(xc_func)


  write(6,'(8f16.10)') epsx, epsc, dexdr, decdr, dexdgr, decdgr, dexdtau, decdtau
  write(6,'(8f16.10)') epsx_lib-epsx, epsc_lib-epsc, dexdr_lib-dexdr,            &
                      decdr_lib-decdr, dexdgr_lib-dexdgr, decdgr_lib-decdgr,    &
                      dexdtau_lib-dexdtau, decdtau_lib-decdtau



  write(6,*)
  write(6,*)  '  testing TASK'
  write(6,*)

!  call xc_mgga_x_task(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau )
  call xc_mgga( 'TASK', rho, grho, tau,                                  &
                   epsx, epsc, dexdr, decdr, dexdgr, decdgr,             &
                   dexdtau, decdtau )


  func_id = 707                             ! TASK exchange

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_mgga_exc_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp,     &
    etmp, vtmp, dedgtmp, vlaptmp, vtautmp)
    epsx_lib = etmp(1)
    dexdr_lib = vtmp(1)
    dexdgr_lib = 2*dedgtmp(1)*grho
    dexdtau_lib = vtautmp(1)
  call xc_f03_func_end(xc_func)


  write(6,'(4f16.10)') epsx, dexdr, dexdgr, dexdtau
  write(6,'(4f16.10)') epsx_lib-epsx, dexdr_lib-dexdr,                      &
                       dexdgr_lib-dexdgr, dexdtau_lib-dexdtau
  write(6,*)



  write(6,*)
  write(6,*)  '  testing Becke-Roussel'
  write(6,*)

!  call xc_mgga_x_task(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau )
  call xc_mgga_vxc('BR', 'PW92', rho, grho, lap_rho, tau,                    &
                     epsx, epsc, vx, vc, 1.09_REAL64 )



  func_id = 586                             ! BR

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_mgga_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp, vtmp, dedgtmp, vlaptmp, vtautmp )
    vx_lib = vtmp(1)
  call xc_f03_func_end(xc_func)

  write(6,'(f16.10)') vx
  write(6,'(f16.10,"    LIBXC")') vx_lib



  write(6,*)
  write(6,*)  '  testing r2SCAN'
  write(6,*)

  call xc_mgga_x_r2SCAN(rho, grho, tau, epsx, dexdr, dexdgr, dexdtau )
!   call xc_mgga( 'r2SC', rho, grho,  lap_rho, tau,                        &
!                    epsx, epsc, dexdr, decdr, dexdgr, decdgr,             &
!                    dexdtau, decdtau )
  call xc_mgga_c_r2SCAN(rho, grho, tau, epsc, decdr, decdgr, decdtau )



  func_id = 497                             ! r2SCAN exchange

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_mgga_exc_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp,     &
    etmp, vtmp, dedgtmp, vlaptmp, vtautmp)
    epsx_lib = etmp(1)
    dexdr_lib = vtmp(1)
    dexdgr_lib = 2*dedgtmp(1)*grho
    dexdtau_lib = vtautmp(1)
  call xc_f03_func_end(xc_func)

  func_id = 498                             ! r2SCAN correlation

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
    call xc_f03_mgga_exc_vxc(xc_func, np, rhotmp, sigtmp, laptmp, tautmp,     &
    etmp, vtmp, dedgtmp, vlaptmp, vtautmp)
    epsc_lib = etmp(1)
    decdr_lib = vtmp(1)
    decdgr_lib = 2*dedgtmp(1)*grho
    decdtau_lib = vtautmp(1)
  call xc_f03_func_end(xc_func)


  write(6,'(8f16.10)') epsx, epsc, dexdr, decdr, dexdgr, decdgr, dexdtau, decdtau
  write(6,'(8f16.10)') epsx_lib-epsx, epsc_lib-epsc, dexdr_lib-dexdr,            &
                      decdr_lib-decdr, dexdgr_lib-dexdgr, decdgr_lib-decdgr,    &
                      dexdtau_lib-dexdtau, decdtau_lib-decdtau


  write(6,*)

  stop

end program test_xc
