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

!>     META-GGA routines for cpw. Inspired by libxc 1.2.0

       subroutine xc_cpw_tran_blaha(rho, sigma, lap_rho, tau, set_ctb09, &
     & vrho, vsigma, vtau, vlapl)

!      Written September 2015.  CLR
!      version 1.0 of xc, September 2015
!      Modified, documentation, December 2019. JLM
!      Copyright Carlos Loia Reis INESC-MN

!      version 4.94


       implicit none

       integer, parameter           :: REAL64 = selected_real_kind(12)

!      input
       
       real(REAL64), intent(in)           :: rho                         !<  charge density  (electrons per unit cell)
       real(REAL64), intent(in)           :: sigma                       !<  gradient of charge density
       real(REAL64), intent(in)           :: lap_rho                     !<  Laplacian of charge density
       real(REAL64), intent(in)           :: tau                         !<  kinetic energy density
       real(REAL64), intent(in)           :: set_ctb09                   !<  TB constant

!      output

       real(real64), intent(out)    :: vrho, vsigma, vtau, vlapl
       

!      local variables
              
       real(REAL64)                 :: rho1D, r_rs, r_x, r_t, r_u
       real(REAL64)                 :: v_BR, dv_BRdbx, dxdQ
       real(REAL64)                 :: r_f, r_dfdrs,r_dfdx
       real(REAL64)                 :: r_dfdt, r_dfdu, c_HEG
       real(REAL64)                 :: lrho,sfact,lsigma,gdm              

!      parameters

       real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter   :: X_FACTOR_C = 0.9305257363491_REAL64  !/* 3/8*cur(3/pi)*4^(2/3) */          
       real(REAL64), parameter   :: MIN_DENS = 5.0D-13
       real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64


       vrho   = ZERO
       vsigma = ZERO
       vtau   = ZERO
       vlapl  = ZERO

       if(max(rho, ZERO) < MIN_DENS) return
       

       call xc_cpw_becke_roussel_core(rho, sigma, lap_rho, tau,          &
     & rho1D, lrho, sfact, lsigma, gdm, r_rs, r_x, r_t, r_u,             &
     & v_BR, dv_BRdbx, dxdQ)
       
       r_f     =  ZERO  
       
       r_dfdx  =  ZERO
       r_dfdt  =  ZERO
       r_dfdu  =  ZERO
        
       r_dfdrs = -set_ctb09*v_BR
       c_HEG  = (3*set_ctb09 - 2*UM)*sqrt((5*UM)/(12*UM)) /              &
     &                                          (X_FACTOR_C*PI)    
       r_dfdrs = r_dfdrs -c_HEG*sqrt(r_t)
       r_dfdrs = r_dfdrs/(-r_rs) !/* due to the definition of dfdrs */
       
       call xc_cpw_out_mgga(rho1D, lrho, sfact, lsigma, gdm, r_rs,       &
     & r_x, r_t, r_u, r_f, r_dfdrs, r_dfdx, r_dfdt, r_dfdu,              &
     & vrho, vsigma, vtau, vlapl)
       
       end subroutine xc_cpw_tran_blaha


       subroutine xc_cpw_becke_roussel(rho, sigma, lap_rho, tau,         &
     & vrho, vsigma, vtau, vlapl)

       implicit none

       integer, parameter           :: REAL64 = selected_real_kind(12)
       
       real(REAL64), intent(in)     :: rho
       real(REAL64), intent(in)     :: sigma
       real(REAL64), intent(in)     :: lap_rho
       real(REAL64), intent(in)     :: tau
       

       real(REAL64), parameter   :: BR_GAMMA = 0.8_REAL64
       real(REAL64), parameter   :: MIN_DENS = 5.0D-13
       real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

       real(REAL64)                 :: lrho,sfact,lsigma, gdm              
       real(REAL64)                 :: rho1D, r_rs, r_x, r_t, r_u
       real(REAL64)                 :: v_BR, dv_BRdbx, dxdQ
       real(REAL64)                 :: r_f, r_dfdrs,r_dfdx
       real(REAL64)                 :: r_dfdt, r_dfdu
       
       real(real64), intent(out)    :: vrho, vsigma, vtau, vlapl
       
       vrho   = ZERO
       vsigma = ZERO
       vtau   = ZERO
       vlapl  = ZERO
        
       if(max(rho, ZERO) < MIN_DENS) return
       

       call xc_cpw_becke_roussel_core(rho, sigma, lap_rho, tau, rho1D,   &
     & lrho, sfact, lsigma, gdm, r_rs, r_x, r_t, r_u,                    &
     & v_BR, dv_BRdbx, dxdQ)
       
       r_f     =  - v_BR / 2
       r_dfdx  =  -r_x*br_gamma*dv_BRdbx*dxdQ/(12*UM)
       r_dfdt  =     2*br_gamma*dv_BRdbx*dxdQ/(12*UM)
       r_dfdu  =               -dv_BRdbx*dxdQ/(12*UM)
       r_dfdrs = ZERO ! is it??

       call xc_cpw_out_mgga(rho1D, lrho, sfact, lsigma, gdm, r_rs,       &
     & r_x, r_t, r_u, r_f, r_dfdrs, r_dfdx, r_dfdt, r_dfdu,              &
     & vrho, vsigma, vtau, vlapl)
       
       end subroutine xc_cpw_becke_roussel
!-----------------------------------------------------------------------                     
       subroutine xc_cpw_out_mgga(rho1D, lrho, sfact, lsigma, gdm,       &
     & r_rs, r_x, r_t, r_u, r_f, r_dfdrs, r_dfdx, r_dfdt, r_dfdu,        &
     & vrho, vsigma, vtau, vlapl)

       implicit none

       integer, parameter           :: REAL64 = selected_real_kind(12)

       real(REAL64), intent(in)     :: lrho,sfact,lsigma,gdm
       real(REAL64), intent(in)     :: rho1D,r_rs,r_x, r_t, r_u
       real(REAL64), intent(in)     :: r_f, r_dfdrs,r_dfdx
       real(REAL64), intent(in)     :: r_dfdt, r_dfdu
       real(real64), intent(out)    :: vrho, vsigma, vtau, vlapl
       

       real(REAL64), parameter   :: X_FACTOR_C = 0.9305257363491_REAL64  !/* 3/8*cur(3/pi)*4^(2/3) */          
       real(REAL64), parameter   :: MIN_GRAD = 5.0D-13
       real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       
       vrho   = ZERO
       vsigma = ZERO
       vtau   = ZERO
       vlapl  = ZERO
       
       
       vrho  = -X_FACTOR_C*rho1D*(-r_rs*r_dfdrs +                        &
     &                   (4*UM)/(3*UM)*(r_f - r_dfdx*r_x) -              &
     &                   (5.0*UM)/(3*UM)*(r_dfdt*r_t + r_dfdu*r_u))
       vtau  = -X_FACTOR_C*r_dfdt/rho1D
       vlapl = -X_FACTOR_C*r_dfdu/rho1D
      
       if(gdm > min_grad) then
        vsigma = -X_FACTOR_C*(rho1D*lrho)*r_dfdx*r_x /                   &
     &           ((2*UM)*sfact*lsigma)
       endif
        
       end subroutine xc_cpw_out_mgga
       
!-----------------------------------------------------------------------                     
       subroutine xc_cpw_becke_roussel_core(rho, sigma, lap_rho, tau,    &
     & rho1D, lrho, sfact, lsigma, gdm, r_rs, r_x, r_t, r_u,             &
     & v_BR, dv_BRdbx, dxdQ)


       implicit none

       integer, parameter           :: REAL64 = selected_real_kind(12)
       
       real(REAL64), intent(in)     :: rho
       real(REAL64), intent(in)     :: sigma
       real(REAL64), intent(in)     :: lap_rho
       real(REAL64), intent(in)     :: tau
       
       real(REAL64)                 :: br_Q
       real(REAL64)                 :: br_x

       real(REAL64), parameter   :: BR_GAMMA = 0.8_REAL64
       real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter   :: cbrtpi = PI**(1.0_REAL64/3.0_REAL64)
       real(REAL64), parameter   :: X_FACTOR_C = 0.9305257363491_REAL64  !/* 3/8*cur(3/pi)*4^(2/3) */          
       real(REAL64), parameter   :: UM = 1.0_REAL64
       real(REAL64), parameter :: CNST_RS = 0.6203504908994000866_REAL64      !  /* = POW(3.0/(4*M_PI), 1.0/3.0)*/       
              
       real(REAL64), parameter     :: MIN_GRAD = 5.0D-13
       real(REAL64), parameter     :: MIN_TAU  = 5.0D-13
       
       real(REAL64)                :: cnst,exp1,exp2
       real(REAL64)                :: ff,dffdx
       real(REAL64), intent(out)   :: rho1D, r_rs, r_x, r_t, r_u
       real(REAL64), intent(out)   :: v_BR, dv_BRdbx, dxdQ 
       
       real(REAL64), intent(out)   :: lrho,sfact,lsigma,gdm              
       
       real(REAL64)                :: sfact2,ltau,lnr2,rho2pD_D

        !sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
        sfact  = 2*UM
        sfact2 = sfact*sfact


!      lsigma= max(sigma[js]/sfact2, p->info->min_grad*p->info->min_grad);
!      gdm   = sqrt(lsigma);
!      lrho  = rho[is]/sfact;
!      rho1D = POW(lrho, 1.0/XC_DIMENSIONS);
!      rho2pD_D = lrho*rho1D*rho1D;
!      r.x   = gdm/(lrho*rho1D);
    
!      ltau  = tau[is]/sfact;
!      r.t   = ltau/rho2pD_D;  /* tau/rho^((2+D)/D) */

!      lnr2  = lapl[is]/sfact; /* this can be negative !*/
!      r.u   = lnr2/rho2pD_D;  /* lapl/rho^((2+D)/D) */

       r_rs = cnst_rs*rho**(-1.0D0/3.0D0)
!       r.rs = cnst_rs*POW(dens, -1.0/XC_DIMENSIONS);
               
       lsigma= max(sigma/sfact2, MIN_GRAD*MIN_GRAD)
       gdm   = sqrt(lsigma)
       lrho  = rho/sfact
       rho1D = lrho**(UM/(3*UM))
       rho2pD_D = lrho*rho1D*rho1D
       r_x   = gdm/(lrho*rho1D) 
       
       ltau  = tau/sfact   
       r_t   = ltau/rho2pD_D                                             !/* tau/rho^((2+D)/D) */
       
       lnr2    = lap_rho/sfact
       r_u     = lnr2/rho2pD_D;                                          !/* lapl/rho^((2+D)/D) */
    
       br_Q  = (r_u - 2*br_gamma*r_t + br_gamma*r_x*r_x/(2*UM))/(6*UM)
     
       call xc_cpw_becke_roussel_x(br_Q,br_x)
       
       cnst = -2*cbrtpi/X_FACTOR_C;
       exp1 = exp(br_x/(3*UM));
       exp2 = exp(-br_x);
      
      
       if (dABS(br_x) > MIN_TAU) then
        v_BR = exp1*(UM - exp2*(UM + br_x/(2*UM)))/br_x
       else
        v_BR = UM/(2*UM) + br_x/(6*UM) - br_x*br_x/(18*UM)
       endif
       
       v_BR = v_BR*cnst
              
       if(ABS(br_x) > MIN_TAU) then
          dv_BRdbx =(3*UM + br_x*(br_x + 2*UM) + (br_x - 3*UM)/exp2) /   &
     &                                (3*exp1*exp1*br_x*br_x)
       else 
          dv_BRdbx = UM/(6*UM) - br_x/(9*UM)
       endif

       dv_BRdbx  = dv_BRdbx * cnst
       ff    = br_x*exp(-2*UM/(3*UM)*br_x)/(br_x - 2*UM)
       dffdx = ff*(-2*UM/(3*UM) + UM/br_x - UM/(br_x - 2*UM))
       dxdQ  = -ff/(br_Q*dffdx)
    
       end subroutine xc_cpw_becke_roussel_core
!-----------------------------------------------------------------------                     
       subroutine xc_cpw_becke_roussel_x(br_Q,br_x)
       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

       real(REAL64),intent(in)     :: br_Q
       real(REAL64),intent(out)    :: br_x
       
       real(REAL64)                :: res
       real(REAL64)                :: rhs
       
       integer                     :: ierr
       
       real(REAL64), parameter   :: TOL = 5.0D-12
       real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter   :: UM = 1.0_REAL64
       
       real(REAL64)                :: xc_cpw_br_newt_raph 
       real(REAL64)                :: xc_cpw_br_bisect 
       
       rhs = ((2*UM/(3*UM))*PI**(2*UM/(3*UM)))/br_Q;
       
       br_x = xc_cpw_br_newt_raph(rhs, TOL, res, ierr);
       if(ierr == 0)then
          br_x = xc_cpw_br_bisect(rhs, TOL, ierr)
          if(ierr == 0) then
            write(6,*) 'Warning: Convergence not reached in',            &
     &          'Becke-Roussel functional'
            write(6,*) 'rhs:      ', rhs
            write(6,*) 'residual: ', res
          endif
       endif
              
       end subroutine       
!-----------------------------------------------------------------------       
       function xc_cpw_br_newt_raph( a, tol, res, ierr)
       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)
       
       real(REAL64)                :: xc_cpw_br_newt_raph
       real(REAL64)                :: a
       real(REAL64)                :: tol
       real(REAL64)                :: res
       integer                     :: ierr
      
       integer                     :: icount
       real(REAL64)                :: x
       real(REAL64)                :: f
       integer, parameter          :: max_iter = 50

       real(REAL64)                ::  arg, eee, xm2, fp

       real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

       ierr = 1
       if(a == ZERO) then 
        xc_cpw_br_newt_raph = ZERO
        return      
       endif

!   /* starting point */
!   x = (a < 0.0) ? -1.0 : 1.0;
       
       if(a < ZERO) then
         x = -UM
       else
         x =  UM
       endif
       
   
       icount = 0
       do icount=1, max_iter
          xm2 = x - 2*UM
          arg = 2*x/(3*UM)
          eee = exp(-arg)/a
          
          f  = x*eee - xm2
          fp = eee*(UM - 2*UM/(3*UM)*x) - UM

          x  = x - f/fp
          x  = abs(x)

          res = abs(f)
!          write(*,*) 'ntrf ', icount,x,res
          
          if(res<=tol) exit

       enddo 
      
       if(icount == max_iter) ierr=0 
       xc_cpw_br_newt_raph = x       

       end function xc_cpw_br_newt_raph
!-----------------------------------------------------------------------       
       function xc_cpw_br_bisect( a, tol, ierr)
       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)
       
       real(REAL64)                :: xc_cpw_br_bisect
       real(REAL64)                :: a
       real(REAL64)                :: tol
       integer                     :: ierr

       integer                     :: icount 
       real(REAL64)                :: f, x, x1, x2
       integer, parameter          :: max_iter = 500 

       real(REAL64)                :: arg, eee, xm2

       real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

       ierr = 1
       if(a == ZERO) then 
        xc_cpw_br_bisect = ZERO
        return      
       endif
        		   
!  /* starting interval */ 
       if(a > ZERO) then 
          x1 = 2*UM + tol
          x2 = UM/a + 2*UM
       else
          x2 = 2*UM - tol 
          x1 = ZERO        
       endif

!  /* bisection */ 
       do icount=1,max_iter
          x   = (x1 + x2)/(2*UM) 
          xm2 = x - 2*UM
          arg = 2*x/(3*UM)
          eee = exp(-arg) 
          f   = x*eee - a*xm2 
 
          if(f > ZERO) x1 = x
          if(f < ZERO) x2 = x 
 
          if (abs(f) > tol) then
          
            exit
            
          endif
       enddo       

       
       if(icount == max_iter) ierr=0
       xc_cpw_br_bisect = x 
       
       end function xc_cpw_br_bisect
!-----------------------------------------------------------------------       





