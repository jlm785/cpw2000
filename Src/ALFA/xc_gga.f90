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

!>     Computes the exchange correlation energy and potential
!>     in the generalized gradient approximation (GGA).
!>     non-spin-polarized 
!>     Lengths in Bohr, energies in Hartrees.
!>     Adapted from a package by L.C. Balbas and J.M. Soler, Dec'96. Version 0.5. 

       subroutine xc_gga( author, rho, grho,                             &
     &                  epsx, epsc, dexdr, decdr, dexdgr, decdgr ) 


!      written 23 february 1999. jlm
!      copyright inesc/uvalladolid/uam jose soler,carlos balbas,
!      jorge pacheco,jose luis martins
!      Written 23 february 1999. jlm
!      Modified by J.L. Martins and J.M. Pacheco.
!      Modified, documentation, December 2019. JLM
!      copyright inesc/uvalladolid/uam Jose Soler, Carlos Balbas,
!      Jorge Pacheco, Jose Luis Martins

!      version 4.94

       implicit none 

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len = *), intent(in)     ::  author                     !<  type of xc wanted (ca=pz , pw92 , vwn, wi)
       real(REAL64), intent(in)           ::  rho                        !<  electron density (1/bohr^3)
       real(REAL64), intent(in)           ::  grho                       !<  gradient of the electron density (1/bohr^4)

!      output

       real(REAL64), intent(out)          ::  epsx                       !<  exchange energy density (hartree/bohr^3) 
       real(REAL64), intent(out)          ::  epsc                       !<  correlation energy density (hartree/bohr^3)
       real(REAL64), intent(out)          ::  dexdr                      !<  derivative with respect to rho of rho*epsx (hartree/bohr^3)
       real(REAL64), intent(out)          ::  decdr                      !<  derivative with respect to rho of rho*epsc (hartree/bohr^3)
       real(REAL64), intent(out)          ::  dexdgr                     !<  derivative with respect to grho of rho*epsx (hartree/bohr^2)
       real(REAL64), intent(out)          ::  decdgr                     !<  derivative with respect to grho of rho*epsc (hartree/bohr^2)

!      local variables

       real(REAL64)     :: rs                                            !  Wigner-Seitz radius
       real(REAL64)     :: a0                                            !  exchange constant
       real(REAL64)     :: sqrs
       real(REAL64)     :: grloc                                         !  local value of grho

       real(REAL64)     ::  ecunif,vcunif,exlda
       real(REAL64)     ::  xkf,xks,s,t,expec,a,at2
       real(REAL64)     ::  qat2,arglog,h,fx
       real(REAL64)     ::  decudrs,drsdr
!       real(REAL64)     ::  dkfdr,dksdr
       real(REAL64)     ::  dtdr,dsdr,dtdgr,dsdgr
       real(REAL64)     ::  dadr,dat2dr,dqat2dat2
       real(REAL64)     ::  dargdr,dargdt,dhdr,dhdt
       real(REAL64)     ::  dfxds

!      parameters

       real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: EPS = 1.0d-24

!      PW92 parameters

       real(REAL64)             ::  b_pw, c_pw, dc_pwdrs, db_pwdrs
       real(REAL64), parameter  ::  BETA1_PW = 7.5957_REAL64
       real(REAL64), parameter  ::  BETA2_PW = 3.5876_REAL64
       real(REAL64), parameter  ::  BETA3_PW = 1.6382_REAL64
       real(REAL64), parameter  ::  BETA4_PW = 0.49294_REAL64
       real(REAL64), parameter  ::  A_PW = 0.062182_REAL64
       real(REAL64), parameter  ::  ALP_PW = 0.21370_REAL64

!      PBE parameters

       real(REAL64), parameter  ::  BETA = 0.066725_REAL64
       real(REAL64), parameter  ::  GAMMA =  0.0310906908696549_REAL64
       real(REAL64), parameter  ::  KAPPA = 0.804_REAL64
       real(REAL64), parameter  ::  MU = BETA * PI*PI / 3

!      initial stuff 

       epsx = ZERO
       epsc = ZERO
       dexdr = ZERO
       decdr = ZERO
       dexdgr = ZERO
       decdgr = ZERO

       if(rho < EPS) return

       grloc = grho
       if(grho < ZERO) grloc = ZERO 

       if (author == 'pbe' .or. author == 'PBE') then 

         rs = (3/(4*PI*rho))**(UM/3)

!        P e r d e w     &     W a n g     1 9 9 2 

!        Implements the perdew-wang'92 local correlation (beyond rpa). 
!        J.P. Perdew and Y. Wang, PRB 45, 13244 (1992) using eq.(10). lcb, jms

         sqrs = sqrt(rs)
         
         b_pw = BETA1_pw*sqrs +   BETA2_pw*rs +  BETA3_pw*rs*sqrs +      &
     &          BETA4_pw*rs*rs 
         db_pwdrs = BETA1_pw / (2*sqrs) + BETA2_pw +                     &
     &              (3*BETA3_pw*sqrs) / 2 + 2*BETA4_PW*rs 
         c_pw = UM + UM / (A_PW*b_pw) 
         dc_pwdrs = - (c_pw - UM)*db_pwdrs / b_pw

         ecunif = - A_PW * (UM + ALP_PW*rs) *  log(c_pw) 
         decudrs = -A_PW * (ALP_PW * log(c_pw) +                         &
     &                         (UM + ALP_PW*rs) * dc_pwdrs / c_pw ) 
         
         vcunif = ecunif - rs*decudrs / 3

!        adds non local correlation 
 
         xkf = (3*PI*PI * rho)**(UM/3) 
         xks = sqrt(4*xkf / PI )  
         t = grloc / (2 * xks * rho)
         s = grloc / (2 * xkf * rho)
         expec = exp(-ecunif / GAMMA)
         a = (BETA / GAMMA) * (UM / (expec - UM)) 
         at2 = a*t*t
         qat2 = (UM + at2)/(UM + at2 + at2*at2)
         arglog = UM + (BETA/GAMMA) * t* t * qat2
         h = GAMMA * log(arglog)

         epsc = ecunif + h 

         drsdr = - rs / (3*rho) 
!         dkfdr =   xkf / (3*rho) 
!         dksdr =  xks / (6*rho)
         dtdr = -7*t / (6*rho)
         dsdr = -4*s / (3*rho)
         dtdgr = UM / (2*xks*rho)
         dsdgr = UM / (2*xkf*rho)

         dadr = +(BETA / (GAMMA*GAMMA))*                                 &
     &    ( expec / ((expec - UM)*(expec - UM)) )*decudrs*drsdr
         dat2dr = t*t*dadr + 2*a*t*dtdr
         dqat2dat2 = - at2 * (2*UM + at2) /                              &
     &                 ((UM + at2 + at2*at2)*(UM +at2 + at2*at2))
         dargdr = 2*(BETA/GAMMA)*t*qat2*dtdr +                           &
     &               (BETA/GAMMA)*t*t*dqat2dat2*dat2dr
         dargdt = 2*(BETA/GAMMA)*t*qat2 +                                &
     &            2*(BETA/GAMMA)*t*t*dqat2dat2*a*t
         dhdr = (GAMMA/arglog)*dargdr
         dhdt = (GAMMA/arglog)*dargdt
         decdr = vcunif + h + rho*dhdr

         decdgr = rho*dhdt*dtdgr

       else 
         write(6,'("   STOPPED in xc_gga:   unknown correlation")')
         write(6,*) '     ',author

         stop 

       endif

!                                 e x c h a n g e 

       a0 = (4 / (9*PI))**(UM/3)
       exlda = -((3*UM) / (4*UM)) / (PI*a0*rs)
       fx = UM + KAPPA - KAPPA / (UM + MU*s*s/KAPPA) 

       epsx = exlda * fx 

       dfxds = 2*MU*s /  ((UM + MU*s*s/KAPPA)*(UM + MU*s*s/KAPPA))

       dexdr = ((4*UM) / (3*UM))*exlda*fx + rho*exlda*dfxds*dsdr
       dexdgr = rho*exlda*dfxds*dsdgr

       return
       end subroutine xc_gga
