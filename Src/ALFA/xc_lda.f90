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
!>     in the local density approximation .
!>     Lengths in Bohr, energies in Hartrees.
!>     Adapted from a package by L.C. Balbas and J.M. Soler, Dec'96. Version 0.5. 

       subroutine xc_lda( author, rho, epsx, epsc, vx, vc ) 

!      Written 23 february 1999. jlm
!      Modified by J.L. Martins and J.M. Pacheco.
!      copyright inesc/uvalladolid/uam Jose Soler, Carlos Balbas
!      Jorge Pacheco, Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len = *), intent(in)     ::  author                     !<  type of xc wanted (ca=pz , pw92 , vwn, wi)
       real(REAL64), intent(in)           ::  rho                        !<  electron density (1/bohr^3)

!      output

       real(REAL64), intent(out)          ::  epsx                       !<  exchange energy density (hartree/bohr^3) 
       real(REAL64), intent(out)          ::  epsc                       !<  correlation energy density (hartree/bohr^3)
       real(REAL64), intent(out)          ::  vx                         !<  exchange potential (hartree/bohr^3)
       real(REAL64), intent(out)          ::  vc                         !<  correlation potential (hartree/bohr^3)

!      local variables

       real(REAL64)     :: rs                                            !  Wigner-Seitz radius
       real(REAL64)     :: a0                                            !  exchange constant
       real(REAL64)     :: sqrs,rslog,rsdecdrs


!      parameters

       real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: EPS = 1.0d-24

!      PZ parameters

       real(REAL64)             ::  be
       real(REAL64), parameter  ::  A1_PZ = 1.0529_REAL64
       real(REAL64), parameter  ::  A2_PZ = 0.3334_REAL64
       real(REAL64), parameter  ::  B_PZ = 0.1423_REAL64
       real(REAL64), parameter  ::  C1_PZ = 0.0311_REAL64
       real(REAL64), parameter  ::  C2_PZ = 0.002_REAL64
       real(REAL64), parameter  ::  C3_PZ = 0.048_REAL64
       real(REAL64), parameter  ::  C4_PZ = 0.0116_REAL64

!      VWN parameters

       real(REAL64)             ::  xd_vwn, q_vwn
       real(REAL64)             ::  fac1, fac2, fac3
       real(REAL64)             ::  b1, b2, b3
       real(REAL64)             ::  arg1, arg2, arg3
       real(REAL64), parameter  ::  X_VWN = 0.10498_REAL64
       real(REAL64), parameter  ::  A_VWN = 0.0310907_REAL64
       real(REAL64), parameter  ::  B_VWN = 3.72744_REAL64
       real(REAL64), parameter  ::  C_VWN = 12.9352_REAL64

!      Wigner parameters

       real(REAL64), parameter  ::  A_WI = -0.4377645_REAL64
       real(REAL64), parameter  ::  B_WI = 7.8_REAL64

!      PW92 parameters

       real(REAL64)             ::  b_pw, c_pw, dc_pwdrs, db_pwdrs
       real(REAL64), parameter  ::  BETA1_PW = 7.5957_REAL64
       real(REAL64), parameter  ::  BETA2_PW = 3.5876_REAL64
       real(REAL64), parameter  ::  BETA3_PW = 1.6382_REAL64
       real(REAL64), parameter  ::  BETA4_PW = 0.49294_REAL64
       real(REAL64), parameter  ::  A_PW = 0.062182_REAL64
       real(REAL64), parameter  ::  ALP_PW = 0.21370_REAL64


!      initial stuff 

       epsx = ZERO 
       epsc = ZERO 
       vx = ZERO 
       vc = ZERO

       if(rho < EPS) return 

       rs = (3/(4*PI*rho))**(UM/3)

!      E x c h a n g e 

       a0 = (4 / (9*PI))**(UM/3)
       vx = -UM / (PI*a0*rs) 
       epsx = ( (3*UM) / 4) * vx

!      C o r r e l a t i o n  

       if ( author == 'ca' .or. author == 'CA' .or.                      &
     &      author == 'pz' .or. author == 'PZ') then 

!        P e r d e w    a n d     Z u n g e r 

!        Perdew-Zunger parameterization of Ceperley-Alder exchange and  
!        correlation. Ref: Perdew and Zunger, Phys. Rev. B 23 5075 (1981). jlm

         if (rs > UM) then   
           sqrs = sqrt(rs)
           be = UM + A1_PZ*sqrs + A2_PZ*rs 
           epsc = -B_PZ / be
           rsdecdrs = B_PZ * (A1_PZ*sqrs/2 + A2_PZ*rs) / (be*be)
           vc = epsc - rsdecdrs/3  
         else 
           rslog = log(rs) 
           epsc = (C1_PZ + C2_PZ*rs)*rslog - C3_PZ - C4_PZ*rs
           rsdecdrs = (C1_PZ + C2_PZ*rs) + C2_PZ*rslog*rs - C4_PZ*rs
           vc = epsc - rsdecdrs/3 
         endif 

       elseif ( author == 'vwn' .or. author == 'VWN' ) then

!        V o s k o     W i l k     and     N u s a i r 

!        Based on the papers of W. Pickett cpr(89)117 
!        (this paper is for LDA and has the right formula, eq. b.11)
!        Wilk and Vosko, J. Phys. C15 (82) 2139. jmp 
 
         sqrs  = sqrt(rs) 
         xd_vwn = rs + B_VWN*sqrs + C_VWN 
         q_vwn  = sqrt(4*C_VWN - B_VWN * B_VWN)
   
         arg1 = rs / xd_vwn 
         arg2 = q_vwn / (2*sqrs + B_VWN) 
         arg3 = (sqrs + X_VWN) * (sqrs + X_VWN) / xd_vwn 
   
         fac1 = 2*B_VWN / q_vwn 
         fac2 = -B_VWN * X_VWN / xd_vwn 
         fac3 = 2*(B_VWN - 2*X_VWN) / q_vwn 
!      
         epsc = log(arg1) + (fac1 - fac2 * fac3) * atan(arg2) -          & 
     &         fac2 * log(arg3) 
         epsc = A_VWN * epsc

         b1 = (B_VWN * X_VWN + C_VWN) / (X_VWN * C_VWN) 
         b2 = (X_VWN + B_VWN) / (X_VWN * C_VWN) 
         b3 = UM / (X_VWN * C_VWN) 

         rsdecdrs = A_VWN * (UM + b1*sqrs) /                             &
     &        (UM + b1*sqrs + b2*rs + b3*rs*sqrs) 
         vc = epsc - rsdecdrs / 3

       elseif ( author == 'wi' .or. author == 'WI') then

!        W i g n e r

!        wigner correlation energy

         epsc = A_WI / (rs + B_WI)
         rsdecdrs =A_WI / ((rs + B_WI)*(rs + B_WI))
         vc = epsc - rsdecdrs / 3

       elseif ( author == 'pw92' .or. author == 'PW92' ) then 

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

         epsc = - A_PW * (UM + ALP_PW*rs) *  log(c_pw) 
         rsdecdrs = -A_PW*rs * (ALP_PW * log(c_pw) +                     &
     &                         (UM + ALP_PW*rs) * dc_pwdrs / c_pw ) 
         vc = epsc - rsdecdrs / 3

       else 
         write(6,'("  STOPPED  in xc_lda:   unknown correlation")')
         write(6,*) '     ',author

         stop 

       endif

       return 
       end subroutine xc_lda
