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

!>     calculates a mgga exchange and correlations
!>     it can interface with the libxc by uncommenting a few lines.

       subroutine xc_mgga(author, rho, sigma, lap_rho, tau,              &
     &                    epsx, epsc, vx, vc, ctb09 ) 

!       Written by Carlos Loia Reis, September 2015.
!      Modified, documentation, December 2019. JLM
!      Copyright Carlos Loia Reis INESC-MN

!       uncomment to use libxc 1.2.0       
!       use XcTlk
!       use MXcAux, only :tblaha_constant

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len = *), intent(in)     ::  author                     !<  type of xc wanted (ca=pz , pw92 , vwn, wi)
       real(REAL64), intent(in)           ::  rho                        !<  electron density (1/bohr^3)
       real(REAL64), intent(in)           ::  sigma                      !<  gradient of charge density
       real(REAL64), intent(in)           ::  lap_rho                    !<  Laplacian of charge density
       real(REAL64), intent(in)           ::  tau                        !<  kinetic energy density
       real(REAL64), intent(in)           ::  ctb09                      !<  TB constant

!      output

       real(REAL64), intent(out)          ::  epsx                       !<  exchange energy density (hartree/bohr^3) 
       real(REAL64), intent(out)          ::  epsc                       !<  correlation energy density (hartree/bohr^3)
       real(REAL64), intent(out)          ::  vx                         !<  exchange potential (hartree/bohr^3)
       real(REAL64), intent(out)          ::  vc                         !<  correlation potential (hartree/bohr^3)

       
       real(REAL64)           :: vsigma,vlap_rho, vtau
       real(REAL64)           :: set_ctb09

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64

!       set_ctb09 = tblaha_constant
        set_ctb09 = ctb09

!       uncomment the following 3 lines to use libxc 1.2.0       
!       call xc_f90_mgga_x_tb09_set_par(x_func,set_ctb09)       
!       call xc_f90_lda_vxc(c_func, 1, rho, vc)
!       call xc_f90_mgga_vxc(x_func, 1, rho, sigma,lap_rho,tau, vx, vsigma, vlap_rho, vtau)

!      comment the following 2 lines to use libxc 1.2.0              
       call xc_lda( 'pz', rho, epsx, epsc, vx, vc )        
       call xc_cpw_tran_blaha(rho, sigma, lap_rho, tau, set_ctb09,       &
     & vx, vsigma, vtau, vlap_rho)       
       
       
       epsx = ZERO
       epsc = ZERO

       end subroutine xc_mgga
