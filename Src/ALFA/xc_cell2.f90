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

!>     Computes the exchange correlation potential and energy (Hartree)
!>     for the Tran-Blaha meta-GGA
!>     from the charge density (in electrons per unit cell), density laplacian 
!>     and "kinetic energy" density.
!>     Adapted from older LDA code.

       subroutine xc_cell2( author, tblaha, id1, id2, n1, n2, n3,        &
     &       chdr, taumsh, lapmsh, vxc, adot, exc, rhovxc, strxc )

!      Written 23 February 1999. jlm
!      Modified for MMGA Tran-Blaha. CLR
!      Modified for f90. 18 December 2016.  JLM
!      Modified, documentation, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins, Carlos Loia Reis

!      version 4.94

       implicit none 

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len = *), intent(in)     ::  author                     !<  type of xc wanted (ca=pz , pw92 , vwn, wi, pbe)
       real(REAL64), intent(in)           ::  tblaha                     !<  Tran-Blaha constant, if negative calculates it...
       integer, intent(in)                ::  id1, id2                   !<  first and second dimension of the fft array 
       integer, intent(in)                ::  n1, n2, n3                 !<  fft dimensions in directions 1,2,3
       real(REAL64), intent(in)           ::  chdr(id1,id2,n3)           !<  charge density (1/bohr^3)
       real(REAL64), intent(in)           ::  taumsh(id1,id2,n3)         !<  kinetic energy density (Hartree/bohr^3)
       real(REAL64), intent(in)           ::  lapmsh(id1,id2,n3)         !<  Laplacian of charge density (1/bohr^5)
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space (covariant components)

!      output

       real(REAL64), intent(out)          ::  vxc(id1,id2,n3)            !<  exchange-correlation potential vxc (Hartree).
       real(REAL64), intent(out)          ::  exc                        !<  the total exchange correlation energy given by the integral of the density times epsilon xc (Hartree).
       real(REAL64), intent(out)          ::  rhovxc                     !<  integral of the density times vxc (Hartree)
       real(REAL64), intent(out)          ::  strxc(3,3)                 !<  contribution of xc to the stress tensor (contravariant,a.u.)

!      local variables

       real(REAL64)        ::  bdot(3,3),vcell                           !  bdot/2*pi**2, cell volume
       real(REAL64)        ::  rho, epsx, epsc, vx, vc
       real(REAL64)        ::  strgga(3,3)                               !  contribution to stress
       real(REAL64)        ::  grho,dexdr,decdr,dexdgr,decdgr

       integer, parameter  ::  mxdnn = 3                                 !  Lagrange interpolation uses 2*mxdnn+1 points                   
       real(REAL64)        ::  dgdm(-mxdnn:mxdnn),drhodm(3),drhocon(3)

       real(REAL64)        ::  coef      

       integer, save       ::  xc_call = 0       
       real(REAL64)        ::  tau, lap, sigma
       real(REAL64)        ::  tb09_integral, tb09_const_c

!     learn to count in portuguese

!      parameters

       real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: EPS = 1.0d-12

!      TB parameters

       real(REAL64), parameter  :: TB09_ALPHA = -0.012_REAL64
       real(REAL64), parameter  :: TB09_BETA = 1.023_REAL64

!      counters

       integer       ::  i1, i2, i3
       integer       ::  i, j
       integer       ::  if1,if2
       integer       ::  in,jn,ip,nn

!      initial stuff

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

       if(n1 > id1 .or. n2 > id2) then
         write(6,'("   STOPPED in xc_cell:   wrong dimensions ",5i5)')   &
     &                n1,n2,n3,id1,id2

         stop

       endif

       call adot_to_bdot(adot,vcell,bdot)

!      correct fot the 2*PI

       do j = 1,3
       do i = 1,3
         bdot(i,j) = bdot(i,j) / (2*PI*2*PI)
       enddo
       enddo

       do i = 1,3
       do j = 1,3
         strgga(i,j) = ZERO
       enddo
       enddo

       if (author == 'pz' .or. author == 'PZ'                        &
     &     .or. author == 'ca' .or. author == 'CA'                   &
     &     .or. author == 'pw92' .or. author == 'PW92'               &
     &     .or. author == 'vwn' .or. author == 'VWN'                 &
     &     .or. author == 'wi' .or. author == 'WI'                   &
     &     .or. xc_call == 0) then

!       LDA exchange and correlation

         if(xc_call == 0) then

           do i3 = 1,n3
           do i2 = 1,n2
           do i1 = 1,n1
             rho = chdr(i1,i2,i3)

             call xc_lda('PZ', rho, epsx, epsc, vx, vc )

             exc = exc + rho * (epsx + epsc)
             vxc(i1,i2,i3) = vx + vc
             rhovxc = rhovxc + rho*(vx + vc)
           enddo
           enddo
           enddo
           xc_call = 1
           tb09_integral = ZERO

         else

           do i3 = 1,n3
           do i2 = 1,n2
           do i1 = 1,n1
             rho = chdr(i1,i2,i3)

             call xc_lda( author, rho, epsx, epsc, vx, vc )

             exc = exc + rho * (epsx + epsc)
             vxc(i1,i2,i3) = vx + vc
             rhovxc = rhovxc + rho*(vx + vc)
           enddo
           enddo
           enddo

         endif

       elseif (author == 'tbl' .or. author == 'TBL') then

!        meta-gga exchange and correlation


!        Weights for lagrange interpolation formula

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
         
!        calculation of tb09 constant c :(
              
         if(tblaha > ZERO) then

           tb09_const_c = tblaha

         else

           tb09_integral = ZERO

           do i3=1,n3
           do i2=1,n2
           do i1=1,n1
             rho   = chdr(i1,i2,i3)

!            calculates gradient of rho

             drhodm(1) = ZERO
             do in = -nn,nn
               ip = i1 + in
               ip = mod(ip+n1-1,n1) + 1
               drhodm(1) = drhodm(1) + dgdm(in)*chdr(ip,i2,i3)
             enddo
             drhodm(1) = n1*drhodm(1)

             drhodm(2) = ZERO
             do in = -nn,nn
               ip = i2 + in
               ip = mod(ip+n2-1,n2) + 1
               drhodm(2) = drhodm(2) + dgdm(in)*chdr(i1,ip,i3)
             enddo
             drhodm(2) = n2*drhodm(2)

             drhodm(3) = ZERO
             do in = -nn,nn
               ip = i3 + in
               ip = mod(ip+n3-1,n3) + 1
               drhodm(3) = drhodm(3) + dgdm(in)*chdr(i1,i2,ip)
             enddo
             drhodm(3) = n3*drhodm(3)

             drhocon(1) = bdot(1,1)*drhodm(1) + bdot(1,2)*drhodm(2) +      &
     &                    bdot(1,3)*drhodm(3)
             drhocon(2) = bdot(2,1)*drhodm(1) + bdot(2,2)*drhodm(2) +      &
     &                    bdot(2,3)*drhodm(3)
             drhocon(3) = bdot(3,1)*drhodm(1) + bdot(3,2)*drhodm(2) +      &
     &                    bdot(3,3)*drhodm(3)
             grho = drhodm(1)*drhocon(1) + drhodm(2)*drhocon(2) +          &
     &              drhodm(3)*drhocon(3)
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

!        end of calculation ot TB constant

         do i3=1,n3
         do i2=1,n2
         do i1=1,n1
           rho = chdr(i1,i2,i3)
           tau = taumsh(i1,i2,i3)
           lap = lapmsh(i1,i2,i3)

!          calculates gradient of rho again

           drhodm(1) = ZERO
           do in = -nn,nn
             ip = i1 + in
             ip = mod(ip+n1-1,n1) + 1
             drhodm(1) = drhodm(1) + dgdm(in)*chdr(ip,i2,i3)
           enddo
           drhodm(1) = n1*drhodm(1)

           drhodm(2) = ZERO
           do in = -nn,nn
             ip = i2 + in
             ip = mod(ip+n2-1,n2) + 1
             drhodm(2) = drhodm(2) + dgdm(in)*chdr(i1,ip,i3)
           enddo
           drhodm(2) = n2*drhodm(2)

           drhodm(3) = ZERO
           do in = -nn,nn
             ip = i3 + in
             ip = mod(ip+n3-1,n3) + 1
             drhodm(3) = drhodm(3) + dgdm(in)*chdr(i1,i2,ip)
           enddo
           drhodm(3) = n3*drhodm(3)

           drhocon(1) = bdot(1,1)*drhodm(1) + bdot(1,2)*drhodm(2) +      &
     &                  bdot(1,3)*drhodm(3)
           drhocon(2) = bdot(2,1)*drhodm(1) + bdot(2,2)*drhodm(2) +      &
     &                  bdot(2,3)*drhodm(3)
           drhocon(3) = bdot(3,1)*drhodm(1) + bdot(3,2)*drhodm(2) +      &
     &                  bdot(3,3)*drhodm(3)
           grho = drhodm(1)*drhocon(1) + drhodm(2)*drhocon(2) +          &
     &            drhodm(3)*drhocon(3)
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
                
           sigma= grho*grho

           call xc_mgga('tbl', rho, sigma, lap, tau,                     &
     &                          epsx, epsc, vx, vc, tb09_const_c ) 
           
           exc = exc + rho * (epsx + epsc)
           vxc(i1,i2,i3) = vx + vc
           rhovxc = rhovxc + rho*(vx + vc)

         enddo
         enddo
         enddo
 
       elseif (author == 'pbe' .or. author == 'PBE') then

!        gga exchange and correlation

         if(n1 < 3 .or. n2 < 3 .or. n3 < 3) then
           write(6,'("   STOPPED in xc_cell:   dimensions too small ",   &
     &            3i5)') n1,n2,n3

           stop

         endif

!        Weights for lagrange interpolation formula

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

         do i3 = 1,n3
         do i2 = 1,n2
         do i1 = 1,n1
           rho = chdr(i1,i2,i3)

!          calculates gradient of rho

           drhodm(1) = ZERO
           do in = -nn,nn
             ip = i1 + in
             ip = mod(ip+n1-1,n1) + 1
             drhodm(1) = drhodm(1) + dgdm(in)*chdr(ip,i2,i3)
           enddo
           drhodm(1) = n1*drhodm(1)

           drhodm(2) = ZERO
           do in = -nn,nn
             ip = i2 + in
             ip = mod(ip+n2-1,n2) + 1
             drhodm(2) = drhodm(2) + dgdm(in)*chdr(i1,ip,i3)
           enddo
           drhodm(2) = n2*drhodm(2)

           drhodm(3) = ZERO
           do in = -nn,nn
             ip = i3 + in
             ip = mod(ip+n3-1,n3) + 1
             drhodm(3) = drhodm(3) + dgdm(in)*chdr(i1,i2,ip)
           enddo
           drhodm(3) = n3*drhodm(3)

           drhocon(1) = bdot(1,1)*drhodm(1) + bdot(1,2)*drhodm(2) +      &
     &                  bdot(1,3)*drhodm(3)
           drhocon(2) = bdot(2,1)*drhodm(1) + bdot(2,2)*drhodm(2) +      &
     &                  bdot(2,3)*drhodm(3)
           drhocon(3) = bdot(3,1)*drhodm(1) + bdot(3,2)*drhodm(2) +      &
     &                  bdot(3,3)*drhodm(3)
           grho = drhodm(1)*drhocon(1) + drhodm(2)*drhocon(2) +          &
     &            drhodm(3)*drhocon(3)
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

           call xc_gga( author, rho, grho,                               &
     &                  epsx, epsc, dexdr, decdr, dexdgr, decdgr ) 

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
             vxc(ip,i2,i3) = vxc(ip,i2,i3) +                             &
     &                n1*(dexdgr + decdgr)*drhocon(1)*dgdm(in)
           enddo
           do in = -nn,nn
             ip = i2 + in
             ip = mod(ip+n2-1,n2) + 1
             vxc(i1,ip,i3) = vxc(i1,ip,i3) +                             &
     &                n2*(dexdgr + decdgr)*drhocon(2)*dgdm(in)
           enddo
           do in = -nn,nn
             ip = i3 + in
             ip = mod(ip+n3-1,n3) + 1
             vxc(i1,i2,ip) = vxc(i1,i2,ip) +                             &
     &                n3*(dexdgr + decdgr)*drhocon(3)*dgdm(in)
           enddo
         enddo
         enddo
         enddo

         do i3=1,n3
         do i2=1,n2
         do i1=1,n1
           rhovxc = rhovxc + chdr(i1,i2,i3)*vxc(i1,i2,i3)           
         enddo
         enddo
         enddo

       else

         write(6,'("    STOPPED in xc_cell:   unknown correlation")')
         write(6,*) '     ',author

         stop

       endif

!      scales by volume factor

       exc  = exc  * vcell / dble(n1*n2*n3)
       rhovxc  = rhovxc * vcell / dble(n1*n2*n3)
       do i=1,3
       do j=1,3
         strgga(i,j) = strgga(i,j) * vcell / dble(n1*n2*n3)
       enddo
       enddo

!      calculates stress

       do i=1,3
       do j=1,3
         strxc(i,j) = (rhovxc - exc) * bdot(j,i) + strgga(i,j)
       enddo
       enddo

       return
       end subroutine xc_cell2

