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

!>     calculates the apropriate sizes for the
!>     most efficient fast fourier transform

       subroutine size_fft(kmscr,nsfft,mfft,mwrk)

!      written april 8 1988. jlm
!      modified june 7, 1993. jlm
!      modified August 6, 2012. jlm
!      Modified, documentation, January 2020. JLM

!      copyright INESC-MN/Jose Luis Martins
!
!      This file is distributed under the terms of the GNU General Public License.
!      See the file COPYING for license details.

!      version 4.94

       implicit none

!      input

       integer, intent(in)                ::  kmscr(3)                   !<  max value of kgv(i,n) used for the potential fft mesh

!      output

       integer, intent(out)               ::  nsfft(3)                   !<  fft dimensions in directions 1,2,3

       integer, intent(out)               ::  mfft                       !<  minimum size of fft data array
       integer, intent(out)               ::  mwrk                       !<  minimum size of fft work array

!      local data

       integer   ::    inprim(11)
       data inprim /2,4,6,8,10,12,16,18,20,24,32/
       integer   ::    insec(6)
       data insec /36,40,48,50,54,64/

       integer   ::  n,j,k

       do n=1,3
         if(2*kmscr(n)+1 < 32) then
           do j=1,11
             nsfft(n) = inprim(j)

             if(2*kmscr(n)+1 <= nsfft(n)) exit

           enddo
         else
           k = int(log(real(2*kmscr(n)+1))/log(real(2)))
           do j=1,6
             nsfft(n) = 2**(k-5) * insec(j)

             if(2*kmscr(n)+1 <= nsfft(n)) exit

           enddo   
         endif
       enddo

!      This values maybe changed depending on the fft package !!!!!!!!!!

       mfft = (nsfft(1)+1)*nsfft(2)*nsfft(3)
       mwrk = 4*max(nsfft(1)*nsfft(2),nsfft(1)*nsfft(3),                 &
     &                nsfft(2)*nsfft(3))

       return

       end subroutine size_fft
