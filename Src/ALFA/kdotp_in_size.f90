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

!>     Reads the file with the data to interface with kdotp
!>     and reports the dimension of the hamiltonian

       subroutine kdotp_in_size(filename,iotape,mxdbnd)


!      Written May 3, 2014, from the writing subroutine. jlm
!      Modified, name, documnetation, 15 June 20120. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.98

       implicit none

!       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len=*), intent(in)       ::  filename                   !<  file to be written
       integer, intent(in)                ::  iotape                     !<  tape number 

!      output

       integer, intent(out)               ::  mxdbnd                     !<  array dimension for the number of bands


       open(unit=iotape,file=filename,form='formatted',status='old')

       read(iotape,*) mxdbnd


       close(unit=iotape)

       return

       end subroutine kdotp_in_size
