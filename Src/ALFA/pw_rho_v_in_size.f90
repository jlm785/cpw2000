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

!>     This subroutines reads the file "filename" (default PW_RHO_V.DAT) and returns the
!>     array dimensions for future allocations.

       subroutine pw_rho_v_in_size(filename,io,                          &
     & mxdtyp,mxdatm,mxdgve,mxdnst,mxdlqp,mxdlao)

!      written January 12, 2014. jlm
!      modified documentation February 4, 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95

       implicit none

!       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len=*), intent(in)       ::  filename                   !<  name of file
       integer, intent(in)                ::  io                         !<  number of tape associated with filename.

!      output
       integer, intent(out)               ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(out)               ::  mxdatm                     !<  array dimension of types of atoms
       integer, intent(out)               ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(out)               ::  mxdnst                     !<  array dimension for g-space stars
       integer, intent(out)               ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(out)               ::  mxdlao                     !<  array dimension of orbital per atom type

!      local allocatable array

       integer,allocatable                ::  natom(:)                   !  number of atoms of type i

!      local variables

       integer      ::  ioerr 

!      counters

       integer    ::  i


       open(unit=io,file=trim(filename),status='old',form='UNFORMATTED')

       read(io,IOSTAT=ioerr) mxdtyp,mxdgve,mxdnst,mxdlqp,mxdlao
      
       if(ioerr /= 0) then
         backspace(io)
         read(io,IOSTAT=ioerr) mxdtyp,mxdgve,mxdnst,mxdlqp
         mxdlao = 4
       endif

       allocate(natom(mxdtyp))

       read(io) (natom(i),i=1,mxdtyp)

       mxdatm = natom(1)
       do i = 1,mxdtyp
         if(mxdatm < natom(i)) mxdatm = natom(i)
       enddo
       
       deallocate(natom)

       close(unit=io)

       return
       end subroutine pw_rho_v_in_size
