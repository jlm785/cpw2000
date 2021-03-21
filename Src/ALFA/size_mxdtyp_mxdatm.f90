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

!>     finds the array sizes for the atomic structure, PW.DAT version

       subroutine size_mxdtyp_mxdatm(filename,iopw,mxdtyp,mxdatm)

!      written 2 May 2014. JLM
!      Modified, documentation, ioerr, January 2020. JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none

!      input


       character(len=*), intent(in)       ::  filename                   !<  file to be written
       integer, intent(in)                ::  iopw                       !<  tape number 

!      output

       integer, intent(out)               ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(out)               ::  mxdatm                     !<  array dimension of types of atoms

!      local variables

       integer     ::  ntype, natom
       integer     ::  ioerr

!      counters

       integer     ::  nt, j



       open(unit=iopw,file=filename,form='formatted',status='old',       &
     &       IOSTAT=ioerr)


       if(ioerr /= 0) then
         write(6,*)
         write(6,*)
         write(6,*) '  stopped in size_mxdtyp_mxdatm'
         write(6,*) '  Unable to open ',filename,' file'

         stop

       endif

!      reads title

       read(iopw,*)
       
!      read scale factor

       read(iopw,*)

!      read the basis vectors

       read(iopw,*)
       read(iopw,*)
       read(iopw,*)

!      read the number of atoms of each type and their names.

       read(iopw,*) ntype
       mxdtyp = ntype
       mxdatm = 0
       do nt=1,ntype
         read(iopw,*) natom
         mxdatm = max(mxdatm,natom)
         do j=1,natom
           read(iopw,*)
         enddo
       enddo

       close(unit=iopw)

       return

       end subroutine size_mxdtyp_mxdatm


