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

!>    chrcap accepts a string of nchar characters and replaces
!>    any lowercase letters by uppercase ones.

      subroutine chrcap(string,nchar)
      
!     Modified documentation August 2019.  JLM

!     version 4.94

      implicit none

      integer, intent(in)                ::  nchar                       !<  string length
      character(len=*), intent(inout)    ::  string                      !<  string to be uppercased

      integer       ::  ncopy, i, itemp

      ncopy = nchar
      if(ncopy < 1) ncopy = len(string)
      do i=1,ncopy

        if(LGE(STRING(I:I),'a') .AND. LLE(STRING(I:I),'z'))then
          ITEMP=ICHAR(STRING(I:I))+ICHAR('A')-ICHAR('a')
          STRING(I:I)=CHAR(ITEMP)
          endif
      enddo
      
      return
      end subroutine chrcap

