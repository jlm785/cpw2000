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

!>     adds a tag to a string

       subroutine add_tag(str,tagid,tag)

!      Written 8 August 2017.
!      Modified documentation august 2019.  JLM
!      Copyright José Luís Martins, INESC-MN.

!      version 4.94

       implicit none

!      input

       character(len=*), intent(in)       ::  tagid                      !<  identifier od the tag
       character(len=*), intent(in)       ::  tag                        !<  tag to be added

!      input and output

       character(len=*), intent(inout)    ::  str                        !<  string to be modified by adding tag

!      local variables

       integer   ::   ls, lslast
       integer   ::   ltagid, ltagidlast
       integer   ::   ltag, ltaglast
       integer   ::   lfree, l1, l2

!      counters
!      parameters

       character(len=1)     ::  tagchar = '#'
       ls = len(str)
       lslast = len(trim(str))

       ltagid = len(trim(adjustl(tagid)))
       ltagidlast = len(trim(tagid))

       ltag = len(trim(adjustl(tag)))
       ltaglast = len(trim(tag))

       lfree = ls - lslast
       if(lfree > ltag + ltagid + 3) then
         str(lslast+3:lslast+3) = tagchar
         l1 = lslast + 3
         l2 = ltagidlast - ltagid
         str(l1+1:l1+ltagid) = tagid(l2+1:l2+ltagid)
         l1 = lslast + ltagid + 4
         l2 = ltaglast - ltag
         str(l1+1:l1+ltag) = tag(l2+1:l2+ltag)
       endif

       return
       end subroutine add_tag
