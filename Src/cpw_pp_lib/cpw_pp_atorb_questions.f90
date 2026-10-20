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

!>  If some atom type has several atomic basis sets,
!>  asks which one should be used.  The default is the first.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         20 October 2026.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_atorb_questions(ioreplay,                              &
     ntype, nameat, n_bsets, norbat, lorb,                               &
     jbset,                                                              &
     mxdtyp, mxdlao, mxdset)

! Written 20 October 2026 from cpw_pp_band_questions. claude

  implicit none

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type
  integer, intent(in)                ::  mxdset                          !<  array dimension of number of atomic basis sets

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  integer, intent(in)                ::  n_bsets(mxdtyp)                 !<  number of atomic basis sets for atom type nt
  integer, intent(in)                ::  norbat(mxdset,mxdtyp)           !<  number of atomic orbitals in basis set ib of atom type nt
  integer, intent(in)                ::  lorb(mxdlao,mxdset,mxdtyp)      !<  angular momentum of orbital n in basis set ib of atom type nt

! output

  integer, intent(out)               ::  jbset(mxdtyp)                   !<  basis set chosen for atom type nt

! local variables

  integer            ::  nsetmax                                         !  largest number of basis sets of an atom type
  integer            ::  ioerr
  character(len=1)   ::  yesno
  character(len=1)   ::  lname(0:3)                                      !  spectroscopic name of angular momentum

! counters

  integer            ::  nt, ib, n

  data lname / 's', 'p', 'd', 'f' /


! default is the first atomic basis set

  do nt = 1,mxdtyp
    jbset(nt) = 1
  enddo

  nsetmax = 1
  do nt = 1,ntype
    nsetmax = max(nsetmax,n_bsets(nt))
  enddo

  if(nsetmax > 1) then

    write(6,*)
    write(6,*) '  Some atoms have several atomic basis sets.'
    write(6,*) '  Do you want to use other than the first (y/n)?'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   choose atomic basis set'

    if(yesno == 'y' .or. yesno == 'Y') then

      do nt = 1,ntype

        if(n_bsets(nt) > 1) then

          write(6,*)
          write(6,'("   The atom ",a2," (type ",i3,") has ",i3,          &
              &     " atomic basis sets:")') nameat(nt), nt, n_bsets(nt)
          write(6,*)
          do ib = 1,n_bsets(nt)
            write(6,'(6x,i3,":",i5," orbitals   ",40(1x,a1))')           &
                 ib, norbat(ib,nt),                                      &
                 (lname(min(max(lorb(n,ib,nt),0),3)),n=1,norbat(ib,nt))
          enddo
          write(6,*)
          write(6,'("   Enter the basis set for ",a2)') nameat(nt)

          read(5,*,iostat=ioerr) jbset(nt)

          if(ioerr /= 0) then
            jbset(nt) = 1
            write(6,*)
            write(6,*) '  Error reading the basis set, using the first'
          endif

          write(ioreplay,*) jbset(nt),'   atomic basis set for ',nameat(nt)

          if(jbset(nt) < 1 .or. jbset(nt) > n_bsets(nt)) then
            jbset(nt) = min(max(jbset(nt),1),n_bsets(nt))
            write(6,*)
            write(6,'("   Basis set not available, will use the ",       &
                &     "nearest, ",i3)') jbset(nt)
          endif

          write(6,*)
          write(6,'("   Will use basis set ",i3," for ",a2)')            &
               jbset(nt), nameat(nt)

        endif

      enddo

    endif

  endif

  return

end subroutine cpw_pp_atorb_questions
