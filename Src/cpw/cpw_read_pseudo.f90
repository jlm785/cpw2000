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

!>  Interface subroutine for read_pseudo
!>
!>  \author       Jose Luis Martins
!>  \version      5.10
!>  \date         19 November 2019. 13 January 2024.
!>  \copyright    GNU Public License v2

subroutine cpw_read_pseudo(iprglob, author,                              &
       crys_, pseudo_, atorb_, dims_)

! Written 19 November 2019. JLM
! Modifierd indentation, author, 13 January 2024. JLM

  use cpw_variables

  implicit none

  type(dims_t)                       ::  dims_                           !<  array dimensions
  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space

  integer,intent(in)                 ::  iprglob                         !<  level of detail of printout
  character(len=*),intent(in)        ::  author                          !<  eXchange-Correlation choice

  integer              ::  ipr

! counters

  integer              ::  nt, l

  ipr = 0
  if(iprglob > 0) ipr = 1

  call size_mxdlqp_lao(crys_%ntype, crys_%nameat,                        &
         dims_%mxdtyp, dims_%mxdlqp, dims_%mxdlao)

! this is for compatibility between old and new version

  dims_%mxdlqp = dims_%mxdlqp + 10

  allocate(pseudo_%nq(dims_%mxdtyp))
  allocate(pseudo_%delq(dims_%mxdtyp))
  allocate(pseudo_%zv(dims_%mxdtyp))
  allocate(pseudo_%nkb(0:3,-1:1,dims_%mxdtyp))
  allocate(pseudo_%vloc(-1:dims_%mxdlqp,dims_%mxdtyp))
  allocate(pseudo_%dcor(-1:dims_%mxdlqp,dims_%mxdtyp))
  allocate(pseudo_%dval(-1:dims_%mxdlqp,dims_%mxdtyp))

  allocate(pseudo_%vkb(-2:dims_%mxdlqp,0:3,-1:1,dims_%mxdtyp))

  allocate(atorb_%norbat(dims_%mxdtyp))
  allocate(atorb_%nqwf(dims_%mxdtyp))
  allocate(atorb_%delqwf(dims_%mxdtyp))
  allocate(atorb_%wvfao(-2:dims_%mxdlqp,dims_%mxdlao,dims_%mxdtyp))
  allocate(atorb_%lorb(dims_%mxdlao,dims_%mxdtyp))

  call read_pseudo(ipr, author,                                          &
       pseudo_%ealraw, PSEUDO_%NQ, PSEUDO_%DELQ, pseudo_%vkb,            &
       pseudo_%nkb,pseudo_%vloc, pseudo_%dcor, pseudo_%dval,             &
       atorb_%norbat, atorb_%nqwf, atorb_%delqwf, atorb_%wvfao,          &
       atorb_%lorb, atorb_%latorb,                                       &
       crys_%ntype, crys_%natom, crys_%nameat,                           &
       pseudo_%zv, pseudo_%ztot,                                         &
       dims_%mxdtyp, dims_%mxdlqp, dims_%mxdlao)


  do nt = 1,crys_%ntype
  do l = 0,3
    pseudo_%nkb(l,-1,nt) = 0
    pseudo_%nkb(l, 1,nt) = 0
  enddo
  enddo


  return

end subroutine cpw_read_pseudo
