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

!>  Sets up the hamiltonian diagonal for a given k-point
!>  for a Kleinman-Bylander type pseudopotential
!>
!>  \author       Jose Luis Martins
!>  \version      4.94
!>  \date         8 February 1990,
!>  \copyright    GNU Public License v2

subroutine hamilt_diag_kb(mtxd, hdiag, qmod, ekpg, veffr1,               &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, adot,                                                &
      mxdtyp, mxdlqp, mxddim)

! written February 8 1990. jlm
! modified 24 March 1999. jlm
! modified June 3, 1999. jlm
! modified February 21, 2000. jlm
! modified 19 November 2013 (veffr1). jlm
! modified (f90) 14 January 2014.  jlm
! Split matrix_diag_kb into ham_struct and ham_diag_kb, 7 February 2014. JLM
! modified, vkb,nkb 16 April 2014.  JLM
! modified, normalization of vkb, 19 September 2015. JLM
! Modified documentation, January 2020. JLM
! Modified, indentation, 16 March 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  real(REAL64), intent(in)           ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  real(REAL64), intent(in)           ::  veffr1                          !<  Average value (veff(1)) of the effective potential.

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

! output

  real(REAL64), intent(out)          ::  hdiag(mxddim)                   !<  hamiltonian diagonal

! local varaibles

  integer         ::  ni
  real(real64)    ::  qi,xni,xi,vqnl
  real(real64)    ::  vq(0:3)
  real(REAL64)    ::  vcell, bdot(3,3),fac

! parameters
  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer   ::  i, nt, l


  call adot_to_bdot(adot, vcell, bdot)
  fac = UM / sqrt(vcell)

! calculates hamiltonian matrix diagonal

  do i = 1,mtxd

!   kinetic energy

    hdiag(i) = ekpg(i)

!   add local potential

    hdiag(i) = hdiag(i) + veffr1

!   add nonlocal potential

    qi = qmod(i)
    do nt = 1,ntype
      xni = qi/delqnl(nt)
      ni  = nint(xni)
      if(ni < 0) ni = 0
      do l = 0,3
        vq(l) = ZERO
        if (ni < nqnl(nt)-1) then
          xi  = xni - UM*ni
          if(nkb(l,0,nt) /= 0) then
            vq(l) =  vkb(ni,l,0,nt) * (UM+xi) * (UM-xi)                  &
                  + (vkb(ni+1,l,0,nt)*(UM+xi)                            &
                  -  vkb(ni-1,l,0,nt)*(UM-xi)) * xi / 2
            vq(l) = fac * vq(l) * qi**l
          endif
        endif
      enddo

      vqnl = ZERO
      do l = 0,3
        vqnl = vqnl + vq(l)*vq(l)*nkb(l,0,nt)
      enddo

      hdiag(i) = hdiag(i) + vqnl*natom(nt)

    enddo
  enddo

  return

end subroutine hamilt_diag_kb
