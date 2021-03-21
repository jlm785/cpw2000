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

!>     interface to the molecular dynamics/minimization steppers.
!>     updates atomic positions and cell metric.

       subroutine move_atom_cell(flgcal,newcel,iconv,istmd,              &
     & rat,vat,adot,vadot,energy,force,stress,                           &
     & tstep,ekin,ekcell,strext,press,celmas,                            &
     & beta,tempk,iseed,                                                 &
     & mkeep,pgtol,dxmax,                                                &
     & rat1,frc1,adot1,frcel1,                                           &
     & ntype,natom,nameat,atmass,                                        &
     & ntrans,mtrx,tnp,                                                  &
     & mxdtyp,mxdatm)

!      Written in January 2017. JLM
!      Modified (newcel) 24 February 2019. JLM
!      added move_epi_langevin, August 2019. JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       character(len=6), intent(in)       ::  flgcal

       real(REAL64), intent(in)           ::  energy                     !<  energy in hartree
       real(REAL64), intent(in)           ::  force(3,mxdatm,mxdtyp)     !<  d energy / d rat,  force on the n-th atom of type i (contravarian components, hartree/bohr)
       real(REAL64), intent(in)           ::  stress(3,3)                !<  d energy / d adot,  stress tensor (contravariant)

       integer,intent(in)                 ::  istmd                      !<  md step. Equal to 1 in first step of molecular dynamics
       real(REAL64), intent(in)           ::  tstep                      !<  molecular dynamics time step (in a.u.)

       real(REAL64), intent(in)           ::  tempk                      !<  ionic temperature (in Kelvin)
       real(REAL64), intent(in)           ::  beta                       !<  friction coefficient/mass (in a.u.)

       integer, intent(in)                ::  mkeep                      !<  number of iterations remembered by lbfgs
       real(REAL64), intent(in)           ::  pgtol                      !<  tolerance for minimization convergence
       real(REAL64), intent(in)           ::  dxmax                      !<  maximum step size

       real(REAL64), intent(in)           ::  strext(3,3)                !<  external applied stress
       real(REAL64), intent(in)           ::  press                      !<  external pressure
       real(real64), intent(in)           ::  celmas                     !<  fictitious cell mass       

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(in)           ::  atmass(mxdtyp)             !<  atomic mass of atoms of type i

       integer, intent(in)                ::  ntrans                     !<  number of symmetry operations in the factor group
       integer, intent(in)                ::  mtrx(3,3,48)               !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
       real(REAL64), intent(in)           ::  tnp(3,48)                  !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

!      input and output
       
       real(REAL64), intent(inout)        ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  vat(3,mxdatm,mxdtyp)       !<  d rat / d t  velocity in lattice coordinates of atom j of type i

       real(REAL64), intent(inout)        ::  adot(3,3)                  !<  metric in real space
       real(REAL64), intent(inout)        ::  vadot(3,3)                 !<  d adot / d t  rate of change of metric

       integer, intent(inout)             ::  iseed                      !<  seed for rndom number generator

       real(REAL64), intent(inout)        ::  rat1(3,mxdatm,mxdtyp)      !<  previous value of lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  frc1(3,mxdatm,mxdtyp)      !<  previous value of force on atom j of type i
       real(REAL64), intent(inout)        ::  adot1(3,3)                 !<  previous value of adot
       real(REAL64), intent(inout)        ::  frcel1(3,3)                !<  previous cell "force" (covariant components)

!      output

       real(REAL64), intent(out)          ::  ekin                       !<  kinetic energy of the nuclei (atomic units)
       real(REAL64), intent(out)          ::  ekcell                     !<  kinetic energy in Hartree of cell (fictitious)C

       logical, intent(out)               ::  newcel                     !<  indicates that adot has changed on output
       integer, intent(out)               ::  iconv                      !<  iconv = 1, lbfgs converged of flgcal = one

!      local variables

       integer               ::  minrat                     !  if =1 minimize with respect to atomic positions
       integer               ::  minstr                     !  if =1 minimize with respect to all adot variables, if =2 minimize with respect to adot(3,3)
       integer, save         ::  ilbfgs = 0                 !  minimization step in lbfgs           





       if(flgcal == 'ONE   ') then
         iconv = 1
       else

         iconv = 0
         ekcell = 0.0

         if(flgcal == 'MICRO ') then
       
           call move_md_micro(rat,vat,force,istmd,tstep,ekin,            &
     &     rat1,frc1,                                                    &
     &     ntype,natom,atmass,adot,                                      &
     &     mxdatm,mxdtyp)

           newcel = .FALSE.

         elseif(flgcal == 'LANG  ') then

           call move_md_langevin(rat,vat,force,istmd,tstep,ekin,         &
     &     beta,tempk,iseed,                                             &
     &     rat1,frc1,                                                    &
     &     ntype,natom,atmass,adot,                                      &
     &     mxdatm,mxdtyp)

           newcel = .FALSE.

         elseif(flgcal == 'LBFSYM' .or. flgcal == 'VCSLBF' .or.          &
     &          flgcal == 'EPILBF') then

           minrat = 1
           if(flgcal == 'LBFSYM') then
             minstr = 0
             newcel = .FALSE.
           elseif(flgcal == 'VCSLBF') then
             minstr = 1
             newcel = .TRUE.
           elseif(flgcal == 'EPILBF') then
             minstr = 2
             newcel = .TRUE.
           endif
           
!           if(istmd == ist0 + 1) ilbfgs = 0 

           call move_lbfgs_call(energy,rat,force,adot,stress,               &
     &     strext,press,ilbfgs,iconv,minrat,minstr,                         &
     &     ntype,natom,                                                     &
     &     mkeep,pgtol,dxmax,                                               &
     &     mxdtyp,mxdatm)
     
           if(iconv < 0) then
             ilbfgs = 0
           else
             ilbfgs = ilbfgs + 1
           endif

!          comment the next call if you want to allow symmetry breaking
!          by noise

           call sym_rat(ntype,natom,rat,ntrans,mtrx,tnp,                 &
     &     mxdtyp,mxdatm)

         elseif(flgcal == 'VCSMIC') then

           call move_vcs_micro(rat,vat,adot,vadot,force,stress,          &
     &     istmd,tstep,ekin,ekcell,strext,press,celmas,                  &
     &     rat1,frc1,adot1,frcel1,                                       &
     &     ntype,natom,atmass,                                           &
     &     mxdatm,mxdtyp)

           newcel = .TRUE.

         elseif(flgcal == 'VCSLNG') then

           call move_vcs_langevin(rat,vat,adot,vadot,force,stress,       &
     &     istmd,tstep,ekin,ekcell,strext,press,celmas,                  &
     &     beta,tempk,iseed,                                             &
     &     rat1,frc1,adot1,frcel1,                                       &
     &     ntype,natom,atmass,                                           &
     &     mxdatm,mxdtyp)

           newcel = .TRUE.

         elseif(flgcal == 'EPILNG') then

           call move_epi_langevin(rat,vat,adot,vadot,force,stress,       &
     &     istmd,tstep,ekin,ekcell,strext,press,celmas,                  &
     &     beta,tempk,iseed,                                             &
     &     rat1,frc1,adot1,frcel1,                                       &
     &     ntype,natom,atmass,                                           &
     &     mxdatm,mxdtyp)

           newcel = .TRUE.

         endif

!        prints velocities

         call move_print_velocity(flgcal,                                &
     &   adot,ntype,natom,nameat,vat,vadot,                              &
     &   mxdtyp,mxdatm)

       endif

       return
       end subroutine move_atom_cell
