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

!>  Calculates the total Kohn-Sham energy,
!>  prints out potentials and charge densities,
!>  and checks for convergence.

subroutine total_ks_energy(ipr, icmplx, iter, itmix, eharrfou,           &
  iconv, errvhxc, epscv,                                                 &
  energy, eband, ektot, exc, ealpha, enerew,                             &
  ng, kgv, ns, mstar, ek,                                                &
  vion, vhxc, vhxcout, den,                                              &
  ztot, adot,                                                            &
  mxdgve, mxdnst)

! Adapted from Sverre Froyen plane wave program.
! written june 8 1987.jlm
! version 4.0. 14 october 93. jlm
! modified 24 march 1999
! modified 6 sept 2002. jlm
! modified 28 sept 2013 (average potential). jlm
! modified, f90, complex*16, 24 September 2015. JLM
! Modified documentation, August 2019. JLM
! Modified, itmix,errvhxc,eharrfou,printing, 2 December 2020. JLM
! copyright INESC-MN/Jose Luis Martins

! Version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  ipr                             !<  print switch
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex
  integer, intent(in)                ::  iter                            !<  iteration number.
  integer, intent(in)                ::  itmix                           !<  iteration for the mixer (itmix = iter in normal cases)

  real(REAL64), intent(in)           ::  eharrfou                        !<  total energy of the Harris-Weinert-Foulkes functional
  real(REAL64), intent(in)           ::  epscv                           !<  convergence criteria

  real(REAL64), intent(in)           ::  eband                           !<  integrated band energy.
  real(REAL64), intent(in)           ::  ektot                           !<  total kinetic energy of electrons
  real(REAL64), intent(in)           ::  exc                             !<  exchange+correlation energy
  real(REAL64), intent(in)           ::  ealpha                          !<  alpha term. (G=0)
  real(REAL64), intent(in)           ::  enerew                          !<  Ewald energy

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  complex(REAL64), intent(in)        ::  vion(mxdnst)                    !<  ionic potential for the prototype G-vector in star j
  complex(REAL64), intent(in)        ::  vhxc(mxdnst)                    !<  input hxc potential (hartree) for the prototype G-vector
  complex(REAL64), intent(in)        ::  vhxcout(mxdnst)                 !<  output hxc potential (hartree) for the prototype G-vector
  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  density for the prototype G-vector

  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

! output

  real(REAL64), intent(out)          ::  energy                          !<  total electronic(+ewald) energy.
  integer, intent(out)               ::  iconv                           !<  set equal to 1 if selfconsistency has been achieved.
  real(REAL64), intent(out)          ::  errvhxc                         !<  maximum value of abs(vhxcout(i) - vhxc(i))

! local variables

  real(REAL64)        ::  evhxc,epsloc,epsnl,ehart
  real(REAL64), save  ::  oevhxc,oepsloc,oepsnl,oehart
  real(REAL64), save  ::  oeband,oektot,oexc
  real(REAL64), save  ::  oenergy, oeharrfou
  real(REAL64)        ::  vcell,bdot(3,3)
  real(REAL64)        ::  dvhxc

  integer             ::  iprint,maxprt,ivec

! counters

  integer       ::  i, j

! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: EPS = 0.00001_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  call adot_to_bdot(adot,vcell,bdot)

! check for convergence.

  iconv = 1
  errvhxc = ZERO
  do i=2,ns
    dvhxc = abs(vhxcout(i) - vhxc(i))
    if(dvhxc > errvhxc) errvhxc = dvhxc
  enddo
  if(errvhxc > epscv) iconv = 0
  if(iter == 1) iconv = 0

! determine what printing should be done.

  iprint = 0
  if (ipr > 0) iprint = 2
  if (ipr == 0 .and. iconv == 1) iprint = 1
  if (ipr < 0 .and. iconv == 0) iprint = 1
  if (ipr < 0 .and. iconv == 1) iprint = 2
  maxprt = min(iabs(ipr),ns)
  if (iprint == 0) maxprt = 0

! printout heading and charge density.

  if (iprint == 2 .and. maxprt > 1) then
    write(6,*)
    write(6,'("   iteration number",i3)') iter
    write(6,*)
    write(6,*)
    write(6,'("    i   k-prot",9x,"Ek",7x,"den",5x,"V(out)",             &
            4x,"V(in)",5x,"delta V",3x,"Vionic")')
    write(6,*)
    write(6,'(1x,i4,1x,3i3,2x,3f10.5)') 1,(kgv(j,1),j=1,3),              &
                    ek(1),real(den(1)),real(vhxc(1))
    if(icmplx /= 0) then
      ivec = 2
      do i = 2,maxprt
        write(6,'(1x,i4,1x,3i3,2x,6f10.5)')                              &
             i,(kgv(j,ivec),j=1,3),ek(i),                                &
             real(den(i)),real(vhxcout(i)),real(vhxc(i)),                &
             real(vhxcout(i)-vhxc(i)),real(vion(i))
        write(6,'(27x,5f10.5)')                                          &
             aimag(den(i)),aimag(vhxcout(i)),aimag(vhxc(i)),             &
             aimag(vhxcout(i)-vhxc(i)),aimag(vion(i))
        ivec = ivec + mstar(i)
      enddo
    else
      ivec = 2
      do i = 2,maxprt
        write(6,'(1x,i4,1x,3i3,2x,6f10.5)')                              &
             i,(kgv(j,ivec),j=1,3),ek(i),                                &
             real(den(i)),real(vhxcout(i)),real(vhxc(i)),                &
             real(vhxcout(i)-vhxc(i)),real(vion(i))
        ivec = ivec + mstar(i)
      enddo
    endif
  endif

! initialize hxc-correction, local pseudo, hartree (g>0) energies
! includes correction for average potential

  evhxc = real(vhxc(1),REAL64)*real(den(1),REAL64)
  epsloc = real(vion(1),REAL64)*real(den(1),REAL64)
  ehart = ZERO

! start loop over stars

  do i=2,ns

!   Compute the hxc correction from the potential of the
!   valence electrons. This will be subtracted from the sum
!   of the eigenvalues.
!   Then compute the contribution to the total energy from
!   the interaction with the ions and the electrostatic
!   (Hartree) interaction between the electrons.

    evhxc = evhxc + mstar(i)*real(vhxc(i)*conjg(den(i)),REAL64)
    epsloc = epsloc + mstar(i)*real(vion(i)*conjg(den(i)),REAL64)
    ehart = ehart + mstar(i)*real(den(i)*conjg(den(i)),REAL64) / (2*ek(i))

  enddo

  ehart = 4*PI*ehart/(2*vcell)

! find the nonlocal energy.

  epsnl = eband - (evhxc + ektot + epsloc)

! the total energy. modified for average potential

  energy = eband - evhxc + ehart + exc + enerew

! printout

  if (iprint == 2) then
    if(itmix == 1) then
      write(6,*)
      write(6,*)
      write(6,*)
      if(iter == 1) then
        write(6,'(11x,"Iteration number",i3," :")') iter
      else
        write(6,'(11x,"Iteration number",i3," , mixer iteration, i3")') iter, itmix
      endif
      write(6,'(11x,"---------------------")')
      write(6,*)
      write(6,'(11x,"Alpha Term        =",2x,f15.6)') ealpha
      write(6,'(11x,"Kinetic  Energy   =",2x,f15.6)') ektot
      write(6,'(11x,"Local PP Energy   =",2x,f15.6)') epsloc
      write(6,'(11x,"Nonlocal Energy   =",2x,f15.6)') epsnl
      write(6,'(11x,"-------------------------------------")')
      write(6,'(11x,"Eigenvalue Sum    =",2x,f15.6)') eband
      write(6,'(11x,"HXC  Correction   =",2x,f15.6)') -evhxc
      write(6,'(11x,"Hartree  Energy   =",2x,f15.6)') ehart
      write(6,'(11x,"XC       Energy   =",2x,f15.6)') exc
      write(6,'(11x,"Ewald    Energy   =",2x,f15.6)') enerew
      write(6,'(11x,"-------------------------------------")')
      write(6,'(11x,"Total    Energy   =",2x,f15.6)') energy
      write(6,*)
      write(6,*)
    else
      write(6,*)
      write(6,*)
      write(6,*)
      write(6,'(11x,"Iteration number",i3,4x,"Energies (Ha)",10x,"Changes")') iter
      write(6,'(11x,"---------------------")')
      write(6,*)
      write(6,'(11x,"Alpha Term        =",2(2x,f15.6))') ealpha
      write(6,'(11x,"Kinetic  Energy   =",2(2x,f15.6))') ektot, ektot-oektot
      write(6,'(11x,"Local PP Energy   =",2(2x,f15.6))') epsloc, epsloc-oepsloc
      write(6,'(11x,"Nonlocal Energy   =",2(2x,f15.6))') epsnl, epsnl-oepsnl
      write(6,'(11x,"-------------------------------------")')
      if(itmix /= 2) then
        write(6,'(11x,"Harris-Foulkes    =",2(2x,f15.6))') eharrfou, eharrfou-oeharrfou
      else
        write(6,'(11x,"Harris-Foulkes    =",2(2x,f15.6))') eharrfou
      endif
      write(6,'(11x,"-------------------------------------")')
      write(6,'(11x,"Eigenvalue Sum    =",2(2x,f15.6))') eband, eband-oeband
      write(6,'(11x,"HXC  Correction   =",2(2x,f15.6))') -evhxc, -(evhxc-oevhxc)
      write(6,'(11x,"Hartree  Energy   =",2(2x,f15.6))') ehart, ehart-oehart
      write(6,'(11x,"XC       Energy   =",2(2x,f15.6))') exc, exc-oexc
      write(6,'(11x,"Ewald    Energy   =",2(2x,f15.6))') enerew
      write(6,'(11x,"-------------------------------------")')
      write(6,'(11x,"Total    Energy   =",2(2x,f15.6))') energy, energy-oenergy
      write(6,*)
      write(6,*)
    endif
  endif

  if (iprint == 1) then
    write(6,*)
    write(6,'(" total energy =",f17.10)') energy
  endif
  if (iprint > 0 .and. iconv == 1) then
    write(6,*)
    write(6,'(2x,f20.10,"    total kohn-sham energy")') energy
    write(6,*)
  endif

! keeps old energies

  oeband = eband
  oevhxc = evhxc
  oektot = ektot
  oepsloc = epsloc
  oepsnl = epsnl
  oexc = exc
  oehart = ehart
  oenergy = energy
  if(itmix /= 1) oeharrfou = eharrfou

! final checks

  if (abs(real(den(1),REAL64) - ztot) > eps) then
    write(6,*)
    write(6,'("     WARNING  total_ks_energy: non-neutral system",  &
         " charge inbalance ",e12.4)') real(den(1)) - ztot
    write(6,*)
  endif
  if (abs(aimag(den(1))) > eps) then
    write(6,*)
    write(6,'("     WARNING  total_ks_energy: strange system ",     &
         "complex charge ",e12.4)') abs(aimag(den(1)))
    write(6,*)
  endif

  return
end subroutine total_ks_energy
