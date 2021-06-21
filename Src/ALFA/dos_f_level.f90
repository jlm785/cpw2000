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

!>  finds the fermi level and intrinsic carrier concentration
!>  and prints them to the standard output

subroutine dos_f_level(tempk, nel, nhist, ehist, dhist, chist, lidos,    &
      efguess, vcell, nvbm, ncbm, xmu, xmdosv, xmdosc, xni)

! written 5 December 2013.
! copyright Jose Luis Martins / INESC-MN
! Modified, documentation, 19 September 2020. JLM
! Modified, Fermi level and effective masses. 20 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)         ::  tempk                             !<  temperature in K
  integer, intent(in)              ::  nel                               !<  number of electrons in unit cell
  integer, intent(in)              ::  nhist                             !<  size of energy array
  real(REAL64), intent(in)         ::  ehist(nhist)                      !<  energy array (Hartree)
  real(REAL64), intent(in)         ::  dhist(nhist)                      !<  density of states array (electron/cell/Hartree)
  real(REAL64), intent(in)         ::  chist(nhist)                      !<  integrated density of states array (electron/cell)
  logical, intent(in)              ::  lidos                             !<  true if integrated density of states is to be computed.
  real(REAL64), intent(in)         ::  efguess                           !<  estimate of the Fermi energy
  real(REAL64), intent(in)         ::  vcell                             !<  unit cell volume (in atomic units)

! output

  integer, intent(out)             ::  nvbm, ncbm                        !<  approximate indices for conduction and valence band
  real(REAL64), intent(out)        ::  xmu                               !<  electron chemical potential, AKA Fermi level
  real(REAL64), intent(out)        ::  xmdosv, xmdosc                    !<  dos effective masses for valence and conduction
  real(REAL64), intent(out)        ::  xni                               !<  intrinsic carrier concentration

! local variables

  real(REAL64)    :: tau                                                 !  temperature in atomic units
  real(REAL64)    :: xn,xp                                               !  concentration of electrons and holes (per unit cell)
  real(REAL64)    :: dxndef,dxpdef                                       !  d xn / d ef    and    d xp / d ef
  real(REAL64)    :: xl,xh                                               !  brackets the Fermi level
  real(REAL64)    :: ef, efold                                           !  Fermi energy
  integer         :: maxit, nave
  real(REAL64)    :: dx, dxold                                           !  "xh-xl"
  real(REAL64)    :: f, df                                               !  xn - xp,    d f / d ef
  real(REAL64)    :: del                                                 !  precision

! counters

  integer               ::  n, m

! constants

  real(REAL64), parameter :: EPS = 0.000001
  real(REAL64), parameter :: ZERO = 0.0_REAL64
  real(REAL64), parameter :: HARTREE = 27.21138386_REAL64
  real(REAL64), parameter :: BOHR = 0.5291772109
  real(REAL64), parameter :: TAUTOK = HARTREE*11604.9_REAL64

  
  tau = tempk / TAUTOK

  xmu = ZERO
  xmdosv = ZERO
  xmdosc = ZERO
  xni = ZERO

  if(lidos) then

!   finds low temperature Fermi level

    m = 1
    do n = 1,nhist
      if(chist(n) > (nel - EPS)) exit
      m = n
    enddo

!   finds index of valence band maximum and conduction band minimum

    do n = m,nhist
      nvbm = n
      if(dhist(n) < EPS) exit
    enddo
    m = nvbm
    do n = nvbm,nhist
      if(chist(n) > (nel + EPS)) exit
      m = n
    enddo
    do n = m,1,-1
      ncbm = n
      if(dhist(n) < EPS) exit
    enddo

  else

!   finds guess in histogram

    m = 1
    do n = 1,nhist
      if(ehist(n) > efguess) exit
      m = n
    enddo

! finds nearest gap
   
    do n = 0,min(nhist-m,m-1)
    
      if(dhist(m-n) < EPS) then
        m = m-n
        exit
      endif

      if(dhist(m+n) < EPS) then
        m = m+n
        exit
      endif
    
    enddo

    ncbm = m
    do n = m,nhist
      if(dhist(n) > EPS) exit
      ncbm = n
    enddo
    nvbm = m
    do n = m,1,-1
      if(dhist(n) > EPS) exit
      nvbm = n
    enddo

  endif

! traps errors if number of electrons falls on a band

  if(ncbm <= nvbm) then
    write(6,*)
    write(6,*) '  ERROR in f_level'
    write(6,*)
    write(6,*) '  nvbm,ncbm,nhist  ',nvbm,ncbm,nhist
    write(6,*) '  ehist(nvbm),ehist(ncbm) ',ehist(nvbm),ehist(ncbm)
    write(6,*) '  dhist(nvbm),dhist(ncbm) ',dhist(nvbm),dhist(ncbm)
    write(6,*) '  chist(nvbm),chist(ncbm) ',chist(nvbm),chist(ncbm)

    stop

  endif

  del = max(ehist(nvbm+1)-ehist(nvbm),ehist(nvbm)-ehist(nvbm-1),        &
            ehist(ncbm+1)-ehist(ncbm),ehist(ncbm)-ehist(ncbm-1) )

! checks if the the calculation is reasonably accurate

  if(tau < 2*del) then

    write(6,*)
    write(6,'("  The energy mesh is too coarse for such low temperature")')
    write(6,*)

    return

  endif

  if(tau < 6*del) then

    write(6,*)
    write(6,'("  The energy mesh is too coarse for such low temperature")')
    write(6,'("  The accuracy of calculated values will be low")')

  endif

  if(tau > (ehist(nhist) - ehist(ncbm))/5 .or.                           &
     tau > (ehist(nvbm) - ehist(1))/5 ) then

    write(6,*)
    write(6,'("  The energy mesh range is too small for such ",          &
        & "high temperature")')
    write(6,'("  The accuracy of calculated values will be low")')
    write(6,*)

  endif

  del = del*HARTREE

! finds the fermi level using code adapted from rtsafe of numerical recipes
! adapted from fermi_level

  maxit = 300

! bracket roots

  nave = (nvbm + ncbm) / 2
  do n = nave,1,-1
    xl = ehist(n)
    m = n

    call dos_carrier_conc(xn,xp,dxndef,dxpdef,tau,xl,nvbm,ncbm,           &
                                 nhist,ehist,dhist)

    if(xp > xn) exit

  enddo

  do n = m,nhist
    xh = ehist(n)

    call dos_carrier_conc(xn,xp,dxndef,dxpdef,tau,xh,nvbm,ncbm,          &
                                  nhist,ehist,dhist)

    if(xp < xn) exit

  enddo

  ef = (xl + xh) / 2
  dx = abs(xh - xl) / 2
  dxold = dx

  call dos_carrier_conc(xn,xp,dxndef,dxpdef,tau,ef,nvbm,ncbm,            &
                                  nhist,ehist,dhist)
  f = xn - xp
  df = dxndef - dxpdef

  do n = 1,maxit
    m = n
    if(((ef - xh)*df - f)*((ef - xl)*df - f) > ZERO .or.                 &
          abs(2*f) > abs(dxold*df)) then

!     bisection

      dxold = dx
      dx = (xh - xl) / 2
      ef = xl + dx

      if (xl == ef) exit

    else

!     Newton

      dxold = dx
      dx = f / df
      efold = ef
      ef = ef - dx

      if (efold == ef) exit

    endif

    if (abs(dx) < EPS*EPS) exit

    call dos_carrier_conc(xn,xp,dxndef,dxpdef,tau,ef,nvbm,ncbm,          &
                                  nhist,ehist,dhist)
    f = xn - xp
    df = dxndef - dxpdef
    if (f < zero) then
      xl = ef
    else
      xh = ef
    endif
  enddo

  if(m == maxit) then
    write(6,*) '  ERROR: maximum number of iterations exceeded',m

    stop

  endif

! Fermi level and effective masses

  xmu = ef

  xmdosv = ( 0.5*(xn/vcell)*exp( (ef-ehist(nvbm)) / tau ) )**(2.0/3.0)
  xmdosv = xmdosv*6.28/tau

  xmdosc = ( 0.5*(xn/vcell)*exp( (ehist(ncbm)-ef) / tau ) )**(2.0/3.0)
  xmdosc = xmdosc*6.28/tau

  xni = sqrt(xn*xp) / vcell

! Prints the results

  write(6,*)
  write(6,*)
  write(6,'("  For a temperature of ",f10.1," K")') tempk
  write(6,*)
  write(6,*)
  write(6,'("  The Fermi level is at ",f10.3," eV")') ef*HARTREE
  write(6,*)
  write(6,'("  It is ",f10.3," eV  above VBM")') (ef-ehist(nvbm))*HARTREE
  write(6,'("  and ",f10.3," eV  below CBM")') (ehist(ncbm)-ef)*HARTREE
  write(6,*)
  if(del > 0.001) then
     write(6,'("  The accuracy of those values is ",f10.3," eV")') del
    write(6,*)
  endif
  write(6,'("  The band gap is between ",f10.3," and ",f10.3, " eV")')   &
           (ehist(ncbm)-ehist(nvbm))*HARTREE,                            &
           (ehist(ncbm+1)-ehist(nvbm-1))*HARTREE
  write(6,*)
  write(6,*)
  write(6,'("  The intrinsic carrier concentration is ",e12.3,           &
           & " per unit cell")') xn
  write(6,*)
  write(6,'("  The intrinsic carrier concentration is ",e12.3,           &
           & " per cubic centimeter")')                                  &
          1.0E24 * xn / (vcell*BOHR*BOHR*BOHR)
  write(6,*)
  write(6,'("  The DOS effective mass is ",f10.3," for the ",            &
       & "valence band and ",f10.3," for the conduction band")')         &
       xmdosv,xmdosc
  write(6,*)
  write(6,*)

  return
end subroutine dos_f_level
