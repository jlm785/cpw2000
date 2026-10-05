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

!>  Writes files for the gnuplot plot of E(V),
!>  calculated points and fitted curves.
!>  Energies and volumes are per formula unit.
!>  The plot can be done with gnuplot (eofv_com.gp) or
!>  python/matplotlib (eofv_plot.py).
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         before 1993, 5 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_gpeofv(xarz, yarz, nptz, ez, vz, bz, bp, nstr,           &
        label, ftype,                                                    &
        mxdstr,  mxdnpt)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! mxdstr (and mxdnpt) passed as arguments. 3 October 2026. JLM+claude
! Intent and description of all arguments. 3 October 2026. JLM+claude
! REAL64 literals, HARTREE instead of EV. 3 October 2026. JLM+claude
! Energies and volumes per formula unit (labels). 4 October 2026. JLM+claude
! Also writes a python (matplotlib) script eofv_plot.py. 4 October 2026. JLM+claude
! 5% padding also at the top of the energy range. 4 October 2026. JLM+claude
! No titles for the fitted curves, no continuation after the last curve,
! trimmed labels. 5 October 2026. JLM+claude

  implicit none
  integer, parameter  ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  mxdstr                               !<  maximum number of structures
  integer, intent(in)           ::  mxdnpt                               !<  maximum number of energy points (per structure)

  real(REAL64), intent(in)      ::  xarz(mxdnpt,mxdstr)                  !<  calculated volumes per formula unit (bohr^3)
  real(REAL64), intent(in)      ::  yarz(mxdnpt,mxdstr)                  !<  calculated energies per formula unit (Rydberg)

  integer, intent(in)           ::  nptz(mxdstr)                         !<  number of calculated points for each structure
  real(REAL64), intent(in)      ::  ez(mxdstr)                           !<  energy at equilibrium per formula unit (Rydberg)
  real(REAL64), intent(in)      ::  vz(mxdstr)                           !<  equilibrium volume per formula unit (bohr^3)
  real(REAL64), intent(in)      ::  bz(mxdstr)                           !<  bulk modulus for each structure (Rydberg/bohr^3)
  real(REAL64), intent(in)      ::  bp(mxdstr)                           !<  pressure derivative of the bulk modulus for each structure
  integer, intent(in)           ::  nstr                                 !<  number of structures
  character(len=10), intent(in) ::  label(mxdstr)                        !<  label of the structure
  character(len=5), intent(in)  ::  ftype                                !<  type of equation of state, MURNA or BIRCH

! local

  real(REAL64)          ::  xminz(mxdstr), xmaxz(mxdstr)
  real(REAL64)          ::  xmin, xmax, emin, emax
  real(REAL64)          ::  erange                                       !  range of calculated energies
  integer               ::  npt
  real(REAL64)          ::  a, b, c, dd
  real(REAL64)          ::  x, y
  real(REAL64)          ::  step

  integer               ::  iodat, iocom
  integer               ::  iopy                                         !  tape for the python script
  character(len=20)     ::  pylabel                                      !  label for python
  character(len=3)      ::  cnext                                        !  continuation of the gnuplot plot command

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   ::  HARTREE = 27.211386246_REAL64            !  Hartree to eV

! counters

  integer                   ::  i, j

  iodat = 9
  iocom = 10
  open(unit=iodat,file='eofv_dat.gp')
  open(unit=iocom,file='eofv_com.gp')
  iopy = 11
  open(unit=iopy,file='eofv_plot.py')

  xmin = xarz(1,1)
  xmax = xarz(1,1)
  emin = yarz(1,1)
  emax = yarz(1,1)
  do j = 1,nstr
    xminz(j) = xarz(1,j)
    xmaxz(j) = xarz(1,j)
    npt = nptz(j)
    do i = 1,npt
      if(xarz(i,j) < xminz(j)) xminz(j) = xarz(i,j)
      if(xarz(i,j) > xmaxz(j)) xmaxz(j) = xarz(i,j)
      if(yarz(i,j) < emin) emin = yarz(i,j)
      if(yarz(i,j) > emax) emax = yarz(i,j)
    enddo
    if(xminz(j) < xmin) xmin = xminz(j)
    if(xmaxz(j) > xmax) xmax = xmaxz(j)
  enddo
  erange = emax-emin
  emin = emin-0.05_REAL64*erange
  emax = emax+0.05_REAL64*erange
  emin = emin*HARTREE/2
  emax = emax*HARTREE/2

  write(iocom,'("set term wxt persist")')
  write(iocom,'("# set output ''eofv.pdf''")')
  write(iocom,'("# set term pdf color font ''Helvetica,12''")')
  write(iocom,'("set title ''Equation of state E(V)''")')

  write(iocom,'("set format x ''%.0f''")')
  write(iocom,'("set format y ''%.1f''")')

  write(iocom,'("set yrange [",f10.2,":",f10.2,"]")')  emin, emax
  write(iocom,'("set xrange [",f8.1,":",f8.1,"]")')  xmin, xmax

  write(iocom,'("set xlabel ''V per formula unit (a.u.)''")')
  write(iocom,'("set ylabel ''E per formula unit (eV)'' ")')
  write(iocom,*)

  write(iocom,'("plot \")')

  call eqst_py_header(iopy, 'eofv_dat.gp', 'eofv_com.gp',                &
      'Equation of state E(V)', 'V per formula unit (a.u.)',             &
      'E per formula unit (eV)')
  write(iopy,'("ax.set_xlim(",f8.1,",",f8.1,")")')  xmin, xmax
  write(iopy,'("ax.set_ylim(",f10.2,",",f10.2,")")')  emin, emax
  write(iopy,'("ax.xaxis.set_major_formatter(FormatStrFormatter(''%.0f''))")')
  write(iopy,'("ax.yaxis.set_major_formatter(FormatStrFormatter(''%.1f''))")')
  write(iopy,*)

  do j = 1,nstr
    npt = nptz(j)
    do i = 1,npt
      write(iodat,'(2f15.8)') xarz(i,j),yarz(i,j)*HARTREE/2
    enddo
    write(iodat,*)
    write(iodat,*)
    write(iocom,'("''eofv_dat.gp'' using 1:2 index ",i3," w p ls ",      &
        &     i3," ti ''",a," calculated'', \")') j-1, j, trim(label(j))
    call eqst_label_py(label(j), pylabel)
    write(iopy,'("ax.plot(blocks[",i3,"][:,0], blocks[",i3,"][:,1], ",   &
        &     "''o'', color=''C",i0,"'', label=r""",a,                   &
        &     " calculated"")")') j-1, j-1, mod(j-1,10), trim(pylabel)
  enddo

  if(ftype == 'murna' .or. ftype == 'MURNA') then
    do j = 1,nstr
      cnext = ', \'
      if(j == nstr) cnext = ' '
      dd = UM-bp(j)
      b = bz(j)/bp(j)
      c = -(vz(j)**bp(j))*b/dd
      a = ez(j)+b*bp(j)*vz(j)/dd
      step = 0.01_REAL64*(xmaxz(j)-xminz(j))
      do i = 1,121
        x = xminz(j) + real(i-11,REAL64)*step
        y = a + b*x + c*x**dd
        write(iodat,'(2f15.8)') x,y*HARTREE/2
      enddo
      write(iocom,'("''eofv_dat.gp'' using 1:2 index ",i3,               &
          &     " w li ls ",i3," notitle",a)') nstr+j-1, j, trim(cnext)
      write(iopy,'("ax.plot(blocks[",i3,"][:,0], blocks[",i3,"][:,1], ", &
          &     "''-'', color=''C",i0,"'')")') nstr+j-1, nstr+j-1, mod(j-1,10)
      write(iodat,*)
      write(iodat,*)
    enddo
  else
    do j = 1,nstr
      cnext = ', \'
      if(j == nstr) cnext = ' '
      step = 0.01_REAL64*(xmaxz(j)-xminz(j))
      do i = 1,121
        x = xminz(j) + real(i-11,REAL64)*step
        call eqst_eofv(y,x,ez(j),vz(j),bz(j),bp(j))
        write(iodat,'(2f15.8)') x,y*HARTREE/2
      enddo
      write(iocom,'("''eofv_dat.gp'' using 1:2 index ",i3,               &
          &     " w li ls ",i3," notitle",a)') nstr+j-1, j, trim(cnext)
      write(iopy,'("ax.plot(blocks[",i3,"][:,0], blocks[",i3,"][:,1], ", &
          &     "''-'', color=''C",i0,"'')")') nstr+j-1, nstr+j-1, mod(j-1,10)
      write(iodat,*)
      write(iodat,*)
    enddo
  endif
  write(iocom,*)

  write(iopy,*)
  write(iopy,'("ax.legend()")')
  write(iopy,'("# fig.savefig(''eofv.pdf'')")')
  write(iopy,'("plt.show()")')

  close(unit=iocom)
  close(unit=iodat)
  close(unit=iopy)

  return

end subroutine eqst_gpeofv
