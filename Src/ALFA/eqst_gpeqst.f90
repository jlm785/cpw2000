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

!>  Writes files for the plot of the equation of state V(p)
!>  of the most stable structures, and finds the transition pressures.
!>  Energies and volumes are per formula unit.
!>  The plot can be done with gnuplot (eqst_com.gp) or
!>  python/matplotlib (eqst_plot.py).
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         before 1993, 4 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_gpeqst(ez, vz, bz, bp, nstr, label, ftype, pmax, npts,   &
        mxdstr)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! mxdstr (and mxdnpt) passed as arguments. 3 October 2026. JLM+claude
! Intent and description of all arguments. 3 October 2026. JLM+claude
! Warnings when the pressure is outside the stable range of the Birch equation. 3 October 2026. JLM+claude
! Murnaghan V(p) for the plot range of a Murnaghan fit. 3 October 2026. JLM+claude
! REAL64 literals, GPa conversion from AUTOGPA. 3 October 2026. JLM+claude
! Also writes a python (matplotlib) script eqst_plot.py. 4 October 2026. JLM+claude
! gnuplot index is the number of the curve, not of the structure; last label. 4 October 2026. JLM+claude

  implicit none
  integer, parameter  ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  mxdstr                               !<  maximum number of structures

  real(REAL64), intent(in)      ::  ez(mxdstr)                           !<  energy at equilibrium per formula unit (Rydberg)
  real(REAL64), intent(in)      ::  vz(mxdstr)                           !<  equilibrium volume per formula unit (bohr^3)
  real(REAL64), intent(in)      ::  bz(mxdstr)                           !<  bulk modulus for each structure (Rydberg/bohr^3)
  real(REAL64), intent(in)      ::  bp(mxdstr)                           !<  pressure derivative of the bulk modulus for each structure
  integer, intent(in)           ::  nstr                                 !<  number of structures
  character(len=10), intent(in) ::  label(mxdstr)                        !<  label of the structure
  character(len=5), intent(in)  ::  ftype                                !<  type of equation of state, MURNA or BIRCH

  real(REAL64), intent(in)      ::  pmax                                 !<  maximum pressure for plot (GPa)
  integer, intent(in)           ::  npts                                 !<  number of points in plot

! local

  real(REAL64)                  ::  h(mxdstr)

  real(REAL64)          ::  emin, vref, hmin
  integer               ::  jmin, jmp
  integer               ::  nseg                                         !  number of curves (blocks) already written
  real(REAL64)          ::  p                                            !  pressure
  real(REAL64)          ::  xx, vr, vmin

  integer               ::  iodat, iocom
  integer               ::  iopy                                         !  tape for the python script
  character(len=20)     ::  pylabel                                      !  label for python

  integer               ::  ierr
  logical               ::  lwarn(mxdstr)                                !  warning already written for structure

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   ::  AUTOGPA = 29421.0158_REAL64              !  Hartree/bohr^3 to GPa
  real(REAL64), parameter   ::  GPA = AUTOGPA/2                          !  Rydberg/bohr^3 to GPa

! counters

  integer                   ::  i, j

  iodat = 9
  iocom = 10
  open(unit=iodat,file='eqst_dat.gp')
  open(unit=iocom,file='eqst_com.gp')
  iopy = 11
  open(unit=iopy,file='eqst_plot.py')

  write(iocom,'("set term wxt persist")')
  write(iocom,'("# set output ''eqst.pdf''")')
  write(iocom,'("# set term pdf color font ''Helvetica,12''")')
  write(iocom,'("set title ''Equation of state p(V)''")')

  write(iocom,'("set format x ''%.0f''")')
  write(iocom,'("set format y ''%.2f''")')

! gets volumes for the range of the plot

  emin = ez(1)
  jmin = 1
  do i = 1,nstr
    if(ez(i) < emin) then
      jmin = i
      emin = ez(i)
    endif
  enddo
  vref = vz(jmin)
  vmin = UM
  do i = 1,nstr
    if(ftype == 'murna' .or. ftype == 'MURNA') then
      xx = UM + bp(i)*pmax/(GPA*bz(i))
      vr = vz(i) / xx**(UM/bp(i))
    else
      call eqst_vofp(vr,pmax/GPA,vz(i),bz(i),bp(i),ierr)
    endif
    if(vr/vref < vmin) vmin = vr/vref
  enddo

  do i = 1,nstr
    lwarn(i) = .FALSE.
  enddo
  vmin = 0.95_REAL64*vmin

  write(iocom,'("set yrange [",f8.2,":",f8.2,"]")')  vmin, UM
  write(iocom,'("set xrange [",f8.1,":",f8.1,"]")')  ZERO, pmax

  write(iocom,'("set xlabel ''p (GPa)''")')
  write(iocom,'("set ylabel ''V / V_0'' ")')
  write(iocom,*)

  write(iocom,'("plot \")')

  call eqst_py_header(iopy, 'eqst_dat.gp', 'eqst_com.gp',                &
      'Equation of state p(V)', 'p (GPa)', '$V / V_0$')
  write(iopy,'("ax.set_xlim(",f8.1,",",f8.1,")")')  ZERO, pmax
  write(iopy,'("ax.set_ylim(",f8.2,",",f8.2,")")')  vmin, UM
  write(iopy,'("ax.xaxis.set_major_formatter(FormatStrFormatter(''%.0f''))")')
  write(iopy,'("ax.yaxis.set_major_formatter(FormatStrFormatter(''%.2f''))")')
  write(iopy,*)

  nseg = 0
  do j = 1,npts
    p = (j-1)*pmax/(GPA*(npts-1))
    jmp = jmin

    if(ftype == 'murna' .or. ftype == 'MURNA') then
      do i = 1,nstr
        xx = UM + bp(i)*p/bz(i)
        xx = xx ** ((bp(i)-UM)/bp(i))
        h(i) = ez(i) + bz(i)*vz(i)*(xx-UM)/(bp(i)-UM)
      enddo
    else
      do i = 1,nstr
        call eqst_hofp(h(i),p,ez(i),vz(i),bz(i),bp(i),ierr)

!       structure is not (mechanically) stable at this pressure

        if(ierr /= 0) then
          h(i) = huge(UM)
          if(.not. lwarn(i)) then
            write(6,*)
            write(6,'("  WARNING: structure ",a10,                       &
                &     " is not stable above",f10.2,                      &
                &     " GPa in the Birch fit")') label(i), p*GPA
            write(6,*)
            lwarn(i) = .TRUE.
          endif
        endif
      enddo
    endif

    hmin = h(1)
    jmin = 1
    do i = 1,nstr
      if(h(i) < hmin) then
        jmin = i
        hmin = h(i)
      endif
    enddo

    if(jmin /= jmp) then
      write(6,'(" TRANSITION PRESSURE P= ",F8.3,"  BETWEEN",2I4)')       &
           p*GPA,jmp,jmin
      write(iocom,'("''eqst_dat.gp'' using 1:2 index ",i3,               &
          &    " w li ls ",i3," ti ''",a10,"'', \")') nseg, jmp,label(jmp)
      call eqst_label_py(label(jmp), pylabel)
      write(iopy,'("ax.plot(blocks[",i3,"][:,0], blocks[",i3,"][:,1], ", &
          &     "''-'', color=''C",i0,"'', label=r""",a,""")")')         &
          nseg, nseg, mod(jmp-1,10), trim(pylabel)
      nseg = nseg + 1
      write(iodat,*)
      write(iodat,*)
   endif

    if(ftype == 'murna' .or. ftype == 'MURNA') then
      xx = UM + bp(jmin)*p/bz(jmin)
      xx = UM / xx ** (UM/bp(jmin))
      vr = xx*vz(jmin)/vref
    else
      call eqst_vofp(vr,p,vz(jmin),bz(jmin),bp(jmin),ierr)
      vr=vr/vref
    endif

    write(iodat,'(2f15.8)') p*GPA, vr

  enddo
  write(iocom,'("''eqst_dat.gp'' using 1:2 index ",i3," w li ls ",       &
      &     i3," ti ''",a10,"''")')  nseg, jmin,label(jmin)
  call eqst_label_py(label(jmin), pylabel)
  write(iopy,'("ax.plot(blocks[",i3,"][:,0], blocks[",i3,"][:,1], ",     &
      &     "''-'', color=''C",i0,"'', label=r""",a,""")")')             &
      nseg, nseg, mod(jmin-1,10), trim(pylabel)

  write(iopy,*)
  write(iopy,'("ax.legend()")')
  write(iopy,'("# fig.savefig(''eqst.pdf'')")')
  write(iopy,'("plt.show()")')

  close(unit=iocom)
  close(unit=iodat)
  close(unit=iopy)

  return

end subroutine eqst_gpeqst
