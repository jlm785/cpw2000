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

!>  Writes the beginning of a python (matplotlib) script that
!>  plots the data in a gnuplot data file.  The data file is read
!>  in blocks, separated by two blank lines, that correspond to
!>  the "index" of gnuplot.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         4 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_py_header(iopy, datafile, comfile, title, xlabel, ylabel)

! Written 4 October 2026. JLM+claude

  implicit none

! input

  integer, intent(in)                   ::  iopy                         !<  tape number of the (open) python file
  character(len=*), intent(in)          ::  datafile                     !<  name of the gnuplot data file
  character(len=*), intent(in)          ::  comfile                      !<  name of the equivalent gnuplot command file
  character(len=*), intent(in)          ::  title                        !<  title of the plot
  character(len=*), intent(in)          ::  xlabel                       !<  label of the x axis
  character(len=*), intent(in)          ::  ylabel                       !<  label of the y axis


  write(iopy,'(a)') '#!/usr/bin/env python3'
  write(iopy,'(a)') '#'
  write(iopy,'(a)') '# Written by eqst.  Plot with matplotlib, equivalent to '//trim(comfile)
  write(iopy,'(a)') '#'
  write(iopy,'(a)') '#   python3 '//trim(comfile(1:index(comfile,'_com.gp')-1))//'_plot.py'
  write(iopy,'(a)')
  write(iopy,'(a)') 'import numpy as np'
  write(iopy,'(a)') 'import matplotlib.pyplot as plt'
  write(iopy,'(a)') 'from matplotlib.ticker import FormatStrFormatter'
  write(iopy,'(a)')
  write(iopy,'(a)')
  write(iopy,'(a)') 'def read_blocks(filename):'
  write(iopy,'(a)') '    """Reads a gnuplot data file.  Two blank lines start a new block (gnuplot index)."""'
  write(iopy,'(a)') '    blocks = [[]]'
  write(iopy,'(a)') '    nblank = 0'
  write(iopy,'(a)') '    with open(filename) as f:'
  write(iopy,'(a)') '        for line in f:'
  write(iopy,'(a)') '            if line.strip():'
  write(iopy,'(a)') '                for i in range(nblank // 2):'
  write(iopy,'(a)') '                    blocks.append([])'
  write(iopy,'(a)') '                nblank = 0'
  write(iopy,'(a)') '                blocks[-1].append([float(x) for x in line.split()[:2]])'
  write(iopy,'(a)') '            else:'
  write(iopy,'(a)') '                nblank = nblank + 1'
  write(iopy,'(a)') '    return [np.array(b).reshape(-1, 2) for b in blocks]'
  write(iopy,'(a)')
  write(iopy,'(a)')
  write(iopy,'(a)') 'blocks = read_blocks("'//trim(datafile)//'")'
  write(iopy,'(a)')
  write(iopy,'(a)') 'fig, ax = plt.subplots()'
  write(iopy,'(a)') 'ax.set_title("'//trim(title)//'")'
  write(iopy,'(a)') 'ax.set_xlabel("'//trim(xlabel)//'")'
  write(iopy,'(a)') 'ax.set_ylabel("'//trim(ylabel)//'")'

  return

end subroutine eqst_py_header
