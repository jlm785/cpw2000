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


!>   Interactive master program that calls independent programs desguised as subroutines.
!>   All "communications" are through external files.
!>
!>   All the input is recorded in a file so that calculations can be reproduced
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         18 january 2022. 1 April 2023. 8 November 2023.
!>  \copyright    GNU Public License v2


program cpw_post_process

! Written August 2, 2014. JLM
! Modified documentation, details, 3 January 2019. JLM
! Added radiative recombination, plot of wave-functions. January 2021. JLM
! minor details.  8 November 2023. JLM

  implicit none

  integer              :: i
  integer              :: ichoice
  integer              :: ioreplay = 35                             !<  hard coded tape unit for replay
  character(len = 15)  :: fname = 'replay_post.dat'                 !<  hard coded filename for replay

  character(len=4)                   ::  vdriv                      !<  version of this program

! Driver program version

  vdriv = '5.09'

  call tpage(vdriv)

  open(unit = ioreplay, file = fname,  status = 'UNKNOWN',          &
                    form = 'FORMATTED')

  if(ioreplay < 7) stop   'ioreplay < 7'

  write(6,*)
  write(6,*) '  Welcome to cpw_analysis'
  write(6,*)

! loop over tasks

  do i=1,100

    write(6,*)
    write(6,*) '  Choose what you want to do:'
    write(6,*)
    write(6,*) '  0)  Exit program'
    write(6,*)
    write(6,*) '  1) Generate band structures,'
    write(6,*) '              density-of-states,'
    write(6,*) '              oscillator strengths,'
    write(6,*) '              k.p files,'
    write(6,*) '              optical response,'
    write(6,*) '              effective masses,'
    write(6,*) '              topological quantities.'
    write(6,*)
    write(6,*) '  2) Use the TB interpolation to generate band '
    write(6,*) '     structures and density-of-states'
    write(6,*)
    write(6,*) '  3) Plot charge densities and potentials'
    write(6,*)
    write(6,*) '  4) Make density of states plots from previously'
    write(6,*) '     generated data'
    write(6,*)
    write(6,*) '  5) Make band structure plots from previously'
    write(6,*) '     generated k.p data'
    write(6,*)
    write(6,*) '  6) Make dielectric function plots from previously'
    write(6,*) '     generated data'
    write(6,*)
    write(6,*) '  7) Estimate the radiative recombination rate'
    write(6,*) '     from previously generated data'
    write(6,*)
    write(6,*) '  8) Plot wave-functions'
    write(6,*)
    write(6,*) '  9) Check crystal structure: lattice,'
    write(6,*) '     symmetry, atom neighbors, write CIF file'
    write(6,*)
    write(6,*) '  Enter your choice (0,1,2,...,9)'
    write(6,*)

    read(5,*) ichoice
    write(ioreplay,*) ichoice,'   main task'

    if(ichoice == 0) then

      exit

    elseif(ichoice == 1) then

      call cpw_pp_band_dos_opt(ioreplay)

    elseif(ichoice == 2) then

      call ao_interpolation_sub(ioreplay)

    elseif(ichoice == 3) then

      call plot_rho_v_sub(ioreplay)

    elseif(ichoice == 4) then

      call dos_sub(ioreplay)

    elseif(ichoice == 5) then

      call kdotp_sub(ioreplay)

    elseif(ichoice == 6) then

      call opt_sub(ioreplay)

    elseif(ichoice == 7) then

      call opt_rad_sub(ioreplay)

    elseif(ichoice == 8) then

      call plot_psi_sub(ioreplay)

    elseif(ichoice == 9) then

      call voronoi_sub(ioreplay)

   endif

  enddo

  close(unit = ioreplay)

  stop

end program cpw_post_process
