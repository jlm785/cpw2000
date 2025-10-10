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

!>  Reads filenames and tape numbers from a file
!>  opened by a previous edsf_init
!>
!>  \author       Jose Luis Martins and Carlos Loia Reis
!>  \version      5.12
!>  \date         10 October 2025.
!>  \copyright    GNU Public License v2

subroutine read_esdf_filenames(ipr,                                      &
    pseudo_path, pseudo_suffix, itape_pseudo,                            &
    save_psi_path, itape_save_psi)

! Written 10 October 2025. JLM

  use esdf

  implicit none


! input

  integer, intent(in)                ::  ipr                             !<  print control: ipr < 3 no printing, except warnings.

! output

  character(len=200), intent(out)    ::  pseudo_path                     !<  path to pseudopotentials
  character(len=50), intent(out)     ::  pseudo_suffix                   !<  suffix for the pseudopotentials
  integer, intent(out)               ::  itape_pseudo                    !<  tape number to read pseudo

  character(len=200), intent(out)    ::  save_psi_path                   !<  path to save psi files to disk
  integer, intent(out)               ::  itape_save_psi                  !<  tape number to read and write psi.  If < 10 (default 0) do not use.

! local variables

  integer          ::  i, n


! files for reading pseudopotentials  40 is hard coded

  do i = 1,200
    pseudo_path(i:i) = ' '
  enddo

  do i = 1,50
    pseudo_suffix(i:i) = ' '
  enddo
 itape_pseudo = 40

  pseudo_path = esdf_string('PathToPseudos',pseudo_path)
  n = len(trim(pseudo_path))
  if(n > 0) then
    if(pseudo_path(n:n) /= '/') pseudo_path(n+1:n+1) = '/'
  endif

  pseudo_suffix = '_POTKB_F.DAT'
  pseudo_suffix = esdf_string('PseudoSuffix',pseudo_suffix)

! files for saving wave-functions to disk

  do i = 1,200
    save_psi_path(i:i) = ' '
  enddo
  itape_save_psi = 0

  save_psi_path = esdf_string('PathToSavePsi',save_psi_path)
  n = len(trim(save_psi_path))
  if(n > 0) then
    if(save_psi_path(n:n) /= '/') save_psi_path(n+1:n+1) = '/'
  endif

  itape_save_psi = esdf_integer('TapeToSavePsi',itape_save_psi)
  if(itape_save_psi < 10) itape_save_psi = 0


! prints the results

  if(ipr > 2) then
    write(6,*)
    write(6,*)
    write(6,*)  '  The values set by read_esdf_filenames are:'
    write(6,*)
    write(6,'(  "  The tape number for pseudopotentials is: ",i6)') itape_pseudo
    n = len(trim(pseudo_path))
    if(n > 0) then
      write(6,*)  '  The path to pseudopotential files is: ', pseudo_path
    endif
    write(6,*)  '  The suffix for pseudopotential files is: ', pseudo_suffix
    if(itape_save_psi > 9) then
      n = len(trim(save_psi_path))
      if(n > 0) then
        write(6,*)  '  The path to save wave-function files is: ', save_psi_path
      endif
      write(6,'(  "  The tape number for saving wave-functions is: ",i6)') itape_save_psi
    endif
    write(6,*)
    write(6,*)
  endif

  return

end subroutine read_esdf_filenames
