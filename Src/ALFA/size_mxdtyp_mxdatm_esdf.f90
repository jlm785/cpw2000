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

!>     finds the array sizes for the atomic structure, esdf version

       subroutine size_mxdtyp_mxdatm_esdf(ipr,fname,                     &
     & lgeom,mxdtyp,mxdatm)

!      written June 2017. JLM
!      Modified, documentation, January 2020. JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94

       use esdf

       implicit none

!      input

       integer, intent(in)                ::  ipr                        !<  print control: ipr < 3 no printing, except warnings.

       character(len=*), intent(in)       ::  fname                      !<  file to be written

!      output

       logical, intent(out)               ::  lgeom                      !<  indicates if data was present
       integer, intent(out)               ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(out)               ::  mxdatm                     !<  array dimension of types of atoms

!      local variables

       integer     ::  nlines


       
       call esdf_init(fname)

       lgeom = .TRUE.

       if(esdf_defined('NumberOfSpecies')) then

         mxdtyp = esdf_integer('NumberOfSpecies',1)

         if(mxdtyp < 1) then
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   Number of elements < 1.'
           write(6,*) '   Will be set by Chemical_Species_Label.'
           if(esdf_block('Chemical_Species_Label',nlines)) then
             mxdtyp = nlines
           else
             lgeom = .FALSE.
           endif
         endif

       else

         if(esdf_block('Chemical_Species_Label',nlines)) then
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   Number of elements not defined.'
           write(6,*) '   Will be set by Chemical_Species_Label.'
           mxdtyp = nlines
         else
           lgeom = .FALSE.
         endif
       endif

       if(esdf_defined('NumberOfAtoms')) then

         mxdatm = esdf_integer('NumberOfAtoms',1)

         if(mxdatm < 1) then
           if(esdf_block('AtomicCoordinatesAndAtomicSpecies',nlines))    &
     &     then
             write(6,*)
             write(6,*) '   WARNING'
             write(6,*)
             write(6,*) '   Number of atoms < 1.'
             write(6,*) '   Will be set by AtomicCoordinates.'
             mxdatm = nlines
           else
             lgeom = .FALSE.
           endif
         endif

       else
         if(esdf_block('AtomicCoordinatesAndAtomicSpecies',nlines))      &
     &   then
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   Number of atoms not defined.'
           write(6,*) '   Will be set by AtomicCoordinates.'
           mxdatm = nlines
         else
           lgeom = .FALSE.
         endif
       endif

       
       call esdf_close

!      prints the results

       if(ipr > 2 .and. lgeom) then
         write(6,*)
         write(6,*)
         write(6,*)  '  The values set by size_mxdtyp_mxdatm_esdf are:'
         write(6,*)
         write(6,'("   The value of mxdtyp is: ",i5)') mxdtyp
         write(6,'("   The value of mxdatm is: ",i8)') mxdatm
         write(6,*)
         write(6,*)
       endif

       return

       end subroutine size_mxdtyp_mxdatm_esdf

