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

!>     reads the plane-wave cutoff and BZ integration mesh from a file
!>     opened by a previous edsf_init

       subroutine read_esdf_bz(ipr,                                      &
     & emax,nbandin,nx,ny,nz,sx,sy,sz,lbz)

!      Written June 16, 2017. jlm
!      Modified, documentation, August 10 2019. JLM
!      copyright inesc-mn/Jose Luis Martins


!      version 4.94

       use esdf

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  ipr                        !<  print control: ipr < 3 no printing, except warnings.


!      output


       real(REAL64), intent(out)          ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).

       integer, intent(out)               ::  nx, ny, nz                 !<  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64), intent(out)          ::  sx, sy, sz                 !<  offset of the integration mesh (usually 0.5)

       integer, intent(out)               ::  nbandin                    !<  target for number of bands      
       logical, intent(out)               ::  lbz                        !<  indicates if Brillouin Zone data was successfully read.
  
!      local variables

       integer                    ::  nlines
       integer                    ::  nmat(3,3)
       real(REAL64)               ::  sxyz(3)
       integer                    ::  ioerr

!      parameters

       real(REAL64), parameter ::  ZERO = 0.0_REAL64

!      counters

       integer   ::  i, j


!      plane-wave cutoff energy

       lbz = .TRUE.

       if(esdf_defined('PWEnergyCutoff')) then

         emax = ZERO
         emax = esdf_physical('PWEnergyCutoff',emax,'hartree')

         if(emax < 1.0) then
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   Plane-wave cut-off looks very small'
           write(6,'("    The value of emax is: ",e12.3)') emax
         endif

         if(emax > 500.0) then
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   Plane-wave cut-off looks looks very large'
           write(6,'("   The value of emax is: ",e12.3)') emax
         endif

       else
         lbz = .FALSE.
         emax = 5.0
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Using default values for emax'
       endif


!      Number of eigenstates

       if(esdf_defined('NumberOfEigenStates')) then

         nbandin = 0
         nbandin = esdf_integer('NumberOfEigenStates',nbandin)

         if(nbandin < 1) then
           nbandin = 1
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   Number of bands was 0 or negative. Set to 1.'
         endif

         if(nbandin > 5000) then
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   Number of bands looks very large'
           write(6,'("    The value of nbandin is: ",i8)') nbandin
         endif

       else
         lbz = .FALSE.
         nbandin = 10
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Using default values for nbandin'
       endif


!      BZ integration


       if(esdf_block('kgrid_Monkhorst_Pack',nlines)) then
         if(nlines /= 3) then
           write(6,*)
           write(6,*) '  esdf_read_bz:  Wrong Number of Lines'
           write(6,*) '  in kgrid_Monkhorst_Pack!'

           stop

         endif

         ioerr = 0
         do i=1,nlines

           read(block_data(i),*,iostat=ioerr) (nmat(i,j),j=1,3),sxyz(i)
         
           if(ioerr /= 0) then
             write(6,*)
             write(6,*) '    Stopped in read_esdf_bz'
             write(6,*) '    error reading kgrid_Monkhorst_Pack'

             stop

           endif

         enddo

         nx = nmat(1,1)
         ny = nmat(2,2)
         nz = nmat(3,3)

         sx = sxyz(1)
         sy = sxyz(2)
         sz = sxyz(3)

         if(nmat(1,2) /= 0 .OR. nmat(1,3) /= 0 .OR. nmat(2,1) /= 0 .OR.  &
     &      nmat(2,3) /= 0 .OR. nmat(3,1) /= 0 .OR. nmat(3,1) /= 0) then
           write(6,*)
           write(6,*) '   WARNING '
           write(6,*) '   in kgrid_Monkhorst_Pack non-zero off-diagonal'
           write(6,*) '   values are not implemented'
           write(6,*)
           write(6,'(3i6)') (nmat(1,j),j=1,3)
           write(6,'(3i6)') (nmat(2,j),j=1,3)
           write(6,'(3i6)') (nmat(3,j),j=1,3)
           write(6,*)
         endif

         if(nx < 1) then
           nx = 1
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   nx was 0 or negative. Set to 1.'
         endif

         if(ny < 1) then
           ny = 1
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   ny was 0 or negative. Set to 1.'
         endif

         if(nz < 1) then
           nz = 1
           write(6,*)
           write(6,*) '   WARNING'
           write(6,*)
           write(6,*) '   nz was 0 or negative. Set to 1.'
         endif

       else
         lbz = .FALSE.

         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Using default values for BZ integration'

         nx = 4
         ny = 4
         nz = 4

         sx = 0.5
         sy = 0.5
         sz = 0.5
       endif


!      prints the results

       if(ipr > 2) then
         write(6,*)
         write(6,*)
         write(6,*)  '  The values set by read_esdf_bz are:'
         write(6,*)
         write(6,'("   The value of emax is: ",e12.3)') emax
         write(6,'("   The value of nbandin is: ",i8)') nbandin
         write(6,'("   The value of nx,ny,nz is: ",3i5)') nx,ny,nz
         write(6,'("   The value of sx,sy,sz is: ",3f10.3)') sx,sy,sz
         write(6,*)
         write(6,*)
       endif

       return
       end subroutine read_esdf_bz
