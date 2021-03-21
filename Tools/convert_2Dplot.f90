       program convert_2Dplot

!      converts the gnuplot file into a more general file

       implicit none

       character(len=100)  ::  filename
       character(len=10)   ::  file_out = 'rho_xy.dat'
       integer             ::  io = 10
       integer             ::  iout = 11

       character(len=10)   ::  cdummy
       integer             ::  ioerr

       character(len=3)    ::  ctype
       integer             ::  nx, ny
       real                ::  dx, dy
       real                ::  dsclx, dscly

       real                ::  rho

       integer             ::  i, j

       write(6,*)
       write(6,*) '  This program converts the 2D plot format'
       write(6,*)
       write(6,*) '  These are the available gnuplot files'
       write(6,*)

       call system('ls *2D.gp')

       write(6,*)
       write(6,*) '  Type the name of the file you want to convert'
       write(6,*)
       
       read(5,*) filename

       open(unit=io,file=adjustl(trim(filename)),form='FORMATTED')

       read(io,*)  cdummy, cdummy, ctype

       if(ctype == 'wxt') then
         backspace(io)
         read(io,*)  cdummy, cdummy, ctype, cdummy, dsclx, dscly
       else
         read(io,*)  cdummy, cdummy, ctype, cdummy, dsclx, dscly
         read(io,*)
       endif

       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)
       read(io,*)

       read(io,*,iostat=ioerr) cdummy,cdummy,cdummy,cdummy,nx,ny

       dx = dsclx / (800.0*(nx-1))
       dy = dsclx / (800.0*(ny-1))

       read(io,*)

       open(unit=iout,  file=file_out, form="formatted")

       do j = 1,ny
         do i = 1,nx
           read(io,*) rho
           write(iout,'(3f14.7)') (i-1)*dx, (j-1)*dy, rho
         enddo
         read(io,*)
         read(io,*)
       enddo

       write(6,*)
       write(6,*) '  You have a new file  ',file_out,'  with the converted data'
       write(6,*)

       close(unit=iout)
       close(unit=io)
     
       stop
       end program convert_2Dplot
