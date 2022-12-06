       program h_pot

!      writes the kleinman and bylander pseudo for hydrogen
!      formats the pseudopotential file for transfer with ftp
!      or unformats a transfered file

       implicit double precision (a-h,o-z)
       parameter ( pi = 3.141592653589793d0)
       parameter(mxdnqp = 101, mxdlqp = 20003)
       parameter(mxdnga = 30)

       character*2 namel,icorr,icorrt
       character*3 irel
       character*4 icore
       character*10 iray(6),ititle(7)

       dimension vqnl(mxdnqp,mxdnqp)
       dimension ngaunl(3),xxnl(mxdnga),wgnl(mxdnga)
       dimension vda(mxdlqp)
       dimension nkb(3)

       isize=32*max(mxdlqp,mxdnqp*mxdnqp)


         it=12
         open(unit=it,file='h_potkb.dat',status='new',                   &
     &        form='unformatted')
       open (unit=10,file='H_POTKB_F.DAT',status='new',form='formatted')

!        title

         namel = 'H '
         icorrt = 'CA'
         irel = '   '
         icore = '    '
         do j=1,6
           iray(j) = '          '
         enddo
         iray(1) = 'True Poten'
         iray(2) = 'tial      '
         do j=1,7
           ititle(j) = '          '
         enddo
         write(it) namel,icorrt,irel,icore,(iray(j),j=1,6),              &
     &   (ititle(j),j=1,7)
         write(10,'(1x,a2,1x,a2,1x,a3,1x,a4,1x,6a10,1x,7a10)') namel,    &
     &      icorrt,irel,icore,(iray(i),i=1,6),(ititle(i),i=1,7)

!        valence and mesh sizes

           izv = 1
           nql = 20000
           delql = 0.005d0
           vql0 = 0.0
           write(it) izv,nql,delql,vql0
      write(10,*) izv,nql,delql,vql0
           norb = 1
           nkb(1) = 0
           write(it) norb,(nkb(j),j=1,norb)
      write(10,*) norb
      write(10,*) 0
      write(10,*) (nkb(i),i=1,norb)
      write(10,*) 0.5d0

!        local potential

         do j=1,nql
           q = real(j)*delql
           vda(j) = -8.0*pi/(q*q)
         enddo

         write(it) (vda(j),j=1,nql)
      do i=1,nql
        write(10,*) vda(i)
      enddo

!        non-local potential

           do lp1 = 1,norb
             do j=1,nql+1
               vda(j) = 0.0
             enddo
             write(it) (vda(j),j=1,nql+1)
        do j=1,nql+1
          write(10,*) vda(j)
        enddo
       enddo

!        core charge

         do j=1,nql
           vda(j) = 0.0
         enddo
         write(it) (vda(j),j=1,nql)
      do j=1,nql
        write(10,*) vda(j)
      enddo

!        valence charge

         do j=1,nql
           q = delql*real(j)
           vda(j) = 16.0/((4.0+q*q)*(4.0+q*q))
         enddo
         write(it) (vda(j),j=1,nql)
      do j=1,nql
        write(10,*) vda(j)
      enddo

         close (unit=it)

      nqwf = nql
      deltaq = delql
      norbat = 1

      write(10,*) nqwf,deltaq,norbat
         do j=1,nql
           q = delql*real(j)
           vda(j) = 16.0/((1.0+q*q)*(1.0+q*q))
         enddo

          write(10,*) 0,0.5d0
          do j=1,nqwf
            write(10,*) vda(j)
          enddo

      close(unit=10)


       stop
       end
