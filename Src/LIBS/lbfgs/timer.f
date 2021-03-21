      subroutine timer(ttime)
      real*8 ttime
c     *********
c
c     Subroutine timer, wrapper to zesec
c
      call zesec(ttime)

      return

      end

