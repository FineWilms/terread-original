      program veg

      include 'newmpar.h'

      parameter ( rtd=180.0/3.14159265 )
      parameter ( dl = 1.0/120.0 )
      parameter ( gslon = -180+dl/2. )
      parameter ( gslat = 90-dl/2. )
      parameter ( mx=43200, my=21600 )

      character*13 gfile
      parameter( gfile='gsib2_0ll.img')


      integer i,j
      byte zza(mx)
      integer iz(mx)

! global 1 km data
      write(*,*) 'gfile=',gfile
      open(20,file=gfile,access='direct',form='unformatted',recl=mx)
      do jg=1,my
! read each row of relevant data from 1 km data set (row=jg)
        read(20,rec=jg) (zza(ig),ig=1,mx)
        iz_max=-9999
        do ig=1,mx
          iz(ig)=zza(ig)
          if(iz(ig).ge.19)iz(ig)=0
          iz_max=max(iz_max,iz(ig))
        enddo ! ig=1,mx
        write(6,*)"jg=",jg," iz_max=",iz_max
      enddo ! jg=jmin,jmax

      stop
      end
