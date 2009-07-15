      program global_landuse
c
      integer nx,ny
      character*14 infile
      character*11 file
      parameter (dl = 1.0/120.0 )
      parameter (mx = 43200, my=21600 )
      parameter(infile='gsib22_0ll.img')

      integer i,j
      byte zza(mx,my)
      real lons,lats,lone,late,dl

      write(*,*) ' infile=',infile
      open(21,file=infile,access='direct',form='unformatted',recl=mx)

      do j=1,my
       read(21,rec=j) (zza(i,j),i=1,mx)
      enddo
      close(21)

      nx=4800
      ny=6000
      do i = 1,9
       lons=-180.+real(i-1)*40.
       lone=lons+real(nx-1)*dl
       ew="W"
       if ( nint(lons).gt.0 ) ew="E"

       do j = 1,3
        lats=90.-real(j-1)*50.
        late=lats-real(ny-1)*dl 
        ns="N"
        if ( nint(lats).lt.0 ) ns="S"

        write(file,'(a1,i3.3,a1,i2.2,".LU")')
     &               ew,abs(nint(lons)),ns,abs(nint(lats))

        write(6,*)"call write1km(file,nx,ny,lons,lats,zza,debug)"
        call write1km(file,nx,ny,lons,lats,zza,debug)

       enddo ! j

      enddo ! i

      nx=7200
      ny=3600
      do i = 1,6
        lons=-180.+real(i-1)*60.
        lone=lons+real(nx-1)*dl
        ew="W"
        if ( nint(lons).ge.0 ) ew="E"
        lats=-60.
        late=lats-real(ny-1)*dl 
        ns="S"

        write(file,'(a1,i3.3,a1,i2.2,".LU")')
     &               ew,abs(nint(lons)),ns,abs(nint(lats))

        write(6,*)"call write1km(file,nx,ny,lons,lats,zza,debug)"
        call write1km(file,nx,ny,lons,lats,zza,debug)

      enddo ! i

      stop
      end
c=======================================================================
      subroutine write1km(file,nx,ny,lons,late,zza,debug)

      parameter ( dl = 1.0/120.0 )
      parameter ( slon=-180.+dl/2., slat=-90.+dl/2.)
      parameter ( mx = 43200, my=21600 )

      character*11 file
      byte zza(mx,my)
      real dl
      real lons,lats,lone,late
      integer jrec,is,ie,js,je,i,j,jj,nx,ny

      logical debug

      write(6,*)"write1km file,nx,ny=",file,nx,ny

      lone=lons+(nx-1)*dl
      lats=late-(ny-1)*dl

      write(6,*)"lons,lone=",lons,lone
      write(6,*)"lats,late=",lats,late

      is = 2+(lons-slon)/dl
      ie = 2+(lone-slon)/dl
      js = 1+(lats-slat)/dl
      je = 1+(late-slat)/dl

      write(6,*)"is,ie=",is,ie
      write(6,*)"js,je=",js,je

      open(22,file=file,access='direct',form='unformatted',recl=nx)
!#######################################################################
      do j=js,je
!#######################################################################
        jj=my+1-j
        jrec=je+1-j
        write(22,rec=jrec) (zza(i,j),i=is,ie)
!#######################################################################
      enddo ! j
!#######################################################################
      close(22)

      return
      end ! write1km
c=======================================================================
