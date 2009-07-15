      program global_terrain

      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'map.h'   ! em
      include 'xyzinfo.h'   ! rlat,rlong
      include 'indices.h'
      include 'parm.h'   ! to pass rlong0,rlat0,schmidt  to setxyz
      logical olam

c     various work arrays
      common/work/rmsk(il,jl),inum(il,jl),inumx(il,jl)
     &,zss(il,jl),almsk(il,jl),tmax(il,jl),tmin(il,jl),tsd(il,jl)
      real grid(il,jl)

      logical debug
      character*60 fileout
      character*9 formout

      integer mx,my
      integer nx,ny
      parameter (mx=4800,my=6000)
      parameter (rtd=180.0/3.14159265 )
      integer i,j
      character*1 ns,ew
      character*11 file
      integer lons,lats

      data debug/.false./

      namelist / topnml / ds, du, tanl, rnml, stl1, stl2, debug
     &  ,luout, fileout, olam, wbd, sbd, dlon, dlat
     &  ,idia, jdia ,rlong0, rlat0, schmidt

      open ( unit=5,file='top.nml',status='unknown' )
      read ( 5,topnml, end=5 )
 5    write ( 6,topnml )

      call setxyz

      gridx=-1.e35
      rlatx=-1.e35
      rlonx=-1.e35
      gridn=+1.e35
      rlatn=+1.e35
      rlonn=+1.e35
      imax=-999
      jmax=-999
      imin=+999
      jmin=+999
      do j=1,jl
       do i=1,il
         n=i+(j-1)*il
         grid(i,j)=ds/em(i,j)
         gridx=max(gridx,grid(i,j))
         gridn=min(gridn,grid(i,j))
         if ( grid(i,j).lt.50.e3 ) then
           rlatx=max(rlatx,rtd*rlat(n))
           rlonx=max(rlonx,rtd*rlong(n))
           rlatn=min(rlatn,rtd*rlat(n))
           rlonn=min(rlonn,rtd*rlong(n))
           jmax=max(jmax,j)
           imax=max(imax,i)
           jmin=min(jmin,j)
           imin=min(imin,i)
           write(6,*)i,j,n,grid(i,j),rtd*rlong(n),rtd*rlat(n)
         endif
       enddo
      enddo

      write(6,*)"grid x,n=",gridx,gridn
      write(6,*)"rlon x,n=",rlonx,rlonn
      write(6,*)"rlat x,n=",rlatx,rlatn
      write(6,*)"i    x,n=",imax,imin
      write(6,*)"j    x,n=",jmax,jmin
      do j=1,288,48
       do i=1,48,47
        n=i+(j-1)*il
        write(6,*)i,j,grid(i,j),rtd*rlong(n),rtd*rlat(n)
       enddo
      enddo
      do j=48,288,48
       do i=1,48,47
        n=i+(j-1)*il
        write(6,*)i,j,grid(i,j),rtd*rlong(n),rtd*rlat(n)
       enddo
      enddo
      stop
      end
      include 'setxyz.f'
      include 'jimcc.f'
      include 'latltoij.f'
