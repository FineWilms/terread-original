      program veg

      include 'newmpar.h'

      parameter ( rtd=180.0/3.14159265 )
      parameter ( dl = 1.0/120.0 )
      parameter ( gslon = -180+dl/2. )
      parameter ( gslat = 90-dl/2. )
      parameter ( mx=43200, my=21600 )

      character*13 gfile
      parameter( gfile='gsib2_0ll.img')

      integer nsum(il,jl,0:20) ! 44 with old ones

      integer i,j
      byte zza(mx)
      integer iz(mx)
      character*9 formout
      character*80 vegin,vegout,topo,topout
      character*80 title

      include 'dates.h'   ! to pass ds from setxyz
      include 'indices.h'
      include 'map.h'     ! em
      include 'parm.h'    ! to pass rlong0,rlat0,schmidt  to setxyz
      include 'xyzinfo.h' ! rlat,rlong
      include 'rwork.h'   ! grid,tsd

      integer iarr(il,jl)
      real arr(ifull)
      real zs(ifull),lsm(ifull)
      real waterfrac(ifull),npoints(ifull),nsummax(ifull)

      namelist / vegnml / topo,vegin,vegout,topout,id,jd,ds1km

      data id/37/,jd/81/,ds1km/60./

      open (85,file='veg.nml',status='old')
      read (85,vegnml)
      write (6,vegnml)

      ijd = id+(jd-1)*il

! global 1 km data
      write(*,*) 'gfile=',gfile
      open(20,file=gfile,access='direct',form='unformatted',recl=mx)
      do jg=jmin,jmax
! read each row of relevant data from 1 km data set (row=jg)
        read(20,rec=jg) (zza(ig),ig=1,mx)
        iz(ig)=zza(ig)
        if(iz(ig).ge.19)iz(ig)=0
        write(6,*)"jg=",jg," iz=",iz
      enddo ! jg=jmin,jmax

      if (jg.ne.0)stop

      write(6,*)"usual topo data"
      open (21,file=topo,status='old')
      read (21,'(i3,i4,2f6.2,f6.3,f8.0,a20)')
     &           ili,jli,rlong0,rlat0,schmidt,ds,title
      write(6,*) ili,jli,rlong0,rlat0,schmidt,ds,title
c     read in g*zs(il,jl) to formatted file
      ilout=min(il,30)
      write (formout,'(''(''i3''f7.0)'')')ilout   !  i.e. (<il>f7.0)
      read(21,formout) zs

c     write out land/sea mask to formatted file
      write (formout,'(''(''i3''f4.1)'')')ilout   !  i.e. (<il>f4.1)
      read(21,formout) lsm

c     write out std.dev of top. to formatted file
      write (formout,'(''(''i3''f6.0)'')')ilout   !  i.e. (<il>f6.0)
      read(21,formout) tsd
      print *,'zss, rmsk and tsd written to unit',luout,' for il=',il

! usual veg data and map proj info
      open (22,file=vegin,status='old')
      read (22,'(i3,i4,2f6.2,f6.3,f8.0,a20)')
     &           ili,jli,rlong0,rlat0,schmidt,ds,title
      write(6,*) ili,jli,rlong0,rlat0,schmidt,ds,title

      write (formout,'(1h(,i3,2hi3,1h))')ili   !  i.e. (<il>i3)
      write(6,*)"formout=",formout
      read (22,formout) iarr     ! read old vegetation type

      write(6,*)"old id,jd,zs,lsm,tsd,veg=",id,jd
     &             ,zs(ijd),lsm(ijd),tsd(id,jd),iarr(id,jd)

! open output file (vegout)
      open (23,file=vegout,status='unknown')
      open (24,file=topout,status='unknown')

! now start good stuff
      call setxyz

!     initialize summing array
      do ntype=0,20
       do iq=1,ifull
        nsum(iq,1,ntype)=0
       enddo ! iq
      enddo ! ntype
! preinit. arr array with old data
       do iq=1,ifull
         arr(iq)=iarr(iq,1)
       enddo ! iq

!      do j=1,jl
!        write(6,*)"il/2,j,iarr(il/2)=",il/2,j,iarr(il/2,j)
!      enddo ! iq

!     is = nint((rlons - gslon)/dl) + 1 
!     js = nint((gslat - rlats)/dl) + 1
!     write(6,*)"is=",is," js=",js

      gridx=-1.e35
      rlatx=-1.e35
      rlonx=-1.e35
      gridn=+1.e35
      rlatn=+1.e35
      rlonn=+1.e35

      do j=1,jl
       do i=1,il
         n=i+(j-1)*il
         grid(i,j)=(ds/em(i,j))/1.e3 ! km
         rlond(i,j)=rlong(n)*rtd
         rlatd(i,j)=rlat (n)*rtd
         if(rlond(i,j).gt.180.)rlond(i,j)=rlond(i,j)-360.
         gridx=max(gridx,grid(i,j))
         gridn=min(gridn,grid(i,j))
         !write(6,*)i,j,rlond(i,j),rlatd(i,j)
         if ( grid(i,j).lt.ds1km ) then    ! use 1km when grid spc < 60 km
           rlatx=max(rlatx,rlatd(i,j))
           rlonx=max(rlonx,rlond(i,j))
           rlatn=min(rlatn,rlatd(i,j))
           rlonn=min(rlonn,rlond(i,j))
         endif ! grid < 60 km
       enddo ! il
       !write(6,*)j,rlond(1,j),rlond(il,j)
       !write(6,*)j,rlatd(1,j),rlatd(il,j)
      enddo ! jl

      imin = nint((rlonn - gslon)/dl) + 1
      imax = nint((rlonx - gslon)/dl) + 1
      jmin = nint((gslat - rlatx)/dl) + 1
      jmax = nint((gslat - rlatn)/dl) + 1

      write(6,*)"grid n,x=",gridn,gridx
      write(6,*)"rlon n,x=",rlonn,rlonx
      write(6,*)"rlat n,x=",rlatn,rlatx
      write(6,*)"i    n,x=",imin,imax
      write(6,*)"j    n,x=",jmin,jmax

      write(6,*)"loop over all valid rows in 1 km data set"
      do jg=jmin,jmax

        glat=gslat-(jg-1)*dl

        if (glat.gt.rlatn .and. glat.lt.rlatx) then
          if(mod(jg,10).eq.0)then
            write(6,*)"jg,glat,rlatn,glat,rlatx=",jg,glat,rlatn,glat,rlatx
          endif! (mod(jg,10).eq.0)then

! read each row of relevant data from 1 km data set (row=jg)
          read(20,rec=jg) (zza(ig),ig=1,mx)

          do ig = 1 , mx
!           if(zza(ig).gt.12) write(6,*)"Big input value ",ig,jg,zza(ig)
          enddo
          !write(6,*) jg,glat,rlatn,rlatx
          !write(6,'(100i2)') (zza(ig),ig=1,mx,500)

          do ig=1,mx

            glon = gslon + (ig-1)*dl
            if (glon.gt.rlonn .and. glon.lt.rlonx) then

              iz(ig)=zza(ig)
              if(iz(ig).ge.19)iz(ig)=0
! compute model grid i,j
              call latltoij(glon,glat,alci,alcj,nface)  ! con-cubic/octagon
              lci = nint(alci)
              lcj = nint(alcj)
! convert to "double" (i,j) notation
              lcj=lcj+nface*il

! check to make sure within grid dimensions
              if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
                nsum(lci,lcj,iz(ig))=nsum(lci,lcj,iz(ig))+1
              endif!(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

              if(lci.eq.id.and.lcj.eq.jd)then
                write(6,'("ig=",i6," jg=",i6," glon=",f7.2," glat=",f6.2
     &           ," i=",i4," j=",i4," iz=",i4," nsum=",i6)')
     &            ig,jg,glon,glat,lci,lcj,iz(ig),nsum(lci,lcj,iz(ig))
!                write(6,*)"ig=",ig," jg=",jg," glon=",glon," glat=",glat
!    &           ," i=",lci," j=",lcj," iz=",iz(ig)," nsum="
!    &           ,nsum(lci,lcj,iz(ig))
              endif

            endif ! glon

          enddo ! ig

        endif ! glat

      enddo ! jg
! finished loop over all valid rows in 1 km data set

      write(6,*)"(nsum(id,jd,nn),nn=1,20)=",(nsum(id,jd,nn),nn=0,20)

! close global datafile
      close(20)

      write(6,*)"old veg for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(iarr(i,j),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"old lsm for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(nint(lsm(i+(j-1)*il)),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"nsum"
      do j=96,1,-1
       write(6,'(49i2)')j,(nsum(i,j,12),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"now update vegie by most popular one ifull=",ifull

!-----------------------------------------------------------------------
! loop over all model grid points
      do iq=1,ifull
!-----------------------------------------------------------------------

       if ( lsm(iq).gt..5 .and. iarr(iq,1).eq.0 ) then
         write(6,*)"Input problem iq=",iq," lsm=",lsm(iq)
     &                           ," iarr=",iarr(iq,1)
         stop
       endif ! ( lsm(iq).lt..5 .and. iarr(i,1).eq.0 ) then

       ntypmax=0
       nsummax(iq)=0
       npoints(iq)=0
       nwater=0
       waterfrac(iq) = 0.
       if ( lsm(iq) .gt. .1 ) waterfrac(iq) = 1.

       do ntype=0,20
        !if(ntype.eq.12) write(6,*)iq,ntype,nsum(iq,1,ntype),nsummax(iq)
! sum up number of water input data points
        if(nsum(iq,1,0).gt.0)nwater=nwater+nsum(iq,1,0)
! pick veg type with largest number of input values
        if(nsum(iq,1,ntype).gt.nsummax(iq))then
          nsummax(iq)=nsum(iq,1,ntype)
          ntypmax=ntype
        endif ! (nsum(iq,1,ntype).gt.nsummax(iq))then
! sum up all 1 km data inputs for this grid box
        npoints(iq)=npoints(iq)+nsum(iq,1,ntype)
       enddo ! ntype

! only for points which have some input data!!!
       if ( npoints(iq) .gt. 1 ) then

         waterfrac(iq) = real(nwater)/real(npoints(iq))

! fix land-sea mask
         if ( ntypmax .lt. 1 ) lsm(iq)=0

         if ( ntypmax .gt. 0 ) then
!          print *,'iq,iarr(iq,1),ntypmax,nsummax(iq) ',
!    &            iq,iarr(iq,1),ntypmax,nsummax(iq)

! here we specify new veg type using high-res data
           iarr(iq,1) = ntypmax+31

           if(iarr(iq,1).gt.43)then
             write(6,*)"BIG IARR iq=",iq," iarr=",iarr(iq,1)
           endif
           arr(iq) = iarr(iq,1)
           lsm(iq)=1
         else ! ntypmax
! water point
           iarr(iq,1)=0
           arr(iq) = iarr(iq,1)
           lsm(iq)=0
         endif ! ntypmax

       endif ! ( npoints(iq) .gt. 0 ) then

       if ( lsm(iq).gt..5 .and. iarr(iq,1).eq.0 ) then
         write(6,*)"Output problem iq=",iq," lsm=",lsm(iq)
     &                           ," iarr=",iarr(iq,1)
       endif ! ( lsm(iq).lt..5 .and. iarr(i,1).eq.0 ) then

!-----------------------------------------------------------------------
      enddo   ! iq loop
!-----------------------------------------------------------------------

      write(6,*)"new id,jd,zs,lsm,tsd,veg=",id,jd
     &             ,zs(ijd),lsm(ijd),tsd(id,jd),iarr(id,jd)

!     do j=1,jl
!        write(6,*)"il/2,j,iarr(il/2)=",il/2,j,iarr(il/2,j)
!     enddo ! iq

      write(6,*)"new veg for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(iarr(i,j),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

      write(6,*)"new lsm for high res panel (j=1:96)"
      do j=96,1,-1
       write(6,'(49i2)')j,(nint(lsm(i+(j-1)*il)),i=1,48)
      enddo
      write(6,'(49i2)')j,(i,i=1,48)

!     write(6,*)"veg for high res panel (j=49:96)"
!     do j=96,1,-1
!      write(6,'(49i2)')j,(iarr(i,j),i=1,48)
!     enddo
!     write(6,'(49i2)')j,(i,i=1,48)

      write(23,'(i3,i4,2f6.2,f6.3,f8.0,''  veg type'')')
     &            il,jl,rlong0,rlat0,schmidt,ds
      write(23,formout) iarr     ! write vegetation type
      close(23)

      write(24,'(i3,i4,2f6.2,f6.3,f8.0,'' orog-mask-var'')')
     &            il,jl,rlong0,rlat0,schmidt,ds
      ilout=min(il,30)
      write (formout,'(''(''i3''f7.0)'')')ilout   !  i.e. (<il>f7.0)
      write(24,formout) zs

c     write out land/sea mask to formatted file
      write (formout,'(''(''i3''f4.1)'')')ilout   !  i.e. (<il>f4.1)
      write(24,formout) lsm

c     write out std.dev of top. to formatted file
      write (formout,'(''(''i3''f6.0)'')')ilout   !  i.e. (<il>f6.0)
      write(24,formout) tsd
      print *,'zss, rmsk and tsd written for il=',il


      call ncdf_setup("veg.nc",idnc,3,il,jl,1,20020101,0000)
      call nc2out(grid ,il,jl,1,1.,idnc,"grid","grid","km",0.,5000.)
      call nc2out(rlond,il,jl,1,1.,idnc,"lon","lon","degrees",0.,360.)
      call nc2out(rlatd,il,jl,1,1.,idnc,"lat","lat","degrees",-90.,90.)
      call nc2out(arr  ,il,jl,1,1.,idnc,"veg","veg","none",-1.,50.)
      call nc2out(zs   ,il,jl,1,1.,idnc,"zs ","zs ","m",-1.,50000.)
      call nc2out(lsm  ,il,jl,1,1.,idnc,"lsm","lsm","none",-1.,1.)
      call nc2out(waterfrac,il,jl,1,1.,idnc,"wfrc","wfrc","none",0.,1.)
      call nc2out(npoints,il,jl,1,1.,idnc,"npts","npts","none",0.,6000.)
      call nc2out(nsummax,il,jl,1,1.,idnc,"nmax","nmax","none",0.,6000.)
      call ncclos(idnc,ier)

      stop
      end
!===========================================================================
!     include 'setxyz.f'      ! for conformal-cubic
!     include 'jimcc.f'       ! for conformal-cubic
!     include 'latltoij.f'    ! for conformal-cubic
c=======================================================================
      subroutine ncdf_setup(ofile,idnc,ndim,il,jl,kl,kdate,ktime)

      include "netcdf.inc"   ! comment out on atmos

      character*(*) ofile
* netCDF id
      integer  idnc
      integer  ier
* dimension ids
      integer  dimil,dimjl,dimkl,dimtim
* variable ids
      integer idil,idjl,idkl,idnt
* input variables
      integer ndim,il,jl,kl
      character timorg*20
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/

      idnc = nccre ( ofile,ncclob,ier )
      print *,'###### create netcdf idnc,ndim,file=',idnc,ndim,ofile

c define dimensions
      print *,'define dimensions'
      dimil = ncddef(idnc,'longitude',il,ier)
      dimjl = ncddef(idnc,'latitude',jl,ier)
      if ( ndim.eq.4 ) then
        dimkl = ncddef(idnc,'dim3',kl,ier)
        dimtim = ncddef(idnc,'time',NCUNLIM,ier)
      elseif ( ndim.eq.3 ) then
        dimtim = ncddef(idnc,'time',NCUNLIM,ier)
      endif ! ( ndim.eq.4 ) then

c define variables
      print *,'define variables'
      idil = ncvdef(idnc,'longitude',NCFLOAT,1,dimil,ier)
      call ncaptc(idnc,idil,'point_spacing',NCCHAR,4,'even',ier)
      call ncaptc(idnc,idil,'units',NCCHAR,12,'degrees_east',ier)
      idjl = ncvdef(idnc,'latitude',NCFLOAT,1,dimjl,ier)
      call ncaptc(idnc,idjl,'point_spacing',NCCHAR,4,'even',ier)
      call ncaptc(idnc,idjl,'units',NCCHAR,13,'degrees_north',ier)

      if ( ndim.ge.3 ) then
        if ( ndim.eq.4 ) then
          idkl = ncvdef(idnc,'dim3',NCFLOAT,1,dimkl,ier)
          call ncaptc(idnc,idkl,'units',NCCHAR,5,'layer',ier)
          call ncaptc(idnc,idkl,'positive',NCCHAR,4,'down',ier)
        endif ! ( ndim.eq.4 ) then
        idnt = ncvdef(idnc,'time',NCFLOAT,1,dimtim,ier)
        call ncaptc(idnc,idnt,'units',NCCHAR,5,'hours',ier)
        call ncaptc(idnc,idnt,'point_spacing',NCCHAR,4,'even',ier)
        icy=kdate/10000
        icm=(kdate-icy*10000)/100
        icd=(kdate-icy*10000-icm*100)
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(6,*) icy,icm,icd,ich,icmi,ics
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &               icd,month(icm),icy,ich,icmi,ics
        print *,'timorg=',timorg
        call ncaptc(idnc,idnt,'time_origin',NCCHAR,20,timorg,ier)
      endif ! ( ndim.ge.3 ) then

      call ncendf(idnc,ier)

      return
      end
c=======================================================================
      subroutine nc2out(data,il,jl,nt,thr,idnc,sname,lname,units,dn,dx)

      logical debug
!     parameter ( debug=.false. )
      parameter ( debug=.true. )
      logical oscale
      parameter ( oscale=.false. )

c this program write out a 2-dim array of data with unlimited time dim.
c in netcdf format

      include "netcdf.inc"        ! comment out on atmos

      character*(*) sname
      character*(*) lname
      character*(*) units
* netCDF id
      integer  idnc
      integer  ier
* dimension ids
      integer  dimil,dimjl,dimkl,dimtim
* variable ids
      integer  idil,idjl,idkl,idnt
      integer  idvar
* variable shapes
      integer idims(3)
* corners and edge lengths
      integer corner(3),edges(3)
* data variables
      real data(il,jl)
* packing variables
      parameter (nx=4800,ny=6000)
      integer*2 ipack(nx*ny)
      integer*2 vrange(2)
      data vrange/-32500,32500/
* dim variables
      real xdim(nx),ydim(ny)

      save vrange

      write(6,'("*** nc2out *** idnc,nt,thr,sname=",2i4,f7.2,a6)')
     &          idnc,nt,thr,sname

      if ( debug ) then
        print *,'il,jl,nt,idnc=',il,jl,nt,idnc
        print *,'lname=',lname
        print *,'units=',units
        print *,'sx,sn=',dx,dn
      endif ! debug

      if ( il*jl .gt. nx*ny ) then
        print *,'il*jl=',il*jl,' too big in nc2out'
        stop
      endif

c***********************************************************************
      if ( nt.eq.1 ) then
c only do the following for the first call to each variable
c***********************************************************************

        call ncredf(idnc,ier)
        if ( debug ) print *,'ncredf ier=',ier

c get dimension ids
        dimil = ncdid(idnc,'longitude',ier)
        dimjl = ncdid(idnc,'latitude',ier)
        dimtim = ncdid(idnc,'time',ier)
        idims(1) = dimil
        idims(2) = dimjl
        idims(3) = dimtim
        if ( debug ) print *,'idims=',idims

c create variable using short name
        write(6,*)"thr=",thr
        if ( thr.ge.0. ) then
          ndims=3
        else
          ndims=2
        endif ! ( thr.ge.0. ) then
          if(oscale)then
            idvar = ncvdef (idnc,sname,ncshort,ndims,idims,ier)
          else
            idvar = ncvdef (idnc,sname,ncfloat,ndims,idims,ier)
          endif
        if(debug)write(6,*)'time idvar,sn,ndims=',idvar,sname,ndims

c give it a long name
        call ncaptc(idnc,idvar,'long_name',NCCHAR
     &             ,lngstr(lname),lname,ier)
        if ( debug ) print *,'lname ier=',ier

c give it units
        call ncaptc(idnc,idvar,'units',NCCHAR
     &             ,lngstr(units),units,ier)
        if ( debug ) print *,'units ier=',ier

        if(oscale)then
c give it a fill value
          call ncapt(idnc,idvar,'_FillValue',NCSHORT,1,vrange(1),ier)
          if ( debug ) print *,'fill ier=',ier
c define valid scaled max/min variables
          scalef= (dx-dn)/float(vrange(2)-vrange(1))
          addoff= dn-scalef*float(vrange(1))
          if ( debug ) then
            write(6,'("sx,sn,a,s=",1p,6e12.3)')dx,dn,addoff,scalef
          endif ! debug
          call ncapt(idnc,idvar,'add_offset',NCFLOAT,1,addoff,ier)
          call ncapt(idnc,idvar,'scale_factor',NCFLOAT,1,scalef,ier)
          call ncaptc(idnc,idvar,'FORTRAN_format',ncchar,5,'G11.4',ier)
          call ncapt(idnc,idvar,'valid_range',ncshort,2,vrange,ier)
        endif!(oscale)then

        if ( debug ) print *,'done attrib'

* leave define mode
        call ncendf(idnc,ier)

* store xdim
        do i=1,il
          xdim(i)=float(i)
        enddo
        idil = ncvid(idnc,'longitude',ier)
c       print *,'idil,xdim=',idil,(xdim(i),i=1,il)
        call ncvpt(idnc,idil,1,il,xdim,ier)

* store ydim
        do j=1,jl
          ydim(j)=float(j)
        enddo
        idjl = ncvid(idnc,'latitude',ier)
c       print *,'idjl,ydim=',idjl,(ydim(i),i=1,jl)
        call ncvpt(idnc,idjl,1,jl,ydim,ier)

c***********************************************************************
      endif ! nt=1
c***********************************************************************

      corner(1) = 1
      corner(2) = 1
      corner(3) = nt
      edges(1) = il
      edges(2) = jl
      edges(3) = 1
      if ( debug ) then
        print *,'corner=',corner
        print *,'edges=',edges
      endif ! debug

      write(6,*)"thr=",thr
      if ( thr.ge.0. ) then
c set time to number of hours since start = thr
        idnt = ncvid(idnc,'time',ier)
        if ( debug ) print *,'ncdiv idnt,ier=',idnt,ier
        call ncvpt1(idnc,idnt,nt,thr,ier)
        if ( debug ) print *,'ncvpt nt,ier=',nt,ier
      endif ! ( thr.ge.0. ) then

c find variable index
      idvar = ncvid(idnc,sname,ier)
      if ( debug ) then
        print *,'ncdiv sname,idvar,ier=',sname,idvar,ier
      endif ! debug

      if(oscale)then
c get scaling factors for this variable
        scalef= (dx-dn)/float(vrange(2)-vrange(1))
        addoff= dn-scalef*float(vrange(1))
        if ( debug ) write(6,'("sx,sn,a,s=",6f10.2)')dx,dn,addoff,scalef
c pack data into ipack
        adx=-1.e29
        adn= 1.e29
        idx= vrange(1)
        idn= vrange(2)
        do j=1,jl
          do i=1,il
           n=i+(j-1)*il
           adx=max(adx,data(i,j))
           adn=min(adn,data(i,j))
           ipack(n)=nint((data(i,j)-addoff)/scalef)
           idx=max(idx,(ipack(n)))
           idn=min(idn,(ipack(n)))
          end do ! i
        end do ! j
        if(adx.gt.dx)print *,'## actual dx .gt. set dx ##',adx,dx
        if(adn.lt.dn)print *,'## actual dn .lt. set dn ##',adn,dn
        if ( debug ) then
          print *,'adx,adn=',adx,adn
          print *,'idx,idn=',idx,idn
        endif ! debug
c put packed data into cdf file
        call ncvpt(idnc,idvar,corner,edges,ipack,ier)
      else!(oscale)then
        call ncvpt(idnc,idvar,corner,edges,data,ier)
      endif!(oscale)then
      if ( debug ) print *,'ier=',ier

      return
      end
c=======================================================================
      function lngstr( string )
      character*(*) string
      ilen = len(string)
c     print*,'string=',string
c     print*,'ilen=',ilen
      do 100 lngstr=ilen,1,-1
        if ( string(lngstr:lngstr) .ne. ' ' ) go to 99
  100 continue
      lngstr = 0
   99 continue
      return
      end
c=======================================================================
