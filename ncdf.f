      subroutine ncdf_setup(ofile,idnc,ndim,il,jl,kl,kdate,ktime,
     &                      rlong0,rlat0,schmidt)

      implicit none

      include "netcdf.inc"   ! comment out on atmos
      
      integer ics,icmi,ich,icd
      integer icy,icm

      character*(*) ofile
* netCDF id
      integer  idnc
      integer  ier
* dimension ids
      integer  dimil,dimjl,dimkl,dimtim
* variable ids
      integer idil,idjl,idkl,idnt
* input variables
      integer kdate,ktime
      integer ndim,il,jl,kl
      real rlong0,rlat0,schmidt
      character timorg*20
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/

#ifdef usenc3
      idnc = nccre ( ofile,ncclob,ier )
#else
      ier=nf_create(ofile,NF_NETCDF4,idnc)
#endif
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

      ier=nf_put_att_real(idnc,nf_global,'lon0',nf_real,1,rlong0)
      ier=nf_put_att_real(idnc,nf_global,'lat0',nf_real,1,rlat0)
      ier=nf_put_att_real(idnc,nf_global,'schmidt',nf_real,1,schmidt)

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

      write(6,'("*** nc2out *** idnc,nt,thr,sname=",2i8,f7.2,a6)')
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
