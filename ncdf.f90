      subroutine ncdf_setup(ofile,idnc,ndim,il,jl,kl,kdate,ktime,rlong0,rlat0,schmidt)

      implicit none

      include "netcdf.inc"   ! comment out on atmos
      
      integer ics,icmi,ich,icd
      integer icy,icm

      character*(*) ofile
! netCDF id
      integer  idnc
      integer  ier
! dimension ids
      integer  dimil,dimjl,dimkl,dimtim
! variable ids
      integer idil,idjl,idkl,idnt
! input variables
      integer kdate,ktime
      integer ndim,il,jl,kl
      integer, dimension(1) :: dimids
      real rlong0,rlat0,schmidt
      character timorg*20
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/

#ifdef usenc3
      !idnc = nccre ( ofile,ncclob,ier )
      ier = nf_create(ofile,nf_64bit_offset,idnc)
#else
      ier=nf_create(ofile,NF_NETCDF4,idnc)
#endif
      print *,'###### create netcdf idnc,ndim,file=',idnc,ndim,ofile

! define dimensions
      print *,'define dimensions'
      ier=nf_def_dim(idnc,'longitude',il,dimil)
      ier=nf_def_dim(idnc,'latitude',jl,dimjl)
      if ( ndim.eq.4 ) then
        ier=nf_def_dim(idnc,'dim3',kl,dimkl)
        ier=nf_def_dim(idnc,'time',nf_unlimited,dimtim)
      elseif ( ndim.eq.3 ) then
        ier=nf_def_dim(idnc,'time',nf_unlimited,dimtim)
      endif ! ( ndim.eq.4 ) then

! define variables
      print *,'define variables'
      dimids=dimil
      ier=nf_def_var(idnc,'longitude',nf_float,1,dimids,idil)
      ier=nf_put_att_text(idnc,idil,'point_spacing',4,'even')
      ier=nf_put_att_text(idnc,idil,'units',12,'degrees_east')
      dimids=dimjl
      ier=nf_def_var(idnc,'latitude',nf_float,1,dimids,idjl)
      ier=nf_put_att_text(idnc,idjl,'point_spacing',4,'even')
      ier=nf_put_att_text(idnc,idjl,'units',13,'degrees_north')

      if ( ndim.ge.3 ) then
        if ( ndim.eq.4 ) then
          dimids=dimkl
          ier=nf_def_var(idnc,'dim3',nf_float,1,dimids,idkl)
          ier=nf_put_att_text(idnc,idkl,'units',5,'layer')
          ier=nf_put_att_text(idnc,idkl,'positive',4,'down')
        endif ! ( ndim.eq.4 ) then
        dimids=dimtim
        ier=nf_def_var(idnc,'time',nf_float,1,dimids,idnt)
        ier=nf_put_att_text(idnc,idnt,'units',5,'hours')
        ier=nf_put_att_text(idnc,idnt,'point_spacing',4,'even')
        icy=kdate/10000
        icm=(kdate-icy*10000)/100
        icd=(kdate-icy*10000-icm*100)
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(6,*) icy,icm,icd,ich,icmi,ics
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
        print *,'timorg=',timorg
        ier=nf_put_att_text(idnc,idnt,'time_origin',20,timorg)
      endif ! ( ndim.ge.3 ) then

      ier=nf_put_att_real(idnc,nf_global,'lon0',nf_real,1,rlong0)
      ier=nf_put_att_real(idnc,nf_global,'lat0',nf_real,1,rlat0)
      ier=nf_put_att_real(idnc,nf_global,'schmidt',nf_real,1,schmidt)

      call ncendf(idnc,ier)

      return
      end
!=======================================================================
      subroutine nc2out(data,il,jl,nt,thr,idnc,sname,lname,units,dn,dx)

      logical debug
!     parameter ( debug=.false. )
      parameter ( debug=.true. )
      logical oscale
      parameter ( oscale=.false. )

! this program write out a 2-dim array of data with unlimited time dim.
! in netcdf format

      include "netcdf.inc"        ! comment out on atmos

      character*(*) sname
      character*(*) lname
      character*(*) units
! netCDF id
      integer  idnc
      integer  ier
! dimension ids
      integer  dimil,dimjl,dimkl,dimtim
! variable ids
      integer  idil,idjl,idkl,idnt
      integer  idvar
! variable shapes
      integer idims(3)
! corners and edge lengths
      integer corner(3),edges(3)
! data variables
      real data(il,jl)
      real rdx,rdn
! packing variables
      integer*2 ipack(il*jl)
      integer*2 vrange(2)
      integer*2 idx, idm
      integer, dimension(1) :: start, ncount
      data vrange/-32500,32500/
! dim variables
      real xdim(il),ydim(jl)
      real, dimension(1) :: rvals

      save vrange

      write(6,'("*** nc2out *** idnc,nt,thr,sname=",2i8,f7.2,a6)') idnc,nt,thr,sname

      if ( debug ) then
        print *,'il,jl,nt,idnc=',il,jl,nt,idnc
        print *,'lname=',lname
        print *,'units=',units
        print *,'sx,sn=',dx,dn
      endif ! debug

!***********************************************************************
      if ( nt.eq.1 ) then
! only do the following for the first call to each variable
!***********************************************************************

        call ncredf(idnc,ier)
        if ( debug ) print *,'ncredf ier=',ier

! get dimension ids
        ier = nf_inq_dimid(idnc,"longitude",dimil)
        ier = nf_inq_dimid(idnc,"latitude",dimjl)
        ier = nf_inq_dimid(idnc,"time",dimtim)
        idims(1) = dimil
        idims(2) = dimjl
        idims(3) = dimtim
        if ( debug ) print *,'idims=',idims

! create variable using short name
        write(6,*)"thr=",thr
        if ( thr.ge.0. ) then
          ndims=3
        else
          ndims=2
        endif ! ( thr.ge.0. ) then
          if(oscale)then
            ier = nf_def_var(idnc,sname,nf_short,ndims,idims,idvar)
          else
            ier = nf_def_var(idnc,sname,nf_float,ndims,idims,idvar)
          endif
        if(debug)write(6,*)'time idvar,sn,ndims=',idvar,sname,ndims

! give it a long name
        ier=nf_put_att_text(idnc,idvar,'long_name',len_trim(lname),lname)
        if ( debug ) print *,'lname ier=',ier

! give it units
        ier=nf_put_att_text(idnc,idvar,'units',len_trim(units),units)
        if ( debug ) print *,'units ier=',ier

        if(oscale)then
! give it a fill value
          ier=nf_put_att_int2(idnc,idvar,'_FillValue',nf_short,1,vrange(1))
          if ( debug ) print *,'fill ier=',ier
! define valid scaled max/min variables
          scalef= (dx-dn)/float(vrange(2)-vrange(1))
          addoff= dn-scalef*float(vrange(1))
          if ( debug ) then
            write(6,'("sx,sn,a,s=",1p,6e12.3)')dx,dn,addoff,scalef
          endif ! debug
          rvals=addoff
          ier=nf_put_att_real(idnc,idvar,'add_offset',nf_float,1,rvals)
          rvals=scalef
          ier=nf_put_att_real(idnc,idvar,'scale_factor',nf_float,1,rvals)
          ier=nf_put_att_text(idnc,idvar,'FORTRAN_format',5,'G11.4')
          ier=nf_put_att_int2(idnc,idvar,'valid_range',nf_short,2,vrange)
        endif!(oscale)then

        if ( debug ) print *,'done attrib'

! leave define mode
        ier=nf_enddef(idnc)

! store xdim
        do i=1,il
          xdim(i)=float(i)
        enddo
        ier = nf_inq_varid(idnc,'longitude',idil)
!       print *,'idil,xdim=',idil,(xdim(i),i=1,il)
        start = 1
        ncount = il
        ier = nf_put_vara_real(idnc,idil,start,ncount,xdim)

! store ydim
        do j=1,jl
          ydim(j)=float(j)
        enddo
        ier = nf_inq_varid(idnc,'latitude',idjl)
!       print *,'idjl,ydim=',idjl,(ydim(i),i=1,jl)
        start = 1
        ncount = jl
        ier = nf_put_vara_real(idnc,idjl,start,ncount,ydim)

!***********************************************************************
      endif ! nt=1
!***********************************************************************

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
! set time to number of hours since start = thr
        ier = nf_inq_varid(idnc,'latitude',idnt)
        if ( debug ) print *,'ncdiv idnt,ier=',idnt,ier
        start = nt
        ncount = 1
        rvals = thr
        ier = nf_put_vara_real(idnc,idnt,start,ncount,rvals)
        if ( debug ) print *,'ncvpt nt,ier=',nt,ier
      endif ! ( thr.ge.0. ) then

! find variable index
      ier = nf_inq_varid(idnc,sname,idvar)
      if ( debug ) then
        print *,'ncdiv sname,idvar,ier=',sname,idvar,ier
      endif ! debug

      if(oscale)then
! get scaling factors for this variable
        scalef= (dx-dn)/float(vrange(2)-vrange(1))
        addoff= dn-scalef*float(vrange(1))
        if ( debug ) write(6,'("sx,sn,a,s=",6f10.2)')dx,dn,addoff,scalef
! pack data into ipack
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
! put packed data into cdf file
        ier = nf_put_vara_int2(idnc,idvar,corner,edges,ipack)
      else!(oscale)then
        ier = nf_put_vara_real(idnc,idvar,corner,edges,data)
      endif!(oscale)then
      if ( debug ) print *,'ier=',ier

      return
      end
!=======================================================================
      function lngstr( string )
      character*(*) string
      ilen = len(string)
!     print*,'string=',string
!     print*,'ilen=',ilen
      do 100 lngstr=ilen,1,-1
        if ( string(lngstr:lngstr) .ne. ' ' ) go to 99
  100 continue
      lngstr = 0
   99 continue
      return
      end
