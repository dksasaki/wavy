!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
!
program umwm_2d_simple

use mod_const
use mod_nondimensional
use mod_precision,only:rk => realkind, ik => intkind
use mod_spectrum,only:spectrum_type
use mod_domain,only:domain_type
use mod_grid,only:grid_type
use mod_spectral_shapes
use mod_source_functions,only:sin_DCCM2012,sds_DCCM2012,sdt_DCCM2012,snl_DCCM2012
use mod_time_integration,only:backward_euler,forward_euler,exact_exponential
use mod_parallel, only: tile_indices, tile_neighbors_2d, num_tiles
use mod_advection
use mod_utility,only:tile,range,ones
use netcdf

implicit none

type(spectrum_type),dimension(:,:),allocatable :: spec
!type(spectrum_type) :: spectrum
type(domain_type) :: domain,tend, domain_parent, domain2, dfdt
type(grid_type) :: spatial_grid
type(domain_type) :: sin_dom,sds_dom,sdt_dom,sbf_dom,snl_dom

real :: t0,t1

real(kind=rk) :: wspd
real(kind=rk) :: fetch
real(kind=rk),parameter :: grav = 9.8_rk
integer :: i,m,n,idm,nfreqs,jp,ip, j

real(kind=rk),dimension(:),allocatable :: freq
real(kind=rk),dimension(:,:),allocatable :: f
real(kind=rk),dimension(:,:),allocatable :: cg
real(kind=rk),dimension(:),allocatable :: dx
real(kind=rk),dimension(:,:),allocatable :: adv

real, dimension(:,:,:), allocatable :: swh

integer :: varid, ncid, xid, yid, ncid2
integer :: stat
integer :: swhid
integer :: xdimid,ydimid,tdimid, imid, dxid

integer(kind=ik),dimension(2) :: halowidth
integer(kind=ik), dimension(2) :: lb,ub, lbp, ubp
integer(kind=ik), dimension(4) :: indices
integer(kind=ik) :: dul1, dul2
real(kind=ik), allocatable, dimension(:,:,:,:) :: f1, f3
real(kind=ik), allocatable, dimension(:,:,:,:), codimension[:] :: f2
real(kind=ik), allocatable, dimension(:,:), codimension[:] :: nimage
character(len=2) :: str4

integer(kind=rk) :: xs, ys

halowidth=1

idm = 20

allocate(f2(43,32,50,20)[*])

!-------------------------------------------------------------------------------
! upper and lower indices of the tiles with (mther grid is the reference)
indices = tile_indices([50, 20])
lb(1) = 1
ub(1) = 50
lb(2) = 1
ub(2) = 20
! size of tiles in each direction
dul1 = abs(indices(2) -indices(1)) +1
dul2 = abs(indices(3) -indices(4)) +1

! upper and lower indices on each tile
lbp(1) = 1
ubp(1) = dul1
lbp(2) = 1
ubp(2) = dul2

allocate(f1(43,32,dul1+2,dul2+2))
!-------------------------------------------------------------------------------

allocate(nimage(50,20)[*])


write(*,*)
write(*,*)'Initialize omnidirectional spectrum with JONSWAP shape;'
write(*,*)

wspd = 1e1_rk
fetch = 1e5_rk

! Initialize spatial grid
spatial_grid = grid_type(lb,ub,&
                 dx = 1e3_rk*tile(ones(50,WAVY_INT),20),&
                 dy = 1e3_rk*tile(ones(50,WAVY_INT),20))

! Initialize domain
domain_parent = domain_type(grid = spatial_grid,&
                spectrum = spectrum_type(0.0313_rk,2._rk,1.1_rk,32,1e3_rk))

domain = domain_type(grid = grid_type([lbp(1),lbp(2)],[ubp(1),ubp(2)],&
                       dx = 1e3_rk*tile(ones(ubp(1),WAVY_INT),ubp(2)),&
                       dy = 1e3_rk*tile(ones(ubp(1),WAVY_INT),ubp(2))),&
                 spectrum = spectrum_type(0.0313_rk,2._rk,1.1_rk,32,1e3_rk))

domain2 = domain_type(grid = grid_type([lbp(1),lbp(2)],[ubp(1)+2,ubp(2)+2],&
                      dx = 1e3_rk*tile(ones(ubp(1),WAVY_INT),ubp(2)),&
                      dy = 1e3_rk*tile(ones(ubp(1),WAVY_INT),ubp(2))),&
                spectrum = spectrum_type(0.0313_rk,2._rk,1.1_rk,32,1e3_rk))


spec = domain_parent % getSpectrum()

spec(25,10) = donelanHamiltonHuiDirectionalSpectrum(f = domain_parent % getFrequency(),&
                            theta = domain_parent % getDirections(),&
                            wspd  = wspd,&
                            fpeak = jonswapPeakFrequency(5._rk,fetch,grav),&
                            theta_mean = 0._rk,&
                            grav  = grav)




sin_dom = domain
sds_dom = domain
snl_dom = domain
sdt_dom = domain
sbf_dom = domain
tend = domain

domain = spec(indices(1):indices(2),indices(3):indices(4))

call cpu_time(t0)
! call domain % writeJSON('domain.json',minify=.true.)
call cpu_time(t1)
write(*,*)'domain.json elapsed',t1-t0,'seconds'

write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss'&
                //'      m0(f)    m1(f)    m2(f)'
write(*,fmt='(a)')'----------------------------------------------------------'&
                //'----------------------'

do n = 1,10

  ip = 1
  jp = 5

  spec = domain % getSpectrum()
  !
  ! ! wind input
  ! sin_dom = sin_DCCM2012(domain % getSpectrum(),wspd,wdir=0._rk,&
  !                        input_height=10._rk,ustar=1e-2_rk,vonkarman=0.4_rk)
  !
  ! ! wave dissipation
  ! sds_dom = sds_DCCM2012(domain % getSpectrum(),sds_coefficient=42._rk,&
  !                        sds_power=2.5_rk,mss_coefficient=120._rk)
  !
  ! ! dissipation from turbulence
  ! sdt_dom = sdt_DCCM2012(domain % getSpectrum(),1e-2_rk,1e-2_rk)
  !
  ! ! non-linear wave-wave transfer
  ! snl_dom = snl_DCCM2012(domain % getSpectrum(),sds_dom % getSpectrum(),5._rk)
  !
  ! domain = exact_exponential(domain,sin_dom-sds_dom-sdt_dom,dt=10._rk)

  write(*,fmt='(9(f8.4,1x),(5(4I3)))')wspd,spec(ip,jp) % significantWaveHeight(),&
    1./spec(ip,jp) % peakFrequency(),spec(ip,jp) % meanPeriod(),&
    spec(ip,jp) % meanPeriodZeroCrossing(),spec(ip,jp) % meanSquareSlope(),&
    spec(ip,jp) % frequencyMoment(0),spec(ip,jp) % frequencyMoment(1),&
    spec(ip,jp) % frequencyMoment(2), shape(spec), this_image(), shape(spec(4,4) % getSpectrum())

  f1(:,:,2:dul1+1,2:dul2+1) = domain % getSpectrumArray()
  ! f1 = domain % getSpectrumArray(halowidth,.false.)
  !

  ! f1(:,:,1,:)=1000
  ! f1(:,:,dul1+1,:)=1000
  ! f1(:,:,:,1)=1000
  ! f1(:,:,:,dul2+1)=1000
  call domain2  % sync_edges(f1, [0,0])
  call domain2 % setSpectrumArray(f1) !
  write(*,*) maxval(domain2 % dx)

  dfdt = domain2 % advect(advectUpwind1stOrder2dRank2,[1,1])
  if(this_image()==1) write(*,*) 'dfdt', shape(dfdt % getSpectrum())

  !
  sync all
  !attention to nan: dx and dy are not updated, which caused nan to appear in the results
  domain2 = forward_euler(domain2, dfdt, dt=10._rk)
  !
  sync all

  f3 = domain2 % getSpectrumArray()
  ! write(*,*)(f3(:,:,dul1+1,:))
  ! write(*,*)(f3(:,:,:,dul2+1))

  if(this_image()==1) write(*,*) 'dom3, f3',shape(domain % getSpectrum())
  if(this_image()==1) write(*,*) 'dom',shape(domain % getSpectrum())
  ! if(this_image()==1) write(*,*)shape(f3(:,:,2:dul1+1, 2:dul2+1))
  if(this_image()==1) write(*,*) indices, dul1, dul2



  f2(:,:,indices(1):indices(2), indices(3):indices(4))[1] = f3(:,:,2:dul1+1, 2:dul2+1)
  nimage(indices(1):indices(2), indices(3):indices(4))[1] = this_image()
  call domain % setSpectrumArray(f3(:,:,2:dul1+1, 2:dul2+1))
  ! call domain % setSpectrumArray(f3)

  ! f3 = dfdt % getSpectrumArray()
  sync all

  ! if (this_image()==1)then
  !
  !   xs = size(f3, dim=3)
  !   ys = size(f3, dim=4)
  !   ! write(*,*) maxval(f2(:,:,:,:)[1]), this_image(), n
  ! ! if (n==1)then
  !   write(str4,'(I0.2)') n
  !   ! call domain_parent % setSpectrumArray(f3)
  !   call nc_check(nf90_create('output2'//str4//'.nc',nf90_clobber,ncid2))
  !
  !   call nc_check(nf90_def_dim(ncid2,'x',xs,xdimid))
  !   call nc_check(nf90_def_dim(ncid2,'y',ys,ydimid))
  !   call nc_check(nf90_def_dim(ncid2,'t',1,tdimid))
  !   stat = nf90_def_var(ncid2,'swh',nf90_float,[xdimid,ydimid,tdimid],swhid)
  !   stat = nf90_def_var(ncid2,'image',nf90_float,[xdimid,ydimid],imid)
  !   stat = nf90_def_var(ncid2,'dx',nf90_float,[xdimid,ydimid],dxid)
  !
  !   call nc_check(nf90_enddef(ncid2))
  !
  !   stat = nf90_put_att(ncid2,swhid,name='description',values='Significant wave height')
  !   stat = nf90_put_att(ncid2,swhid,name='units',values='m')
  !
  !   stat = nf90_put_att(ncid2,swhid,name='description',values='dx')
  !   stat = nf90_put_att(ncid2,swhid,name='units',values='m')
  !
  !
  !   stat = nf90_put_var(ncid2,swhid,reshape(dfdt % significantWaveHeight(),[xs,ys,1]),start=[1,1,1],count=[xs,ys,1])
  !   stat = nf90_put_var(ncid2,imid, nimage,start=[1,1],count=[xs,ys])
  !   stat = nf90_close(ncid2)
  ! ! end if
  ! end if

if (this_image()==1)then
  ! write(*,*) maxval(f2(:,:,:,:)[1]), this_image(), n
! if (n==1)then
  write(str4,'(I2.2)') n
  call domain_parent % setSpectrumArray(f2(:,:,:,:)[1])
  call nc_check(nf90_create('output2'//str4//'.nc',nf90_clobber,ncid2))

  call nc_check(nf90_def_dim(ncid2,'x',50,xdimid))
  call nc_check(nf90_def_dim(ncid2,'y',20,ydimid))
  call nc_check(nf90_def_dim(ncid2,'t',1,tdimid))
  stat = nf90_def_var(ncid2,'swh',nf90_float,[xdimid,ydimid,tdimid],swhid)
  stat = nf90_def_var(ncid2,'image',nf90_float,[xdimid,ydimid,tdimid],imid)
  stat = nf90_def_var(ncid2,'dx',nf90_float,[xdimid,ydimid,tdimid],dxid)

  call nc_check(nf90_enddef(ncid2))

  stat = nf90_put_att(ncid2,swhid,name='description',values='Significant wave height')
  stat = nf90_put_att(ncid2,swhid,name='units',values='m')
  stat = nf90_put_att(ncid2,dxid,name='description',values='dx')

  write(*,*) shape(domain_parent % getGridSpacingYWithHalo([1,1],.false.))
  stat = nf90_put_var(ncid2,dxid,reshape(domain_parent % getGridSpacingXWithHalo([0,0],.false.),[50,20,1]),&
          start=[1,1,1],count=[50,20,1])
  stat = nf90_put_var(ncid2,swhid,reshape(domain_parent % significantWaveHeight(), [50,20,1]),start=[1,1,1],count=[50,20,1])
  stat = nf90_put_var(ncid2,imid, nimage,start=[1,1,1],count=[50,20,1])
  stat = nf90_close(ncid2)
! end if
end if
end do


contains
  SUBROUTINE nc_check(stat,varb)
!======================================================================!
!                                                                      !
! DESCRIPTION: Checks for NetCDF errors and if any, print and abort    !
!                                                                      !
!======================================================================!

INTEGER,INTENT(IN) :: stat
CHARACTER,INTENT(IN), optional :: VARB
!======================================================================!

IF(stat /= nf90_noerr)THEN
if (present(varb))then
  WRITE(*,*)(VARB)
endif
WRITE(*,*)'Error in NetCDF I/O'
WRITE(*,*)TRIM(nf90_strerror(stat))
STOP
ENDIF

ENDSUBROUTINE nc_check

endprogram umwm_2d_simple
