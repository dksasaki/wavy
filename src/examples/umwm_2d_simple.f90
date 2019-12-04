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
use mod_precision,only:rk => realkind
use mod_spectrum,only:spectrum_type
use mod_domain,only:domain_type
use mod_grid,only:grid_type
use mod_spectral_shapes
use mod_source_functions,only:sin_DCCM2012,sds_DCCM2012,sdt_DCCM2012,snl_DCCM2012
use mod_time_integration,only:backward_euler,forward_euler,exact_exponential
use mod_advection
use mod_utility,only:tile,range,ones
use netcdf

implicit none

type(spectrum_type),dimension(:,:),allocatable :: spec
!type(spectrum_type) :: spectrum
type(domain_type) :: domain,tend
!type(grid_type) :: grid
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
integer :: xdimid,ydimid,tdimid


idm = 20

write(*,*)
write(*,*)'Initialize omnidirectional spectrum with JONSWAP shape;'
write(*,*)

wspd = 3e1_rk
fetch = 1e5_rk

! Initialize domain
domain = domain_type(grid = grid_type([1,1],[50,20],&
  dx = 1e3_rk*tile(ones(50,WAVY_INT),20),&
  dy = 1e3_rk*tile(ones(50,WAVY_INT),20)),&
  spectrum = spectrum_type(0.0313_rk,2._rk,1.1_rk,32,1e3_rk))

sin_dom = domain
sds_dom = domain
snl_dom = domain
sdt_dom = domain
sbf_dom = domain


spec = domain % getSpectrum()
tend = domain

! Initialize spectrum arrays

!do concurrent(i = 1:idm)
  spec(5,5) = donelanHamiltonHuiDirectionalSpectrum(f = domain % getFrequency(),&
                              theta = domain % getDirections(),&
                              wspd  = wspd,&
                              fpeak = jonswapPeakFrequency(5._rk,fetch,grav),&
                              theta_mean = -1.17_rk,&
                              grav  = grav)

!enddo
domain = spec

call cpu_time(t0)
! call domain % writeJSON('domain.json',minify=.true.)
call cpu_time(t1)
write(*,*)'domain.json elapsed',t1-t0,'seconds'

write(*,fmt='(a)')'   wspd      Hs      Tp       Tm1      Tm2      mss'&
                //'      m0(f)    m1(f)    m2(f)'
write(*,fmt='(a)')'----------------------------------------------------------'&
                //'----------------------'

do n = 1,1000

  ip = 1
  jp = 5

  spec = domain % getSpectrum()

  ! wind input
  sin_dom = sin_DCCM2012(domain % getSpectrum(),wspd,wdir=1.17_rk,&
                         input_height=10._rk,ustar=1e-2_rk,vonkarman=0.4_rk)

  ! wave dissipation
  sds_dom = sds_DCCM2012(domain % getSpectrum(),sds_coefficient=42._rk,&
                         sds_power=2.5_rk,mss_coefficient=120._rk)

  ! dissipation from turbulence
  sdt_dom = sdt_DCCM2012(domain % getSpectrum(),1e-2_rk,1e-2_rk)

  ! non-linear wave-wave transfer
  snl_dom = snl_DCCM2012(domain % getSpectrum(),sds_dom % getSpectrum(),5._rk)

  !domain = forward_euler(domain,(sin_dom-sds_dom-sdt_dom)*domain+snl_dom,dt=1._rk)
  domain = exact_exponential(domain,sin_dom-sds_dom-sdt_dom,dt=10._rk)
  !domain = forward_euler(domain,snl_dom,dt=10._rk)

  write(*,fmt='(9(f8.4,1x))')wspd,spec(ip,jp) % significantWaveHeight(),&
    1./spec(ip,jp) % peakFrequency(),spec(ip,jp) % meanPeriod(),&
    spec(ip,jp) % meanPeriodZeroCrossing(),spec(ip,jp) % meanSquareSlope(),&
    spec(ip,jp) % frequencyMoment(0),spec(ip,jp) % frequencyMoment(1),&
    spec(ip,jp) % frequencyMoment(2)

  !domain = forward_euler(domain,&
  !  domain % advect(advectUpwind1stOrder1dRank1,1,&
  !  WAVY_OMNIDIRECTIONAL),dt=30._rk)

  !write(*,fmt='(20(f6.4,1x))')domain % significantWaveHeight()
  !write(*,fmt='(20(f7.4,1x))')domain % meanPeriod()

  domain = forward_euler(domain,&
    domain % advect(advectUpwind1stOrder2dRank2,[1,1]),dt=10._rk)


!call domain % writeJSON('domain.json',minify=.false.)
!
!
call nc_check(nf90_create('output.nc',nf90_clobber,ncid2))

call nc_check(nf90_def_dim(ncid2,'x',idm,xdimid))
call nc_check(nf90_def_dim(ncid2,'y',idm,ydimid))
stat = nf90_def_var(ncid2,'swh',nf90_float,[xdimid,ydimid],swhid)
call nc_check(nf90_enddef(ncid2))


stat = nf90_put_att(ncid2,swhid,name='description',values='Significant wave height')
stat = nf90_put_att(ncid2,swhid,name='units',values='m')

stat = nf90_put_var(ncid2,swhid,domain % significantWaveHeight(),start=[1,1],count=[idm,idm])
stat = nf90_close(ncid2)
enddo

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
