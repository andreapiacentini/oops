! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_obsdb_mod

use config_mod
use datetime_mod
use duration_mod
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use netcdf
use qg_constants_mod
use qg_locs_mod
use qg_obsvec_mod
use qg_projection_mod
use qg_tools_mod
use random_mod
use string_f_c_mod

implicit none

private
public :: qg_obsdb
public :: qg_obsdb_registry
public :: qg_obsdb_setup,qg_obsdb_delete,qg_obsdb_get,qg_obsdb_put,qg_obsdb_has,qg_obsdb_locations,qg_obsdb_generate,qg_obsdb_nobs
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 1 !< Random seed (for reproducibility)

type column_data
  character(len=50) :: colname                !< Column name
  type(column_data),pointer :: next => null() !< Next column
  integer :: nlev                             !< Number of levels
  real(kind_real),allocatable :: values(:,:)  !< Values
end type column_data

type group_data
  character(len=50) :: grpname                   !< Group name
  type(group_data),pointer :: next => null()     !< Next group
  integer :: nobs                                !< Number of observations
  type(datetime),allocatable :: times(:)         !< Time-slots
  type(column_data),pointer :: colhead => null() !< Head column
end type group_data

type qg_obsdb
  integer :: ngrp                               !< Number of groups
  character(len=1024) :: filein                 !< Input filename
  character(len=1024) :: fileout                !< Output filename
  type(group_data),pointer :: grphead => null() !< Head group
end type qg_obsdb

#define LISTED_TYPE qg_obsdb

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_obsdb_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup observation data
subroutine qg_obsdb_setup(self,conf)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data
type(c_ptr),intent(in) :: conf       !< Configuration

! Local variables
character(len=1024) :: fin,fout

! Input file
if (config_element_exists(conf,'ObsDataIn')) then
  fin = config_get_string(conf,1024,'ObsDataIn.obsfile')
else
  fin = ''
endif
call fckit_log%info('qg_obsdb_setup: file in = '//trim(fin))

! Output file
fout = config_get_string(conf,1024,'ObsDataOut.obsfile')
call fckit_log%info('qg_obsdb_setup: file out = '//trim(fout))

! Set attributes
self%ngrp = 0
self%filein = fin
self%fileout = fout

! Read observation data
if (self%filein/='') call qg_obsdb_read(self)

end subroutine qg_obsdb_setup
! ------------------------------------------------------------------------------
!> Delete observation data
subroutine qg_obsdb_delete(self)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs

! Write observation data
if (self%fileout/='') call qg_obsdb_write(self)

! Release memory
do while (associated(self%grphead))
  jgrp => self%grphead
  self%grphead => jgrp%next
  do jobs=1,jgrp%nobs
    call datetime_delete(jgrp%times(jobs))
  enddo
  deallocate(jgrp%times)
  do while (associated(jgrp%colhead))
    jcol => jgrp%colhead
    jgrp%colhead => jcol%next
    deallocate(jcol%values)
    deallocate(jcol)
  enddo
  deallocate(jgrp)
enddo

end subroutine qg_obsdb_delete
! ------------------------------------------------------------------------------
!> Get observation data
subroutine qg_obsdb_get(self,grp,col,ovec)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self     !< Observation data
character(len=*),intent(in) :: grp    !< Group
character(len=*),intent(in) :: col    !< Column
type(qg_obsvec),intent(inout) :: ovec !< Observation vector

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs,jlev

! Print Group and column
call fckit_log%info('qg_obsdb_get: grp = '//trim(grp))
call fckit_log%info('qg_obsdb_get: col = '//trim(col))

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) then
  jgrp => self%grphead
  do while (associated(jgrp))
    call fckit_log%info('qg_obsdb_get: group '//trim(jgrp%grpname)//' exists')
    jgrp => jgrp%next
  enddo
  call fckit_log%error('qg_obsdb_get: cannot find '//trim(grp))
  call abor1_ftn('qg_obsdb_get: obs group not found')
endif

! Find observation column
call qg_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) call abor1_ftn('qg_obsdb_get: obs column not found')

! Get observation data
if (allocated(ovec%values)) deallocate(ovec%values)
ovec%nobs = jgrp%nobs
ovec%nlev = jcol%nlev
allocate(ovec%values(ovec%nlev,ovec%nobs))
do jobs=1,jgrp%nobs
  do jlev=1,jcol%nlev
    ovec%values(jlev,jobs) = jcol%values(jlev,jobs)
  enddo
enddo

end subroutine qg_obsdb_get
! ------------------------------------------------------------------------------
!> Put observations data
subroutine qg_obsdb_put(self,grp,col,ovec)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data
character(len=*),intent(in) :: grp   !< Group
character(len=*),intent(in) :: col   !< Column
type(qg_obsvec),intent(in) :: ovec   !< Observation vector

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs,jlev

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) then
  jgrp => self%grphead
  do while (associated(jgrp))
    call fckit_log%info('qg_obsdb_put: group '//trim(jgrp%grpname)//' exists')
    jgrp => jgrp%next
  enddo
  call fckit_log%error('qg_obsdb_put: cannot find '//trim(grp))
  call abor1_ftn('qg_obsdb_put: obs group not found')
endif

! Find observation column (and add it if not there)
call qg_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) then
  if (.not.associated(jgrp%colhead)) call abor1_ftn('qg_obsdb_put: no locations')
  jcol => jgrp%colhead
  do while (associated(jcol%next))
    jcol => jcol%next
  enddo
  allocate(jcol%next)
  jcol => jcol%next
  jcol%colname = col
  jcol%nlev = ovec%nlev
  allocate(jcol%values(jcol%nlev,jgrp%nobs))
endif

! Put observation data
if (ovec%nobs/=jgrp%nobs) call abor1_ftn('qg_obsdb_put: error obs number')
if (ovec%nlev/=jcol%nlev) call abor1_ftn('qg_obsdb_put: error col number')
do jobs=1,jgrp%nobs
  do jlev=1,jcol%nlev
    jcol%values(jlev,jobs) = ovec%values(jlev,jobs)
  enddo
enddo

end subroutine qg_obsdb_put
! ------------------------------------------------------------------------------
!> Test observation data existence
subroutine qg_obsdb_has(self,grp,col,has)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self  !< Observation data
character(len=*),intent(in) :: grp !< Group
character(len=*),intent(in) :: col !< Column
integer,intent(out) :: has         !< Test flag

! Passed variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol

! Initialization
has = 0

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (associated(jgrp)) then
  ! Find observation column
  call qg_obsdb_find_column(jgrp,col,jcol)
  if (associated(jcol)) has = 1
endif

end subroutine qg_obsdb_has
! ------------------------------------------------------------------------------
!> Get locations from observation data
subroutine qg_obsdb_locations(self,grp,t1,t2,locs) 

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self   !< Observation data
character(len=*),intent(in) :: grp  !< Group
type(datetime),intent(in) :: t1     !< Time 1
type(datetime),intent(in) :: t2     !< Time 2
type(qg_locs),intent(inout) :: locs !< Locations

! Local variables
integer :: nobs
integer,allocatable :: mobs(:)
character(len=8),parameter :: col = 'Location'
type(qg_obsvec) :: ovec

! Count observations
call qg_obsdb_count_time(self,grp,t1,t2,nobs)
allocate(mobs(nobs))
call qg_obsdb_count_indx(self,grp,t1,t2,mobs)
call qg_obsdb_time_get(self,grp,col,t1,t2,ovec)

! Get locations
call qg_locs_from_obsvec(locs,ovec,mobs)

deallocate(ovec%values)
deallocate(mobs)

end subroutine qg_obsdb_locations
! ------------------------------------------------------------------------------
!> Generate observation data
subroutine qg_obsdb_generate(self,grp,conf,bgn,step,ktimes,kobs)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data
character(len=*),intent(in) :: grp   !< Group
type(c_ptr),intent(in) :: conf       !< Configuration
type(datetime),intent(in) :: bgn     !< Start time
type(duration),intent(in) :: step    !< Time-step
integer,intent(in) :: ktimes         !< Number of time-slots
integer,intent(inout) :: kobs        !< Number of observations

! Local variables
integer :: nlev,nlocs
real(kind_real) :: err
type(datetime),allocatable :: times(:)
type(qg_obsvec) :: obsloc,obserr

! Get number of observations
nlocs = config_get_int(conf,'obs_density');
kobs = nlocs*ktimes

! Allocation
allocate(times(kobs))

! Generate locations
call qg_obsdb_generate_locations(nlocs,ktimes,bgn,step,times,obsloc)

! Create observations data
call qg_obsdb_create(self,trim(grp),times,obsloc)

! Create observation error
err = config_get_real(conf,'obs_error')
nlev = config_get_int(conf,'nval')
call qg_obsvec_setup(obserr,nlev,kobs)
obserr%values(:,:) = err
call qg_obsdb_put(self,trim(grp),'ObsError',obserr)

! Release memory
deallocate(times)
deallocate(obsloc%values)
deallocate(obserr%values)

end subroutine qg_obsdb_generate
! ------------------------------------------------------------------------------
!> Get observation data size
subroutine qg_obsdb_nobs(self,grp,kobs)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self  !< Observation data
character(len=*),intent(in) :: grp !< Group
integer,intent(inout) :: kobs      !< Number of observations

! Local variables
type(group_data),pointer :: jgrp

! Find group
call qg_obsdb_find_group(self,grp,jgrp)

! Get observation data size
if (associated(jgrp)) then
  kobs = jgrp%nobs
else
  kobs = 0
endif

end subroutine qg_obsdb_nobs
! ------------------------------------------------------------------------------
!  Private
! ------------------------------------------------------------------------------
!> Read observation data
subroutine qg_obsdb_read(self)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data

! Local variables
integer :: igrp,icol,iobs,ncol
integer :: ncid,grpname_id,ngrp_id,nobs_id,ncol_id,times_id,nlev_id,colname_id,values_id
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
character(len=6) :: igrpchar
character(len=50) :: stime
character(len=1024) :: record

! Open NetCDF file
call ncerr(nf90_open(trim(self%filein)//'.nc',nf90_nowrite,ncid))

! Get dimensions ids
call ncerr(nf90_inq_dimid(ncid,'ngrp',ngrp_id))

! Get dimensions
call ncerr(nf90_inquire_dimension(ncid,ngrp_id,len=self%ngrp))

! Get variables ids
call ncerr(nf90_inq_varid(ncid,'grpname',grpname_id))

do igrp=1,self%ngrp 
  ! Allocation
  if (igrp==1) then
    allocate(self%grphead)
    jgrp => self%grphead
  else
    allocate(jgrp%next)
    jgrp => jgrp%next
  endif
  write(igrpchar,'(i6.6)') igrp

  ! Get variables
  call ncerr(nf90_get_var(ncid,grpname_id,jgrp%grpname,(/1,igrp/),(/50,1/)))

  ! Get dimensions ids
  call ncerr(nf90_inq_dimid(ncid,'nobs_'//igrpchar,nobs_id))
  call ncerr(nf90_inq_dimid(ncid,'ncol_'//igrpchar,ncol_id))

  ! Get dimensions
  call ncerr(nf90_inquire_dimension(ncid,nobs_id,len=jgrp%nobs))
  call ncerr(nf90_inquire_dimension(ncid,ncol_id,len=ncol))
  write(record,*) 'qg_obsdb_read: reading ',jgrp%nobs,' ',jgrp%grpname,' observations'
  call fckit_log%info(record)

  ! Get variables ids
  call ncerr(nf90_inq_varid(ncid,'times_'//igrpchar,times_id))
  call ncerr(nf90_inq_varid(ncid,'nlev_'//igrpchar,nlev_id))
  call ncerr(nf90_inq_varid(ncid,'colname_'//igrpchar,colname_id))
  call ncerr(nf90_inq_varid(ncid,'values_'//igrpchar,values_id))

  ! Allocation
  allocate(jgrp%times(jgrp%nobs))

  ! Get variables
  call ncerr(nf90_get_var(ncid,grpname_id,jgrp%grpname,(/1,igrp/),(/50,1/)))
  do iobs=1,jgrp%nobs
    call ncerr(nf90_get_var(ncid,times_id,stime,(/1,iobs/),(/50,1/)))
    call datetime_create(stime,jgrp%times(iobs))
  end do

  ! Loop over columns
  do icol=1,ncol
    ! Allocation
    if (icol==1) then
      allocate(jgrp%colhead)
      jcol => jgrp%colhead
    else
      allocate(jcol%next)
      jcol => jcol%next
    endif

    ! Get variables
    call ncerr(nf90_get_var(ncid,nlev_id,jcol%nlev,(/icol/)))
    call ncerr(nf90_get_var(ncid,colname_id,jcol%colname,(/1,icol/),(/50,1/)))

    ! Allocation
    allocate(jcol%values(jcol%nlev,jgrp%nobs))

    ! Get variables
    call ncerr(nf90_get_var(ncid,values_id,jcol%values(1:jcol%nlev,:),(/1,icol,1/),(/jcol%nlev,1,jgrp%nobs/)))
  enddo
enddo

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_obsdb_read
! ------------------------------------------------------------------------------
!> Write observation data
subroutine qg_obsdb_write(self)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self !< Observation data

! Local variables
integer :: igrp,icol,iobs,ncol,nlevmax
integer :: ncid,nstrmax_id,grpname_id,ngrp_id,nobs_id,ncol_id,nlevmax_id,times_id,nlev_id,colname_id,values_id
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
character(len=6) :: igrpchar
character(len=50) :: stime

! Create NetCDF file
call ncerr(nf90_create(trim(self%fileout)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nstrmax',50,nstrmax_id))
call ncerr(nf90_def_dim(ncid,'ngrp',self%ngrp,ngrp_id))

! Define variable
call ncerr(nf90_def_var(ncid,'grpname',nf90_char,(/nstrmax_id,ngrp_id/),grpname_id))

! End definitions
call ncerr(nf90_enddef(ncid))

! Loop over groups
igrp = 0
jgrp => self%grphead
do while (associated(jgrp))
  igrp = igrp+1
  write(igrpchar,'(i6.6)') igrp

  ! Enter definitions mode
  call ncerr(nf90_redef(ncid))

  ! Compute dimensions
  ncol = 0
  nlevmax = 0
  jcol => jgrp%colhead
  do while (associated(jcol))
    ncol = ncol+1
    nlevmax = max(jcol%nlev,nlevmax)
    jcol => jcol%next
  enddo

  ! Define dimensions
  call ncerr(nf90_def_dim(ncid,'nobs_'//igrpchar,jgrp%nobs,nobs_id))
  call ncerr(nf90_def_dim(ncid,'ncol_'//igrpchar,ncol,ncol_id))
  call ncerr(nf90_def_dim(ncid,'nlevmax_'//igrpchar,nlevmax,nlevmax_id))

  ! Define variable
  call ncerr(nf90_def_var(ncid,'times_'//igrpchar,nf90_char,(/nstrmax_id,nobs_id/),times_id))
  call ncerr(nf90_def_var(ncid,'nlev_'//igrpchar,nf90_int,(/ncol_id/),nlev_id))
  call ncerr(nf90_def_var(ncid,'colname_'//igrpchar,nf90_char,(/nstrmax_id,ncol_id/),colname_id))
  call ncerr(nf90_def_var(ncid,'values_'//igrpchar,nf90_double,(/nlevmax_id,ncol_id,nobs_id/),values_id))

  ! End definitions
  call ncerr(nf90_enddef(ncid))

  ! Put variables
  call ncerr(nf90_put_var(ncid,grpname_id,jgrp%grpname,(/1,igrp/),(/50,1/)))
  do iobs=1,jgrp%nobs
    call datetime_to_string(jgrp%times(iobs),stime)
    call ncerr(nf90_put_var(ncid,times_id,stime,(/1,iobs/),(/50,1/)))
  end do

  ! Loop over columns
  icol = 0
  jcol => jgrp%colhead
  do while (associated(jcol))
    icol = icol+1

    ! Put variables
    call ncerr(nf90_put_var(ncid,nlev_id,jcol%nlev,(/icol/)))
    call ncerr(nf90_put_var(ncid,colname_id,jcol%colname,(/1,icol/),(/50,1/)))
    call ncerr(nf90_put_var(ncid,values_id,jcol%values(1:jcol%nlev,:),(/1,icol,1/),(/jcol%nlev,1,jgrp%nobs/)))

    ! Update
    jcol => jcol%next
  enddo

  ! Update
  jgrp=>jgrp%next
end do

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_obsdb_write
! ------------------------------------------------------------------------------
!> Find observation data group
subroutine qg_obsdb_find_group(self,grp,find)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self              !< Observation data
character(len=*),intent(in) :: grp             !< Group
type(group_data),pointer,intent(inout) :: find !< Result

! Initialization
find => self%grphead

! Loop
do while (associated(find))
  if (find%grpname==grp) exit
  find => find%next
enddo

end subroutine qg_obsdb_find_group
! ------------------------------------------------------------------------------
!> Find observation data column
subroutine qg_obsdb_find_column(grp,col,find)

implicit none

! Passed variables
type(group_data),intent(in) :: grp              !< Observation data
character(len=*),intent(in) :: col              !< Column
type(column_data),pointer,intent(inout) :: find !< Result

! Initialization
find=>grp%colhead

! Loop
do while (associated(find))
  if (find%colname==col) exit
  find => find%next
enddo

end subroutine qg_obsdb_find_column
! ------------------------------------------------------------------------------
!> Generate random locations
subroutine qg_obsdb_generate_locations(nlocs,ntimes,bgn,step,times,obsloc)

implicit none

! Passed variables
integer,intent(in) :: nlocs                         !< Number of locations
integer,intent(in) :: ntimes                        !< Number of time-slots
type(datetime),intent(in) :: bgn                    !< Start time
type(duration),intent(in) :: step                   !< Time-step
type(datetime),intent(inout) :: times(nlocs*ntimes) !< Time-slots
type(qg_obsvec),intent(inout) :: obsloc             !< Observation locations

! Local variables
integer :: jobs,iobs,jstep
real(kind_real) :: x(nlocs),y(nlocs),z(nlocs),lon(nlocs),lat(nlocs)
type(datetime) :: now

! Backward variables - BACKWARD
integer :: ijk, intinv, ii, jj, kk, ij ! - BACKWARD
real(kind_real) :: goldinv ! - BACKWARD

! Generate random locations
call uniform_distribution(x,0.0_kind_real,domain_zonal,rseed)
call uniform_distribution(y,0.0_kind_real,domain_meridional,rseed)
call uniform_distribution(z,0.0_kind_real,domain_depth,rseed)

goldinv = 0.5_kind_real*(sqrt(5.0_kind_real)-1.0_kind_real) ! 1/(golden ratio) ! - BACKWARD
intinv = nint(goldinv*real(40*20*2,kind_real)) ! - BACKWARD
ijk=0 ! - BACKWARD
do jobs=1,nlocs ! - BACKWARD
  ijk=ijk+intinv ! - BACKWARD
  ijk=modulo(ijk,2*40*(20-2)) ! - BACKWARD
  kk=ijk/(40*(20-2)) ! - BACKWARD
  ij=ijk-kk*40*(20-2) ! - BACKWARD
  jj=ij/40 ! - BACKWARD
  ii=ij-jj*40 ! - BACKWARD
  jj=jj+2 ! - BACKWARD
  x(jobs)=(real(ii,kind_real)+0.5)*300000.0 ! - BACKWARD
  y(jobs)=real(jj,kind_real)*300000.0 ! - BACKWARD
  if (kk==0) then ! - BACKWARD
    z(jobs) = 4500.0+0.5*5500.0 ! - BACKWARD
  elseif (kk==1) then ! - BACKWARD
    z(jobs) = 0.5*4500.0 ! - BACKWARD
  end if ! - BACKWARD
enddo ! - BACKWARD

! Convert to lon/lat
do jobs=1,nlocs
  call xy_to_lonlat(x(jobs),y(jobs),lon(jobs),lat(jobs))
enddo

! Setup observation vector
call qg_obsvec_setup(obsloc,3,nlocs*ntimes)

! Set observation locations
now = bgn
iobs=0
do jstep=1,ntimes
  do jobs=1,nlocs
    iobs = iobs+1
    times(iobs) = now
    obsloc%values(:,iobs) = (/lon(jobs),lat(jobs),z(jobs)/)
  enddo
  call datetime_update(now,step)
enddo

! Release memory
call datetime_delete(now)

end subroutine qg_obsdb_generate_locations
! ------------------------------------------------------------------------------
!> Create observation data
subroutine qg_obsdb_create(self,grp,times,locs)

implicit none

! Passed varaibles
type(qg_obsdb),intent(inout) :: self  !< Observation data
character(len=*),intent(in) :: grp    !< Group
type(datetime),intent(in) :: times(:) !< Time-slots
type(qg_obsvec),intent(in) :: locs    !< Locations

! Local variables
type(group_data),pointer :: igrp
integer :: jobs,jlev

! Find observation group
call qg_obsdb_find_group(self,grp,igrp)
if (associated(igrp)) call abor1_ftn('qg_obsdb_create: obs group already exists')
if (associated(self%grphead)) then
  igrp => self%grphead
  do while (associated(igrp%next))
    igrp => igrp%next
  enddo
  allocate(igrp%next)
  igrp => igrp%next
else
  allocate(self%grphead)
  igrp => self%grphead
endif

! Create observation data
igrp%grpname = grp
igrp%nobs = size(times)
allocate(igrp%times(igrp%nobs))
igrp%times(:) = times(:)
allocate(igrp%colhead)
igrp%colhead%colname = 'Location'
igrp%colhead%nlev = 3
allocate(igrp%colhead%values(3,igrp%nobs))
if (locs%nlev/=3) call abor1_ftn('qg_obsdb_create: error locations not 3D')
if (locs%nobs/=igrp%nobs) call abor1_ftn('qg_obsdb_create: error locations number')
do jobs=1,igrp%nobs
  do jlev=1,3
    igrp%colhead%values(jlev,jobs) = locs%values(jlev,jobs)
  enddo
enddo
self%ngrp = self%ngrp+1

end subroutine qg_obsdb_create
! ------------------------------------------------------------------------------
!> Get observation data time
subroutine qg_obsdb_time_get(self,grp,col,t1,t2,ovec)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self     !< Observation data
character(len=*),intent(in) :: grp    !< Group
character(len=*),intent(in) :: col    !< Column
type(datetime),intent(in) :: t1       !< Time 1
type(datetime),intent(in) :: t2       !< Time 2
type(qg_obsvec),intent(inout) :: ovec !< Observation vector

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs,jlev,iobs

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) call abor1_ftn('qg_obsdb_time_get: obs group not found')

! Find observation column
call qg_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) call abor1_ftn('qg_obsdb_time_get: obs column not found')

! Time selection
iobs = 0
do jobs=1,jgrp%nobs
  if ((t1<jgrp%times(jobs)).and.(jgrp%times(jobs)<=t2)) iobs = iobs+1
enddo

! Get data
if ((ovec%nobs/=iobs).or.(ovec%nlev/=jcol%nlev)) then
  if (allocated(ovec%values)) deallocate(ovec%values)
  ovec%nobs = iobs
  ovec%nlev = jcol%nlev
  allocate(ovec%values(ovec%nlev,ovec%nobs))
endif
iobs = 0
do jobs=1,jgrp%nobs
  if ((t1<jgrp%times(jobs)).and.(jgrp%times(jobs)<=t2)) then
    iobs = iobs+1
    do jlev=1,jcol%nlev
      ovec%values(jlev,iobs) = jcol%values(jlev,jobs)
    enddo
  endif
enddo

end subroutine qg_obsdb_time_get
! ------------------------------------------------------------------------------
!> Get observation data size between times
subroutine qg_obsdb_count_time(self,grp,t1,t2,kobs)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self  !< Observation data
character(len=*),intent(in) :: grp !< Group
type(datetime),intent(in) :: t1    !< Time 1
type(datetime),intent(in) :: t2    !< Time 2
integer,intent(inout) :: kobs      !< Number of observations

! Local variables
type(group_data),pointer :: jgrp
integer :: jobs

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) call abor1_ftn('qg_obsdb_count_time: obs group not found')

! Time selection
kobs = 0
do jobs=1,jgrp%nobs
  if ((t1<jgrp%times(jobs)).and.(jgrp%times(jobs)<=t2)) kobs = kobs+1
enddo

end subroutine qg_obsdb_count_time
! ------------------------------------------------------------------------------
!> Get observation data index
subroutine qg_obsdb_count_indx(self,grp,t1,t2,indx)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self  !< Observation data
character(len=*),intent(in) :: grp !< Group
type(datetime),intent(in) :: t1    !< Time 1
type(datetime),intent(in) :: t2    !< Time 2
integer,intent(inout) :: indx(:)   !< Observation index

! Local variables
type(group_data),pointer :: jgrp
integer :: jobs,io

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) call abor1_ftn('qg_obsdb_count_indx: obs group not found')

! Time selection
io = 0
do jobs=1,jgrp%nobs
  if ((t1<jgrp%times(jobs)).and.(jgrp%times(jobs)<=t2)) then
    io = io+1
    indx(io)=jobs
  endif
enddo

end subroutine qg_obsdb_count_indx
! ------------------------------------------------------------------------------
end module qg_obsdb_mod