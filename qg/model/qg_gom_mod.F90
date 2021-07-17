! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_gom_mod

use atlas_module, only: atlas_field
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use netcdf
use oops_variables_mod
use qg_constants_mod
use qg_geom_mod
use qg_locs_mod
! use qg_projection_mod
! use qg_tools_mod
use random_mod

implicit none
private
public :: qg_gom
public :: qg_gom_registry
public :: qg_gom_setup,qg_gom_delete,qg_gom_copy,qg_gom_zero,qg_gom_abs,qg_gom_random,qg_gom_mult, &
        & qg_gom_add,qg_gom_diff,qg_gom_schurmult,qg_gom_divide,qg_gom_rms,qg_gom_dotprod,qg_gom_stats,qg_gom_maxloc, &
        & qg_gom_read_file, qg_gom_write_file,qg_gom_analytic_init
! ------------------------------------------------------------------------------
type :: qg_gom
  integer :: nobs                      !< Number of observations
  real(kind_real), allocatable :: x(:) !< Chemical observations values
  logical :: lalloc = .false.          !< Allocation flag
  type(oops_variables) :: vars         !< Variables
end type qg_gom

#define LISTED_TYPE qg_gom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_gom_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup GOM
subroutine qg_gom_setup(self,nobs)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self      !< GOM
integer, intent(in) :: nobs             !< Number of observations

! Set attributes
self%nobs = nobs

! Allocation
allocate(self%x(self%nobs))

self%lalloc = .true.

end subroutine qg_gom_setup
! ------------------------------------------------------------------------------
!> Delete GOM
subroutine qg_gom_delete(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Release memory
if (allocated(self%x)) deallocate(self%x)

self%lalloc = .false.

end subroutine qg_gom_delete
! ------------------------------------------------------------------------------
!> Copy GOM
subroutine qg_gom_copy(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self            !< GOM
type(qg_gom),intent(in) :: other              !< Other GOM

! Copy attributes
self%nobs = other%nobs

! Allocation
if (.not.self%lalloc) then
  allocate(self%x(self%nobs))
  self%lalloc = .true.
endif

! Copy
self%x = other%x

end subroutine qg_gom_copy
! ------------------------------------------------------------------------------
!> Set GOM to zero
subroutine qg_gom_zero(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Set to zero
self%x = 0.0

end subroutine qg_gom_zero
! ------------------------------------------------------------------------------
!> Get GOM absolute value
subroutine qg_gom_abs(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Get absolute value
self%x = abs(self%x)

end subroutine qg_gom_abs
! ------------------------------------------------------------------------------
!> Generate random GOM values
subroutine qg_gom_random(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Local variables
integer :: nv
real(kind_real),allocatable :: values(:,:)

! TODO(Benjamin): change that in a following PR
nv = 0
nv = nv+1
allocate(values(nv,self%nobs))

! Generate random GOM values
call normal_distribution(values,0.0_kind_real,1.0_kind_real,seed=14871,reset=.true.)

! Split random values
nv = 0
nv = nv+1
self%x = values(nv,:)

end subroutine qg_gom_random
! ------------------------------------------------------------------------------
!> Multiply GOM with a scalar
subroutine qg_gom_mult(self,zz)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
real(kind_real),intent(in) :: zz   !< Multiplier

! Multiply GOM with a scalar
self%x = zz*self%x

end subroutine qg_gom_mult
! ------------------------------------------------------------------------------
!> Add GOM
subroutine qg_gom_add(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Add GOM
self%x = self%x+other%x

end subroutine qg_gom_add
! ------------------------------------------------------------------------------
!> Subtract GOM
subroutine qg_gom_diff(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Subtract GOM
self%x = self%x-other%x

end subroutine qg_gom_diff
! ------------------------------------------------------------------------------
!> Schur product for GOM
subroutine qg_gom_schurmult(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Multiply GOM
self%x = self%x*other%x

end subroutine qg_gom_schurmult
! ------------------------------------------------------------------------------
!> Schur division for GOM
subroutine qg_gom_divide(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Local variables
real(kind_real) :: tol
integer :: jloc

! Set tolerance
tol = epsilon(tol)

! Conditional division
do jloc=1,self%nobs
  if (abs(other%x(jloc))>tol) then
    self%x(jloc) = self%x(jloc)/other%x(jloc)
  else
    self%x(jloc) = 0.0
  endif
enddo

end subroutine qg_gom_divide
! ------------------------------------------------------------------------------
!> Compute GOM RMS
subroutine qg_gom_rms(self,rms)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self   !< GOM
real(kind_real),intent(inout) :: rms !< RMS

! Local variables
integer :: nv

! Initialization
rms = 0.0
nv = 0

! Loop over values
rms = rms+sum(self%x**2)
nv = nv+1

! Normalize and take square-root
rms = sqrt(rms/real(self%nobs*nv,kind_real))

end subroutine qg_gom_rms
! ------------------------------------------------------------------------------
!> GOM dot product
subroutine qg_gom_dotprod(gom1,gom2,prod)

implicit none

! Passed variables
type(qg_gom),intent(in) :: gom1       !< GOM 1
type(qg_gom),intent(in) :: gom2       !< GOM 2
real(kind_real),intent(inout) :: prod !< Dot product

! Local variables
integer :: jo,jv

! Check
if (gom1%nobs/=gom2%nobs) call abor1_ftn('qg_gom_dotprod: inconsistent GOM sizes')

! Initialization
prod = 0.0

! Dot product
prod = prod+sum(gom1%x*gom2%x)

end subroutine qg_gom_dotprod
! ------------------------------------------------------------------------------
!> Compute GOM stats
subroutine qg_gom_stats(self,kobs,pmin,pmax,prms)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self       !< GOM
integer,intent(inout) :: kobs            !< Number of observations
real(kind_real),intent(inout) :: pmin    !< Minimum value
real(kind_real),intent(inout) :: pmax    !< Maximum value
real(kind_real),intent(inout) :: prms    !< RMS

! Local variables
integer :: nv

! Compute GOM stats
kobs = self%nobs
if (self%nobs>0) then
  pmin = huge(1.0)
  pmax = -huge(1.0)
  prms = 0.0
  nv = 0
  pmin = min(pmin,minval(self%x))
  pmax = max(pmax,maxval(self%x))
  prms = prms+sum(self%x**2)
  nv = nv+1
  prms = sqrt(prms/real(self%nobs*nv,kind_real))
else
  pmin = 0.0
  pmax = 0.0
  prms = 0.0
end if

end subroutine qg_gom_stats
! ------------------------------------------------------------------------------
!> Find and locate GOM max. value
subroutine qg_gom_maxloc(self,mxval,mxloc,mxvar)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self          !< GOM
real(kind_real),intent(inout) :: mxval      !< Maximum value
integer,intent(inout) :: mxloc              !< Location of maximum value
type(oops_variables),intent(inout) :: mxvar !< Variable of maximum value

! Local variables
integer :: mxloc_arr(1),mxval_tmp
character(len=1) :: var

! Initialization
mxval = -huge(1.0)

! Find GOM max. value
mxval_tmp = maxval(self%x)
if (mxval_tmp>mxval) then
  mxval = mxval
  mxloc_arr = maxloc(self%x)
  var = 'x'
endif

! Locate GOM max. value
mxloc = mxloc_arr(1)

! Set GOM max. variable
call mxvar%push_back(var)

end subroutine qg_gom_maxloc
! ------------------------------------------------------------------------------
!> Read GOM from file
subroutine qg_gom_read_file(self,f_conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self             !< GOM
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ncid,nobs_id,nobs,x_id,q_id,u_id,v_id
character(len=1024) :: filename
character(len=:),allocatable :: str

! Get filename
call f_conf%get_or_die("filename",str)
filename = str
call fckit_log%info('qg_gom_read_file: reading '//trim(filename))

! Open NetCDF file
! call ncerr(nf90_open(trim(filename)//'.nc',nf90_nowrite,ncid))

! ! Get dimension id
! call ncerr(nf90_inq_dimid(ncid,'nobs',nobs_id))

! ! Get dimension
! call ncerr(nf90_inquire_dimension(ncid,nobs_id,len=nobs))

! ! GOM setup
! call qg_gom_setup(self,nobs)

! ! Get variables ids
! if (self%vars%has('x')) call ncerr(nf90_inq_varid(ncid,'x',x_id))
! if (self%vars%has('q')) call ncerr(nf90_inq_varid(ncid,'q',q_id))
! if (self%vars%has('u')) call ncerr(nf90_inq_varid(ncid,'u',u_id))
! if (self%vars%has('v')) call ncerr(nf90_inq_varid(ncid,'v',v_id))

! ! Get variables
! if (self%vars%has('x')) call ncerr(nf90_get_var(ncid,x_id,self%x))
! if (self%vars%has('q')) call ncerr(nf90_get_var(ncid,q_id,self%q))
! if (self%vars%has('u')) call ncerr(nf90_get_var(ncid,u_id,self%u))
! if (self%vars%has('v')) call ncerr(nf90_get_var(ncid,v_id,self%v))

! ! Close NetCDF file
! call ncerr(nf90_close(ncid))

end subroutine qg_gom_read_file
! ------------------------------------------------------------------------------
!> Write GOM to file
subroutine qg_gom_write_file(self,f_conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ncid,nobs_id,x_id,q_id,u_id,v_id
character(len=1024) :: filename
character(len=:),allocatable :: str

! Check allocation
if (.not.self%lalloc) call abor1_ftn('qg_gom_write_file: gom not allocated')

! Set filename
call f_conf%get_or_die("filename",str)
filename = str
call fckit_log%info('qg_gom_write_file: writing '//trim(filename))

! Create NetCDF file
! call ncerr(nf90_create(trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! ! Define dimensions
! call ncerr(nf90_def_dim(ncid,'nobs',self%nobs,nobs_id))

! ! Define variables
! if (self%vars%has('x')) call ncerr(nf90_def_var(ncid,'x',nf90_double,(/nobs_id/),x_id))
! if (self%vars%has('q')) call ncerr(nf90_def_var(ncid,'q',nf90_double,(/nobs_id/),q_id))
! if (self%vars%has('u')) call ncerr(nf90_def_var(ncid,'u',nf90_double,(/nobs_id/),u_id))
! if (self%vars%has('v')) call ncerr(nf90_def_var(ncid,'v',nf90_double,(/nobs_id/),v_id))

! ! End definitions
! call ncerr(nf90_enddef(ncid))

! ! Put variables
! if (self%vars%has('x')) call ncerr(nf90_put_var(ncid,x_id,self%x))
! if (self%vars%has('q')) call ncerr(nf90_put_var(ncid,q_id,self%q))
! if (self%vars%has('u')) call ncerr(nf90_put_var(ncid,u_id,self%u))
! if (self%vars%has('v')) call ncerr(nf90_put_var(ncid,v_id,self%v))

! ! Close NetCDF file
! call ncerr(nf90_close(ncid))

end subroutine qg_gom_write_file
! ------------------------------------------------------------------------------
!> GOM analytic initialization
subroutine qg_gom_analytic_init(self,locs,f_conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self             !< GOM
type(qg_locs),intent(inout) :: locs            !< Locations
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: iloc
real(kind_real) :: x,y
character(len=30) :: ic
character(len=:),allocatable :: str
real(kind_real), pointer :: lonlat(:,:), z(:)
type(atlas_field) :: lonlat_field, z_field

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

z_field = locs%altitude()
call z_field%data(z)

! Check allocation
if (.not.self%lalloc) call abor1_ftn('qg_gom_analytic init: gom not allocated')

! Get analytic configuration
call f_conf%get_or_die("analytic_init",str)
ic = str
call fckit_log%info('qg_gom_analytic_init: ic = '//trim(ic))
! do iloc=1,locs%nlocs()
!   select case (trim(ic))
!   case ('baroclinic-instability')
!     ! Go to cartesian coordinates
!     call lonlat_to_xy(lonlat(1,iloc),lonlat(2,iloc),x,y)

!     ! Compute values for baroclinic instability
!     if (self%vars%has('x')) call baroclinic_instability(x,y,z(iloc),'x',self%x(iloc))
!     if (self%vars%has('q')) call baroclinic_instability(x,y,z(iloc),'q',self%q(iloc))
!     if (self%vars%has('u')) call baroclinic_instability(x,y,z(iloc),'u',self%u(iloc))
!     if (self%vars%has('v')) call baroclinic_instability(x,y,z(iloc),'v',self%v(iloc))
!   case ('large-vortices')
!     ! Go to cartesian coordinates
!     call lonlat_to_xy(lonlat(1,iloc),lonlat(2,iloc),x,y)

!     ! Compute values for large vortices
!     if (self%vars%has('x')) call large_vortices(x,y,z(iloc),'x',self%x(iloc))
!     if (self%vars%has('q')) call large_vortices(x,y,z(iloc),'q',self%q(iloc))
!     if (self%vars%has('u')) call large_vortices(x,y,z(iloc),'u',self%u(iloc))
!     if (self%vars%has('v')) call large_vortices(x,y,z(iloc),'v',self%v(iloc))
!   case default
!     call abor1_ftn('qg_gom_analytic_init: unknown initialization')
!   endselect
! enddo

call lonlat_field%final()
call z_field%final()

end subroutine qg_gom_analytic_init
! ------------------------------------------------------------------------------
end module qg_gom_mod
