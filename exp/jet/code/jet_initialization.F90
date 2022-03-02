!> A template of a user to code up customized initial conditions.
module jet_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N
use MOM_open_boundary, only : OBC_DIRECTION_S
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>


character(len=40) :: mdl = "jet_initialization" !< This module's name.

! The following routines are visible to the outside world
public jet_set_OBC_data

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A module variable that should not be used.
!! \todo Move this module variable into a control structure.
logical :: first_call = .true.

contains

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine jet_set_OBC_data(OBC, G, GV, US, param_file)
  type(ocean_OBC_type),       pointer    :: OBC   !< This open boundary condition type specifies
                                                  !! whether, where, and what open boundary
                                                  !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G     !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV    !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in) :: US    !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

  ! local variables to set up the inflow
  integer :: i, j, k, n, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment => NULL()
  real    :: km_to_L_scale  ! A scaling factor from longitudes in km to L [L km-1 ~> 1e3]
  real    :: f_0                                    !! Coriolis parameter [T-1 ~> s-1]
  real    :: delta                                  !! asymmetry parameter
  real    :: F                                      !! Inverse Burger number
  real    :: Ro                                     !! Rossby number
  real    :: L                                      !! Length scale [L ~> km]
  real    :: Ls, Ln                                 !! Jet decay scales to the south and north [L -> km]
  real    :: U0                                     !! Peak velocity of jet
  real    :: y

  if (first_call) call write_user_log(param_file)
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) call MOM_error(FATAL, 'jet_initialization.F90: '// &
        'jet_initialization() was called but OBC type was not initialized!')
        
  call get_param(param_file, "MOM_shared_initialization", "F_0", f_0, &
                 "The reference value of the Coriolis parameter with the "//&
                 "betaplane option.", units="s-1", default=0.0, scale=US%T_to_s)
  call get_param(param_file, mdl,"JET_L", L, &
                 "Jet length scale.",&
                 units="km", default=31.0, scale=1.0e3*US%m_to_L)
  call get_param(param_file, mdl, "JET_DELTA", delta, &
                 "Jet asymmetry parameter.",  &
                 units="nondim", default=0.3)
  call get_param(param_file, mdl, "JET_ROSSBY", Ro, &
                 "Jet Rossby number.",  &
                 units="nondim", default=0.6)
  call get_param(param_file, mdl, "JET_F", F, &
                 "Jet inverse Burger number.",  &
                 units="nondim", default=1.0989)
                 
  km_to_L_scale = 1000.0*US%m_to_L
  ! Unravel the dimensionless parameters into dimensional parameters
  Ls = L*(1.0 + delta)
  Ln = L*(1.0 - delta)
  U0 = Ro*f_0*L*km_to_L_scale
  
  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    ! Don't bother to do anything unless there are OB points in this tile
    if (.not. segment%on_pe) cycle
    
    ! Apply values on inflow end only
    if (segment%direction /= OBC_DIRECTION_W) cycle
    
    do k=1,nz
      IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
      jsd = segment%HI%jsd ; jed = segment%HI%jed
      do j=jsd,jed ; do I=IsdB,IedB
        y = G%geoLatCu(I,j)
      
        if (k == 1) then
          segment%normal_vel(I,j,k) = U0*exp(-abs(y)/(L*(1.0 - sign(delta, y))))
        endif
      enddo ; enddo
    enddo
  enddo
  
end subroutine jet_set_OBC_data

!> Write output about the parameter values being used.
subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "jet_initialization" ! This module's name.

  call log_version(param_file, mdl, version)
  first_call = .false.

end subroutine write_user_log

!> \namespace user_initialization
!!
!!  This subroutine initializes the fields for the simulations.
!!  The one argument passed to initialize, Time, is set to the
!!  current time of the simulation.  The fields which might be initialized
!!  here are:
!!  - u - Zonal velocity [Z T-1 ~> m s-1].
!!  - v - Meridional velocity [Z T-1 ~> m s-1].
!!  - h - Layer thickness [H ~> m or kg m-2].  (Must be positive.)
!!  - G%bathyT - Basin depth [Z ~> m].
!!  - G%CoriolisBu - The Coriolis parameter [T-1 ~> s-1].
!!  - GV%g_prime - The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
!!  - GV%Rlay - Layer potential density (coordinate variable) [R ~> kg m-3].
!!  If ENABLE_THERMODYNAMICS is defined:
!!  - T - Temperature [degC].
!!  - S - Salinity [ppt].
!!  If BULKMIXEDLAYER is defined:
!!  - Rml - Mixed layer and buffer layer potential densities [R ~> kg m-3].
!!  If SPONGE is defined:
!!  - A series of subroutine calls are made to set up the damping
!!    rates and reference profiles for all variables that are damped
!!    in the sponge.
!!
!!  Any user provided tracer code is also first linked through this
!!  subroutine.
!!
!!  These variables are all set in the set of subroutines (in this
!!  file) USER_initialize_bottom_depth, USER_initialize_thickness,
!!  USER_initialize_velocity,  USER_initialize_temperature_salinity,
!!  USER_initialize_mixed_layer_density, USER_initialize_sponges,
!!  USER_set_coord, and USER_set_ref_profile.
!!
!!  The names of these subroutines should be self-explanatory. They
!!  start with "USER_" to indicate that they will likely have to be
!!  modified for each simulation to set the initial conditions and
!!  boundary conditions.  Most of these take two arguments: an integer
!!  argument specifying whether the fields are to be calculated
!!  internally or read from a NetCDF file; and a string giving the
!!  path to that file.  If the field is initialized internally, the
!!  path is ignored.

end module jet_initialization
