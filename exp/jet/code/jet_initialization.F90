!> A template of a user to code up customized initial conditions.
module user_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N
use MOM_open_boundary, only : OBC_DIRECTION_S
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
implicit none ; private

#include <MOM_memory.h>


character(len=40) :: mdl = "jet_initialization" !< This module's name.

! The following routines are visible to the outside world
public jet_initialize_topography
public jet_set_OBC_data
public register_jet_OBC, jet_OBC_end

!> Control structure for jet open boundaries.
type, public :: jet_OBC_CS ; private
  real :: delta = 0.3       !< asymmetry parameter
  real :: Ro = 0.6          !< Rossby nubmer
  real :: F = 1.0989        !< Inverse Burger number
  real :: L = 31.0          !< Length scale [L ~> km]
  real :: Ls, Ln            !< Jet decay scales to the south and north [L -> km]
  real :: U0                !< Peak velocity of jet
  real :: F_0               !< Coriolis parameter [T-1 ~> s-1]
end type jet_OBC_CS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Add jet to OBC registry.
function register_jet_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(jet_OBC_CS),         pointer    :: CS         !< jet control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_jet_OBC
  ! Local variables
  real    :: km_to_L_scale  ! A scaling factor from longitudes in km to L [L km-1 ~> 1e3]

  character(len=32)  :: casename = "jet"       !< This case's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "register_jet_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Register the tracer for horizontal advection & diffusion.
  call register_OBC(casename, param_file, OBC_Reg)
  call get_param(param_file, mdl, "F_0", CS%F_0, &
                 default=0.0, units="s-1", scale=US%T_to_s, do_not_log=.true.)
  call get_param(param_file, mdl,"JET_L", CS%L, &
                 "Jet length scale.",&
                 units="km", default=31.0, scale=1.0e3*US%m_to_L)
  call get_param(param_file, mdl, "JET_DELTA", CS%delta, &
                 "Jet asymmetry parameter.",  &
                 units="nondim", default=0.3)
  call get_param(param_file, mdl, "JET_ROSSBY", CS%Ro, &
                 "Jet Rossby number.",  &
                 units="nondim", default=0.6)
  call get_param(param_file, mdl, "JET_F", CS%F, &
                 "Jet inverse Burger number.",  &
                 units="nondim", default=1.0989)
                 
  km_to_L_scale = 1000.0*US%m_to_L
  ! Unravel the dimensionless parameters into dimensional parameters
  CS%Ls = CS%L*(1.0 + CS%delta)
  CS%Ln = CS%L*(1.0 - CS%delta)
  CS%U0 = CS%Ro*CS%F_0*CS%L*km_to_L_scale
                 
  register_jet_OBC = .true.

end function register_jet_OBC

!> Clean up the jet OBC from registry.
subroutine jet_OBC_end(CS)
  type(jet_OBC_CS), pointer    :: CS         !< jet control structure.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine jet_OBC_end

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine jet_set_OBC_data(OBC, CS, G, GV, US, PF)
  type(ocean_OBC_type),       pointer    :: OBC   !< This open boundary condition type specifies
                                                  !! whether, where, and what open boundary
                                                  !! conditions are used.
  type(jet_OBC_CS),           pointer    :: CS    !< tidal bay control structure.
  type(ocean_grid_type),      intent(in) :: G     !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV    !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in) :: US    !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

  ! local variables to set up the inflow
  integer :: i, j, k, n, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  integer :: n
  type(OBC_segment_type), pointer :: segment => NULL()
  real    :: km_to_L_scale  ! A scaling factor from longitudes in km to L [L km-1 ~> 1e3]
  real    :: delta                                  !! asymmetry parameter
  real    :: F                                      !! Inverse Burger number
  real    :: Ro                                     !! Rossby number
  real    :: L                                      !! Length scale [L ~> km]
  real    :: Ls, Ln                                 !! Jet decay scales to the south and north [L -> km]
  real    :: U0                                     !! Peak velocity of jet

  if (first_call) call write_user_log(param_file)
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) call MOM_error(FATAL, 'jet_initialization.F90: '// &
        'jet_initialization() was called but OBC type was not initialized!')
  
  if (OBC%number_of_segments /= 2) then
    call MOM_error(WARNING, 'Error in jet OBC segment setup', .true.)
    return   !!! Need a better error message here
  endif

  
  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    ! Don't bother to do anything unless there are OB points in this tile
    if (.not. segment%on_pe) cycle
    
    ! Apply values on inflow end only
    if (segment%direction /= OBC_DIRECTION_W) cycle
  enddo
  
end subroutine USER_set_OBC_data

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

end module user_initialization
