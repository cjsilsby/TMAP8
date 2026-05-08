# Verification Problem #1e from TMAP4/TMAP7 V&V document
# Permeation problem in a composite layer
# No Soret effect, trapping, or solubility included.

# Numerical parameters
num_Cu = 3000 # - 1000
num_W = 1500
simulation_time = '${units 0.2 s}' # previous 1
step = '${units 0.00001 s}' # previous 0.000001

# System properties
T_Cu = '${units 0.05 mm -> m}' # Thickness of Copper substrate. # T_PyC
T_W = '${units 0.025 mm -> m}' # Thickness of Tungsten film  # SiC
D_ver = '${units 0.02 mm -> m}' # Distance past T_Cu to calculate theoretical diffusion in tungsten layer # D_ver
D_ver_Cu = '${units 0.045 mm -> m}' # Distance past 0 to calculate theoretical diffusion on copper layer # D_ver_PyC
#length_Cu = ${T_Cu} # '${units 1 mm}'


# Initial Conditions
initial_concentration = '${units 0.1 mol/m^3}'
temperature_i = '${units 1500 K}'
initial_Cu = 0
initial_W = 0

# Constants
k_B = '${units 8.61733e-5 eV/K}'
ideal_R = '${units 8.31446261815324 J/mol/K}'

# Material Properties
Do_Cu = '${units 1.74e-6 m^2/s}' # Diffusivity of H in Cu (m^2/s)
E_D_Cu = '${units 42000 J/mol}'
# Diffusivity_Cu = '${units ${fparse Do_Cu * exp(-E_D_Cu/ideal_R/temperature_i)} m^2/s}'
Do_W = '${units 7.44e-8 m^2/s}' # Diffusivity of H in W (m^2/s)
E_D_W = '${units 0.13 eV}'
# Diffusivity_W = '${units ${fparse Do_W * exp(-E_D_W/k_B/temperature_i)} m^2/s}'


[Mesh]
  [total_geo]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${T_Cu} ${T_W}'
    ix = '${num_Cu} ${num_Cu}'
    subdomain_id = '1 2'
   # allow_renumbering = false
  []
  [Cu_inter]
    type = SideSetsBetweenSubdomainsGenerator
    input = total_geo
    new_boundary = 'inter_Cu_to_W'
    primary_block = '1' # Cu
    paired_block = '2' # W
  []
  [W_inter]
    type = SideSetsBetweenSubdomainsGenerator
    input = Cu_inter
    new_boundary = 'inter_W_to_Cu'
    primary_block = '2' # W
    paired_block = '1' # Cu
  []
[]

[Variables]
  [u_Cu]
    block = 1
    order = FIRST
    family = LAGRANGE
    initial_condition = '${initial_Cu}'
  []
  [u_W]
    block = 2
    order = FIRST
    family = LAGRANGE
    initial_condition = '${initial_W}'
  []
[]

[AuxVariables]
  [T] # Temperature
    initial_condition = '${temperature_i}'
  []
[]

# [Kernels]
#   [diff]
#     type = FunctionDiffusion
#     variable = u
#     function = diffusivity_value
#   []
#   [time]
#     type = TimeDerivative
#     variable = u
#   []
# []
[Kernels]
  # Transient diffusion of deuterium in BeO
  [time_Cu]
    type = TimeDerivative
    variable = u_Cu
    block = 1
  []
  [diffusion_Cu]
    type = ADMatDiffusion
    variable = u_Cu
    diffusivity = Diffusivity_Cu
    block = 1
  []
  # Transient diffusion of deuterium in Be
  [time_W]
    type = TimeDerivative
    variable = u_W
    block = 2
  []
  [diffusion_W]
    type = ADMatDiffusion
    variable = u_W
    diffusivity = Diffusivity_W
    block = 2
  []
[]

[InterfaceKernels]
  # Penalized continuity with solubility jump at BeO/Be interface
  [tied]
    type = ADPenaltyInterfaceDiffusion
    variable = u_Cu
    neighbor_var = u_W
    penalty = ${jump_penalty}
    jump_prop_name = solubility_ratio
    boundary = 'Cu_inter'
  []
[]

[AuxKernels]
  [constant_temperature]
    type = FunctionAux
    variable = T
    function = temp_func
    execute_on = 'INITIAL LINEAR'
  []
[]

[BCs]
  [left]
    type = EquilibriumBC
    Ko =
    enclosure_var =

    variable = u
    boundary = left
    value = ${initial_concentration}
  []
  [right]
    type = ADNeumannBC
    variable = u_W
    boundary = right
    value = 0
  []
[]

[Functions]
  [temp_func]
    type = ParsedFunction
    expression = '${temperature_i} + 0'
  []
  [Diffusivity_Cu]
    type = ParsedFunction
    symbol_names = 'T'
    symbol_values = 'temp_func'
    expression = '${Do_Cu} * exp(-${E_D_Cu} / ${ideal_R} / T)'
  []
  [Diffusivity_W]
    type = ParsedFunction
    symbol_names = 'T'
    symbol_values = 'temp_func'
    expression = '${Do_W} * exp(-${E_D_W} / ${k_B} / T)'
  []
  # Diffusivity assign based on different material domains
  [diffusivity_value]
    type = ParsedFunction
    expression = 'if(x < ${T_Cu}, Diffusivity_Cu, Diffusivity_W)'
  []
[]

# Used while obtaining steady-state solution
#
# [VectorPostprocessors]
#   [line]
#     type = LineValueSampler
#     start_point = '0 0 0'
#     end_point = '${Mesh/xmax} 0 0'
#     num_points = ${Mesh/nx}
#     sort_by = 'x'
#     variable = u
#     outputs = vector_postproc
#   []
# []

[Postprocessors]
  # Used to obtain varying concentration with time at a
  # point in W layer 'x' um from Interface Cu/W boundary
  # x = 8 um for TMAP4 verification case,
  # x = 15.75 um for TMAP7 verification case
  [concentration_at_x_W] #_x_SiC
    type = PointValue
    variable = u
    point = '${fparse ${T_Cu} + ${D_ver}} 0 0'
    outputs = 'csv'
  []
  [concentration_at_x_Cu] # _x_PyC
    type = PointValue
    variable = u
    point = '${D_ver_Cu} 0 0'
    outputs = 'csv'
  []
[]

[Executioner]
  type = Transient
  end_time = ${simulation_time}
  dtmax = 0.010
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  scheme = 'bdf2'
  nl_rel_tol = 1e-50 # Make this really tight so that our absolute tolerance criterion is the one
  # we must meet
  nl_abs_tol = 1e-12
  abort_on_solve_fail = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = '${step}'
    optimal_iterations = 4
    growth_factor = 1.1
    cutback_factor = 0.9
  []
[]

[Outputs]
  [exodus]
    type = Exodus
  []
  [csv]
    type = CSV
  []
  [vector_postproc]
    type = CSV
    sync_times = ${simulation_time}
    sync_only = true
  []
[]
