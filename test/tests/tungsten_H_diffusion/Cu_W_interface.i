# Verification Problem #1e from TMAP4/TMAP7 V&V document
# Permeation problem in a composite layer
# No Soret effect, trapping, or solubility included.

# Numerical parameters
num_Cu = 3000 # - 1000
num_W = 1500
simulation_time = '${units 0.2 s}' # previous 0.2
step = '${units 0.00001 s}' # previous 0.000001

# System properties
T_Cu = '${units 0.05 mm -> m}' # Thickness of Copper substrate.
T_W = '${units 0.025 mm -> m}' # Thickness of Tungsten film
D_ver = '${units 0.02 mm -> m}' # Distance past T_Cu to calculate theoretical diffusion in tungsten layer
D_ver_Cu = '${units 0.045 mm -> m}' # Distance past 0 to calculate theoretical diffusion on copper layer

# Initial Conditions
#initial_concentration = '${units 0.1 mol/m^3}'
temperature_i = '${units 1500 K}'
H_partialPressure = '${units 3.35466447 Pa}' # Calcualted using EquilibriumBC equation # 1247.1 Pa = 0.1 mol H2/m^3 by ideal gas law
initial_Cu = 0.0
initial_W = 0.0

# Constants
#k_B = '${units 8.6173303e-5 eV/K}'
ideal_R = '${units 8.31446261815324 J/mol/K}'

# Material Properties
Do_Cu = '${units 1.74e-6 m^2/s}' # Diffusivity of H in Cu (m^2/s)
E_D_Cu = '${units 42000 J/mol}'
So_Cu = '${units ${fparse 0.8125*2} mol/m^3/Pa^(1/2)}' # H2 -> H
E_S_Cu = '${units 42320 J/mol}'
# Diffusivity_Cu = '${units ${fparse Do_Cu * exp(-E_D_Cu/ideal_R/temperature_i)} m^2/s}'
Do_W = '${units 7.44e-8 m^2/s}' # Diffusivity of H in W (m^2/s)
E_D_W = '${units 1.2543097e+04 J/mol}' # = eV*R/k_B to cancel out use of R in BC, Kernals, and Funtions, 0.13 eV
So_W = '${units 3.2651 mol/m^3/Pa^(1/2)}' # H
E_S_W = '${units 1.2060671e+05 J/mol}' # = eV*R/k_B to cancel out use of R in BC, Kernals, and Funtions, 1.25 eV
# Diffusivity_W = '${units ${fparse Do_W * exp(-E_D_W/k_B/temperature_i)} m^2/s}'
jump_penalty = 1e0


[Mesh]
  [total_geo]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${T_Cu} ${T_W}'
    ix = '${num_Cu} ${num_W}'
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
    boundary = 'inter_Cu_to_W'
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
    variable = u_Cu
    boundary = left
    enclosure_var = ${H_partialPressure}
    Ko = ${So_Cu}
    activation_energy = ${E_S_Cu}
    temperature = T
    p = 0.5 # Sieverts law
  []
  [right]
    type = DirichletBC
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
  [Diffusivity_Cu_func]
    type = ParsedFunction
    symbol_names = 'T'
    symbol_values = 'temp_func'
    expression = '${Do_Cu} * exp(-${E_D_Cu} / ${ideal_R} / T)'
  []
  [Diffusivity_W_func]
    type = ParsedFunction
    symbol_names = 'T'
    symbol_values = 'temp_func'
    expression = '${Do_W} * exp(-${E_D_W} / ${ideal_R} / T)'
  []
  [Solubility_Cu_func]
    type = ParsedFunction
    symbol_names = 'T'
    symbol_values = 'temp_func'
    expression = '${So_Cu} * exp(-${E_S_Cu} / ${ideal_R} / T)'
  []
  [Solubility_W_func]
    type = ParsedFunction
    symbol_names = 'T'
    symbol_values = 'temp_func'
    expression = '${So_W} * exp(-${E_S_W} / ${ideal_R} / T)'
  []
[]

[Materials]
  # Temperature-dependent diffusivity and solubility materials plus interface jump ratio
  [diffusion_solubility]
    type = ADGenericFunctionMaterial
    prop_names = 'Diffusivity_Cu Diffusivity_W Solubility_Cu Solubility_W'
    prop_values = 'Diffusivity_Cu_func Diffusivity_W_func Solubility_Cu_func Solubility_W_func'
    outputs = 'exodus'
  []
  # Used in val-2b to check flux at boundaries. May be needed later.
  # [converter_to_nonAD]
  #   type = MaterialADConverter
  #   ad_props_in = 'diffusivity_Be diffusivity_BeO'
  #   reg_props_out = 'diffusivity_Be_nonAD diffusivity_BeO_nonAD'
  #   outputs = 'exodus'
  # []
  [interface_jump]
    type = SolubilityRatioMaterial
    solubility_primary = Solubility_Cu
    solubility_secondary = Solubility_W
    boundary = inter_Cu_to_W
    concentration_primary = u_Cu
    concentration_secondary = u_W
    outputs = 'exodus'
  []
[]

# Used while obtaining steady-state solution
#
[VectorPostprocessors]
  [line_Cu]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${T_Cu} 0 0'
    num_points = ${num_Cu}
    sort_by = 'x'
    variable = u_Cu
    outputs = vector_postproc_Cu
  []
  [line_W]
    type = LineValueSampler
    start_point = '${T_Cu} 0 0'
    end_point = '${fparse ${T_Cu} + ${T_W}} 0 0'
    num_points = ${num_W}
    sort_by = 'x'
    variable = u_W
    outputs = vector_postproc_W
  []
[]

[Postprocessors]
  # Used to obtain varying concentration with time at a
  # point in W layer 'x' um from Interface Cu/W boundary
  [concentration_at_x_W] #_x_SiC
    type = PointValue
    variable = u_W
    point = '${fparse ${T_Cu} + ${D_ver}} 0 0'
    outputs = 'csv'
  []
  [concentration_at_x_Cu] # _x_PyC
    type = PointValue
    variable = u_Cu
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
  [vector_postproc_Cu]
    type = CSV
    sync_times = ${simulation_time}
    sync_only = true
  []
  [vector_postproc_W]
    type = CSV
    sync_times = ${simulation_time}
    sync_only = true
  []
[]
