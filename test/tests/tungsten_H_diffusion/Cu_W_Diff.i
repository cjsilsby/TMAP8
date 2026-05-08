# Verification Problem #1e from TMAP4/TMAP7 V&V document
# Permeation problem in a composite layer
# No Soret effect, trapping, or solubility included.

# Numerical parameters
nx_num = 6000 # - 1000
simulation_time = '${units 0.2 s}' # previous 1
step = '${units 0.00001 s}' # previous 0.000001

# System properties
T_Cu = '${units 0.05 mm -> m}' # Thickness of Copper substrate. # T_PyC
T_W = '${units 0.025 mm -> m}' # Thickness of Tungsten film  # SiC
D_ver = '${units 0.02 mm -> m}' # Distance past T_Cu to calculate theoretical diffusion in tungsten layer # D_ver
D_ver_Cu = '${units 0.045 mm -> m}' # Distance past 0 to calculate theoretical diffusion on copper layer # D_ver_PyC
length_Cu = ${T_Cu} # '${units 1 mm}'

# Initial Conditions
initial_concentration = '${units 0.1 mol/m^3}'
temperature = '${units 1500 K}'

# Constants
k_B = '${units 8.61733e-5 eV/K}'
ideal_R = '${units 8.31446261815324 J/mol/K}'

# Material Propertoes
Do_Cu = '${units 1.74e-6 m^2/s}' # Diffusivity of H in Cu (m^2/s)
E_D_Cu = '${units 42000 J/mol}'
Diffusivity_Cu = '${units ${fparse Do_Cu * exp(-E_D_Cu/ideal_R/temperature)} m^2/s}'
Do_W = '${units 7.44e-8 m^2/s}' # Diffusivity of H in W (m^2/s)
E_D_W = '${units 0.13 eV}'
Diffusivity_W = '${units ${fparse Do_W * exp(-E_D_W/k_B/temperature)} m^2/s}'

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${nx_num}
  xmax = '${fparse ${T_Cu} + ${T_W} }'
  allow_renumbering = false
[]

[Variables]
  [u]
  []
[]

[Functions]
  # Diffusivity assign based on different material domains
  [diffusivity_value]
    type = ParsedFunction
    expression = 'if(x < ${length_Cu}, ${Diffusivity_Cu}, ${Diffusivity_W} )'
  []
[]

[Kernels]
  [diff]
    type = FunctionDiffusion
    variable = u
    function = diffusivity_value
  []
  [time]
    type = TimeDerivative
    variable = u
  []
[]

[AuxVariables]
  [T] # Temperature
    initial_condition = ${temperature}
  []
[]

[AuxKernels]
  [constant_temperature]
    type = ConstantAux
    variable = T
    value = '${temperature}'
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = ${initial_concentration}
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 0
  []
[]

# Used while obtaining steady-state solution
#
[VectorPostprocessors]
  [line]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${Mesh/xmax} 0 0'
    num_points = ${Mesh/nx}
    sort_by = 'x'
    variable = u
    outputs = vector_postproc
  []
[]

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
