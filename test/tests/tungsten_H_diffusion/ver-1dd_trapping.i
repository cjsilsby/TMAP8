# Verification Problem #1dd from TMAP7 V&V document
# Permeation Problem without Trapping sites
# No Soret effect or solubility included.

# modeling parameters
nx_num = 5000 # (-)
diffusivity = '${units 1.632 mm^2/s}' # ----- 1.632 mm^2/min (H2 in W)-- input as mm^2/s
simulation_time = '${units 0.0005 s}' # actually min
interval_time_min = '${units 0.000001 s}' # actually min # 0.00000001 s needed for mass balance
interval_time = '${units 0.000001 s}' # actually min
thickness = '${units 0.025 mm}'
temperature = '${units 1500 K}' # Trapping has limited effect at 1500 K

# Trapping parameters
cl = '${units 3.1622e18 at/m^3 -> at/mm^3}'
N = '${units 3.1622e22 at/m^3 -> at/mm^3}'
trapping_prefactor = '${fparse ${units 1e15 1/s} / time_scaling}' # assumed from TMAP8 documentation. No alternative rate in literature
release_prefactor = '${fparse ${units 1e13 1/s} / time_scaling}' # assumed from TMAP8 documentation. No alternative rate in literature
epsilon = ${units 9864 K} # epsilon/k, (trapping energy)/(Boltzman constant) = 0.85 eV / 8.61733326e-5 eV/K = 9980 K
trapping_fraction = 0.01 # fraction of surface sites that are trapping sites
trap_per_free = 1e3 # C_T multiplier for numerical purposes

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${nx_num}
  xmax = ${thickness}
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Variables]
  [mobile]
  []
  [trapped]
  []
[]

[AuxVariables]
  [empty_sites]
  []
  [scaled_empty_sites]
  []
  [trapped_sites]
  []
  [total_sites]
  []
[]

[AuxKernels]
  [empty_sites]
    variable = empty_sites
    type = EmptySitesAux
    N = '${fparse N / cl}'
    Ct0 = ${trapping_fraction}
    trap_per_free = ${trap_per_free}
    trapped_concentration_variables = trapped
  []
  [scaled_empty]
    variable = scaled_empty_sites
    type = NormalizationAux
    normal_factor = ${cl}
    source_variable = empty_sites
  []
  [trapped_sites]
    variable = trapped_sites
    type = NormalizationAux
    normal_factor = ${trap_per_free}
    source_variable = trapped
  []
  [total_sites]
    variable = total_sites
    type = ParsedAux
    expression = 'trapped_sites + empty_sites'
    coupled_variables = 'trapped_sites empty_sites'
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    variable = mobile
    diffusivity = ${diffusivity}
    extra_vector_tags = ref
  []
  [time]
    type = TimeDerivative
    variable = mobile
    extra_vector_tags = ref
  []
  [coupled_time]
    type = ScaledCoupledTimeDerivative
    variable = mobile
    v = trapped
    factor = ${trap_per_free}
    extra_vector_tags = ref
  []
[]

[NodalKernels]
  [time]
    type = TimeDerivativeNodalKernel
    variable = trapped
  []
  [trapping]
    type = TrappingNodalKernel
    variable = trapped
    alpha_t = ${trapping_prefactor}
    N = '${fparse N / cl}'
    Ct0 = ${trapping_fraction}
    mobile_concentration = 'mobile'
    temperature = ${temperature}
    trap_per_free = ${trap_per_free}
    extra_vector_tags = ref
  []
  [release]
    type = ReleasingNodalKernel
    alpha_r = ${release_prefactor}
    temperature = ${temperature}
    detrapping_energy = ${epsilon}
    variable = trapped
  []
[]

[BCs]
  [left]
    type = FunctionDirichletBC
    variable = mobile
    function = 'BC_func'
    boundary = left
  []
  [right]
    type = DirichletBC
    variable = mobile
    value = 0
    boundary = right
  []
[]
[Functions]
  [BC_func]
    type = ParsedFunction
    expression = '${fparse cl / cl}*tanh( 3 * t )'
  []
[]

[Postprocessors]
  ############
  [mass_in_domain] # Calculates total mass in the domain
    type = ElementIntegralVariablePostprocessor
    variable = mobile
    # outputs = csv_data
  []

  [influx]
    type = ADSideDiffusiveFluxIntegral
    boundary = left
    variable = mobile
    diffusivity = ${diffusivity}
    # outputs = csv_data
  []
 
  [outflux]
    type = ADSideDiffusiveFluxIntegral
    boundary = right
    variable = mobile
    diffusivity = ${diffusivity}
    # outputs = csv_data
  []

  [scaled_outflux]
    type = ScalePostprocessor
    value = outflux
    scaling_factor = ${cl}
  []

  [min_trapped]
    type = NodalExtremeValue
    value_type = MIN
    variable = trapped
  []

  [flux_difference] # Ensure that we are accounting for atomic vs molecular hydrogen
    type = ParsedPostprocessor
    # expression = 'scaled_outflux - scaled_influx'
    # pp_names = 'scaled_influx scaled_outflux'
    # expression = 'outflux - influx'
    expression = '-influx - outflux'
    pp_names = 'influx outflux'
    # outputs = csv_data
  []

  [time_integrated_flux]
    type = TimeIntegratedPostprocessor
    value = flux_difference
    # outputs = csv_data
  [] 
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = ${simulation_time}
  dt = ${interval_time}
  dtmin = ${interval_time_min}
  solve_type = NEWTON
  scheme = BDF2
  nl_abs_tol = 1e-13
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  #automatic_scaling = true
  #verbose = true
  #compute_scaling_once = false
  line_search = 'none'
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-6
    optimal_iterations = 9
    growth_factor = 1.1
    cutback_factor = 0.909
  []
[]

[Outputs]
  exodus = true
  csv = true
  [dof]
    type = DOFMap
    execute_on = initial
  []
  perf_graph = true
[]
