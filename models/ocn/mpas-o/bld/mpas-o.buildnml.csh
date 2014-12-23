#!/bin/csh

# For now, manually build the namelist. Soon this will call the standard CESM
# build-namelist script.

mkdir -p $CASEROOT/CaseDocs

/bin/cp $CODEROOT/ocn/mpas-o/bld/m120/mpas-o.graph.info.part.* $CASEROOT/Buildconf/

# The script refences ocean120km.nc from the input_data location. So,
# ensure it's in the proper location:
# #DIN_LOC_ROOT/ocn/mpas-o/mpas120/ocean120km.nc


set MPAS_NML = $CASEROOT/CaseDocs/mpaso.in
touch $MPAS_NML
chmod 644 $MPAS_NML

set MPAS_STREAMS = $CASEROOT/CaseDocs/streams.ocean_forward
touch $MPAS_STREAMS
chmod 644 $MPAS_STREAMS

if ($CONTINUE_RUN == 'TRUE') then
    set config_do_restart = .true.
    set config_start_time = 'file'
	set config_set_restingThickness_to_IC = .false.
	set config_alter_ICs_for_pbcs = .false.
    #TODO - config_start_time must not be read in - but obtained from the coupler
else
    set config_do_restart = .false.
    set config_start_time = '0001-01-01_00:00:00' 
	set config_set_restingThickness_to_IC = .true.
	set config_alter_ICs_for_pbcs = .true.
endif

cat >! $MPAS_NML << EOF
&time_management
    config_do_restart = $config_do_restart
    config_Restart_timestamp_name = 'rpointer.ocn'
    config_start_time = $config_start_time
    config_stop_time = 'none'
    config_run_duration = '0001-00-00_00:00:00'
    config_calendar_type = 'gregorian_noleap'
/

&io
    config_stats_interval = '0001_00:00:00'
    config_write_stats_on_startup = .true.
    config_write_output_on_startup = .true.
    config_pio_num_iotasks = 0
    config_pio_stride = 1
/

&time_integration
    config_dt = '00:30:00'
    config_time_integrator = 'split_explicit'
/

&ALE_vertical_grid
    config_vert_coord_movement = 'uniform_stretching'
    config_use_min_max_thickness = .false.
    config_min_thickness = 1.0
    config_max_thickness_factor = 6.0
    config_set_restingThickness_to_IC = $config_set_restingThickness_to_IC
    config_dzdk_positive = .false.
/

&ALE_frequency_filtered_thickness
    config_use_freq_filtered_thickness = .false.
    config_thickness_filter_timescale = 5.0
    config_use_highFreqThick_restore = .false.
    config_highFreqThick_restore_time = 30.0
    config_use_highFreqThick_del2 = .false.
    config_highFreqThick_del2 = 100.0
/

&partial_bottom_cells
    config_alter_ICs_for_pbcs = $config_alter_ICs_for_pbcs
    config_pbc_alteration_type = 'partial_cell'
    config_min_pbc_fraction = 0.10
    config_check_ssh_consistency = .true.
/

&decomposition
    config_num_halos = 3
    config_block_decomp_file_prefix = '$CASEROOT/Buildconf/mpas-o.graph.info.part.'
    config_number_of_blocks = 0
    config_explicit_proc_decomp = .false.
    config_proc_decomp_file_prefix = '$CASEROOT/Buildconf/mpas-o.graph.info.part.'
/

&hmix
    config_hmix_scaleWithMesh = .false.
    config_maxMeshDensity = -1.0
    config_apvm_scale_factor = 0.0
/

&hmix_del2
    config_use_mom_del2 = .false.
    config_use_tracer_del2 = .false.
    config_mom_del2 = 10.0
    config_tracer_del2 = 10.0
/

&hmix_del4
    config_use_mom_del4 = .true.
    config_use_tracer_del4 = .false.
    config_mom_del4 = 2.6e13
    config_tracer_del4 = 0.0
/

&hmix_Leith
    config_use_Leith_del2 = .false.
    config_Leith_parameter = 1.0
    config_Leith_dx = 15000.0
    config_Leith_visc2_max = 2.5e3
/

&mesoscale_eddy_parameterization
    config_use_standardGM = .false.
    config_standardGM_tracer_kappa = 0.0
    config_Redi_kappa = 0.0
    config_gravWaveSpeed_trunc = 0.3
    config_max_relative_slope = 0.01
/

&hmix_del2_tensor
    config_use_mom_del2_tensor = .false.
    config_mom_del2_tensor = 10.0
/

&hmix_del4_tensor
    config_use_mom_del4_tensor = .false.
    config_mom_del4_tensor = 5.0e13
/

&Rayleigh_damping
    config_Rayleigh_friction = .false.
    config_Rayleigh_damping_coeff = 0.0
/

&vmix
    config_convective_visc = 1.0
    config_convective_diff = 1.0
/

&vmix_const
    config_use_const_visc = .false.
    config_use_const_diff = .false.
    config_vert_visc = 1.0e-4
    config_vert_diff = 1.0e-5
/

&vmix_rich
    config_use_rich_visc = .true.
    config_use_rich_diff = .true.
    config_bkrd_vert_visc = 1.0e-4
    config_bkrd_vert_diff = 1.0e-5
    config_rich_mix = 0.005
/

&vmix_tanh
    config_use_tanh_visc = .false.
    config_use_tanh_diff = .false.
    config_max_visc_tanh = 2.5e-1
    config_min_visc_tanh = 1.0e-4
    config_max_diff_tanh = 2.5e-2
    config_min_diff_tanh = 1.0e-5
    config_zMid_tanh = -100
    config_zWidth_tanh = 100
/

&cvmix
    config_use_cvmix = .false.
    config_cvmix_prandtl_number = 1.0
    config_use_cvmix_background = .true.
    config_cvmix_background_diffusion = 1.0e-5
    config_cvmix_background_viscosity = 1.0e-4
    config_use_cvmix_convection = .false.
    config_cvmix_convective_diffusion = 1.0
    config_cvmix_convective_viscosity = 1.0
    config_cvmix_convective_basedOnBVF = .true.
    config_cvmix_convective_triggerBVF = 0.0
    config_use_cvmix_shear = .false.
    config_cvmix_shear_mixing_scheme = 'PP'
    config_cvmix_shear_PP_nu_zero = 0.005
    config_cvmix_shear_PP_alpha = 5.0
    config_cvmix_shear_PP_exp = 2.0
    config_cvmix_shear_KPP_nu_zero = 0.005
    config_cvmix_shear_KPP_Ri_zero = 0.7
    config_cvmix_shear_KPP_exp = 3
    config_use_cvmix_tidal_mixing = .false.
    config_use_cvmix_double_diffusion = .false.
    config_use_cvmix_kpp = .false.
    config_cvmix_kpp_niterate = 2
    config_cvmix_kpp_criticalBulkRichardsonNumber = 0.25
    config_cvmix_kpp_matching = 'SimpleShapes'
    config_cvmix_kpp_EkmanOBL = .false.
    config_cvmix_kpp_MonObOBL = .false.
    config_cvmix_kpp_interpolationOMLType = 'quadratic'
    config_cvmix_kpp_surface_layer_extent = 0.1
/

&forcing
    config_forcing_type = 'bulk'
    config_restoreT_timescale = 30.0
    config_restoreS_timescale = 30.0
    config_restoreT_lengthscale = 50.0
    config_restoreS_lengthscale = 50.0
    config_flux_attenuation_coefficient = 0.001
    config_frazil_ice_formation = .true.
    config_sw_absorption_type = 'jerlov'
    config_jerlov_water_type = 3
    config_fixed_jerlov_weights = .true.
/

&advection
    config_vert_tracer_adv = 'stencil'
    config_vert_tracer_adv_order = 3
    config_horiz_tracer_adv_order = 3
    config_coef_3rd_order = 0.25
    config_monotonic = .true.
/

&bottom_drag
    config_bottom_drag_coeff = 1.0e-3
/

&pressure_gradient
    config_pressure_gradient_type = 'pressure_and_zmid'
    config_density0 = 1014.65
    config_common_level_weight = 0.5
/

&eos
    config_eos_type = 'jm'
/

&eos_linear
    config_eos_linear_alpha = 0.2
    config_eos_linear_beta = 0.8
    config_eos_linear_Tref = 5.0
    config_eos_linear_Sref = 35.0
    config_eos_linear_densityref = 1000.0
/

&split_explicit_ts
    config_n_ts_iter = 2
    config_n_bcl_iter_beg = 1
    config_n_bcl_iter_mid = 2
    config_n_bcl_iter_end = 2
    config_n_btr_subcycles = 20
    config_n_btr_cor_iter = 2
    config_vel_correction = .true.
    config_btr_subcycle_loop_factor = 2
    config_btr_gam1_velWt1 = 0.5
    config_btr_gam2_SSHWt1 = 1.0
    config_btr_gam3_velWt2 = 1.0
    config_btr_solve_SSH2 = .false.
/

&testing
    config_conduct_tests = .false.
    config_test_tensors = .false.
    config_tensor_test_function = 'sph_uCosCos'
/

&debug
    config_disable_redi_k33 = .false.
    config_disable_redi_horizontal_term1 = .false.
    config_disable_redi_horizontal_term2 = .false.
    config_disable_redi_horizontal_term3 = .false.
    config_check_zlevel_consistency = .false.
    config_filter_btr_mode = .false.
    config_prescribe_velocity = .false.
    config_prescribe_thickness = .false.
    config_include_KE_vertex = .false.
    config_check_tracer_monotonicity = .false.
    config_disable_thick_all_tend = .false.
    config_disable_thick_hadv = .false.
    config_disable_thick_vadv = .false.
    config_disable_thick_sflux = .true.
    config_disable_vel_all_tend = .false.
    config_disable_vel_coriolis = .false.
    config_disable_vel_pgrad = .false.
    config_disable_vel_hmix = .false.
    config_disable_vel_windstress = .false.
    config_disable_vel_vmix = .false.
    config_disable_vel_vadv = .false.
    config_disable_tr_all_tend = .false.
    config_disable_tr_adv = .false.
    config_disable_tr_hmix = .false.
    config_disable_tr_vmix = .false.
    config_disable_tr_sflux = .false.
    config_disable_tr_nonlocalflux = .false.
/

&global_stats
    config_use_global_stats = .true.
    config_global_stats_compute_interval = 'same_as_output'
    config_global_stats_compute_startup = .true.
/

&zonal_mean
    config_use_zonal_mean = .false.
    config_zonal_mean_compute_interval = 'same_as_output'
    config_zonal_mean_compute_startup = .true.
    config_number_zonal_mean_bins = 180
    config_min_zonal_mean_bin = -1.0e34
    config_max_zonal_mean_bin = -1.0e34
/
EOF

cat >! $MPAS_STREAMS << 'EOF'
<streams>

<immutable_stream name="mesh"
                  type="none"
                  filename_template="mesh_variables.nc"
/>

<immutable_stream name="input"
                  type="input"
'EOF'

# Breaking file to insert input file location.
cat >>! $MPAS_STREAMS << EOF
                  filename_template="$DIN_LOC_ROOT/ocn/mpas-o/mpas120/ocean120km.nc"
EOF

cat >>! $MPAS_STREAMS << 'EOF'
                  input_interval="initial_only"/>

<!--
The restart stream is actually controlled via the coupler.
Changing output_interval here will not have any affect on
the frequency restart files are written.

The output_interval is set to 1 second to ensure each restart frame has a
unique file.
-->
<immutable_stream name="restart"
                  type="input;output"
                  filename_template="rst/rst.ocn.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="output_interval"
                  reference_time="0000-01-01_00:00:00"
                  clobber_mode="truncate"
                  input_interval="initial_only"
                  output_interval="00-00-00_00:00:01"/>

<!--
output is the main history output stream. You can add auxiliary streams to
this stream to include more fields.
-->

<stream name="output"
        type="output"
        filename_template="hist/hist.ocn.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="00-01-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
        output_interval="00-00-01_00:00:00">

    <stream name="mesh"/>
    <stream name="real_world"/>
	<stream name="forcing"/>
    <var_array name="tracers"/>
    <var name="layerThickness"/>
    <var name="ssh"/>
    <var name="maxLevelEdgeTop"/>
    <var name="vertCoordMovementWeights"/>
    <var name="edgeMask"/>
    <var name="vertexMask"/>
    <var name="cellMask"/>
    <var name="refZMid"/>
    <var name="refLayerThickness"/>
    <var name="xtime"/>
    <var name="zMid"/>
    <var name="zTop"/>
    <var name="kineticEnergyCell"/>
    <var name="relativeVorticityCell"/>
    <var name="areaCellGlobal"/>
    <var name="areaEdgeGlobal"/>
    <var name="areaTriangleGlobal"/>
    <var name="volumeCellGlobal"/>
    <var name="volumeEdgeGlobal"/>
    <var name="CFLNumberGlobal"/>

</stream>

<!--
Streams between this line and the auxiliary stream line below are analysis member streams.
They can be used to perform online analysis of the simulation and control the output of
the analysis data.
-->

<stream name="globalStatsOutput"
        type="output"
        filename_template="analysis_members/ocn.globalStats.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        packages="amGlobalStats"
        clobber_mode="truncate"
        output_interval="0010_00:00:00">

    <var_array name="minGlobalStats"/>
    <var_array name="maxGlobalStats"/>
    <var_array name="sumGlobalStats"/>
    <var_array name="rmsGlobalStats"/>
    <var_array name="avgGlobalStats"/>
    <var_array name="vertSumMinGlobalStats"/>
    <var_array name="vertSumMaxGlobalStats"/>
    <var name="xtime"/>

</stream>

<stream name="zonalMeanOutput"
        type="output"
        filename_template="analysis_members/ocn.zonalMeans.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        packages="amZonalMean"
        clobber_mode="truncate"
        output_interval="0000_12:00:00">

    <var_array name="tracersZonalMean"/>
    <var name="xtime"/>
    <var name="binCenterZonalMean"/>
    <var name="binBoundaryZonalMean"/>
    <var name="velocityZonalZonalMean"/>
    <var name="velocityMeridionalZonalMean"/>

</stream>

<!--
All streams below this line are auxiliary streams. They are provided as
groupings of fields that one might be interested in. You can either enable the
stream to write a file for the fileds, or add the stream to another stream that
will already be written.  
-->

<stream name="additional_output"
        type="none"
        filename_template="hist/ocn.additional_output.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
>

    <var name="normalVelocity"/>
    <var name="density"/>
    <var name="pressure"/>
    <var name="divergence"/>
    <var name="viscosity"/>
    <var name="vertViscTopOfEdge"/>
    <var name="vertViscTopOfCell"/>
    <var name="vertDiffTopOfCell"/>
    <var name="BruntVaisalaFreqTop"/>
    <var name="RiTopOfCell"/>
    <var name="bulkRichardsonNumber"/>
    <var name="vertAleTransportTop"/>
    <var name="vertVelocityTop"/>

</stream>

<stream name="real_world"
        type="none"
        filename_template="hist/ocn.real_world_variables.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
>

    <stream name="mesh"/>
    <var name="normalVelocity"/>
    <var name="velocityZonal"/>
    <var name="velocityMeridional"/>
    <var name="displacedDensity"/>
    <var name="potentialDensity"/>
    <var name="boundaryLayerDepth"/>
    <var name="boundaryLayerDepthEdge"/>
    <var name="indexBoundaryLayerDepth"/>
    <var name="indexSurfaceLayerDepth"/>
    <var name="surfaceFrictionVelocity"/>
    <var name="windStressZonalDiag"/>
    <var name="windStressMeridionalDiag"/>
    <var name="surfaceBuoyancyForcing"/>
    <var name="seaSurfacePressure"/>

</stream>

<stream name="averages"
        type="none"
        filename_template="hist/ocn.average_variables.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
>

    <stream name="mesh"/>
    <var_array name="avgTracersSurfaceValue"/>
    <var_array name="avgSurfaceVelocity"/>
    <var_array name="avgSSHGradient"/>
    <var name="nAverage"/>
    <var name="avgSSH"/>
    <var name="avgNormalVelocity"/>
    <var name="avgVelocityZonal"/>
    <var name="avgVelocityMeridional"/>
    <var name="avgVertVelocityTop"/>
    <var name="avgNormalTransportVelocity"/>
    <var name="avgTransportVelocityZonal"/>
    <var name="avgTransportVelocityMeridional"/>
    <var name="avgVertTransportVelocityTop"/>
    <var name="varSSH"/>
    <var name="varNormalVelocity"/>
    <var name="varVelocityZonal"/>
    <var name="varVelocityMeridional"/>

</stream>

<stream name="Cartesian"
        type="none"
        filename_template="hist/ocn.Cartesian_variables.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
>

    <stream name="mesh"/>
    <var name="velocityX"/>
    <var name="velocityY"/>

</stream>

<stream name="forcing"
        type="none"
        filename_template="hist/ocn.forcing_variables.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
>

    <stream name="mesh"/>
    <var_array name="tracersSurfaceValue"/>
    <var_array name="surfaceVelocity"/>
    <var_array name="SSHGradient"/>
    <var_array name="surfaceTracerFlux"/>
    <var_array name="vertNonLocalFlux"/>
    <var name="surfaceWindStressMagnitude"/>
    <var name="surfaceWindStress"/>
    <var name="surfaceMassFlux"/>
    <var name="seaIceEnergy"/>
    <var name="penetrativeTemperatureFlux"/>
    <var name="transmissionCoefficients"/>
    <var name="windStressZonal"/>
    <var name="windStressMeridional"/>
    <var name="latentHeatFlux"/>
    <var name="sensibleHeatFlux"/>
    <var name="longWaveHeatFluxUp"/>
    <var name="longWaveHeatFluxDown"/>
    <var name="seaIceHeatFlux"/>
    <var name="shortWaveHeatFlux"/>
    <var name="evaporationFlux"/>
    <var name="seaIceSalinityFlux"/>
    <var name="seaIceFreshWaterFlux"/>
    <var name="riverRunoffFlux"/>
    <var name="iceRunoffFlux"/>
    <var name="rainFlux"/>
    <var name="snowFlux"/>
    <var name="iceFraction"/>
    <var name="prognosticCO2"/>
    <var name="diagnosticCO2"/>
    <var name="squaredWindSpeed10Meter"/>
    <var name="CO2Flux"/>
    <var name="DMSFlux"/>
    <var name="nAccumulatedCoupled"/>
    <var name="thermalExpansionCoeff"/>
    <var name="salineContractionCoeff"/>

</stream>

<stream name="Gent_McWilliams_spherical"
        type="none"
        filename_template="hist/ocn.Gent_McWilliams_spherical_variables.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
>

    <stream name="mesh"/>
    <var name="relativeSlopeTopOfCell"/>
    <var name="relativeSlopeTaperingCell"/>
    <var name="relativeSlopeTopOfCellZonal"/>
    <var name="relativeSlopeTopOfCellMeridional"/>
    <var name="k33"/>
    <var name="GMBolusVelocityZonal"/>
    <var name="GMBolusVelocityMeridional"/>
    <var name="normalGMBolusVelocity"/>
    <var name="vertGMBolusVelocityTop"/>
    <var name="gmStreamFuncTopOfEdge"/>
    <var name="avgNormalGMBolusVelocity"/>
    <var name="avgGMBolusVelocityZonal"/>
    <var name="avgGMBolusVelocityMeridional"/>
    <var name="avgVertGMBolusVelocityTop"/>

</stream>

<stream name="Gent_McWilliams_Cartesian"
        type="none"
        filename_template="hist/ocn.Gent_McWilliams_Cartesian_variables.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="01-00-00_00:00:00"
        reference_time="0000-01-01_00:00:00"
        clobber_mode="truncate"
>

    <stream name="mesh"/>
    <var name="relativeSlopeTopOfCell"/>
    <var name="relativeSlopeTaperingCell"/>
    <var name="relativeSlopeTopOfCellX"/>
    <var name="relativeSlopeTopOfCellY"/>
    <var name="relativeSlopeTopOfCellZ"/>
    <var name="k33"/>
    <var name="GMBolusVelocityX"/>
    <var name="GMBolusVelocityY"/>
    <var name="normalGMBolusVelocity"/>
    <var name="vertGMBolusVelocityTop"/>
    <var name="gmStreamFuncTopOfEdge"/>
    <var name="avgNormalGMBolusVelocity"/>
    <var name="GMStreamFuncX"/>
    <var name="GMStreamFuncY"/>

</stream>



</streams>
'EOF'

/bin/cp $MPAS_NML $RUNDIR
/bin/cp $MPAS_STREAMS $RUNDIR
