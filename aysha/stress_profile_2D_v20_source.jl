# ============================================================================
# Gate-free v20 source-shape extension at the selected structural stress cell.
#
# T3, T9, and T10 remain diagnostics. Candidate selection uses only the three
# side thermocouples and T2 through the existing v20 plant objective.
# This is a structural stress test, not a source-parameter calibration.
# ============================================================================

using Statistics

include("profile_2D_v20_ua_plant.jl")

const STRESS_SOURCE_IDS_2D_v20 = ("E72", "E69", "E81")
const STRESS_SOURCE_EXPONENT_2D_v20 = 0.50
const STRESS_SOURCE_NU50_RATIO_2D_v20 = 2.00
const STRESS_SOURCE_FELT_K_2D_v20 = 0.40
const STRESS_SOURCE_FELT_CP_2D_v20 = 0.55
const STRESS_SOURCE_POWERS_2D_v20 = (1.05, 1.23, 0.70)

const STRESS_SOURCE_CANDIDATES_2D_v20 = (
    (
        label="beer_lambert",
        model=:beer_lambert,
        deep_fraction=0.0,
        deep_length_m=0.08,
    ),
    (
        label="near_deep_f0.9_L120mm",
        model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
    ),
    (
        label="near_deep_f0.9_L200mm",
        model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.20,
    ),
    (
        label="near_deep_f1.0_L120mm",
        model=:near_deep,
        deep_fraction=1.00,
        deep_length_m=0.12,
    ),
    (
        label="near_deep_f1.0_L200mm",
        model=:near_deep,
        deep_fraction=1.00,
        deep_length_m=0.20,
    ),
)

function _stress_source_parameters_2D_v20(candidate)
    nu50 = STRESS_SOURCE_NU50_RATIO_2D_v20 *
           UA_PROFILE_NU50_REFERENCE_2D_v20
    return parameters_2D_v20(
        mesh=:nominal,
        source_model=candidate.model,
        deep_fraction=candidate.deep_fraction,
        deep_length_m=candidate.deep_length_m,
        nu_prefactor=nu50 / 50.0^STRESS_SOURCE_EXPONENT_2D_v20,
        reynolds_exponent=STRESS_SOURCE_EXPONENT_2D_v20,
        distributed_tube_flange_h_W_m2_K=0.0,
        probe_capacity_areal_J_m2_K=1200.0,
        probe_stem_conductance_areal_W_m2_K=0.0,
        felt_conductivity_scale=STRESS_SOURCE_FELT_K_2D_v20,
        felt_heat_capacity_scale=STRESS_SOURCE_FELT_CP_2D_v20,
        felt_contact_scale=0.30,
        power_scales=STRESS_SOURCE_POWERS_2D_v20,
        t3_location=:receiver_136,
    )
end

function stress_profile_2D_v20_source(
    ; max_points=31, suffix="",
)
    mkpath(OUTPUT_DIR_2D_v20)
    inputs = [
        case_inputs_2D_v20(id; max_points=max_points)
        for id in STRESS_SOURCE_IDS_2D_v20
    ]
    rows = NamedTuple[]
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_source_profile$(suffix)_2D_v20.csv",
    )
    for (index, candidate) in
        enumerate(STRESS_SOURCE_CANDIDATES_2D_v20)
        println(
            "v20 source stress $index/" *
            "$(length(STRESS_SOURCE_CANDIDATES_2D_v20)) " *
            candidate.label,
        )
        flush(stdout)
        cases = [
            simulate_case_2D_v20(
                input, _stress_source_parameters_2D_v20(candidate);
                reltol=5e-4, abstol=1e-4, dtmax=120.0,
            ) for input in inputs
        ]
        score = _aggregate_plant_score_2D_v20(
            _plant_case_score_2D_v20.(cases),
        )
        flow = all(
            all(case.result.flow_solver_converged)
            for case in cases
        )
        push!(rows, merge((
            candidate=candidate.label,
            source_model=String(candidate.model),
            deep_fraction=candidate.deep_fraction,
            deep_length_m=candidate.deep_length_m,
            reynolds_exponent=STRESS_SOURCE_EXPONENT_2D_v20,
            Nu50_ratio=STRESS_SOURCE_NU50_RATIO_2D_v20,
            felt_conductivity_scale=STRESS_SOURCE_FELT_K_2D_v20,
            felt_heat_capacity_scale=STRESS_SOURCE_FELT_CP_2D_v20,
            power_scale_456=STRESS_SOURCE_POWERS_2D_v20[1],
            power_scale_304=STRESS_SOURCE_POWERS_2D_v20[2],
            power_scale_256=STRESS_SOURCE_POWERS_2D_v20[3],
            fraction_at_upper_boundary=(
                candidate.model === :near_deep &&
                candidate.deep_fraction == 1.0
            ),
            length_at_upper_boundary=(
                candidate.model === :near_deep &&
                candidate.deep_length_m == 0.20
            ),
            flow_converged=flow,
        ), score))
        _write_namedtuples_csv_2D_v20(output, rows)
    end
    sort!(rows; by=row -> (
        row.flow_converged ? 0 : 1, row.objective,
    ))
    _write_namedtuples_csv_2D_v20(output, rows)
    println("v20 source stress best: ", first(rows))
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    points = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 31
    suffix = length(ARGS) >= 2 ? ARGS[2] : ""
    stress_profile_2D_v20_source(
        ; max_points=points, suffix=suffix,
    )
end
