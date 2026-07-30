# ============================================================================
# Gate-free v20 structural stress profile.
#
# This intentionally extends beyond the identification bounds. It is a
# diagnostic of best achievable fit, not a coefficient-calibration script.
# ============================================================================

using Statistics

include("profile_2D_v20_ua_plant.jl")

const STRESS_EXPONENTS_2D_v20 = (
    0.00, 0.25, 0.50, 0.75, 1.00, 1.25,
)
const STRESS_NU50_RATIOS_2D_v20 = (
    1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00,
)

function stress_profile_2D_v20_ua(
    ; exponent_indices=1:length(STRESS_EXPONENTS_2D_v20),
    max_points=31, suffix="",
    felt_k=0.70, felt_cp=0.75,
    powers=(1.05, 1.23, 0.84),
)
    mkpath(OUTPUT_DIR_2D_v20)
    inputs = [
        case_inputs_2D_v20(id; max_points=max_points)
        for id in UA_PROFILE_IDS_2D_v20
    ]
    rows = NamedTuple[]
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_ua_profile$(suffix)_2D_v20.csv",
    )
    total = length(exponent_indices) *
            length(STRESS_NU50_RATIOS_2D_v20)
    counter = 0
    for exponent_index in exponent_indices
        exponent = STRESS_EXPONENTS_2D_v20[exponent_index]
        for ratio in STRESS_NU50_RATIOS_2D_v20
            counter += 1
            nu50 = ratio * UA_PROFILE_NU50_REFERENCE_2D_v20
            prefactor = nu50 / 50.0^exponent
            println(
                "v20 gate-free UA $counter/$total " *
                "n=$exponent ratio=$ratio",
            )
            flush(stdout)
            p = parameters_2D_v20(
                mesh=:nominal,
                source_model=:near_deep,
                deep_fraction=0.90,
                deep_length_m=0.12,
                nu_prefactor=prefactor,
                reynolds_exponent=exponent,
                distributed_tube_flange_h_W_m2_K=0.0,
                probe_capacity_areal_J_m2_K=1200.0,
                probe_stem_conductance_areal_W_m2_K=0.0,
                felt_conductivity_scale=felt_k,
                felt_heat_capacity_scale=felt_cp,
                felt_contact_scale=0.30,
                power_scales=powers,
                t3_location=:receiver_136,
            )
            cases = [
                simulate_case_2D_v20(
                    input, p; reltol=5e-4, abstol=1e-4,
                    dtmax=120.0,
                ) for input in inputs
            ]
            scores = _plant_case_score_2D_v20.(cases)
            aggregate = _aggregate_plant_score_2D_v20(scores)
            push!(rows, merge((
                reynolds_exponent=exponent,
                Nu50_ratio=ratio,
                Nu50=nu50,
                nu_prefactor=prefactor,
                felt_conductivity_scale=felt_k,
                felt_heat_capacity_scale=felt_cp,
                power_scale_456=powers[1],
                power_scale_304=powers[2],
                power_scale_256=powers[3],
                flow_converged=all(
                    all(case.result.flow_solver_converged)
                    for case in cases
                ),
            ), aggregate))
            _write_namedtuples_csv_2D_v20(output, rows)
        end
    end
    sort!(rows; by=row -> row.objective)
    _write_namedtuples_csv_2D_v20(output, rows)
    println("v20 gate-free UA best: ", first(rows))
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    first_index = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    last_index = length(ARGS) >= 2 ? parse(Int, ARGS[2]) :
                 length(STRESS_EXPONENTS_2D_v20)
    points = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 31
    suffix = length(ARGS) >= 4 ? ARGS[4] : ""
    stress_profile_2D_v20_ua(
        ; exponent_indices=first_index:last_index,
        max_points=points, suffix=suffix,
    )
end
