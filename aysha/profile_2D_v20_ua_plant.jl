# ============================================================================
# V20 T3-free plant profile for the integrated UA law.
#
# A and n are parameterized as Nu50 and n to expose their covariance. T3,
# T9, and T10 are diagnostics only and never enter the objective.
# ============================================================================

using Statistics

include("run_2D_v20.jl")

const UA_PROFILE_IDS_2D_v20 = ("E72", "E69", "E81")
const UA_PROFILE_EXPONENTS_2D_v20 = (
    1.25, 1.35, 1.45, 1.55, 1.65,
)
const UA_PROFILE_RATIOS_2D_v20 = (0.75, 1.35, 1.95)
const UA_PROFILE_NU50_REFERENCE_2D_v20 =
    3.07786e-4 * 50.0^1.44346

_profile_rho_2D_v20(value) =
    2.0 * (sqrt(1.0 + Float64(value)^2) - 1.0)

function _plant_case_score_2D_v20(case)
    rows = 2:size(case.model, 1)
    final_start = findfirst(
        time -> time >= first(case.inputs.times) +
                0.9 * (
                    last(case.inputs.times) -
                    first(case.inputs.times)
                ),
        case.inputs.times,
    )
    final = max(2, something(
        final_start, size(case.model, 1),
    )):size(case.model, 1)
    side_error = case.model[rows, 1:3] .-
                 case.observed[rows, 1:3]
    T2_error = case.model[rows, 7] .-
               case.observed[rows, 7]
    side_curve = mean(
        _profile_rho_2D_v20(value / 35.0)
        for value in side_error
    )
    side_level = mean(
        _profile_rho_2D_v20((
            mean(case.model[final, column]) -
            mean(case.observed[final, column])
        ) / 30.0) for column in 1:3
    )
    T2_curve = mean(
        _profile_rho_2D_v20(value / 15.0)
        for value in T2_error
    )
    T2_level = _profile_rho_2D_v20((
        mean(case.model[final, 7]) -
        mean(case.observed[final, 7])
    ) / 10.0)
    t3_error = case.model[rows, 6] .-
               case.observed[rows, 6]
    return (
        objective=0.55side_curve + 0.20side_level +
                  0.15T2_curve + 0.10T2_level,
        side_rmse_K=sqrt(mean(abs2, side_error)),
        side_final_bias_K=mean(
            case.model[final, 1:3] .-
            case.observed[final, 1:3],
        ),
        T2_rmse_K=sqrt(mean(abs2, T2_error)),
        T2_final_bias_K=mean(
            case.model[final, 7] .-
            case.observed[final, 7],
        ),
        T3_diagnostic_rmse_K=sqrt(mean(abs2, t3_error)),
    )
end

function _aggregate_plant_score_2D_v20(scores)
    return NamedTuple(
        name => mean(
            getproperty(score, name) for score in scores
        )
        for name in propertynames(first(scores))
    )
end

function profile_2D_v20_ua_plant(
    ; mesh=:nominal, max_points=31,
    exponent_indices=1:length(UA_PROFILE_EXPONENTS_2D_v20),
    suffix="",
)
    mkpath(OUTPUT_DIR_2D_v20)
    inputs = [
        case_inputs_2D_v20(id; max_points=max_points)
        for id in UA_PROFILE_IDS_2D_v20
    ]
    rows = NamedTuple[]
    total = length(exponent_indices) *
            length(UA_PROFILE_RATIOS_2D_v20)
    counter = 0
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "ua_plant_profile$(suffix)_2D_v20.csv",
    )
    for exponent_index in exponent_indices
        exponent = UA_PROFILE_EXPONENTS_2D_v20[exponent_index]
        for ratio in UA_PROFILE_RATIOS_2D_v20
            counter += 1
            nu50 = ratio * UA_PROFILE_NU50_REFERENCE_2D_v20
            prefactor = nu50 / 50.0^exponent
            println(
                "v20 UA plant $counter/$total n=$exponent " *
                "Nu50 ratio=$ratio",
            )
            flush(stdout)
            p = parameters_2D_v20(
                mesh=mesh,
                source_model=:near_deep,
                deep_fraction=0.90,
                deep_length_m=0.12,
                nu_prefactor=prefactor,
                reynolds_exponent=exponent,
                distributed_tube_flange_h_W_m2_K=0.0,
                probe_capacity_areal_J_m2_K=3000.0,
                probe_stem_conductance_areal_W_m2_K=0.0,
                felt_conductivity_scale=0.70,
                felt_heat_capacity_scale=0.75,
                felt_contact_scale=0.30,
                power_scales=(1.05, 1.23, 0.84),
                t3_location=:receiver_136,
            )
            cases = [
                simulate_case_2D_v20(
                    input, p; reltol=5e-4, abstol=1e-4,
                    dtmax=120.0,
                ) for input in inputs
            ]
            score = _aggregate_plant_score_2D_v20(
                _plant_case_score_2D_v20.(cases),
            )
            flow_converged = all(
                all(case.result.flow_solver_converged)
                for case in cases
            )
            push!(rows, merge((
                mesh=String(mesh),
                source_model="near_deep",
                deep_fraction=0.90,
                deep_length_m=0.12,
                nu_prefactor=prefactor,
                reynolds_exponent=exponent,
                Nu50=nu50,
                Nu50_ratio=ratio,
                flow_converged=flow_converged,
            ), score))
            _write_namedtuples_csv_2D_v20(output, rows)
        end
    end
    sort!(rows; by=row -> row.objective)
    _write_namedtuples_csv_2D_v20(output, rows)
    println("v20 UA plant best: ", first(rows))
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    mesh = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :nominal
    points = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 31
    first_index = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1
    last_index = length(ARGS) >= 4 ? parse(Int, ARGS[4]) :
                 length(UA_PROFILE_EXPONENTS_2D_v20)
    suffix = length(ARGS) >= 5 ? ARGS[5] : ""
    profile_2D_v20_ua_plant(
        ; mesh=mesh, max_points=points,
        exponent_indices=first_index:last_index,
        suffix=suffix,
    )
end
