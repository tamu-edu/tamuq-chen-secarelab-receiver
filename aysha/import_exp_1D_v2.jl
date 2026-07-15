# External experimental-data import for 1D_v2.jl.
#
# Contract intentionally follows import_exp_0D.jl:
#   measurements, measurements_cooling
#   simulation_conditions, simulation_conditions_cooling
#   IDs, IDs_cooling
#
# Required symbols supplied by 1D_v2.jl:
#   Io, qlpm, Tinit, T_in, T_amb

dec = 50
raw_data_path = joinpath(@__DIR__, "SolarSimulator", "RAW")

filenames = [
    "Data_FPT0067_231125_161757",
    "Data_FPT0068_231126_115725",
    "Data_FPT0069_231126_140153",
    "Data_FPT0070_231127_090339",
    "Data_FPT0071_231128_102707",
    "Data_FPT0072_231129_104140",
    "Data_FPT0073_231129_132744",
    "Data_FPT0074_231130_123228",
    "Data_FPT0075_231201_162138",
    "Data_FPT0076_231203_120521",
    "Data_FPT0077_231203_161315",
    "Data_FPT0078_231204_132252",
    "Data_FPT0079_231204_172244",
    "Data_FPT0080_231205_095122",
    "Data_FPT0081_231205_135354",
]

IDs = ["E67", "E68", "E69", "E70", "E71",
       "E72", "E73", "E74", "E75", "E76",
       "E77", "E78", "E79", "E80", "E81"]

filenames_cooling = [
    "Data_FPT0069-Cooling_231126_153148",
    "Data_FPT0080-cooling_231205_112837",
    "Data_FPT0081-cooling_231205_153409",
]

IDs_cooling = ["C69", "C80", "C81"]

# Nominal aperture irradiances reported in manuscript Table 2. No unexplained
# group multipliers are applied here.
ArIo = vcat(fill(456000.0, 5), fill(304000.0, 5), fill(256000.0, 5))
ArIo_cooling = zeros(3)

const MeasurementRowV2 = NamedTuple{
    (:simulation_id, :obs_id, :time, :temperatures),
    Tuple{String,String,Vector{Float64},Vector{Float64}},
}

measurements = MeasurementRowV2[]
measurements_cooling = MeasurementRowV2[]

experiment_metadata = NamedTuple{
    (:simulation_id, :phase, :irradiance, :mean_flow_lpm,
     :mean_inlet_temperature, :mean_ambient_temperature),
    Tuple{String,Symbol,Float64,Float64,Float64,Float64},
}[]

function decimate_vector(values, step=dec)
    indices = unique(vcat(collect(1:step:length(values)), length(values)))
    return Float64.(values[indices])
end

function read_receiver_file(filename)
    path = joinpath(raw_data_path, filename * ".csv")
    isfile(path) || error("Experimental file not found: $path")
    selected = (2, 7, 8, 9, 10, 36, 37, 42, 43, 44, 45, 46, 49, 50)
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        first_line = true
        for line in eachline(io)
            if first_line
                first_line = false
                continue
            end
            fields = split(line, ',')
            length(fields) >= maximum(selected) || continue
            try
                push!(rows, [parse(Float64, fields[i]) for i in selected])
            catch error
                error isa ArgumentError || rethrow()
            end
        end
    end
    isempty(rows) && error("No numeric experimental rows found in $path")
    return reduce(vcat, permutedims.(rows))
end

function import_receiver_experiment!(target, id, filename, irradiance, phase)
    data = read_receiver_file(filename)

    time = decimate_vector(data[:, 1])
    time .-= time[1]
    T2 = decimate_vector(data[:, 6]) .+ 273.15
    T3 = decimate_vector(data[:, 7]) .+ 273.15
    T8 = decimate_vector(data[:, 8]) .+ 273.15
    T9 = decimate_vector(data[:, 9]) .+ 273.15
    T10 = decimate_vector(data[:, 10]) .+ 273.15
    T11 = decimate_vector(data[:, 11]) .+ 273.15
    T12 = decimate_vector(data[:, 12]) .+ 273.15
    T15 = decimate_vector(data[:, 13]) .+ 273.15
    T16 = decimate_vector(data[:, 14]) .+ 273.15

    # The four controllers all carry air despite their inherited gas labels.
    flow_full = data[:, 2] .+ data[:, 3] .+ data[:, 4] .+ data[:, 5]
    flow = decimate_vector(flow_full)
    ambient = 0.5 .* (T15 .+ T16)
    inlet = copy(ambient)  # TI-1 is strongly radiation-biased in these tests.

    T_avg = (T8 .+ T9 .+ T10) ./ 3.0
    T_avg_v4 = 0.248 .* T8 .+ 0.365 .* T9 .+ 0.387 .* T10

    for (obs_id, values) in (
        ("_Tavg", T_avg), ("_Tavg_v4", T_avg_v4), ("_Tf", T3),
        ("_T2", T2), ("_T3", T3), ("_T8", T8), ("_T9", T9),
        ("_T10", T10), ("_T11", T11), ("_T12", T12),
        ("_flow", flow), ("_Tin", inlet), ("_Tamb", ambient),
    )
        push!(target, (simulation_id=id, obs_id=obs_id,
                       time=copy(time), temperatures=Float64.(values)))
    end

    push!(experiment_metadata, (
        simulation_id=id, phase=phase, irradiance=Float64(irradiance),
        mean_flow_lpm=mean(flow), mean_inlet_temperature=mean(inlet),
        mean_ambient_temperature=mean(ambient),
    ))
    return nothing
end

for i in eachindex(IDs)
    import_receiver_experiment!(measurements, IDs[i], filenames[i], ArIo[i], :heating)
end

for i in eachindex(IDs_cooling)
    import_receiver_experiment!(measurements_cooling, IDs_cooling[i],
                                filenames_cooling[i], ArIo_cooling[i], :cooling)
end

function observation(dataframe, simulation_id, obs_id)
    rows = filter(row -> row.simulation_id == simulation_id && row.obs_id == obs_id,
                  dataframe)
    length(rows) == 1 || error("Expected one row for $simulation_id/$obs_id; found $(length(rows))")
    return rows[1].temperatures
end

function observation_time(dataframe, simulation_id)
    rows = filter(row -> row.simulation_id == simulation_id && row.obs_id == "_Tf",
                  dataframe)
    length(rows) == 1 || error("Expected one time row for $simulation_id")
    return rows[1].time
end

simulation_conditions = Dict{String,Dict{Symbol,Float64}}()
for (i, id) in pairs(IDs)
    meta = only(filter(row -> row.simulation_id == id, experiment_metadata))
    simulation_conditions[id] = Dict(
        Io => ArIo[i],
        qlpm => meta.mean_flow_lpm,
        Tinit => observation(measurements, id, "_Tf")[1],
        T_in => meta.mean_inlet_temperature,
        T_amb => meta.mean_ambient_temperature,
    )
end

simulation_conditions_cooling = Dict{String,Dict{Symbol,Float64}}()
for (i, id) in pairs(IDs_cooling)
    meta = only(filter(row -> row.simulation_id == id, experiment_metadata))
    simulation_conditions_cooling[id] = Dict(
        Io => ArIo_cooling[i],
        qlpm => meta.mean_flow_lpm,
        Tinit => observation(measurements_cooling, id, "_Tf")[1],
        T_in => meta.mean_inlet_temperature,
        T_amb => meta.mean_ambient_temperature,
    )
end
