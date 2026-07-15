using Plots, XLSX, DataFrames, CSV, Statistics, Dates, LinearAlgebra, DelimitedFiles

dec = 50
function decimate!(x, step)
    x = [x[i] for i in 1:step:length(x)]
end
function getT(df, col, dec)
    temp = df[:, col]
    temp = decimate!(temp, dec)
    return temp
end
path = "./SolarSimulator/RAW/"

    #Exp. data to extract temp.
    measurements = DataFrame(
        simulation_id=String[],
        obs_id=String[],
        time=Array[],
        temperatures=Array[]
    )
    measurements_cooling = DataFrame(
        simulation_id=String[],
        obs_id=String[],
        time=Array[],
        temperatures=Array[]
    )
    filenames=["Data_FPT0067_231125_161757",
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
        "Data_FPT0081_231205_135354"]
    IDs=["E67", "E68", "E69", "E70", "E71", 
        "E72", "E73", "E74", "E75", "E76", 
        "E77", "E78", "E79", "E80", "E81"]

    filenames_cooling=["Data_FPT0069-Cooling_231126_153148",
        "Data_FPT0080-cooling_231205_112837",
        "Data_FPT0081-cooling_231205_153409"]
    IDs_cooling=["C69", "C80", "C81"]

    IoA = 456000.0 * 1.
    IoB = 304000.0 * 1.15
    IoC = 256000.0 * 0.7 # arbitrary fix KK

    ArIo=[IoA, IoA, IoA, IoA, IoA,
        IoB, IoB, IoB, IoB, IoB,
        IoC, IoC, IoC, IoC, IoC] * 1.0
    Arqplm=[15.27, 12.50, 10.50, 9.10, 7.12,
        18.34, 13.16, 9.03, 6.95, 4.53,
        13.85, 10.02, 8.04, 6.62, 4.53]
    Arqplm_cooling=[10.50, 6.62, 4.53]

    # Arqplm=[1.31, 1.19, 1.12, 1.064, 0.93,
    #      1.58, 1.22, 1.06, 0.92, 0.61,
    #      1.24, 1.11, 1.00, 0.88, 0.61]

    function rd_data(df, iID, IDs)
        f = CSV.File(path * filenames[iID] * ".csv"; skipto=2, delim=",", header=false) |> DataFrame
        t = getT(f, 2, dec) #float.(Z[:, 2]) #xz_data
        T2 = getT(f, 36, dec) .+ 273.15 #insulation
        T3 = getT(f, 37, dec) .+ 273.15 #y1z_data around 136 mm in
        T8 = getT(f, 42, dec) .+ 273.15 #around 5 mm in
        T9 = getT(f, 43, dec) .+ 273.15 #58mm (wall)
        T10 = getT(f, 44, dec) .+ 273.15 #107mm (wall)
        T11 = getT(f, 45, dec) .+ 273.15 #107mm (in)
        T12 = getT(f, 46, dec) .+ 273.15 #58mm (in)
  
        Tf = T3
        #T2 = (T9 .+ T11) ./ 2
        #T3 = (T12 .+ T10) ./ 2
        # Average of T8, T9, and T10
        #T_avg = (T9 .+ T10 .+ T11 .+ T12) ./ 4
        T_avg = (T8 .+ T9 .+ T10) ./ 3
        #T_avg = (T8 .+ T9 .+ T10 .+ T11 .+ T12 ./ 5)
        #T_avg = (T9 .+ T10) ./ 2

        scatter(t, Tf, ylim=(200, 1200))
        plot!(t, T_avg, label="T_avg")
        #plot!(t, T2)
        #plot!(t, T3)
        push!(df, (IDs[iID], "_Tavg", t, T_avg)) #first and second position for the 0D model
        push!(df, (IDs[iID], "_Tf", t, Tf))
        push!(df, (IDs[iID], "_T2", t, T2))
        push!(df, (IDs[iID], "_T3", t, T3))
        push!(df, (IDs[iID], "_T8", t, T8))
        push!(df, (IDs[iID], "_T9", t, T9))
        push!(df, (IDs[iID], "_T10", t, T10))
        push!(df, (IDs[iID], "_T11", t, T11))
        push!(df, (IDs[iID], "_T12", t, T12))

    end

    #rd_data(measurements, i) for i=1:length(filenames)
    map(x -> rd_data(measurements, x, IDs), 1:length(filenames))

# create a dataframe with the last temperature (steady state) reported for each thermocouple [Tx] and simulation ID [EXX]
# each row is a simulation ID and each column is a thermocouple
    # Create a DataFrame with the last temperature for each thermocouple
    last_temps = DataFrame(simulation_id=String[], T2=Float64[], T3=Float64[], T8=Float64[], T9=Float64[], T10=Float64[], T11=Float64[], T12=Float64[])
    for id in unique(measurements.simulation_id)
        row = (simulation_id=id,
               T2=measurements[(measurements.simulation_id .== id) .& (measurements.obs_id .== "_T2"), :temperatures][1][end],
               T3=measurements[(measurements.simulation_id .== id) .& (measurements.obs_id .== "_T3"), :temperatures][1][end],
               T8=measurements[(measurements.simulation_id .== id) .& (measurements.obs_id .== "_T8"), :temperatures][1][end],
               T9=measurements[(measurements.simulation_id .== id) .& (measurements.obs_id .== "_T9"), :temperatures][1][end],
               T10=measurements[(measurements.simulation_id .== id) .& (measurements.obs_id .== "_T10"), :temperatures][1][end],
               T11=measurements[(measurements.simulation_id .== id) .& (measurements.obs_id .== "_T11"), :temperatures][1][end],
               T12=measurements[(measurements.simulation_id .== id) .& (measurements.obs_id .== "_T12"), :temperatures][1][end])
        push!(last_temps, row)
    end
    # Save the measurements DataFrame to a CSV file
    CSV.write("measurements.csv", last_temps)
           
#Defining simulation conditions
simulation_conditions = Dict(
   IDs[i] => Dict(Io => ArIo[i], qlpm => Arqplm[i], Tinit => measurements[(measurements.simulation_id .== IDs[i]) .& (measurements.obs_id .== "_Tf"), :temperatures][1][1]) 
   for i=1:length(IDs)
)

# process the cooling measurements
map(x -> rd_data(measurements_cooling, x, IDs_cooling), 1:length(filenames_cooling))
last_temps = DataFrame(simulation_id=String[], T2=Float64[], T3=Float64[], T8=Float64[], T9=Float64[], T10=Float64[], T11=Float64[], T12=Float64[])
    for id in unique(measurements_cooling.simulation_id)
        row = (simulation_id=id,
               T2=measurements_cooling[(measurements_cooling.simulation_id .== id) .& (measurements_cooling.obs_id .== "_T2"), :temperatures][1][end],
               T3=measurements_cooling[(measurements_cooling.simulation_id .== id) .& (measurements_cooling.obs_id .== "_T3"), :temperatures][1][end],
               T8=measurements_cooling[(measurements_cooling.simulation_id .== id) .& (measurements_cooling.obs_id .== "_T8"), :temperatures][1][end],
               T9=measurements_cooling[(measurements_cooling.simulation_id .== id) .& (measurements_cooling.obs_id .== "_T9"), :temperatures][1][end],
               T10=measurements_cooling[(measurements_cooling.simulation_id .== id) .& (measurements_cooling.obs_id .== "_T10"), :temperatures][1][end],
               T11=measurements_cooling[(measurements_cooling.simulation_id .== id) .& (measurements_cooling.obs_id .== "_T11"), :temperatures][1][end],
               T12=measurements_cooling[(measurements_cooling.simulation_id .== id) .& (measurements_cooling.obs_id .== "_T12"), :temperatures][1][end])
        push!(last_temps, row)
    end
# Save the measurements DataFrame to a CSV file
CSV.write("measurements_cooling.csv", last_temps)

# Defining cooling simulation conditions
simulation_conditions_cooling = Dict(
   IDs_cooling[i] => Dict(qlpm => Arqplm_cooling[i], Tinit => measurements_cooling[(measurements_cooling.simulation_id .== IDs_cooling[i]) .& (measurements_cooling.obs_id .== "_Tf"), :temperatures][1][1]) 
   for i=1:length(IDs_cooling)
)

