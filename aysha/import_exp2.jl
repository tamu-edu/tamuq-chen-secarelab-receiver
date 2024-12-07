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
    
    IoA = 456000.0 * 1.
    IoB = 304000.0 * 1.
    IoC = 256000.0 * 1.

    ArIo=[IoA, IoA, IoA, IoA, IoA,
        IoB, IoB, IoB, IoB, IoB,
        IoC, IoC, IoC, IoC, IoC]

    Arqplm=[15.27, 12.50, 10.50, 9.10, 7.12,
        18.34, 13.16, 9.03, 6.95, 4.53,
        13.85, 10.02, 8.04, 6.62, 4.53]

    function rd_data(df, iID)
        f = CSV.File(path * filenames[iID] * ".csv"; skipto=2, delim=",", header=false) |> DataFrame
        t = getT(f, 2, dec) #float.(Z[:, 2]) #xz_data
        T3 = getT(f, 37, dec) .+ 273.15 #y1z_data around 136 mm in
        T8 = getT(f, 42, dec) .+ 273.15 #around 5 mm in
        T9 = getT(f, 43, dec) .+ 273.15 #58mm (wall)
        T10 = getT(f, 44, dec) .+ 273.15 #107mm (wall)
        T11 = getT(f, 45, dec) .+ 273.15 #107mm (in)
        T12 = getT(f, 46, dec) .+ 273.15 #58mm (in)
  
        Tf = T3
        T2 = (T9 .+ T11) ./ 2
        T3 = (T12 .+ T10) ./ 2
        scatter(t, Tf, ylim=(200, 1200))
        plot!(t, T8)
        plot!(t, T2)
        plot!(t, T3)
        push!(df, (IDs[iID], "_T8", t, T8))
        push!(df, (IDs[iID], "_T2", t, T2))
        push!(df, (IDs[iID], "_T3", t, T3))
        push!(df, (IDs[iID], "_Tf", t, Tf))

    end

    #rd_data(measurements, i) for i=1:length(filenames)
    map(x -> rd_data(measurements, x), 1:length(filenames))   
    
#Defining simulation conditions
    simulation_conditions = Dict(
        IDs[i] => Dict(Io => ArIo[i], qlpm => Arqplm[i], Tinit => measurements[(measurements.simulation_id .== IDs[i]) .& (measurements.obs_id .== "_Tf"), :temperatures][1][1]) 
        for i=1:length(IDs)
    )