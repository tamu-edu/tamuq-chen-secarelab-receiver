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
    begin
        #Exp 67 - T3, T8
        #Z = XLSX.readxlsx("./SolarSimulator/RAW/Data_FPT0067_231125_161757.xlsx")["Sheet 1 - Data_FPT0067_231125_1"]["A3:C3932"]
        Z = CSV.File(path * "Data_FPT0067_231125_161757.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E67t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E67T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E67T8 = getT(Z, 42, dec) .+ 273.15
        E67T9 = getT(Z, 43, dec) .+ 273.15
        E67T10 = getT(Z, 44, dec) .+ 273.15
        E67T11 = getT(Z, 45, dec) .+ 273.15
        E67T12 = getT(Z, 46, dec) .+ 273.15
  
        E67Tf = E67T3
        E67Ts2 = (E67T9 .+ E67T11) ./ 2
        E67Ts3 = (E67T12 .+ E67T10) ./ 2
        scatter(E67t, E67Tf, ylim=(200, 1200))
        plot!(E67t, E67T8)
        plot!(E67t, E67Ts2)
        plot!(E67t, E67Ts3)
    end
    begin
        #Exp 68 - T3, T8
        #A1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0068_231126_115725.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:C5365"]
        Z = CSV.File(path * "Data_FPT0068_231126_115725.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E68t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E68T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E68T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E68T8 = getT(Z, 42, dec) .+ 273.15
        E68T9 = getT(Z, 43, dec) .+ 273.15
        E68T10 = getT(Z, 44, dec) .+ 273.15
        E68T11 = getT(Z, 45, dec) .+ 273.15
        E68T12 = getT(Z, 46, dec) .+ 273.15
        
        E68Tf = E68T3
        E68Ts2 = (E68T9 .+ E68T11) ./ 2
        E68Ts3 = (E68T12 .+ E68T10) ./ 2
        scatter(E68t, E68Tf, ylim=(200, 1200))
        plot!(E68t, E68T8)
        plot!(E68t, E68Ts2)
        plot!(E68t, E68Ts3)
    end

    begin
        #Exp 69 - T3, T8
        #B1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0069_231126_140153.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:C5366"]
        Z = CSV.File(path * "Data_FPT0069_231126_140153.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E69t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E69T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E69T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E69T8 = getT(Z, 42, dec) .+ 273.15
        E69T9 = getT(Z, 43, dec) .+ 273.15
        E69T10 = getT(Z, 44, dec) .+ 273.15
        E69T11 = getT(Z, 45, dec) .+ 273.15
        E69T12 = getT(Z, 46, dec) .+ 273.15
        
        E69Tf = E69T3
        E69Ts2 = (E69T9 .+ E69T11) ./ 2
        E69Ts3 = (E69T12 .+ E69T10) ./ 2
        scatter(E69t, E69Tf, ylim=(200, 1200))
        plot!(E69t, E69T8)
        plot!(E69t, E69Ts2)
        plot!(E69t, E69Ts3)
    end


    begin
        #Exp 70 - T3, T8
        #C1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0070_231127_090339.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:C6705"]
        Z = CSV.File(path * "Data_FPT0070_231127_090339.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E70t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E70T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E70T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E70T8 = getT(Z, 42, dec) .+ 273.15
        E70T9 = getT(Z, 43, dec) .+ 273.15
        E70T10 = getT(Z, 44, dec) .+ 273.15
        E70T11 = getT(Z, 45, dec) .+ 273.15
        E70T12 = getT(Z, 46, dec) .+ 273.15
        
        E70Tf = E70T3
        E70Ts2 = (E70T9 .+ E70T11) ./ 2
        E70Ts3 = (E70T12 .+ E70T10) ./ 2
        scatter(E70t, E70Tf, ylim=(200, 1200))
        plot!(E70t, E70T8)
        plot!(E70t, E70Ts2)
        plot!(E70t, E70Ts3)
    end
    begin
        #Exp 71 - T3, T8
        #D1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_231128_102707.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:C7087"]
        Z = CSV.File(path * "Data_FPT0071_231128_102707.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E71t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E71T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E71T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E71T8 = getT(Z, 42, dec) .+ 273.15
        E71T9 = getT(Z, 43, dec) .+ 273.15
        E71T10 = getT(Z, 44, dec) .+ 273.15
        E71T11 = getT(Z, 45, dec) .+ 273.15
        E71T12 = getT(Z, 46, dec) .+ 273.15
        
        E71Tf = E71T3
        E71Ts2 = (E71T9 .+ E71T11) ./ 2
        E71Ts3 = (E71T12 .+ E71T10) ./ 2
        scatter(E71t, E71Tf, ylim=(200, 1200))
        plot!(E71t, E71T8)
        plot!(E71t, E71Ts2)
        plot!(E71t, E71Ts3)
    end

    begin
        #Exp 72 - T3, T8
        #E1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0072_231129_104140.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:C3217"]
        Z = CSV.File(path * "Data_FPT0072_231129_104140.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E72t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E72T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E72T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E72T8 = getT(Z, 42, dec) .+ 273.15
        E72T9 = getT(Z, 43, dec) .+ 273.15
        E72T10 = getT(Z, 44, dec) .+ 273.15
        E72T11 = getT(Z, 45, dec) .+ 273.15
        E72T12 = getT(Z, 46, dec) .+ 273.15
        
        E72Tf = E72T3
        E72Ts2 = (E72T9 .+ E72T11) ./ 2
        E72Ts3 = (E72T12 .+ E72T10) ./ 2
        scatter(E72t, E72Tf, ylim=(200, 1200))
        plot!(E72t, E72T8)
        plot!(E72t, E72Ts2)
        plot!(E72t, E72Ts3)
    end

    begin
        #Exp 73 - T3, T8
        #F1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0073_231129_132744.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:C4575"]
        Z = CSV.File(path * "Data_FPT0073_231129_132744.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E73t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E73T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E73T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E73T8 = getT(Z, 42, dec) .+ 273.15
        E73T9 = getT(Z, 43, dec) .+ 273.15
        E73T10 = getT(Z, 44, dec) .+ 273.15
        E73T11 = getT(Z, 45, dec) .+ 273.15
        E73T12 = getT(Z, 46, dec) .+ 273.15
        
        E73Tf = E73T3
        E73Ts2 = (E73T9 .+ E73T11) ./ 2
        E73Ts3 = (E73T12 .+ E73T10) ./ 2
        scatter(E73t, E73Tf, ylim=(200, 1200))
        plot!(E73t, E73T8)
        plot!(E73t, E73Ts2)
        plot!(E73t, E73Ts3)

    end

    begin
        #Exp 74 - T3, T8
        #G = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0074_231130_123228.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:C6018"]
        Z = CSV.File(path * "Data_FPT0074_231130_123228.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E74t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E74T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E74T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E74T8 = getT(Z, 42, dec) .+ 273.15
        E74T9 = getT(Z, 43, dec) .+ 273.15
        E74T10 = getT(Z, 44, dec) .+ 273.15
        E74T11 = getT(Z, 45, dec) .+ 273.15
        E74T12 = getT(Z, 46, dec) .+ 273.15
        
        E74Tf = E74T3
        E74Ts2 = (E74T9 .+ E74T11) ./ 2
        E74Ts3 = (E74T12 .+ E74T10) ./ 2
        scatter(E74t, E74Tf, ylim=(200, 1200))
        plot!(E74t, E74T8)
        plot!(E74t, E74Ts2)
        plot!(E74t, E74Ts3)
    end

    begin
        #Exp 75 - T3, T8
        #H = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0075_231201_162138.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:C6354"]
        Z = CSV.File(path * "Data_FPT0075_231201_162138.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E75t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E75T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E75T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E75T8 = getT(Z, 42, dec) .+ 273.15
        E75T9 = getT(Z, 43, dec) .+ 273.15
        E75T10 = getT(Z, 44, dec) .+ 273.15
        E75T11 = getT(Z, 45, dec) .+ 273.15
        E75T12 = getT(Z, 46, dec) .+ 273.15
        
        E75Tf = E75T3
        E75Ts2 = (E75T9 .+ E75T11) ./ 2
        E75Ts3 = (E75T12 .+ E75T10) ./ 2
        scatter(E75t, E75Tf, ylim=(200, 1200))
        plot!(E75t, E75T8)
        plot!(E75t, E75Ts2)
        plot!(E75t, E75Ts3)

    end

    begin
        #Exp 76 - T3, T8
        #I = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0076_231203_120521.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:C7147"]
        Z = CSV.File(path * "Data_FPT0076_231203_120521.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E76t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E76T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E76T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E76T8 = getT(Z, 42, dec) .+ 273.15
        E76T9 = getT(Z, 43, dec) .+ 273.15
        E76T10 = getT(Z, 44, dec) .+ 273.15
        E76T11 = getT(Z, 45, dec) .+ 273.15
        E76T12 = getT(Z, 46, dec) .+ 273.15
        
        E76Tf = E76T3
        E76Ts2 = (E76T9 .+ E76T11) ./ 2
        E76Ts3 = (E76T12 .+ E76T10) ./ 2
        scatter(E76t, E76Tf, ylim=(200, 1200))
        plot!(E76t, E76T8)
        plot!(E76t, E76Ts2)
        plot!(E76t, E76Ts3)

    end

    begin
        #Exp 77 - T3, T8
        #J = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0077_231203_161315.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:C3044"]
        Z = CSV.File(path * "Data_FPT0077_231203_161315.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E77t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E77T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E77T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E77T8 = getT(Z, 42, dec) .+ 273.15
        E77T9 = getT(Z, 43, dec) .+ 273.15
        E77T10 = getT(Z, 44, dec) .+ 273.15
        E77T11 = getT(Z, 45, dec) .+ 273.15
        E77T12 = getT(Z, 46, dec) .+ 273.15
        
        E77Tf = E77T3
        E77Ts2 = (E77T9 .+ E77T11) ./ 2
        E77Ts3 = (E77T12 .+ E77T10) ./ 2
        scatter(E77t, E77Tf, ylim=(200, 1200))
        plot!(E77t, E77T8)
        plot!(E77t, E77Ts2)
        plot!(E77t, E77Ts3)
    end

    begin
        #Exp 78 - T3, T8
        #K = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0078_231204_132252.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:C5384"]
        Z = CSV.File(path * "Data_FPT0078_231204_132252.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E78t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E78T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E78T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E78T8 = getT(Z, 42, dec) .+ 273.15
        E78T9 = getT(Z, 43, dec) .+ 273.15
        E78T10 = getT(Z, 44, dec) .+ 273.15
        E78T11 = getT(Z, 45, dec) .+ 273.15
        E78T12 = getT(Z, 46, dec) .+ 273.15
        
        E78Tf = E78T3
        E78Ts2 = (E78T9 .+ E78T11) ./ 2
        E78Ts3 = (E78T12 .+ E78T10) ./ 2
        scatter(E78t, E78Tf, ylim=(200, 1200))
        plot!(E78t, E78T8)
        plot!(E78t, E78Ts2)
        plot!(E78t, E78Ts3)

    begin
        #Exp 79 - T3, T8
        #L0 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0079_231204_172244.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:C5233"]
        Z = CSV.File(path * "Data_FPT0079_231204_172244.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E79t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E79T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E79T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E79T8 = getT(Z, 42, dec) .+ 273.15
        E79T9 = getT(Z, 43, dec) .+ 273.15
        E79T10 = getT(Z, 44, dec) .+ 273.15
        E79T11 = getT(Z, 45, dec) .+ 273.15
        E79T12 = getT(Z, 46, dec) .+ 273.15
        
        E79Tf = E79T3
        E79Ts2 = (E79T9 .+ E79T11) ./ 2
        E79Ts3 = (E79T12 .+ E79T10) ./ 2
        scatter(E79t, E79Tf, ylim=(200, 1200))
        plot!(E79t, E79T8)
        plot!(E79t, E79Ts2)
        plot!(E79t, E79Ts3)
    end

    begin
        #Exp 80 - T3, T8
        #M = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0080_231205_095122.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:C5814"]
        Z = CSV.File(path * "Data_FPT0080_231205_095122.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E80t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E80T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E80T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E80T8 = getT(Z, 42, dec) .+ 273.15
        E80T9 = getT(Z, 43, dec) .+ 273.15
        E80T10 = getT(Z, 44, dec) .+ 273.15
        E80T11 = getT(Z, 45, dec) .+ 273.15
        E80T12 = getT(Z, 46, dec) .+ 273.15
        
        E80Tf = E80T3
        E80Ts2 = (E80T9 .+ E80T11) ./ 2
        E80Ts3 = (E80T12 .+ E80T10) ./ 2
        scatter(E80t, E80Tf, ylim=(200, 1200))
        plot!(E80t, E80T8)
        plot!(E80t, E80Ts2)
        plot!(E80t, E80Ts3)
    end

    begin
        #Exp 81 - T3, T8
        #N = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0081_231205_135354.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:C5989"]
        Z = CSV.File(path * "Data_FPT0081_231205_135354.csv"; skipto=2, delim=",", header=false) |> DataFrame
        E81t = getT(Z, 2, dec) #float.(Z[:, 2]) #xz_data
        E81T3 = getT(Z, 37, dec) .+ 273.15 #y1z_data
        E81T5 = getT(Z, 39, dec) .+ 273.15 #y1z_data
        E81T8 = getT(Z, 42, dec) .+ 273.15
        E81T9 = getT(Z, 43, dec) .+ 273.15
        E81T10 = getT(Z, 44, dec) .+ 273.15
        E81T11 = getT(Z, 45, dec) .+ 273.15
        E81T12 = getT(Z, 46, dec) .+ 273.15
        
        E81Tf = E81T3
        E81Ts2 = (E81T9 .+ E81T11) ./ 2
        E81Ts3 = (E81T12 .+ E81T10) ./ 2
        scatter(E81t, E81Tf, ylim=(200, 1200))
        plot!(E81t, E81T8)
        plot!(E81t, E81Ts2)
        plot!(E81t, E81Ts3)

    end


#Defining simulation conditions


    IoA = 456000.0 * 1.
    IoB = 304000.0 * 1.
    IoC = 256000.0 * 1.

    # IoA = 530000.0
    # IoB = 388000.0
    # IoC = 230000.0

    # condition_E67 = Dict(Io => 442320.0, qlpm => 15.27) #1
    # condition_E68 = Dict(Io => 442320.0, qlpm => 12.50) #2
    # condition_E69 = Dict(Io => 442320.0, qlpm => 10.50) #3
    # condition_E70 = Dict(Io => 442320.0, qlpm => 9.10) #4
    # condition_E71 = Dict(Io => 442320.0, qlpm => 7.12) #5 
    # condition_E72 = Dict(Io => 400000.0, qlpm => 18.34) #6
    # condition_E73 = Dict(Io => 450000.0, qlpm => 13.16) #7
    # condition_E74 = Dict(Io => 550000.0, qlpm => 9.03) #8
    # condition_E75 = Dict(Io => 450000.0, qlpm => 6.95) #9
    # condition_E76 = Dict(Io => 470000.0, qlpm => 4.53) #10
    # condition_E77 = Dict(Io => 240000.0, qlpm => 13.85) #11
    # condition_E78 = Dict(Io => 240000.0, qlpm => 10.02) #12
    # condition_E79 = Dict(Io => 300000.0, qlpm => 8.04) #13
    # condition_E80 = Dict(Io => 248320.0, qlpm => 6.62) #14
    # condition_E81 = Dict(Io => 300000.0, qlpm => 4.53) #15

    # condition_E67 = Dict(Io => 442320.0, qlpm => 15.27)
    # condition_E68 = Dict(Io => 442320.0, qlpm => 12.50)
    # condition_E69 = Dict(Io => 442320.0, qlpm => 10.50)
    # condition_E70 = Dict(Io => 442320.0, qlpm => 9.10)
    # condition_E71 = Dict(Io => 442320.0, qlpm => 7.12)
    # condition_E72 = Dict(Io => 371792.0, qlpm => 18.34)
    # condition_E73 = Dict(Io => 371792.0, qlpm => 13.16)
    # condition_E74 = Dict(Io => 371792.0, qlpm => 9.03)
    # condition_E75 = Dict(Io => 371792.0, qlpm => 6.95)
    # condition_E76 = Dict(Io => 371792.0, qlpm => 4.53)
    # condition_E77 = Dict(Io => 313088.0, qlpm => 13.85)
    # condition_E78 = Dict(Io => 313088.0, qlpm => 10.02)
    # condition_E79 = Dict(Io => 313088.0, qlpm => 8.04)
    # condition_E80 = Dict(Io => 313088.0, qlpm => 6.62)
    # condition_E81 = Dict(Io => 313088.0, qlpm => 4.53)

    condition_E67 = Dict(Io => IoA, qlpm => 15.27, Tinit => E67Tf[1])
    condition_E68 = Dict(Io => IoA, qlpm => 12.50, Tinit => E68Tf[1])
    condition_E69 = Dict(Io => IoA, qlpm => 10.50, Tinit => E69Tf[1])
    condition_E70 = Dict(Io => IoA, qlpm => 9.10, Tinit => E70Tf[1])
    condition_E71 = Dict(Io => IoA, qlpm => 7.12, Tinit => E71Tf[1])
    condition_E72 = Dict(Io => IoB, qlpm => 18.34, Tinit => E72Tf[1])
    condition_E73 = Dict(Io => IoB, qlpm => 13.16, Tinit => E73Tf[1])
    condition_E74 = Dict(Io => IoB, qlpm => 9.03, Tinit => E74Tf[1])
    condition_E75 = Dict(Io => IoB, qlpm => 6.95, Tinit => E75Tf[1])
    condition_E76 = Dict(Io => IoB, qlpm => 4.53, Tinit => E76Tf[1])
    condition_E77 = Dict(Io => IoC, qlpm => 13.85, Tinit => E77Tf[1])
    condition_E78 = Dict(Io => IoC, qlpm => 10.02, Tinit => E78Tf[1])
    condition_E79 = Dict(Io => IoC, qlpm => 8.04, Tinit => E79Tf[1])
    condition_E80 = Dict(Io => IoC, qlpm => 6.62, Tinit => E80Tf[1])
    condition_E81 = Dict(Io => IoC, qlpm => 4.53, Tinit => E81Tf[1])


    simulation_conditions = Dict("E67" => condition_E67, "E68" => condition_E68,
        "E69" => condition_E69, "E70" => condition_E70,
        "E71" => condition_E71, "E72" => condition_E72,
        "E73" => condition_E73, "E74" => condition_E74,
        "E75" => condition_E75, "E76" => condition_E76,
        "E77" => condition_E77, "E78" => condition_E78,
        "E79" => condition_E79, "E80" => condition_E80,
        "E81" => condition_E81)
    #Defining measurement data
    # measurements = DataFrame(
    #     simulation_id=repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=6, outer=1),
    #     obs_id=repeat(["_T8", "_T9", "_T10", "_T11", "_T12", "_T3"], inner=1, outer=15),
    #     time=repeat([3929, 5362, 5363, 6702, 7084, 3214, 4572, 6015, 6351, 7144, 3041, 5381, 5230, 5811, 5986], inner=6, outer=1),
    #     temperatures=[965.407, 975.144, 825.592, 880.867, 1004.165, 763.859,
    #         1031.574, 1023.115, 850.099, 898.754, 1050.691, 773.207,
    #         1070.803, 1045.898, 852.837, 896.788, 1072.727, 769.76,
    #         1167.978, 1093.849, 871.496, 912.173, 1120.42, 779.53,
    #         1210.945, 1095.322, 847.476, 882.823, 1120.417, 753.56,
    #         742.125, 778.592, 684.246, 736.626, 807.125, 652.955,
    #         844.257, 870.26, 747.444, 791.958, 898.081, 694.626,
    #         962.113, 938.106, 767.803, 804.082, 965.702, 697.672,
    #         1015.081, 954.214, 757.393, 788.678, 979.795, 681.066,
    #         1069.567, 947.372, 726.308, 751.159, 970.498, 647.019,
    #         569.248, 604.984, 543.727, 574.299, 627.16, 525.356,
    #         634.731, 664.296, 583.704, 612.092, 686.936, 554.455,
    #         677.817, 694.156, 595.766, 622.325, 716.314, 560.033,
    #         711.537, 713.686, 601.77, 626.485, 735.984, 561.254,
    #         763.299, 729.461, 597.766, 618.975, 751.15, 550.499])
    # measurements = DataFrame( #T3 only
    #             simulation_id=repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=1, outer=1),
    #             obs_id=repeat(["_T3"], inner=1, outer=15),
    #             time=repeat([3929, 5362, 5363, 6702, 7084, 3214, 4572, 6015, 6351, 7144, 3041, 5381, 5230, 5811, 5986], inner=1, outer=1),
    #             temperatures=[763.859, 773.207, 769.76, 779.53, 753.56, 652.955, 694.626, 697.672, 681.066,
    #                 647.019, 525.356, 554.455, 560.033, 561.254, 550.499])
    # measurements = DataFrame(
    #     simulation_id = repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=6, outer=1),
    #     obs_id=repeat(["_T8", "_T9", "_T10", "_T11", "_T12", "_T3"], inner=1, outer=15),
    #     time= repeat([E67t, E68t, E69t, E70t, E71t, E72t, E73t, E74t, E75t, E76t, E77t, E78t, E79t, E80t, E81t], inner=6, outer=1),
    #         temperatures= [E67T8, E67T9, E67T10, E67T11, E67T12, E67Tf,
    #         E68T8, E68T9, E68T10, E68T11, E68T12, E68Tf,
    #         E69T8, E69T9, E69T10, E69T11, E69T12, E69Tf,
    #         E70T8, E70T9, E70T10, E70T11, E70T12, E70Tf,
    #         E71T8, E71T9, E71T10, E71T11, E71T12, E71Tf,
    #         E72T8, E72T9, E72T10, E72T11, E72T12, E72Tf,
    #         E73T8, E73T9, E73T10, E73T11, E73T12, E73Tf,
    #         E74T8, E74T9, E74T10, E74T11, E74T12, E74Tf,
    #         E75T8, E75T9, E75T10, E75T11, E75T12, E75Tf,
    #         E76T8, E76T9, E76T10, E76T11, E76T12, E76Tf,
    #         E77T8, E77T9, E77T10, E77T11, E77T12, E77Tf,
    #         E78T8, E78T9, E78T10, E78T11, E78T12, E78Tf,
    #         E79T8, E79T9, E79T10, E79T11, E79T12, E79Tf,
    #         E80T8, E80T9, E80T10, E80T11, E80T12, E80Tf,
    #         E81T8, E81T9, E81T10, E81T11, E81T12, E81Tf])   
    measurements = DataFrame(
        simulation_id=repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], 
            inner=1, outer=1),
        obs_id=repeat(["_T3"], inner=1, outer=15),
        time=repeat([E67t, E68t, E69t, E70t, E71t, E72t, E73t, E74t, E75t, E76t, E77t, E78t, E79t, E80t, E81t], inner=1, outer=1),
        temperatures=[E67Tf, E68Tf, E69Tf, E70Tf,
            E71Tf, E72Tf, E73Tf, E74Tf, E75Tf, E76Tf,
            E77Tf, E78Tf, E79Tf, E80Tf, E81Tf]
            )
end