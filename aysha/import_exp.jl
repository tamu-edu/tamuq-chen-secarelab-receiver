using Plots, XLSX, DataFrames, CSV, Statistics, Dates, LinearAlgebra, DelimitedFiles

dec = 100
function decimate!(x, step)
    x = [x[i] for i in 1:step:length(x)]
end
begin
    #Exp. data to extract temp.
    begin
        #Exp 67 - T3, T8
        Z = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0067_231125_161757.xlsx")["Sheet 1 - Data_FPT0067_231125_1"]["A3:C3932"]
        E67t = float.(Z[:, 1]) #xz_data
        E67Tf = Z[:, 2] .+ 273.15 #y1z_data
        E67T8 = Z[:, 3] .+ 273.15
        #decimate experimental data E67 for T3
        E67t = decimate!(E67t, dec)
        E67Tf = decimate!(E67Tf, dec)
        scatter(E67t, E67Tf)

        #Exp 67 - T9, T10, T11, T12
        Z1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0067_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0067_T9"]["A3:E3932"]
        E67T9 = Z1[:, 2] .+ 273.15
        E67T10 = Z1[:, 3] .+ 273.15
        E67T11 = Z1[:, 4] .+ 273.15
        E67T12 = Z1[:, 5] .+ 273.15
        #decimate experimental data E67 for T8
        E67T8 = decimate!(E67T8, dec)
        scatter(E67t, E67T8)
        #decimate experimental data E67 for T9
        E67T9 = decimate!(E67T9, dec)
        scatter(E67t, E67T9)
        #decimate experimental data E67 for T10
        E67T10 = decimate!(E67T10, dec)
        scatter(E67t, E67T10)
        #decimate experimental data E67 for T11
        E67T11 = decimate!(E67T11, dec)
        scatter(E67t, E67T11)
        #decimate experimental data E67 for T12
        E67T12 = decimate!(E67T12, dec)
        scatter(E67t, E67T12)
    end

    begin
        #Exp 68 - T3, T8
        A1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0068_231126_115725.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:C5365"]
        E68t = float.(A1[:, 1]) #xa1_data
        E68Tf = A1[:, 2] .+ 273.15 #y1a1_data
        E68T8 = A1[:, 3] .+ 273.15
        #decimate experimental data E68 for T3
        E68t = decimate!(E68t, dec)
        E68Tf = decimate!(E68Tf, dec)
        scatter(E68t, E68Tf)
        #decimate experimental data E68 for T8
        E68T8 = decimate!(E68T8, dec)
        scatter(E68t, E68T8)


        #Exp 68 - T9, T10, T11, T12
        A11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0068_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:E5365"]
        E68T9 = A11[:, 2] .+ 273.15
        E68T10 = A11[:, 3] .+ 273.15
        E68T11 = A11[:, 4] .+ 273.15
        E68T12 = A11[:, 5] .+ 273.15
        #decimate experimental data E68 for T9
        E68T9 = decimate!(E68T9, dec)
        scatter(E68t, E68T9)
        #decimate experimental data E68 for T10
        E68T10 = decimate!(E68T10, dec)
        scatter(E68t, E68T10)
        #decimate experimental data E68 for T11
        E68T11 = decimate!(E68T11, dec)
        scatter(E68t, E68T11)
        #decimate experimental data E68 for T12
        E68T12 = decimate!(E68T12, dec)
        scatter(E68t, E68T12)
    end

    begin
        #Exp 69 - T3, T8
        B1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0069_231126_140153.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:C5366"]
        E69t = float.(B1[:, 1])
        E69Tf = B1[:, 2] .+ 273.15
        E69T8 = B1[:, 3] .+ 273.15
        #decimate experimental data E69 for T3
        E69t = decimate!(E69t, dec)
        E69Tf = decimate!(E69Tf, dec)
        scatter(E69t, E69Tf)
        #decimate experimental data E69 for T8
        E69T8 = decimate!(E69T8, dec)
        scatter(E69t, E69T8)

        #Exp 69 - T9, T10, T11, T12
        B11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0069_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:E5366"]
        E69T9 = B11[:, 2] .+ 273.15
        E69T10 = B11[:, 3] .+ 273.15
        E69T11 = B11[:, 4] .+ 273.15
        E69T12 = B11[:, 5] .+ 273.15
        #decimate experimental data E69 for T9
        E69T9 = decimate!(E69T9, dec)
        scatter(E69t, E69T9)
        #decimate experimental data E69 for T10
        E69T10 = decimate!(E69T10, dec)
        scatter(E69t, E69T10)
        #decimate experimental data E69 for T11
        E69T11 = decimate!(E69T11, dec)
        scatter(E69t, E69T11)
        #decimate experimental data E69 for T12
        E69T12 = decimate!(E69T12, dec)
        scatter(E69t, E69T12)
    end


    begin
        #Exp 70 - T3, T8
        C1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0070_231127_090339.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:C6705"]
        E70t = float.(C1[:, 1]) #xc1_data 
        E70Tf = C1[:, 2] .+ 273.15 #y1c1_data
        E70T8 = C1[:, 3] .+ 273.15
        #decimate experimental data E70 for T3
        E70t = decimate!(E70t, dec)
        E70Tf = decimate!(E70Tf, dec)
        scatter(E70t, E70Tf)
        #decimate experimental data E70 for T8
        E70T8 = decimate!(E70T8, dec)
        scatter(E70t, E70T8)

        #Exp 70 - T9, T10, T11, T12
        C11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0070_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:E6705"]
        E70T9 = C11[:, 2] .+ 273.15
        E70T10 = C11[:, 3] .+ 273.15
        E70T11 = C11[:, 4] .+ 273.15
        E70T12 = C11[:, 5] .+ 273.15
        #decimate experimental data E70 for T9
        E70T9 = decimate!(E70T9, dec)
        scatter(E70t, E70T9)
        #decimate experimental data E70 for T10
        E70T10 = decimate!(E70T10, dec)
        scatter(E70t, E70T10)
        #decimate experimental data E70 for T11
        E70T11 = decimate!(E70T11, dec)
        scatter(E70t, E70T11)
        #decimate experimental data E70 for T12
        E70T12 = decimate!(E70T12, dec)
        scatter(E70t, E70T12)
    end
    begin
        #Exp 71 - T3, T8
        D1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_231128_102707.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:C7087"]
        E71t = float.(D1[:, 1]) #xd1_data  
        E71Tf = D1[:, 2] .+ 273.15 #y1d1_data
        E71T8 = D1[:, 3] .+ 273.15
        #decimate experimental data E71 for T3
        E71t = decimate!(E71t, dec)
        E71Tf = decimate!(E71Tf, dec)
        scatter(E71t, E71Tf)
        #decimate experimental data E71 for T8
        E71T8 = decimate!(E71T8, dec)
        scatter(E71t, E71T8)
        #Exp 71 - T9, T10, T11, T12
        D11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:E7087"]
        E71T9 = D11[:, 2] .+ 273.15
        E71T10 = D11[:, 3] .+ 273.15
        E71T11 = D11[:, 4] .+ 273.15
        E71T12 = D11[:, 5] .+ 273.15
        #decimate experimental data E71 for T9
        E71T9 = decimate!(E71T9, dec)
        scatter(E71t, E71T9)
        #decimate experimental data E71 for T10
        E71T10 = decimate!(E71T10, dec)
        scatter(E71t, E71T10)
        #decimate experimental data E71 for T11
        E71T11 = decimate!(E71T11, dec)
        scatter(E71t, E71T11)
        #decimate experimental data E71 for T12
        E71T12 = decimate!(E71T12, dec)
        scatter(E71t, E71T12)
    end

    begin
        #Exp 72 - T3, T8
        E1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0072_231129_104140.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:C3217"]
        E72t = float.(E1[:, 1]) #xe1_data
        E72Tf = E1[:, 2] .+ 273.15 #y1e1_data
        E72T8 = E1[:, 3] .+ 273.15
        #decimate experimental data E72 for T3
        E72t = decimate!(E72t, dec)
        E72Tf = decimate!(E72Tf, dec)
        scatter(E72t, E72Tf)
        #decimate experimental data E72 for T8
        E72T8 = decimate!(E72T8, dec)
        scatter(E72t, E72T8)

        #Exp 72 - T9, T10, T11, T12
        E11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0072_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:E3217"]
        E72T9 = E11[:, 2] .+ 273.15
        E72T10 = E11[:, 3] .+ 273.15
        E72T11 = E11[:, 4] .+ 273.15
        E72T12 = E11[:, 5] .+ 273.15
        #decimate experimental data E72 for T9
        E72T9 = decimate!(E72T9, dec)
        scatter(E72t, E72T9)
        #decimate experimental data E72 for T10
        E72T10 = decimate!(E72T10, dec)
        scatter(E72t, E72T10)
        #decimate experimental data E72 for T11
        E72T11 = decimate!(E72T11, dec)
        scatter(E72t, E72T11)
        #decimate experimental data E72 for T12
        E72T12 = decimate!(E72T12, dec)
        scatter(E72t, E72T12)
    end

    begin
        #Exp 73 - T3, T8
        F1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0073_231129_132744.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:C4575"]
        E73t = float.(F1[:, 1])#xf1_data
        E73Tf = F1[:, 2] .+ 273.15 #y1f1_data
        E73T8 = F1[:, 3] .+ 273.15
        #decimate experimental data E73 for T3
        E73t = decimate!(E73t, dec)
        E73Tf = decimate!(E73Tf, dec)
        scatter(E73t, E73Tf)
        #decimate experimental data E73 for T8
        E73T8 = decimate!(E73T8, dec)
        scatter(E73t, E73T8)

        #Exp 73 - T9, T10, T11, T12
        F11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0073_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:E4575"]
        E73T9 = F11[:, 2] .+ 273.15
        E73T10 = F11[:, 3] .+ 273.15
        E73T11 = F11[:, 4] .+ 273.15
        E73T12 = F11[:, 5] .+ 273.15
        #decimate experimental data E73 for T9
        E73T9 = decimate!(E73T9, dec)
        scatter(E73t, E73T9)
        #decimate experimental data E73 for T10
        E73T10 = decimate!(E73T10, dec)
        scatter(E73t, E73T10)
        #decimate experimental data E73 for T11
        E73T11 = decimate!(E73T11, dec)
        scatter(E73t, E73T11)
        #decimate experimental data E73 for T12
        E73T12 = decimate!(E73T12, dec)
        scatter(E73t, E73T12)

    end

    begin
        #Exp 74 - T3, T8
        G = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0074_231130_123228.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:C6018"]
        E74t = float.(G[:, 1]) #xg_data 
        E74Tf = G[:, 2] .+ 273.15 #y1g_data 
        E74T8 = G[:, 3] .+ 273.15
        #decimate experimental data E74 for T3
        E74t = decimate!(E74t, dec)
        E74Tf = decimate!(E74Tf, dec)
        scatter(E74t, E74Tf)
        #decimate experimental data E74 for T8
        E74T8 = decimate!(E74T8, dec)
        scatter(E74t, E74T8)

        #Exp 74 - T9, T10, T11, T12
        G1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0074_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:E6018"]
        E74T9 = G1[:, 2] .+ 273.15
        E74T10 = G1[:, 3] .+ 273.15
        E74T11 = G1[:, 4] .+ 273.15
        E74T12 = G1[:, 5] .+ 273.15
        #decimate experimental data E74 for T9
        E74T9 = decimate!(E74T9, dec)
        scatter(E74t, E74T9)
        #decimate experimental data E74 for T10
        E74T10 = decimate!(E74T10, dec)
        scatter(E74t, E74T10)
        #decimate experimental data E74 for T11
        E74T11 = decimate!(E74T11, dec)
        scatter(E74t, E74T11)
        #decimate experimental data E74 for T12
        E74T12 = decimate!(E74T12, dec)
        scatter(E74t, E74T12)
    end

    begin
        #Exp 75 - T3, T8
        H = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0075_231201_162138.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:C6354"]
        E75t = float.(H[:, 1]) #xh_data
        E75Tf = H[:, 2] .+ 273.15 #y1h_data
        E75T8 = H[:, 3] .+ 273.15
        #decimate experimental data E75 for T3
        E75t = decimate!(E75t, dec)
        E75Tf = decimate!(E75Tf, dec)
        scatter(E75t, E75Tf)
        #decimate experimental data E75 for T8
        E75T8 = decimate!(E75T8, dec)
        scatter(E75t, E75T8)


        #Exp 75 - T9, T10, T11, T12
        H1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0075_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:E6354"]
        E75T9 = H1[:, 2] .+ 273.15
        E75T10 = H1[:, 3] .+ 273.15
        E75T11 = H1[:, 4] .+ 273.15
        E75T12 = H1[:, 5] .+ 273.15
        #decimate experimental data E75 for T9
        E75T9 = decimate!(E75T9, dec)
        scatter(E75t, E75T9)
        #decimate experimental data E75 for T10
        E75T10 = decimate!(E75T10, dec)
        scatter(E75t, E75T10)
        #decimate experimental data E75 for T11
        E75T11 = decimate!(E75T11, dec)
        scatter(E75t, E75T11)
        #decimate experimental data E75 for T12
        E75T12 = decimate!(E75T12, dec)
        scatter(E75t, E75T12)

    end

    begin
        #Exp 76 - T3, T8
        I = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0076_231203_120521.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:C7147"]
        E76t = float.(I[:, 1]) #xi_data
        E76Tf = I[:, 2] .+ 273.15 #y1i_data
        E76T8 = I[:, 3] .+ 273.15
        #decimate experimental data E76 for T3
        E76t = decimate!(E76t, dec)
        E76Tf = decimate!(E76Tf, dec)
        scatter(E76t, E76Tf)
        #decimate experimental data E76 for T8
        E76T8 = decimate!(E76T8, dec)
        scatter(E76t, E76T8)

        #Exp 76 - T9, T10, T11, T12
        I1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0076_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:E7147"]
        E76T9 = I1[:, 2] .+ 273.15
        E76T10 = I1[:, 3] .+ 273.15
        E76T11 = I1[:, 4] .+ 273.15
        E76T12 = I1[:, 5] .+ 273.15
        #decimate experimental data E76 for T9
        E76T9 = decimate!(E76T9, dec)
        scatter(E76t, E76T9)
        #decimate experimental data E76 for T10
        E76T10 = decimate!(E76T10, dec)
        scatter(E76t, E76T10)
        #decimate experimental data E76 for T11
        E76T11 = decimate!(E76T11, dec)
        scatter(E76t, E76T11)
        #decimate experimental data E76 for T12
        E76T12 = decimate!(E76T12, dec)
        scatter(E76t, E76T12)

    end

    begin
        #Exp 77 - T3, T8
        J = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0077_231203_161315.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:C3044"]
        E77t = float.(J[:, 1]) #xj_data
        E77Tf = J[:, 2] .+ 273.15 #y1j_data
        E77T8 = J[:, 3] .+ 273.15
        #decimate experimental data E77 for T3
        E77t = decimate!(E77t, dec)
        E77Tf = decimate!(E77Tf, dec)
        scatter(E77t, E77Tf)
        #decimate experimental data E77 for T8
        E77T8 = decimate!(E77T8, dec)
        scatter(E77t, E77T8)

        #Exp 77 - T9, T10, T11, T12
        J1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0077_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:E3044"]
        E77T9 = J1[:, 2] .+ 273.15
        E77T10 = J1[:, 3] .+ 273.15
        E77T11 = J1[:, 4] .+ 273.15
        E77T12 = J1[:, 5] .+ 273.15
        #decimate experimental data E77 for T9
        E77T9 = decimate!(E77T9, dec)
        scatter(E77t, E77T9)
        #decimate experimental data E77 for T10
        E77T10 = decimate!(E77T10, dec)
        scatter(E77t, E77T10)
        #decimate experimental data E77 for T11
        E77T11 = decimate!(E77T11, dec)
        scatter(E77t, E77T11)
        #decimate experimental data E77 for T12
        E77T12 = decimate!(E77T12, dec)
        scatter(E77t, E77T12)
    end

    begin
        #Exp 78 - T3, T8
        K = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0078_231204_132252.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:C5384"]
        E78t = float.(K[:, 1]) #xk_data
        E78Tf = K[:, 2] .+ 273.15 #y1k_data
        E78T8 = K[:, 3] .+ 273.15
        #decimate experimental data E78 for T3
        E78t = decimate!(E78t, dec)
        E78Tf = decimate!(E78Tf, dec)
        scatter(E78t, E78Tf)
        #decimate experimental data E78 for T8
        E78T8 = decimate!(E78T8, dec)
        scatter(E78t, E78T8)

        #Exp 78 - T9, T10, T11, T12
        K1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0078_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:E5384"]
        E78T9 = K1[:, 2] .+ 273.15
        E78T10 = K1[:, 3] .+ 273.15
        E78T11 = K1[:, 4] .+ 273.15
        E78T12 = K1[:, 5] .+ 273.15
        #decimate experimental data E78 for T9
        E78T9 = decimate!(E78T9, dec)
        scatter(E78t, E78T9)
        #decimate experimental data E78 for T10
        E78T10 = decimate!(E78T10, dec)
        scatter(E78t, E78T10)
        #decimate experimental data E78 for T11
        E78T11 = decimate!(E78T11, dec)
        scatter(E78t, E78T11)
        #decimate experimental data E78 for T12
        E78T12 = decimate!(E78T12, dec)
        scatter(E78t, E78T12)
    end

    begin
        #Exp 79 - T3, T8
        L0 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0079_231204_172244.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:C5233"]
        E79t = float.(L0[:, 1]) #xl_data 
        E79Tf = L0[:, 2] .+ 273.15 #y1l_data
        E79T8 = L0[:, 3] .+ 273.15
        #decimate experimental data E79 for T3
        E79t = decimate!(E79t, dec)
        E79Tf = decimate!(E79Tf, dec)
        scatter(E79t, E79Tf)
        #decimate experimental data E79 for T8
        E79T8 = decimate!(E79T8, dec)
        scatter(E79t, E79T8)

        #Exp 79 - T9, T10, T11, T12
        L1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0079_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:E5233"]
        E79T9 = L1[:, 2] .+ 273.15
        E79T10 = L1[:, 3] .+ 273.15
        E79T11 = L1[:, 4] .+ 273.15
        E79T12 = L1[:, 5] .+ 273.15
        #decimate experimental data E79 for T9
        E79T9 = decimate!(E79T9, dec)
        scatter(E79t, E79T9)
        #decimate experimental data E79 for T10
        E79T10 = decimate!(E79T10, dec)
        scatter(E79t, E79T10)
        #decimate experimental data E79 for T11
        E79T11 = decimate!(E79T11, dec)
        scatter(E79t, E79T11)
        #decimate experimental data E79 for T12
        E79T12 = decimate!(E79T12, dec)
        scatter(E79t, E79T12)
    end

    begin
        #Exp 80 - T3, T8
        M = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0080_231205_095122.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:C5814"]
        E80t = float.(M[:, 1]) # xm_data
        E80Tf = M[:, 2] .+ 273.15 #y1m_data
        E80T8 = M[:, 3] .+ 273.15
        #decimate experimental data E80 for T3
        E80t = decimate!(E80t, dec)
        E80Tf = decimate!(E80Tf, dec)
        scatter(E80t, E80Tf)

        #decimate experimental data E80 for T8
        E80T8 = decimate!(E80T8, dec)
        scatter(E80t, E80T8)
    end
    begin
        #Exp 80 - T9, T10, T11, T12
        M1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0080_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:E5814"]
        E80T9 = M1[:, 2] .+ 273.15
        E80T10 = M1[:, 3] .+ 273.15
        E80T11 = M1[:, 4] .+ 273.15
        E80T12 = M1[:, 5] .+ 273.15
        #decimate experimental data E80 for T9
        E80T9 = decimate!(E80T9, dec)
        scatter(E80t, E80T9)
        #decimate experimental data E80 for T10
        E80T10 = decimate!(E80T10, dec)
        scatter(E80t, E80T10)
        #decimate experimental data E80 for T11
        E80T11 = decimate!(E80T11, dec)
        scatter(E80t, E80T11)
        #decimate experimental data E80 for T12
        E80T12 = decimate!(E80T12, dec)
        scatter(E80t, E80T12)
    end

    begin
        #Exp 81 - T3, T8
        N = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0081_231205_135354.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:C5989"]
        E81t = float.(N[:, 1]) #xn_data 
        E81Tf = N[:, 2] .+ 273.15 #y1n_data 
        E81T8 = N[:, 3] .+ 273.15
        #decimate experimental data E81 for T3
        E81t = decimate!(E81t, dec)
        E81Tf = decimate!(E81Tf, dec)
        scatter(E81t, E81Tf)
        #decimate experimental data E81 for T8
        E81T8 = decimate!(E81T8, dec)
        scatter(E81t, E81T8)

        #Exp 81 - T9, T10, T11, T12
        N1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0081_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:E5989"]
        E81T9 = N1[:, 2] .+ 273.15
        E81T10 = N1[:, 3] .+ 273.15
        E81T11 = N1[:, 4] .+ 273.15
        E81T12 = N1[:, 5] .+ 273.15
        #decimate experimental data E81 for T9
        E81T9 = decimate!(E81T9, dec)
        scatter(E81t, E81T9)
        #decimate experimental data E81 for T10
        E81T10 = decimate!(E81T10, dec)
        scatter(E81t, E81T10)
        #decimate experimental data E81 for T11
        E81T11 = decimate!(E81T11, dec)
        scatter(E81t, E81T11)
        #decimate experimental data E81 for T12
        E81T12 = decimate!(E81T12, dec)
        scatter(E81t, E81T12)

    end
end

#Defining simulation conditions
begin

    IoA = 456000.0 * 1.3
    IoB = 304000.0 * 0.8
    IoC = 256000.0 * 0.5

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

    condition_E67 = Dict(Io => IoA, qlpm => 15.27)
    condition_E68 = Dict(Io => IoA, qlpm => 12.50)
    condition_E69 = Dict(Io => IoA, qlpm => 10.50)
    condition_E70 = Dict(Io => IoA, qlpm => 9.10)
    condition_E71 = Dict(Io => IoA, qlpm => 7.12)
    condition_E72 = Dict(Io => IoB, qlpm => 18.34)
    condition_E73 = Dict(Io => IoB, qlpm => 13.16)
    condition_E74 = Dict(Io => IoB, qlpm => 9.03)
    condition_E75 = Dict(Io => IoB, qlpm => 6.95)
    condition_E76 = Dict(Io => IoB, qlpm => 4.53)
    condition_E77 = Dict(Io => IoC, qlpm => 13.85)
    condition_E78 = Dict(Io => IoC, qlpm => 10.02)
    condition_E79 = Dict(Io => IoC, qlpm => 8.04)
    condition_E80 = Dict(Io => IoC, qlpm => 6.62)
    condition_E81 = Dict(Io => IoC, qlpm => 4.53)


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
        simulation_id=repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=1, outer=1),
        obs_id=repeat(["_T3"], inner=1, outer=15),
        time=repeat([E67t, E68t, E69t, E70t, E71t, E72t, E73t, E74t, E75t, E76t, E77t, E78t, E79t, E80t, E81t], inner=1, outer=1),
        temperatures=[E67Tf, E68Tf, E69Tf, E70Tf,
            E71Tf, E72Tf, E73Tf, E74Tf, E75Tf, E76Tf,
            E77Tf, E78Tf, E79Tf, E80Tf, E81Tf])
end