### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3e50e470-0251-11ef-3b14-7d58b00ae519
begin
	import Pkg
	Pkg.activate("energy")
end

# ╔═╡ 86c7df59-4f03-4730-9cbf-72d8fc7c34bd
begin
	Pkg.resolve()
	using ModelingToolkit, DifferentialEquations, Plots
	#using ModelingToolkit: t_nounits as t, D_nounits as D
	using PEtab, XLSX, Statistics, DataFrames
	using PlutoUI, StatsPlots
	using Optimization, Optim, Ipopt, OptimizationNLopt
end

# ╔═╡ 8c28c8d1-dd57-4a04-9121-cbfed404d824
html"""<style>
main {
    max-width: 70%;
    margin-left: 1%;
    margin-right: 30% !important;
}
"""

# ╔═╡ 232f4256-2095-470e-a669-216e1642c374
Pkg.status()

# ╔═╡ 3f7ecdc3-b35e-4cc7-8a37-4576cd0f02b9
plotly()

# ╔═╡ 2eb83c61-2ea3-4120-9a55-3e0b66e279bb
begin
	@parameters aCp A B C al aIo #to be fitted
	@parameters qlpm Io Tins #varying conditions
	@variables t Tf(t) Ts(t)
	D = Differential(t)
end

# ╔═╡ 2c50a54e-96bd-4f48-ae7f-33cbface0754
#Exp. data to extract temp.
begin
	starttime = 1
	path = "./aysha/SolarSimulator/EXCEL/"
    #Exp 67 - T3, T8
    Z = XLSX.readxlsx(path*"Data_FPT0067_231125_161757.xlsx")["Sheet 1 - Data_FPT0067_231125_1"]["A3:C3932"]
    E67t = Float64.(Z[starttime:end, 1])
    y1z_data = Z[starttime:end, 2] .+ 273.15
    y2z_data = Z[starttime:end, 3] .+ 273.15

    #Exp 67 - T9, T10, T11, T12
    Z1 = XLSX.readxlsx(path*"Data_FPT0067_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0067_T9"]["A3:E3932"]
    y3z_data = Z1[starttime:end, 2] .+ 273.15
    y4z_data = Z1[starttime:end, 3] .+ 273.15
    y5z_data = Z1[starttime:end, 4] .+ 273.15
    y6z_data = Z1[starttime:end, 5] .+ 273.15

	E67Tf = y1z_data
	#E67Ts = mean([y2z_data, y3z_data, y4z_data, y5z_data, y6z_data])
    E67Ts = mean([y3z_data, y4z_data, y5z_data, y6z_data])
    #Exp 68 - T3, T8
    A1 = XLSX.readxlsx(path*"Data_FPT0068_231126_115725.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:C5365"]
    E68t = Float64.(A1[starttime:end, 1])
    y1a1_data = A1[starttime:end, 2] .+ 273.15
    y2a1_data = A1[starttime:end, 3] .+ 273.15

    #Exp 68 - T9, T10, T11, T12
    A11 = XLSX.readxlsx(path*"Data_FPT0068_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:E5365"]
    y3a_data = A11[starttime:end, 2] .+ 273.15
    y4a_data = A11[starttime:end, 3] .+ 273.15
    y5a_data = A11[starttime:end, 4] .+ 273.15
    y6a_data = A11[starttime:end, 5] .+ 273.15

	E68Tf = y1a1_data
	#E68Ts = mean([y2a1_data, y3a_data, y4a_data, y5a_data, y6a_data])
	E68Ts = mean([y3a_data, y4a_data, y5a_data, y6a_data])
	
    #Exp 69 - T3, T8
    B1 = XLSX.readxlsx(path*"Data_FPT0069_231126_140153.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:C5366"]
    E69t = Float64.(B1[starttime:end, 1])
    y1b1_data = B1[starttime:end, 2] .+ 273.15
    y2b1_data = B1[starttime:end, 3] .+ 273.15

    #Exp 69 - T9, T10, T11, T12
    B11 = XLSX.readxlsx(path*"Data_FPT0069_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:E5366"]
    y3b_data = B11[starttime:end, 2] .+ 273.15
    y4b_data = B11[starttime:end, 3] .+ 273.15
    y5b_data = B11[starttime:end, 4] .+ 273.15
    y6b_data = B11[starttime:end, 5] .+ 273.15

	E69Tf = y1b1_data
	#E69Ts = mean([y2b1_data, y3b_data, y4b_data, y5b_data, y6b_data])
	E69Ts = mean([y3b_data, y4b_data, y5b_data, y6b_data])
	
    #Exp 70 - T3, T8
    C1 = XLSX.readxlsx(path*"Data_FPT0070_231127_090339.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:C6705"]
    E70t = Float64.(C1[starttime:end, 1])
    y1c1_data = C1[starttime:end, 2] .+ 273.15
    y2c1_data = C1[starttime:end, 3] .+ 273.15

    #Exp 70 - T9, T10, T11, T12

    C11 = XLSX.readxlsx(path*"Data_FPT0070_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:E6705"]
    y3c_data = C11[starttime:end, 2] .+ 273.15
    y4c_data = C11[starttime:end, 3] .+ 273.15
    y5c_data = C11[starttime:end, 4] .+ 273.15
    y6c_data = C11[starttime:end, 5] .+ 273.15

	E70Tf = y1c1_data
	#E70Ts = mean([y2c1_data, y3c_data, y4c_data, y5c_data, y6c_data])
	E70Ts = mean([y3c_data, y4c_data, y5c_data, y6c_data])
	
    #Exp 71 - T3, T8
    D1 = XLSX.readxlsx(path*"Data_FPT0071_231128_102707.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:C7087"]
    E71t = Float64.(D1[starttime:end, 1])
    y1d1_data = D1[starttime:end, 2] .+ 273.15
    y2d1_data = D1[starttime:end, 3] .+ 273.15

    #Exp 71 - T9, T10, T11, T12

    D11 = XLSX.readxlsx(path*"Data_FPT0071_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:E7087"]
    y3d_data = D11[starttime:end, 2] .+ 273.15
    y4d_data = D11[starttime:end, 3] .+ 273.15
    y5d_data = D11[starttime:end, 4] .+ 273.15
    y6d_data = D11[starttime:end, 5] .+ 273.15

	E71Tf = y1d1_data
	#E71Ts = mean([y2d1_data, y3d_data, y4d_data, y5d_data, y6d_data])
	E71Ts = mean([y3d_data, y4d_data, y5d_data, y6d_data])
	
    #Exp 72 - T3, T8
    E1 = XLSX.readxlsx(path*"Data_FPT0072_231129_104140.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:C3217"]
    E72t = Float64.(E1[starttime:end, 1])
    y1e1_data = E1[starttime:end, 2] .+ 273.15
    y2e1_data = E1[starttime:end, 3] .+ 273.15

    #Exp 72 - T9, T10, T11, T12

    E11 = XLSX.readxlsx(path*"Data_FPT0072_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:E3217"]
    y3e_data = E11[starttime:end, 2] .+ 273.15
    y4e_data = E11[starttime:end, 3] .+ 273.15
    y5e_data = E11[starttime:end, 4] .+ 273.15
    y6e_data = E11[starttime:end, 5] .+ 273.15

	E72Tf = y1e1_data
	#E72Ts = mean([y2e1_data, y3e_data, y4e_data, y5e_data, y6e_data])
	E72Ts = mean([y3e_data, y4e_data, y5e_data, y6e_data])
	
    #Exp 73 - T3, T8
    F1 = XLSX.readxlsx(path*"Data_FPT0073_231129_132744.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:C4575"]
    E73t = Float64.(F1[starttime:end, 1])
    y1f1_data = F1[starttime:end, 2] .+ 273.15
    y2f1_data = F1[starttime:end, 3] .+ 273.15

    #Exp 73 - T9, T10, T11, T12

    F11 = XLSX.readxlsx(path*"Data_FPT0073_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:E4575"]
    y3f_data = F11[starttime:end, 2] .+ 273.15
    y4f_data = F11[starttime:end, 3] .+ 273.15
    y5f_data = F11[starttime:end, 4] .+ 273.15
    y6f_data = F11[starttime:end, 5] .+ 273.15

	E73Tf = y1f1_data
	#E73Ts = mean([y2f1_data, y3f_data, y4f_data, y5f_data, y6f_data])
    E73Ts = mean([y3f_data, y4f_data, y5f_data, y6f_data])
	
	#Exp 74 - T3, T8
    G = XLSX.readxlsx(path*"Data_FPT0074_231130_123228.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:C6018"]
    E74t = Float64.(G[starttime:end, 1])
    y1g_data = G[starttime:end, 2] .+ 273.15
    y2g_data = G[starttime:end, 3] .+ 273.15

	#Exp 74 - T9, T10, T11, T12

    G1 = XLSX.readxlsx(path*"Data_FPT0074_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:E6018"]
    y3g1_data = G1[starttime:end, 2] .+ 273.15
    y4g1_data = G1[starttime:end, 3] .+ 273.15
    y5g1_data = G1[starttime:end, 4] .+ 273.15
    y6g1_data = G1[starttime:end, 5] .+ 273.15

	E74Tf = y1g_data
	#E74Ts = mean([y2g_data, y3g1_data, y4g1_data, y5g1_data, y6g1_data])
	E74Ts = mean([y3g1_data, y4g1_data, y5g1_data, y6g1_data])
	
	#Exp 75 - T3, T8
    H = XLSX.readxlsx(path*"Data_FPT0075_231201_162138.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:C6354"]
    E75t = Float64.(H[starttime:end, 1])
    y1h_data = H[starttime:end, 2] .+ 273.15
    y2h_data = H[starttime:end, 3] .+ 273.15
	
	#Exp 75 - T9, T10, T11, T12

    H1 = XLSX.readxlsx(path*"Data_FPT0075_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:E6354"]
    y3h1_data = H1[starttime:end, 2] .+ 273.15
    y4h1_data = H1[starttime:end, 3] .+ 273.15
    y5h1_data = H1[starttime:end, 4] .+ 273.15
    y6h1_data = H1[starttime:end, 5] .+ 273.15
	
    E75Tf = y1h_data
	#E75Ts = mean([y2h_data, y3h1_data, y4h1_data, y5h1_data, y6h1_data])
	E75Ts = mean([y3h1_data, y4h1_data, y5h1_data, y6h1_data])

	#Exp 76 - T3, T8
    I = XLSX.readxlsx(path*"Data_FPT0076_231203_120521.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:C7147"]
    E76t = Float64.(I[starttime:end, 1])
    y1i_data = I[starttime:end, 2] .+ 273.15
    y2i_data = I[starttime:end, 3] .+ 273.15

	#Exp 76 - T9, T10, T11, T12

    I1 = XLSX.readxlsx(path*"Data_FPT0076_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:E7147"]
    y3i1_data = I1[starttime:end, 2] .+ 273.15
    y4i1_data = I1[starttime:end, 3] .+ 273.15
    y5i1_data = I1[starttime:end, 4] .+ 273.15
    y6i1_data = I1[starttime:end, 5] .+ 273.15
	
	E76Tf = y1i_data
	#E76Ts = mean([y2i_data, y3i1_data, y4i1_data, y5i1_data, y6i1_data])
    E76Ts = mean([y3i1_data, y4i1_data, y5i1_data, y6i1_data])

	 #Exp 77 - T3, T8
    J = XLSX.readxlsx(path*"Data_FPT0077_231203_161315.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:C3044"]
    E77t = Float64.(J[starttime:end, 1])
    y1j_data = J[starttime:end, 2] .+ 273.15
    y2j_data = J[starttime:end, 3] .+ 273.15

    #Exp 77 - T9, T10, T11, T12

    J1 = XLSX.readxlsx(path*"Data_FPT0077_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:E3044"]
    y3j1_data = J1[starttime:end, 2] .+ 273.15
    y4j1_data = J1[starttime:end, 3] .+ 273.15
    y5j1_data = J1[starttime:end, 4] .+ 273.15
    y6j1_data = J1[starttime:end, 5] .+ 273.15

	E77Tf = y1j_data
	#E77Ts = mean([y2j_data, y3j1_data, y4j1_data, y5j1_data, y6j1_data])
    E77Ts = mean([y3j1_data, y4j1_data, y5j1_data, y6j1_data])
	
    #Exp 78 - T3, T8
    K = XLSX.readxlsx(path*"Data_FPT0078_231204_132252.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:C5384"]
    E78t = Float64.(K[starttime:end, 1])
    y1k_data = K[starttime:end, 2] .+ 273.15
    y2k_data = K[starttime:end, 3] .+ 273.15

	 #Exp 78 - T9, T10, T11, T12

    K1 = XLSX.readxlsx(path*"Data_FPT0078_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:E5384"]
    y3k1_data = K1[starttime:end, 2] .+ 273.15
    y4k1_data = K1[starttime:end, 3] .+ 273.15
    y5k1_data = K1[starttime:end, 4] .+ 273.15
    y6k1_data = K1[starttime:end, 5] .+ 273.15

	E78Tf = y1k_data 
	#E78Ts = mean([y2k_data, y3k1_data, y4k1_data, y5k1_data, y6k1_data])
	E78Ts = mean([y3k1_data, y4k1_data, y5k1_data, y6k1_data])

	#Exp 79 - T3, T8
    L1 = XLSX.readxlsx(path*"Data_FPT0079_231204_172244.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:C5233"]
    E79t = Float64.(L1[starttime:end, 1])
    y1l_data = L1[starttime:end, 2] .+ 273.15
    y2l_data = L1[starttime:end, 3] .+ 273.15  

	 #Exp 79 - T9, T10, T11, T12

    L11 = XLSX.readxlsx(path*"Data_FPT0079_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:E5233"]
    y3l1_data = L11[starttime:end, 2] .+ 273.15
    y4l1_data = L11[starttime:end, 3] .+ 273.15
    y5l1_data = L11[starttime:end, 4] .+ 273.15
    y6l1_data = L11[starttime:end, 5] .+ 273.15
	
    E79Tf = y1l_data
	#E79Ts = mean([y2l_data, y3l1_data, y4l1_data, y5l1_data, y6l1_data])
	E79Ts = mean([y3l1_data, y4l1_data, y5l1_data, y6l1_data])

	#Exp 80 - T3, T8
    M = XLSX.readxlsx(path*"Data_FPT0080_231205_095122.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:C5814"]
    E80t = Float64.(M[starttime:end, 1])
    y1m_data = M[starttime:end, 2] .+ 273.15
    y2m_data = M[starttime:end, 3] .+ 273.15   
	
    #Exp 80 - T9, T10, T11, T12

    M1 = XLSX.readxlsx(path*"Data_FPT0080_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:E5814"]
    y3m1_data = M1[starttime:end, 2] .+ 273.15
    y4m1_data = M1[starttime:end, 3] .+ 273.15
    y5m1_data = M1[starttime:end, 4] .+ 273.15
    y6m1_data = M1[starttime:end, 5] .+ 273.15

	E80Tf = y1m_data
	#E80Ts = mean([y2m_data, y3m1_data, y4m1_data, y5m1_data, y6m1_data])
    E80Ts = mean([y3m1_data, y4m1_data, y5m1_data, y6m1_data])
	 #Exp 81 - T3, T8
    N = XLSX.readxlsx(path*"Data_FPT0081_231205_135354.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:C5989"]
    E81t = Float64.(N[starttime:end, 1])
    y1n_data = N[starttime:end, 2] .+ 273.15
    y2n_data = N[starttime:end, 3] .+ 273.15   

	#Exp 81 - T9, T10, T11, T12

    N1 = XLSX.readxlsx(path*"Data_FPT0081_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:E5989"]
    y3n1_data = N1[starttime:end, 2] .+ 273.15
    y4n1_data = N1[starttime:end, 3] .+ 273.15
    y5n1_data = N1[starttime:end, 4] .+ 273.15
    y6n1_data = N1[starttime:end, 5] .+ 273.15

	E81Tf = y1n_data
	#E81Ts = mean([y2n_data, y3n1_data, y4n1_data, y5n1_data, y6n1_data])
	E81Ts = mean([y3n1_data, y4n1_data, y5n1_data, y6n1_data])
	
	md"Exp. data to extract temperature"

end

# ╔═╡ 5d04df82-9b4f-43c4-83c7-13dda03f4cbb
# ╠═╡ disabled = true
#=╠═╡
begin
	plot(E67t, E67Ts)
	plot!(E67t, E67Tf)
end
  ╠═╡ =#

# ╔═╡ ce81f78a-6722-497f-b031-9ac6f933b7fe
begin
	L = 137e-3 #m
    Tamb = 296. #K (same for all exp) 
	th_s = 0.4e-3 #m
    th_f = 0.7e-3 #m
    w_t = 19.e-3 #m
	Lc = 4 * (w_t * w_t) / (4 * w_t) #m
	w_chnl = 1.36e-3 #m
	h_ext = 10 #W/m2.K natural convection
    kins = 0.078 #W/m*K thermal conductivity of insulation
    r_H_chnl = 4 * (w_chnl * w_chnl) / (4 * w_chnl) #hydraulic channel diameter
    r_H = 4 * (w_t * w_t) / (4 * w_t) #hydraulic receiver diameter
    r_ins1 = 23e-3
	r_ins2 = 42e-3 #m location of insulation thermocouple
	r_ins3 = 65e-3
    A_frt = w_t * w_t #m2 - for the whole receiver (19x19mm2)
    A_s_p = w_t * L * 4 #total area solid periphery m2
	A_ins = π*(r_ins1)^2
    A_chnl_p = w_chnl * L * 4  #m2 channel periphery
	A_chnl_frt = w_chnl * w_chnl #m2 channel front
	n_chnl = 10*10
    # A_exchange = A_chnl_p * n_chnl #m2 - contact area between fluid and solid
	A_chnl_frt_all = A_chnl_frt * n_chnl #m2 all frontal area of channels
    Vs =  A_frt * L #m3
	Vf = n_chnl * A_chnl_frt * L #m3
	Af = n_chnl * A_chnl_frt #m2
	#Af =  (w_chnl - (th_s * 2))^2 * n_chnl #m2 - fluid area
	#Vs = (w_chnl + 2 * th_s)^2 * L #m3
	#As = A_chnl_frt_all - Af #m2 - solid area
	#ε = Vf / Vs  #porosity
	ε = ((w_chnl)^2 * L) /  ((w_chnl+(2*th_s))^2 * L)
	Av = (4 * w_chnl * L) / ((w_chnl + th_s)^2 * L) #1/m
	Vi = w_chnl^2 * L * n_chnl
	A_exchange = Av * Vi #m2
	dh = 4 * ε / Av
    q = qlpm / 1000 / 60 #m3/s
    ρf = 1. #kg/m3 at lab conditions
	mf = q * ρf #kg/s
	ρf2 = 0.5 #kg/m3 # at high temp.
	#kf = 0.056 #W/m.K
    em = 0.8
	ab = 0.8
    σ = 5.17e-8 #W/m2.K^4 Stefan-Boltzmann constant
    kf = 0.056 #W/m.K
	#k_eff = 1.97 * (1-ε)^1.5 * ks #Archie's Law
    Cpf = 1090 #J/kg.K
	ρCp_s_0 = 3200 * 1290  #J/kg.K * kg/m3
	ks_0 = 0.1165 * 4184 #W/m.K
	hf_0 = 500 #W/m2.K
    μ = 2.0921e-5 #Pa.s
    
end;

# ╔═╡ 974dac64-8fc4-4eec-ba5e-77bdc65fad92
Vf/Vs

# ╔═╡ a94a8bd4-00f5-47db-816e-77a8e75b8ee0
A_exchange

# ╔═╡ c253cc70-ce35-4ba6-a22f-c6de75b826db
Av

# ╔═╡ 07422568-a64d-4d91-ae07-291a655587e3
Vi

# ╔═╡ 644d36c5-70e9-4f1d-88f3-52859a24e5d0
A_ins

# ╔═╡ 3b6fa5a1-f5ba-4e48-ba8d-0e11315fae5d
ε

# ╔═╡ 01553e7d-e33a-4bac-98d0-7669f410d637
dh

# ╔═╡ 77c3135c-6534-49a6-8ad3-4c22a2495f71
r_H*1000

# ╔═╡ d71c355c-891e-4195-a482-652749bbb5cc
begin
	ρCp_sf(T) = 3200. * (0.27 + 0.135E-4 * T -9720.0 * T^-2 + 0.204E-7 * T^2)  * 4187
	# Cp_sf(T) = (0.27+(0.135e-4)*(T)-9720*((T)^-2)+(0.204e-7)*(T)^2) * 4169
	@register_symbolic ρCp_sf(T)  #kg/m3 * J/kg.K 
	ρf_f(T)= 3.018 * exp(-0.00574*T) + 0.8063*exp(-0.0008381*T) #kg/m3 (Roldan et al.)
	@register_symbolic ρf_f(T)  #kg/m3 * J/kg.K 

	#ρCp_s = 3290. * 46. /40. * 1000.
	#Tx=Tamb:1000
	#plot(Tx, ρCp_sf.(Tx))
	#plot(Tx, ρf_f.(Tx))
end

# ╔═╡ 763c37a2-1169-4ada-8a87-3dd62f0f5ff4
Cps = ρCp_sf(2200)/3200

# ╔═╡ 8f9afcbf-5a98-4709-b763-844b058d155e
#measurements and conditions#Defining simulation conditions
begin
    condition_E67 = Dict(Io => 456000.0, qlpm => 15.27, aIo => :g1_aIo, Tins => 326.437)	
    condition_E68 = Dict(Io => 456000.0, qlpm => 12.50, aIo => :g1_aIo, Tins => 338.52)
    condition_E69 = Dict(Io => 456000.0, qlpm => 10.50, aIo => :g1_aIo, Tins => 344.308)
    condition_E70 = Dict(Io => 456000.0, qlpm => 9.10, aIo => :g1_aIo, Tins => 352.422)
    condition_E71 = Dict(Io => 456000.0, qlpm => 7.12, aIo => :g1_aIo, Tins => 356.004)	
    condition_E72 = Dict(Io => 304000.0, qlpm => 18.34, aIo => :g2_aIo, Tins => 309.928)
    condition_E73 = Dict(Io => 304000.0, qlpm => 13.16, aIo => :g2_aIo, Tins => 325.12)
    condition_E74 = Dict(Io => 304000.0, qlpm => 9.03, aIo => :g2_aIo, Tins => 333.964)	
    condition_E75 = Dict(Io => 304000.0, qlpm => 6.95, aIo => :g2_aIo, Tins=> 336.517)	
    condition_E76 = Dict(Io => 304000.0, qlpm => 4.53, aIo => :g2_aIo, Tins => 338.123)	
    condition_E77 = Dict(Io => 256000.0, qlpm => 13.85, aIo => :g3_aIo, Tins => 308.37)	
    condition_E78 = Dict(Io => 256000.0, qlpm => 10.02, aIo => :g3_aIo, Tins => 312.959)	
    condition_E79 = Dict(Io => 256000.0, qlpm => 8.04, aIo => :g3_aIo, Tins => 314.96)	
    condition_E80 = Dict(Io => 256000.0, qlpm => 6.62, aIo => :g3_aIo, Tins => 316.119)	
    condition_E81 = Dict(Io => 256000.0, qlpm => 4.53, aIo => :g3_aIo, Tins => 319.315)	

	# condition_E67 = Dict(Io => 456000.0, qlpm => 1.22*1000*60*A_chnl_frt_all, aIo => :g1_aIo)
 #    condition_E68 = Dict(Io => 456000.0, qlpm => 1.00*1000*60*A_chnl_frt_all, aIo => :g1_aIo)
 #    condition_E69 = Dict(Io => 456000.0, qlpm => 0.84*1000*60*A_chnl_frt_all, aIo => :g1_aIo)
 #    condition_E70 = Dict(Io => 456000.0, qlpm => 0.73*1000*60*A_chnl_frt_all, aIo => :g1_aIo)
 #    condition_E71 = Dict(Io => 456000.0, qlpm => 0.57*1000*60*A_chnl_frt_all, aIo => :g1_aIo)
 #    condition_E72 = Dict(Io => 304000.0, qlpm => 1.46*1000*60*A_chnl_frt_all, aIo => :g2_aIo)
 #    condition_E73 = Dict(Io => 304000.0, qlpm => 1.05*1000*60*A_chnl_frt_all, aIo => :g2_aIo)
 #    condition_E74 = Dict(Io => 304000.0, qlpm => 0.72*1000*60*A_chnl_frt_all, aIo => :g2_aIo)
 #    condition_E75 = Dict(Io => 304000.0, qlpm => 0.55*1000*60*A_chnl_frt_all, aIo => :g2_aIo)
 #    condition_E76 = Dict(Io => 304000.0, qlpm => 0.36*1000*60*A_chnl_frt_all, aIo => :g2_aIo)
 #    condition_E77 = Dict(Io => 256000.0, qlpm => 1.10*1000*60*A_chnl_frt_all, aIo => :g3_aIo)
 #    condition_E78 = Dict(Io => 256000.0, qlpm => 0.80*1000*60*A_chnl_frt_all, aIo => :g3_aIo)
 #    condition_E79 = Dict(Io => 256000.0, qlpm => 0.64*1000*60*A_chnl_frt_all, aIo => :g3_aIo)
 #    condition_E80 = Dict(Io => 256000.0, qlpm => 0.53*1000*60*A_chnl_frt_all, aIo => :g3_aIo)
 #    condition_E81 = Dict(Io => 256000.0, qlpm => 0.36*1000*60*A_chnl_frt_all, aIo => :g3_aIo)

    simulation_conditions = Dict("E67" => condition_E67, "E68" => condition_E68,
        "E69" => condition_E69, 
		"E70" => condition_E70,
        "E71" => condition_E71, "E72" => condition_E72,
        "E73" => condition_E73, "E74" => condition_E74,
        "E75" => condition_E75, "E76" => condition_E76,
        "E77" => condition_E77, "E78" => condition_E78,
        "E79" => condition_E79, "E80" => condition_E80,
        "E81" => condition_E81)

	measurements = DataFrame(simulation_id = "E67", obs_id="obs_Tf", time = E67t, measurement = E67Tf)
		meas = DataFrame(simulation_id = "E67", obs_id="obs_Ts", time = E67t, measurement = E67Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E68", obs_id="obs_Tf", time = E68t, measurement = E68Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E68", obs_id="obs_Ts", time = E68t, measurement = E68Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E69", obs_id="obs_Tf", time = E69t, measurement = E69Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E69", obs_id="obs_Ts", time = E69t, measurement = E69Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E70", obs_id="obs_Tf", time = E70t, measurement = E70Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E70", obs_id="obs_Ts", time = E70t, measurement = E70Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E71", obs_id="obs_Tf", time = E71t, measurement = E71Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E71", obs_id="obs_Ts", time = E71t, measurement = E71Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E72", obs_id="obs_Tf", time = E72t, measurement = E72Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E72", obs_id="obs_Ts", time = E72t, measurement = E72Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E73", obs_id="obs_Tf", time = E73t, measurement = E73Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E73", obs_id="obs_Ts", time = E73t, measurement = E73Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E74", obs_id="obs_Tf", time = E74t, measurement = E74Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E74", obs_id="obs_Ts", time = E74t, measurement = E74Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E75", obs_id="obs_Tf", time = E75t, measurement = E75Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E75", obs_id="obs_Ts", time = E75t, measurement = E75Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E76", obs_id="obs_Tf", time = E76t, measurement = E76Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E76", obs_id="obs_Ts", time = E76t, measurement = E76Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E77", obs_id="obs_Tf", time = E77t, measurement = E77Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E77", obs_id="obs_Ts", time = E77t, measurement = E77Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E78", obs_id="obs_Tf", time = E78t, measurement = E78Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E78", obs_id="obs_Ts", time = E78t, measurement = E78Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E79", obs_id="obs_Tf", time = E79t, measurement = E79Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E79", obs_id="obs_Ts", time = E79t, measurement = E79Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E80", obs_id="obs_Tf", time = E80t, measurement = E80Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E80", obs_id="obs_Ts", time = E80t, measurement = E80Ts)
		measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E81", obs_id="obs_Tf", time = E81t, measurement = E81Tf)
		measurements = vcat(measurements, meas)
		meas = DataFrame(simulation_id = "E81", obs_id="obs_Ts", time = E81t, measurement = E81Ts)
		measurements = vcat(measurements, meas)
end

# ╔═╡ be19295e-483d-4655-8e8e-55eb349c9958
1.22*1000*60*A_chnl_frt_all

# ╔═╡ 3e44b5ad-ca49-4799-aa9d-af1df817430c
qlpm/(1000*60*A_chnl_frt_all)

# ╔═╡ bc5bf598-46e5-4beb-a9b9-e23c75137fa1
begin
	# #hf = hfa * qlpm^hfn
	
	# Re = ρf * (qlpm/(1000*60*A_chnl_frt_all)) * dh / μ
	# Pr = (Cpf * μ) / kf
	# Gz = Lc * Re * Pr/ L #last position point L=x since model is lumped
	# #Nu = A * (1+(B*((Gz)^n)*exp(-C/Gz)))
	# #Nu = hfa * (Re^hfb) * (Pr^hfn)
	# Nu = A * (Re^B) * (Pr^C)
	# hf = Nu * kf/ dh
	# eq1 = [D(Ts) ~ (1/((1-ε) * aCp * (3200. * (0.27 + 0.135E-4 * Ts -9720.0 * Ts^-2 + 0.204E-7 * Ts^2)  * 4187) * Vs)) * (aIo*Io * A_frt - al * (kins * (r_ins2/r_ins1) * (Ts - Tins) * A_ins / (r_ins2 - r_ins1)) - h_ext * A_frt * (Ts - Tamb) - em * σ * A_frt * (Ts^4 - Tamb^4) - hf * A_exchange * (Ts - Tf)),
	# D(Tf) ~ (1/(ε * (3.018 * exp(-0.00574*Tf) + 0.8063*exp(-0.0008381*Tf)) * Cpf * Vf)) * (hf * A_exchange * (Ts - Tf) - mf * Cpf * (Tf - Tamb))
 #    ]

	# u0 = [Ts => Tamb, Tf => Tamb]

	# state_param = [qlpm => 15.27, Io => 456. *1e3, Tins=>(40. + 273.15)]
	# fit_param = [aCp => 1., aIo => 1., A => 69.9, B => 0.352, C => 6.5, al => 2.]
	# p = vcat(state_param, fit_param)
	# tspan = (0, 3600.)

	#hf = hfa * qlpm^hfn
	
	Re = (ρf * (qlpm/(1000*60*A_chnl_frt_all)) * Lc) / μ
	Pr = (Cpf * μ) / kf
	Gz = Lc * Re * Pr/ L #last position point L=x since model is lumped
	#Nu = A * (1+(B*((Gz)^n)*exp(-C/Gz)))
	#Nu = hfa * (Re^hfb) * (Pr^hfn)
	#Nu = A * (Re^B) * (Pr^C)
	Nu = A * Re
	hf = Nu * kf/ Lc
	eq1 = [D(Ts) ~ (1/((1-ε) * aCp * (3200. * (0.27 + 0.135E-4 * Ts -9720.0 * Ts^-2 + 0.204E-7 * Ts^2)  * 4187) * Vs)) * (aIo*Io * A_frt - al * (kins * (r_ins2/r_ins1) * (Ts - Tins) * A_s_p / (r_ins2 - r_ins1)) - h_ext * A_frt * (Ts - Tamb) - em * σ * A_frt * (Ts^4 - Tamb^4) - hf * A_exchange * (Ts - Tf)),
	D(Tf) ~ (1/(ε * ρf2 * Cpf * Vf)) * (hf * A_exchange * (Ts - Tf) - mf * Cpf * (Tf - Tamb))
    ]

	u0 = [Ts => Tamb, Tf => Tamb]

	state_param = [qlpm => 15.27, Io => 456. *1e3, Tins=>(40. + 273.15)]
	fit_param = [aCp => 1., aIo => 1., A => 69.9, B => 0.352, C => 6.5, al => 2.]
	p = vcat(state_param, fit_param)
	tspan = (0, 3600.)
end

# ╔═╡ 50b80f28-226a-4681-8f39-a721b7aa6d67
begin
	@mtkbuild odes = ODESystem(eq1, t)
	prob = ODEProblem(odes, u0, tspan, p)
	odes
end

# ╔═╡ a1127074-1338-4e12-8ded-bf43c4f32515
sol = solve(prob);

# ╔═╡ 7a997a1e-9731-4824-9f59-731ada97c83e
plot(sol)

# ╔═╡ 8063f54c-0256-4dd3-a0bf-7e2ea10cb671
@bind slaIo Slider(0.8:0.05:1.5, show_value=true)

# ╔═╡ 1485b13d-ab99-467f-9eb0-57db2a7f6c07
@bind slaCp Slider(0.8:0.05:5, show_value=true)

# ╔═╡ 37d17dae-3718-4a8a-bad4-407e6f4b4978
@bind slha Slider(0.1:0.1:2., show_value=true)

# ╔═╡ 1180068b-287b-4038-ac87-b689d42c8c98
# rmp = ModelingToolkit.varmap_to_vars([aCp => slaCp, aIo => slaIo , hfa => slha, hfb => 0., hfn => 0., Io => 456000, Tins => 313., qlpm => 16.47], parameters(odes))
rmp = ModelingToolkit.varmap_to_vars([aCp => slaCp, aIo => slaIo , A => 8, 0*B => 0.352, 0*C => 6.5, al => 1., Io => 456000, Tins => 313., qlpm => 16.47], parameters(odes))

# ╔═╡ b50df44c-f6ef-44b7-b477-0e4516540e6f
begin
	rmprob = remake(prob; p = rmp)
	rmsol = solve(rmprob)
	plot(rmsol, ylim=(270, 1000.), lw=2.)
	plot!(E67t, E67Ts, label="ETs", ls=:dash)
	plot!(E67t, E67Tf, label="ETf", ls=:dash)
end

# ╔═╡ 634082dc-6050-42fc-a208-1a68802fe5da
begin
	_g1_aIo = PEtabParameter(:g1_aIo, lb=1., ub=3., scale=:lin)
	_g2_aIo = PEtabParameter(:g2_aIo, lb=1., ub=3., scale=:lin)
	_g3_aIo = PEtabParameter(:g3_aIo, lb=1., ub=3., scale=:lin)
	#_g1_al = PEtabParameter(:g1_al, lb=0.221, ub=3., scale=:lin)
	#_g2_al = PEtabParameter(:g2_al, lb=0.269, ub=3., scale=:lin)
    #_g3_al = PEtabParameter(:g3_al, lb=0.206, ub=3., scale=:lin)	
	# _hfa = PEtabParameter(hfa, lb=0.01, ub=20., scale=:lin)
	# _hfb = PEtabParameter(hfb, lb=0.1, ub=10., scale=:lin)
	# _hfn = PEtabParameter(hfn, lb=0.01, ub=10., scale=:lin)
	_A = PEtabParameter(A, lb=0.001, ub=5., scale=:lin)
	#_B = PEtabParameter(B, lb=0.01, ub=1., scale=:lin)
	# _n = PEtabParameter(n, lb=0.01, ub=5., scale=:lin)
	#_C = PEtabParameter(C, lb=0.01, ub=20., scale=:lin)
	#_al = PEtabParameter(al, lb=0.01, ub=6., scale=:lin)
	_aCp = PEtabParameter(aCp, lb=0.05, ub=4., scale=:lin)
	params = [_g1_aIo, _g2_aIo, _g3_aIo, _A, _aCp]
	obs_Ts = PEtabObservable(Ts, 0.5)
	obs_Tf = PEtabObservable(Tf, 0.5)
	observables = Dict("obs_Tf" => obs_Tf, "obs_Ts" => obs_Ts)
	petab_model = PEtabModel(odes, simulation_conditions, observables, measurements, params, state_map = [:Tf => 295., :Ts => 295.], parameter_map = [Tins => 313.15], verbose=false)
	petab_problem = PEtabODEProblem(petab_model, verbose=false)
end

# ╔═╡ 4cb05ead-4ff7-4049-a5d5-c6618650d60c
begin
	#p0 = generate_startguesses(petab_problem, 1)
	#res = calibrate_model(petab_problem, [1., 1., 1., 1., 1., 1., 1., 1.], Optim.LBFGS(), options=Optim.Options(iterations = 1000, time_limit=90))
#alg = IpoptOptimiser(true)
#res = calibrate_model_multistart(petab_problem, alg, 10, [])
end

# ╔═╡ b118d8a4-6de8-42b0-9811-53b6dfe82ac9
# ╠═╡ disabled = true
#=╠═╡
begin
	optprob = PEtab.OptimizationProblem(petab_problem)
	optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method = NLopt.LD_LBFGS(), maxtime = 10.0, local_maxiters = 10)
	optsol.u
end
  ╠═╡ =#

# ╔═╡ 7477df3b-f6e0-4c52-8279-63e9a836ff2b
begin
	optprob = PEtab.OptimizationProblem(petab_problem)
	optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method = NLopt.LD_LBFGS(), maxtime = 500.0, local_maxiters = 10000)
	optsol.u
end

# ╔═╡ fa801a0d-3342-4647-9f5d-de18fac34a7a
res = calibrate_model(petab_problem, optsol.u, Optim.LBFGS(), options=Optim.Options(iterations = 1, time_limit=10))

# ╔═╡ add4c024-8c97-4ff8-9fe2-43d22b24c5e6
Pr

# ╔═╡ 8f110eaa-232f-4739-b1e1-d740832e919c
Re

# ╔═╡ 3cb41c98-5494-4276-8147-f3b2e361debd
Nu__ = 67.272*0.00309208*10.5

# ╔═╡ 3f090dba-5a3d-483e-b40b-e51d5002e9a0
h_f = (Nu__ *kf)/ Lc

# ╔═╡ de7b5502-6578-4fb4-9490-a0896bb379b4
begin
	all_cond = [k for (k,v) in simulation_conditions]
	@bind casesim Select(all_cond)
end

# ╔═╡ af38678f-b6b8-42ad-be75-9455ea41edc3
pnew = get_ps(res, petab_problem, condition_id = casesim)

# ╔═╡ e2613b68-6e5c-44dd-a631-65a367726753
plot(res, petab_problem; observable_ids=["obs_Tf"], condition_id=casesim, ylim=(300, 800), ylabel = "Temperature (K)", xlabel = "Time (s)")

# ╔═╡ fef4139f-caad-4cc8-ab05-dcbad38f49f0
plot(res, petab_problem; observable_ids=["obs_Ts"], condition_id=casesim, ylim=(0, 1500), ylabel = "Temperature (K)", xlabel = "Time (s)")

# ╔═╡ b78293e0-a118-4575-be93-a3a1c7b9b0fd
begin
	Tend_m = []
	Tend_e = []
	cond = []
	df_bar = DataFrame(simcond=[], me=[], T=[])
	for (k,v) in simulation_conditions
		simcond = k
		sim_sol = get_odesol(res, petab_problem; condition_id=simcond)
		Temp_m = sim_sol[Tf][end]
		push!(Tend_m, Temp_m)
		Temp_e = subset(measurements, :simulation_id => ByRow(==(simcond)))[!, :measurement][end]
		push!(Tend_e, Temp_e)
		push!(cond, simcond)
		append!(df_bar, DataFrame(simcond=simcond, me="m", T=Temp_m))
		append!(df_bar, DataFrame(simcond=simcond, me="e", T=Temp_e))
	end
end

# ╔═╡ 78264642-e4db-477a-ae15-80a9a086c31b
Tend_e

# ╔═╡ 9c2359ca-0956-4a9a-8706-ab96824ba35e
df_bar

# ╔═╡ af72a10b-bf37-45b3-8593-b57e4bf2ef24
Tend_e

# ╔═╡ e884b556-79a6-4cf4-a09b-0d3afe2489f8
begin

# Tend_model = [726.56, 760.764, 769.829, 766.411, 737.15, 609.797, 677.812, 705.514, 685.29, 613.815, 508.378, 556.195, 564.398, 563.821, 531.047]
# Tend_exp = [763.414, 773.207, 769.76, 779.518, 753.56, 652.955, 694.626, 
# 697.672, 681.066, 647.019, 525.356, 554.455, 560.033, 561.254, 550.499]

# Define the group labels
labels = ["E67", "E68", "E69","E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]


# Plot the grouped bar chart
groupedbar([Tend_e Tend_m], bar_position = :dodge, bar_width=0.5, xticks=(1:15, labels), label=["Experimental T" "Model T"], xlabel="Experimental Runs", ylabel="Temperature (K)", title="Air Outlet Temperature")
end

# ╔═╡ 4fd3723a-531f-4f09-8806-458e5b33f629
begin
	# Initialize an empty array to store relative errors
	rel_errors = []
	
	# Loop through each pair of corresponding values in Tend_model and Tend_exp
	for i in 1:length(Tend_m)
	    rel_error = (Tend_m[i] - Tend_e[i]) * 100/ Tend_e[i]
	    push!(rel_errors, rel_error)
	end
	
	# Print the relative errors
	println("Relative Errors:")
	for (i, rel_error) in enumerate(rel_errors)
	    println("E$i: ", rel_error)
	end
end

# ╔═╡ ead9f712-5a2e-42b3-9cb2-3ef7656472b8
# ╠═╡ disabled = true
#=╠═╡
rel_error = (Tend_model[15] - Tend_exp[15]) / Tend_exp[15]
  ╠═╡ =#

# ╔═╡ 0eb4b057-48c2-4f77-9be9-adace62a9544
Tend_e

# ╔═╡ 8ad40098-80a4-4a44-811d-d94b68c6980e
# ╠═╡ disabled = true
#=╠═╡
begin
	using StatsPlots
	bar(cond, [Tend_m, Tend_e])
	groupedbar(cond, [Tend_m, Tend_e], bar_position = :dodge, bar_width=0.7)

end
  ╠═╡ =#

# ╔═╡ 8fdb601b-55f4-477c-8160-0faa678e9ed1
begin
	scatter(Tend_e, Tend_m, label = "Experimental and Model Temperatures ", legend = :bottomright, xlabel = "Experimental Temperature (K)", ylabel = "Model Temperature (K)", title = "Air Outlet Temperature")
	plot!([560, 1200], [560, 1200], label = "Best Fit")
end

# ╔═╡ dbf6e835-239d-4c71-9a0a-fbb36f0a3f48
# begin
# 	# Nu = A * (1+(B*((Gz)^n)*exp(-C/Gz)))
# 	fpA, fpB, fpC, fpn = Tuple(ModelingToolkit.varmap_to_vars(pnew, [A, B, C, n]))
# 	fNu(Gx) = fpA * (1+(fpB*((Gx)^fpn)*exp(-fpC/Gx)))
# end

# ╔═╡ 51650a85-38c5-4000-b7f7-8a36b85a40d9
Tuple(ModelingToolkit.varmap_to_vars(pnew, [A, B, C]))

# ╔═╡ 6c30cb8d-5756-4173-b867-c99992d666f8
begin
# Nu = A * (Re^B) * (Pr^C)
	fpA, fpB, fpC = Tuple(ModelingToolkit.varmap_to_vars(pnew, [A, B, C]))
	fNu = fpA * (Re^fpB) * (Pr^fpC)
end

# ╔═╡ 04a208ba-3918-408e-877b-e2467b6e53f5
Nu_ = 0.014728*(15.27^0.86106)

# ╔═╡ 1c47435c-1f91-43c2-86eb-6875a240067c
h_ = 0.15399*kf/L

# ╔═╡ 351b2777-c1ae-4fb3-a89e-7742a1acfabe
begin
	Gzx = 10.:50.
	plot(1 ./Gzx, fNu.(Gzx))
end

# ╔═╡ Cell order:
# ╠═8c28c8d1-dd57-4a04-9121-cbfed404d824
# ╠═3e50e470-0251-11ef-3b14-7d58b00ae519
# ╠═232f4256-2095-470e-a669-216e1642c374
# ╠═86c7df59-4f03-4730-9cbf-72d8fc7c34bd
# ╠═3f7ecdc3-b35e-4cc7-8a37-4576cd0f02b9
# ╠═2eb83c61-2ea3-4120-9a55-3e0b66e279bb
# ╠═2c50a54e-96bd-4f48-ae7f-33cbface0754
# ╠═5d04df82-9b4f-43c4-83c7-13dda03f4cbb
# ╠═ce81f78a-6722-497f-b031-9ac6f933b7fe
# ╠═974dac64-8fc4-4eec-ba5e-77bdc65fad92
# ╠═a94a8bd4-00f5-47db-816e-77a8e75b8ee0
# ╠═c253cc70-ce35-4ba6-a22f-c6de75b826db
# ╠═07422568-a64d-4d91-ae07-291a655587e3
# ╠═644d36c5-70e9-4f1d-88f3-52859a24e5d0
# ╠═3b6fa5a1-f5ba-4e48-ba8d-0e11315fae5d
# ╠═01553e7d-e33a-4bac-98d0-7669f410d637
# ╠═77c3135c-6534-49a6-8ad3-4c22a2495f71
# ╠═d71c355c-891e-4195-a482-652749bbb5cc
# ╠═763c37a2-1169-4ada-8a87-3dd62f0f5ff4
# ╠═8f9afcbf-5a98-4709-b763-844b058d155e
# ╠═be19295e-483d-4655-8e8e-55eb349c9958
# ╠═3e44b5ad-ca49-4799-aa9d-af1df817430c
# ╠═bc5bf598-46e5-4beb-a9b9-e23c75137fa1
# ╠═50b80f28-226a-4681-8f39-a721b7aa6d67
# ╠═a1127074-1338-4e12-8ded-bf43c4f32515
# ╠═7a997a1e-9731-4824-9f59-731ada97c83e
# ╠═8063f54c-0256-4dd3-a0bf-7e2ea10cb671
# ╠═1485b13d-ab99-467f-9eb0-57db2a7f6c07
# ╠═37d17dae-3718-4a8a-bad4-407e6f4b4978
# ╠═1180068b-287b-4038-ac87-b689d42c8c98
# ╠═b50df44c-f6ef-44b7-b477-0e4516540e6f
# ╠═634082dc-6050-42fc-a208-1a68802fe5da
# ╠═4cb05ead-4ff7-4049-a5d5-c6618650d60c
# ╠═b118d8a4-6de8-42b0-9811-53b6dfe82ac9
# ╠═7477df3b-f6e0-4c52-8279-63e9a836ff2b
# ╠═fa801a0d-3342-4647-9f5d-de18fac34a7a
# ╠═add4c024-8c97-4ff8-9fe2-43d22b24c5e6
# ╠═8f110eaa-232f-4739-b1e1-d740832e919c
# ╠═3cb41c98-5494-4276-8147-f3b2e361debd
# ╠═3f090dba-5a3d-483e-b40b-e51d5002e9a0
# ╠═de7b5502-6578-4fb4-9490-a0896bb379b4
# ╠═af38678f-b6b8-42ad-be75-9455ea41edc3
# ╠═e2613b68-6e5c-44dd-a631-65a367726753
# ╠═fef4139f-caad-4cc8-ab05-dcbad38f49f0
# ╠═b78293e0-a118-4575-be93-a3a1c7b9b0fd
# ╠═78264642-e4db-477a-ae15-80a9a086c31b
# ╠═9c2359ca-0956-4a9a-8706-ab96824ba35e
# ╠═af72a10b-bf37-45b3-8593-b57e4bf2ef24
# ╠═e884b556-79a6-4cf4-a09b-0d3afe2489f8
# ╠═4fd3723a-531f-4f09-8806-458e5b33f629
# ╠═ead9f712-5a2e-42b3-9cb2-3ef7656472b8
# ╠═0eb4b057-48c2-4f77-9be9-adace62a9544
# ╠═8ad40098-80a4-4a44-811d-d94b68c6980e
# ╠═8fdb601b-55f4-477c-8160-0faa678e9ed1
# ╠═dbf6e835-239d-4c71-9a0a-fbb36f0a3f48
# ╠═51650a85-38c5-4000-b7f7-8a36b85a40d9
# ╠═6c30cb8d-5756-4173-b867-c99992d666f8
# ╠═04a208ba-3918-408e-877b-e2467b6e53f5
# ╠═1c47435c-1f91-43c2-86eb-6875a240067c
# ╠═351b2777-c1ae-4fb3-a89e-7742a1acfabe
