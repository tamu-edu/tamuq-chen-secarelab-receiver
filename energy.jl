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

# ╔═╡ 1ba0ed08-6224-4def-8272-39a31c9cad63
begin
	import Pkg
	Pkg.activate("energy")
end

# ╔═╡ bbad2d96-14cf-4f88-83e7-ccd203486ab4
begin
	Pkg.resolve()
	using ModelingToolkit, DifferentialEquations, Plots
	#using ModelingToolkit: t_nounits as t, D_nounits as D
	using PEtab, XLSX, Statistics, DataFrames
	using PlutoUI, StatsPlots
	using Optimization, Optim, Ipopt, OptimizationNLopt
end

# ╔═╡ a0be0e54-32e8-43d3-b626-6e1cfb6a15fa
Pkg.status()

# ╔═╡ 9be816ec-db63-4737-a966-1b505cd67f94
begin
	@parameters aCp A B aIo #to be fitted
	@parameters qlpm Io Tins #varying conditions
	@variables t Tf(t) Ts(t)
	D = Differential(t)
end

# ╔═╡ 6a720c4e-4bf6-42b7-b447-ad32a6ea44a1
#Exp. data to extract temp.
begin
	starttime = 1
	path = "/Users/aishamelhim/Documents/GitHub/tamuq-chen-secarelab-receiver/aysha/SolarSimulator/EXCEL/"
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

# ╔═╡ 30f79ea4-8c87-41ae-ba41-64ac28d93021
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

# ╔═╡ 19a35210-757c-4298-a1f4-1b6bf92cbd62
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

# ╔═╡ ed9ac6d6-9ae1-4d51-bd44-0278bfd767ac
#measurements and conditions#Defining simulation conditions
begin
    condition_E67 = Dict(Io => 456000.0, qlpm => 15.27, aIo => :g1_aIo, Tins => 326.437) #, al => 0.553)	
    condition_E68 = Dict(Io => 456000.0, qlpm => 12.50, aIo => :g1_aIo, Tins => 338.52) #, al => 0.553) 
    condition_E69 = Dict(Io => 456000.0, qlpm => 10.50, aIo => :g1_aIo, Tins => 344.308) #, al => 0.553)
    condition_E70 = Dict(Io => 456000.0, qlpm => 9.10, aIo => :g1_aIo, Tins => 352.422) #, al => 0.553)
    condition_E71 = Dict(Io => 456000.0, qlpm => 7.12, aIo => :g1_aIo, Tins => 356.004) #, al => 0.553)
    condition_E72 = Dict(Io => 304000.0, qlpm => 18.34, aIo => :g2_aIo, Tins => 309.928) #, al => 0.672)
    condition_E73 = Dict(Io => 304000.0, qlpm => 13.16, aIo => :g2_aIo, Tins => 325.12) #, al => 0.672)
    condition_E74 = Dict(Io => 304000.0, qlpm => 9.03, aIo => :g2_aIo, Tins => 333.964) #, al => 0.672)
    condition_E75 = Dict(Io => 304000.0, qlpm => 6.95, aIo => :g2_aIo, Tins=> 336.517) #, al => 0.672)
    condition_E76 = Dict(Io => 304000.0, qlpm => 4.53, aIo => :g2_aIo, Tins => 338.123) #, al => 0.672)	
    condition_E77 = Dict(Io => 256000.0, qlpm => 13.85, aIo => :g3_aIo, Tins => 308.37) #, al => 0.515)	
    condition_E78 = Dict(Io => 256000.0, qlpm => 10.02, aIo => :g3_aIo, Tins => 312.959) #, al => 0.515)	
    condition_E79 = Dict(Io => 256000.0, qlpm => 8.04, aIo => :g3_aIo, Tins => 314.96) #, al => 0.515)
    condition_E80 = Dict(Io => 256000.0, qlpm => 6.62, aIo => :g3_aIo, Tins => 316.119) #, al => 0.515)
    condition_E81 = Dict(Io => 256000.0, qlpm => 4.53, aIo => :g3_aIo, Tins => 319.315) #, al => 0.515)

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
		#meas = DataFrame(simulation_id = "E67", obs_id="obs_Ts", time = E67t, measurement = E67Ts)
		#measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E68", obs_id="obs_Tf", time = E68t, measurement = E68Tf)
		measurements = vcat(measurements, meas)
		#meas = DataFrame(simulation_id = "E68", obs_id="obs_Ts", time = E68t, measurement = E68Ts)
		#measurements = vcat(measurements, meas)
	#meas = DataFrame(simulation_id = "E69", obs_id="obs_Tf", time = E69t, measurement = E69Tf)
	#	measurements = vcat(measurements, meas)
		#meas = DataFrame(simulation_id = "E69", obs_id="obs_Ts", time = E69t, measurement = E69Ts)
		#measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E70", obs_id="obs_Tf", time = E70t, measurement = E70Tf)
		measurements = vcat(measurements, meas)
		#meas = DataFrame(simulation_id = "E70", obs_id="obs_Ts", time = E70t, measurement = E70Ts)
		#measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E71", obs_id="obs_Tf", time = E71t, measurement = E71Tf)
		measurements = vcat(measurements, meas)
		#meas = DataFrame(simulation_id = "E71", obs_id="obs_Ts", time = E71t, measurement = E71Ts)
		#measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E72", obs_id="obs_Tf", time = E72t, measurement = E72Tf)
		measurements = vcat(measurements, meas)
		#meas = DataFrame(simulation_id = "E72", obs_id="obs_Ts", time = E72t, measurement = E72Ts)
		#measurements = vcat(measurements, meas)
	meas = DataFrame(simulation_id = "E73", obs_id="obs_Tf", time = E73t, measurement = E73Tf)
		measurements = vcat(measurements, meas)
		#meas = DataFrame(simulation_id = "E73", obs_id="obs_Ts", time = E73t, measurement = E73Ts)
		#measurements = vcat(measurements, meas)
end

# ╔═╡ 7f2ac27f-5a77-489f-970c-d3add0b03201
begin
	Re = (ρf * (qlpm/(1000*60*A_chnl_frt_all)) * Lc) / μ
	Pr = (Cpf * μ) / kf
	Gz = Lc * Re * Pr/ L #last position point L=x since model is lumped
	#Nu = A * (1+(B*((Gz)^n)*exp(-C/Gz)))
	#Nu = hfa * (Re^hfb) * (Pr^hfn)
	#Nu = A * (Re^B) * (Pr^C)
	Nu = A * (Re^B)
	hf = Nu * kf/ Lc
	eq1 = [D(Ts) ~ (1/((1-ε) * aCp * (3200. * (0.27 + 0.135E-4 * Ts -9720.0 * Ts^-2 + 0.204E-7 * Ts^2)  * 4187) * Vs)) * (aIo*Io * A_frt - (kins * (r_ins2/r_ins1) * (Ts - Tins) * A_s_p / (r_ins2 - r_ins1)) - h_ext * A_frt * (Ts - Tamb) - em * σ * A_frt * (Ts^4 - Tamb^4) - hf * A_exchange * (Ts - Tf)),
	D(Tf) ~ (1/(ε * ρf2 * Cpf * Vf)) * (hf * A_exchange * (Ts - Tf) - mf * Cpf * (Tf - Tamb))
    ]

	u0 = [Ts => Tamb, Tf => Tamb]

	state_param = [qlpm => 15.27, Io => 456. *1e3, Tins=>(40. + 273.15)]
	fit_param = [aCp => 1., aIo => 1., A => 2., B => 0.5]
	p = vcat(state_param, fit_param)
	tspan = (0, 3600.)
end

# ╔═╡ 038e82bd-63cc-448c-9eb9-4f63fb7a7314
begin
	@mtkbuild odes = ODESystem(eq1, t)
	prob = ODEProblem(odes, u0, tspan, p)
	odes
end

# ╔═╡ b69c4adc-5438-4cf5-b5c8-454b2da41e3e
sol = solve(prob);

# ╔═╡ 4731d76c-75fc-4a54-9902-c6e3587f5dee
@bind slaIo Slider(0.8:0.05:1.5, show_value=true)

# ╔═╡ 7ab499a8-0b7e-4418-9ea7-092ee539c252
@bind slaCp Slider(0.8:0.05:5, show_value=true)

# ╔═╡ c2568acb-72af-49e3-aa88-e99cfcc86592
@bind slha Slider(0.1:0.1:2., show_value=true)

# ╔═╡ 6da2f0d3-dacf-4d4d-bf90-8e36ab292e4c
rmp = ModelingToolkit.varmap_to_vars([aCp => slaCp, aIo => slaIo, A => 2., B => 0.5, Io => 456000, Tins => 313., qlpm => 16.47], parameters(odes))

# ╔═╡ 1a204d41-2c49-4189-84f0-dbfd2aa21b54
begin
	rmprob = remake(prob; p = rmp)
	rmsol = solve(rmprob)
	plot(rmsol, ylim=(270, 1000.), lw=2.)
	plot!(E67t, E67Ts, label="ETs", ls=:dash)
	plot!(E67t, E67Tf, label="ETf", ls=:dash)
end

# ╔═╡ afa2ecef-5874-44dc-a05e-f8b2cf5692e2
begin
	_g1_aIo = PEtabParameter(:g1_aIo, lb=0.7, ub=2., scale=:lin)
	_g2_aIo = PEtabParameter(:g2_aIo, lb=0.7, ub=2., scale=:lin)
	_g3_aIo = PEtabParameter(:g3_aIo, lb=0.7, ub=2., scale=:lin)
	#_g1_al = PEtabParameter(:g1_al, lb=0.221, ub=3., scale=:lin)
	#_g2_al = PEtabParameter(:g2_al, lb=0.269, ub=3., scale=:lin)
    #_g3_al = PEtabParameter(:g3_al, lb=0.206, ub=3., scale=:lin)	
	#_hfa = PEtabParameter(hfa, lb=0.001, ub=5., scale=:lin)
	# _hfb = PEtabParameter(hfb, lb=0.1, ub=10., scale=:lin)
	#_hfn = PEtabParameter(hfn, lb=0.001, ub=10., scale=:lin)
	_A = PEtabParameter(A, lb=0.01, ub=20., scale=:lin)
	_B = PEtabParameter(B, lb=0.01, ub=1., scale=:lin)
	#_n = PEtabParameter(n, lb=0.01, ub=10., scale=:lin)
	#_C = PEtabParameter(C, lb=0.001, ub=60., scale=:lin)
	#_al = PEtabParameter(al, lb=0.01, ub=1., scale=:lin)
	_aCp = PEtabParameter(aCp, lb=0.05, ub=5., scale=:lin)
	params = [_g1_aIo, _g2_aIo, _g3_aIo, _A, _B, _aCp]
	#obs_Ts = PEtabObservable(Ts, 0.5)
	obs_Tf = PEtabObservable(Tf, 0.5)
	#observables = Dict("obs_Tf" => obs_Tf, "obs_Ts" => obs_Ts)
	observables = Dict("obs_Tf" => obs_Tf)
	petab_model = PEtabModel(odes, simulation_conditions, observables, measurements, params, state_map = [:Tf => 295.], parameter_map = [Tins => 313.15], verbose=false)
	petab_problem = PEtabODEProblem(petab_model, verbose=false)
end

# ╔═╡ 4471e831-3b6e-4ffa-96a1-528ad3c00be6
# begin
# 	optprob = PEtab.OptimizationProblem(petab_problem)
# 	optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method = NLopt.LD_LBFGS(), maxtime = 10.0, local_maxiters = 10)
# 	optsol.u
# end

# ╔═╡ c3644645-44d7-4b6a-b719-b6df82b25787
begin
	p0 = generate_startguesses(petab_problem, 1)
	res = calibrate_model(petab_problem, [1., 1., 1., 1., 1., 1.], Optim.LBFGS(), options=Optim.Options(iterations = 1000, time_limit=90))
	#res = calibrate_model_multistart(petab_problem, IpoptOptimiser(false), 10)
end

# ╔═╡ a7357e86-69ed-404b-89dd-a2dfd7e0708b
begin
	all_cond = [k for (k,v) in simulation_conditions]
	@bind casesim Select(all_cond)
end


# ╔═╡ 5dee22dc-b9c7-4299-9e67-2ea440e16251
pnew = get_ps(res, petab_problem, condition_id = casesim)

# ╔═╡ ab7b647a-18fd-4110-ad00-62168e1758c3
plot(res, petab_problem; observable_ids=["obs_Tf"], condition_id=casesim, ylim=(300, 1000))

# ╔═╡ Cell order:
# ╠═a0be0e54-32e8-43d3-b626-6e1cfb6a15fa
# ╠═bbad2d96-14cf-4f88-83e7-ccd203486ab4
# ╠═1ba0ed08-6224-4def-8272-39a31c9cad63
# ╠═9be816ec-db63-4737-a966-1b505cd67f94
# ╠═6a720c4e-4bf6-42b7-b447-ad32a6ea44a1
# ╠═30f79ea4-8c87-41ae-ba41-64ac28d93021
# ╠═19a35210-757c-4298-a1f4-1b6bf92cbd62
# ╠═ed9ac6d6-9ae1-4d51-bd44-0278bfd767ac
# ╠═7f2ac27f-5a77-489f-970c-d3add0b03201
# ╠═038e82bd-63cc-448c-9eb9-4f63fb7a7314
# ╠═b69c4adc-5438-4cf5-b5c8-454b2da41e3e
# ╠═4731d76c-75fc-4a54-9902-c6e3587f5dee
# ╠═7ab499a8-0b7e-4418-9ea7-092ee539c252
# ╠═c2568acb-72af-49e3-aa88-e99cfcc86592
# ╠═6da2f0d3-dacf-4d4d-bf90-8e36ab292e4c
# ╠═1a204d41-2c49-4189-84f0-dbfd2aa21b54
# ╠═afa2ecef-5874-44dc-a05e-f8b2cf5692e2
# ╠═4471e831-3b6e-4ffa-96a1-528ad3c00be6
# ╠═c3644645-44d7-4b6a-b719-b6df82b25787
# ╠═a7357e86-69ed-404b-89dd-a2dfd7e0708b
# ╠═5dee22dc-b9c7-4299-9e67-2ea440e16251
# ╠═ab7b647a-18fd-4110-ad00-62168e1758c3
