begin

	using CSV
	using DataFrames
	using Plots
	using XLSX
	using Loess
	using PythonPlot
	using Plotly

end

begin
#Exp 26
A = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0026_230302_141014.xlsx")["Sheet 1 - Data_FPT0026_230302_1"]["A3:C4398"]
x_data = A[:,1]
y1_data = A[:,2]
y2_data = A[:,3]
end

# ╔═╡ 745def32-1888-414a-9160-9ae790428057
begin
#Exp 27
B = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0027_230302_131154.xlsx")["Sheet 1 - Data_FPT0027_230302_1"]["A3:C3197"]
xb_data = B[:,1]
y1b_data = B[:,2]
y2b_data = B[:,3]
end

# ╔═╡ 9ef2b698-3b4f-4ed9-92a6-29e7cf020cd4
begin
#Exp 28
C = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0028_230302_111716.xlsx")["Sheet 1 - Data_FPT0028_230302_1"]["A3:C6376"]
xc_data = C[:,1]
y1c_data = C[:,2]
y2c_data = C[:,3]
end

# ╔═╡ c5b65f1e-433d-4491-aa3d-1b1a0cfaba63
begin
#Exp 25
D = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0025_230305_104314.xlsx")["Sheet 1 - Data_FPT0025_230305_1"]["A3:C6690"]
xd_data = D[:,1]
y1d_data = D[:,2]
y2d_data = D[:,3]
end

begin
#Exp 31 - after moving T3 **
E = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0031_230308_105238.xlsx")["Sheet 1 - Data_FPT0031_230308_1"]["A3:C6843"]
xe_data = E[:,1]
y1e_data = E[:,2]
y2e_data = E[:,3]
end


begin
#Exp 29
F = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0029_230305_130623.xlsx")["Sheet 1 - Data_FPT0029_230305_1"]["A3:C3595"]
xf_data = F[:,1]
y1f_data = F[:,2]
y2f_data = F[:,3]
end

begin
#Exp 30
G = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0030_230306_123042.xlsx")["Sheet 1 - Data_FPT0030_230306_1"]["A3:C3595"]
xg_data = G[:,1]
y1g_data = G[:,2]
y2g_data = G[:,3]
end

begin
#Exp 32
H = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0032_230312_104843.xlsx")["Sheet 1 - Data_FPT0032_230312_1"]["A3:C8860"]
xh_data = H[:,1]
y1h_data = H[:,2]
y2h_data = H[:,3]
end

# ╔═╡ 081244ea-8dba-486e-a499-bf45f4d38aa3
begin
#Exp 33
I = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0033_230314_095046.xlsx")["Sheet 1 - Data_FPT0033_230314_0"]["A3:C8064"]
xi_data = I[:,1]
y1i_data = I[:,2]
y2i_data = I[:,3]
end

# ╔═╡ 6cf77b9d-ce89-4433-a6dd-6e5eddd385a0
begin
#Exp 34
J = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0034_230315_104549.xlsx")["Sheet 1 - Data_FPT0034_230315_1"]["A3:C8268"]
xj_data = J[:,1]
y1j_data = J[:,2]
y2j_data = J[:,3]
end

# ╔═╡ 24457e6f-1262-458e-a5e6-ef3b8a85cb2a
begin
#Exp 35
K = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0035_230316_102600.xlsx")["Sheet 1 - Data_FPT0035_230316_1"]["A3:C5194"]
xk_data = K[:,1]
y1k_data = K[:,2]
y2k_data = K[:,3]
end

# ╔═╡ 6303de5b-d281-4412-828d-cc6d3e573566
begin
#Exp 36
L = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0036_230319_095339.xlsx")["Sheet 1 - Data_FPT0036_230319_0"]["A3:C7571"]
xl_data = L[:,1]
y1l_data = L[:,2]
y2l_data = L[:,3]
end

# ╔═╡ 4447eaec-3102-4d0d-8065-a190dec0fd29
begin
#Exp 37
M = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0037_230321_102008.xlsx")["Sheet 1 - Data_FPT0037_230321_1"]["A3:C8212"]
xm_data = M[:,1]
y1m_data = M[:,2]
y2m_data = M[:,3]
end

# ╔═╡ 91cec167-4c6e-4756-9f26-8b314ae08d9e
begin
#Exp 38
N = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0038_230322_102200.xlsx")["Sheet 1 - Data_FPT0038_230322_1"]["A3:C6746"]
xn_data = N[:,1]
y1n_data = N[:,2]
y2n_data = N[:,3]
end

# ╔═╡ 906dd1f7-580f-4c07-92f2-22c5c26367e9
begin
#Exp 39
O = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0039_230323_095708.xlsx")["Sheet 1 - Data_FPT0039_230323_0"]["A3:C7946"]
xo_data = O[:,1]
y1o_data = O[:,2]
y2o_data = O[:,3]
end

# ╔═╡ f684efc9-9155-4e1a-b060-6758d62b4c06
begin
#Exp 40
P = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0040_230323_125951.xlsx")["Sheet 1 - Data_FPT0040_230323_1"]["A3:C5441"]
xp_data = P[:,1]
y1p_data = P[:,2]
y2p_data = P[:,3]
end

# ╔═╡ ca65ac5d-8629-4478-aacc-0160da78b025
begin
#Exp 41
Q = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0041_230326_125415.xlsx")["Sheet 1 - Data_FPT0041_230326_1"]["A3:C7822"]
xq_data = Q[:,1]
y1q_data = Q[:,2]
y2q_data = Q[:,3]
end

# ╔═╡ 9f9232c5-c164-4519-9cbc-72c7ad6fc9c8
begin
#Exp 42
R = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0042_230327_112345.xlsx")["Sheet 1 - Data_FPT0042_230327_1"]["A3:C6428"]
xr_data = R[:,1]
y1r_data = R[:,2]
y2r_data = R[:,3]
end

# ╔═╡ 9449ee55-e2e2-400c-afc0-eff86f3b9f4c
begin
#Exp 43
S = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0043_230327_133110.xlsx")["Sheet 1 - Data_FPT0043_230327_1"]["A3:C6434"]
xs_data = S[:,1]
y1s_data = S[:,2]
y2s_data = S[:,3]
end

# ╔═╡ 5a61fe6e-0c44-4284-bc52-ea865dbb5add
begin
#Exp 44
T = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0044_230530_122524.xlsx")["Sheet 1 - Data_FPT0044_230530_1"]["A3:C5303"]
xt_data = T[:,1]
y1t_data = T[:,2]
y2t_data = T[:,3]
end

# ╔═╡ 5f19e797-94e1-4982-810a-f868a27545c2
begin
#Exp 45
U = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0045_230603_133128.xlsx")["Sheet 1 - Data_FPT0045_230603_1"]["A3:C6412"]
xu_data = U[:,1]
y1u_data = U[:,2]
y2u_data = U[:,3]
end

# ╔═╡ 9a3ddbf7-b3e1-4084-b1dc-7475a2bcb977
begin
#Exp 46
V = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0046_230604_105602.xlsx")["Sheet 1 - Data_FPT0046_230604_1"]["A3:C4741"]
xv_data = V[:,1]
y1v_data = V[:,2]
y2v_data = V[:,3]
end

# ╔═╡ 71d9c864-b501-4404-9bd4-eb1c60b0142f
begin
#Exp 47
W = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0047_230604_132234.xlsx")["Sheet 1 - Data_FPT0047_230604_1"]["A3:C6479"]
xw_data = W[:,1]
y1w_data = W[:,2]
y2w_data = W[:,3]
end

# ╔═╡ 857ec2a7-1ead-4f54-ae62-5cd9649d0330
begin
#Exp 48
X = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0048_230605_114237.xlsx")["Sheet 1 - Data_FPT0048_230605_1"]["A3:C4398"]
xx_data = X[:,1]
y1x_data = X[:,2]
y2x_data = X[:,3]
end

# ╔═╡ 01bd43c3-c9c0-4183-bf59-bdf1633e52c5
begin 
#Exp 59
Y = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0059_230615_094340.xlsx")["Sheet 1 - Data_FPT0059_230615_0"]["A3:C4398"]
xy_data = Y[:,1]
y1y_data = Y[:,2]
y2y_data = Y[:,3]
end

# ╔═╡ e5138587-707f-4a89-b5c1-534abb7582ca
begin
	plot(xc_data,
			y1c_data, label="4.52 lpm")
	
	plot!(x_data,
			y1_data, label= "5.95 lpm", title="T3 - Gas Temp. @ 296 kW/m2 Before Moving T3", ylabel= "Temperature (°C)", xlabel=" Time (s)") 
		plot!(xd_data,
			y1d_data, label="9.73 lpm")
	
		plot!(xb_data, 
			y1b_data, label= "13.5 lpm")
	
end

# ╔═╡ b068c166-a476-4131-b451-bbc0714b6127
begin
	plot(xc_data,
			y2c_data, label="4.52 lpm")

		plot!(x_data,
			y2_data, label= "5.95 lpm", title= "T8 - Solid Temp. @ 296 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
	plot!(xd_data,
			y2d_data, label="9.73 lpm")
		
		plot!(xb_data, 
			y2b_data, label= "13.5 lpm")
	
end

# ╔═╡ 6cd2d1e7-3e0a-4f23-ad85-c615b41519c6
begin
plot(xe_data,
			y1e_data, label= "10.41 lpm **", title= "T3 - Gas Temp. @ 127 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
		#plot!(xf_data,
			#y1f_data, label="15.21 lpm")
		#plot!(xg_data,
			#y1g_data, label="10.41 lpm")
		plot!(xh_data,
			y1h_data, label="7.63 lpm **")
		plot!(xi_data,
			y1i_data, label="5.93 lpm **")
		plot!(xj_data,
			y1j_data, label="4.52 lpm **")
		plot!(xk_data,
			y1k_data, label="15.21 lpm **")

end

# ╔═╡ f6988d49-77df-4a6a-93a7-01e8ae9b01d4
begin
plot(xe_data,
			y2e_data, label= "10.41 lpm **", title= "T8 - Solid Temp. @ 127 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
		plot!(xf_data,
			y2f_data, label="15.21 lpm")
		plot!(xg_data,
			y2g_data, label="10.41 lpm")
			plot!(xh_data,
			y2h_data, label="7.63 lpm **")
		plot!(xi_data,
			y2i_data, label="5.93 lpm **")
		plot!(xj_data,
			y2j_data, label="4.52 lpm **")
		plot!(xk_data,
			y2k_data, label="15.21 lpm **")

end

# ╔═╡ 80ce5b38-0f7b-423b-8dfb-7ffddba1545f
begin
	plot(xn_data,
			y1n_data, label="13.50 lpm **")
plot!(xl_data,
			y1l_data, label= "9.73 lpm **", title= "T3 - Gas Temp. @ 308 kW/m2 After Moving T3",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
		plot!(xm_data,
			y1m_data, label="5.95 lpm **")
		plot!(xo_data,
			y1o_data, label="4.52 lpm **")
end

# ╔═╡ 62b82d90-98ed-4c91-8463-0e8b287a3298
begin
	plot(xn_data,
			y1n_data, label="13.50 lpm")
plot!(xl_data,
			y1l_data, label= "9.73 lpm", title= "T3 - Gas Temperature @ 308 kW/m^2",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
		plot!(xm_data,
			y1m_data, label="5.95 lpm")
		plot!(xo_data,
			y1o_data, label="4.52 lpm")
end

# ╔═╡ 98e4a48e-c360-4801-9547-3ed349f91e4e
begin
	plot(xn_data,
			y2n_data, label="13.50 lpm **")
plot!(xl_data,
			y2l_data, label= "9.73 lpm **", title= "T8 - Solid Temp. @ 308 kW/m2 After Moving T3",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
		plot!(xm_data,
			y2m_data, label="5.95 lpm **")
		plot!(xo_data,
			y2o_data, label="4.52 lpm **")
end

# ╔═╡ 349302df-b49a-4372-80d4-ebd0e130928a
begin
	plot(xn_data,
			y2n_data, label="13.50 lpm")
plot!(xl_data,
			y2l_data, label= "9.73 lpm", title= "T8 - Solid Temperature @ 308 kW/m^2",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
		plot!(xm_data,
			y2m_data, label="5.95 lpm")
		plot!(xo_data,
			y2o_data, label="4.52 lpm")
end

# ╔═╡ a0f81172-2765-4ec0-ad15-ef954d35cc9b
begin
plot(xp_data,
			y1p_data, label= "17.01 lpm **", title= "T3 - Gas Temp. @ 601 kW/m2 After Moving T3",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
	plot!(xq_data,
			y1q_data, label="10.01 lpm **")
	plot!(xr_data,
			y1r_data, label="13.15 lpm **")
	plot!(xs_data,
			y1s_data, label="8.03 lpm **")
end

# ╔═╡ 75e98afc-8f4b-44a8-8b5d-6c8b6c4d91ab
begin
plot(xp_data,
			y2p_data, label= "17.01 lpm **", title= "T8 - Solid Temp. @ 601 kW/m2 After Moving T3",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
	plot!(xq_data,
			y2q_data, label="10.01 lpm **")
	plot!(xr_data,
			y2r_data, label="13.15 lpm **")
	plot!(xs_data,
			y2s_data, label="8.03 lpm **")
	
end

# ╔═╡ 3a73811e-101e-4964-bc9c-d811c946d171
begin 
plot(xt_data,
			y1t_data, label= "9.73 lpm ** - positive defocus ", title= "T3 - Gas Temp. @ 294 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
plot!(xl_data,
			y1l_data, label= "9.73 lpm ** - negative defocus")

end


# ╔═╡ a5bfe045-24e7-452e-aa5b-18988747c1a7
begin 
plot(xx_data, y1x_data, label= "17.01 lpm")
	
plot!(xv_data, y1v_data, label= "13.50 lpm")
	
plot!(xt_data,
			y1t_data, label= "9.73 lpm", title= "T3 - Gas Temp. @ 294 kW/m2 (positive defocus bits)",  ylabel= "Temperature (°C)", xlabel=" Time (s)")
plot!(xu_data, y1u_data, label= "5.95 lpm")
	
plot!(xw_data, y1w_data, label= "4.53 lpm")
end

# ╔═╡ f29410b7-887d-4032-8f62-2cfe69077e70
begin 
plot(xt_data,
			y2t_data, label= "9.73 lpm ** - positive defocus", title= "T8 - Solid Temp. @ 294 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
plot!(xl_data,
			y2l_data, label= "9.73 lpm ** - negative defocus")

end


# ╔═╡ 5852cef1-3e71-4931-89bf-a5d7f2ca398e
begin 
plot(xx_data, y2x_data, label= "17.01 lpm")
	
plot!(xv_data, y2v_data, label= "13.50 lpm")
	
plot!(xt_data,
			y2t_data, label= "9.73 lpm", title= "T8 - Solid Temp. @ 294 kW/m2 (positive defocus bits)",  ylabel= "Temperature (°C)", xlabel=" Time (s)")
plot!(xu_data, y2u_data, label= "5.95 lpm")
	
plot!(xw_data, y2w_data, label= "4.53 lpm")
end

# ╔═╡ a773c2e9-53be-4e89-b16e-6dcc9da65aaa
@bind step Slider(1:1:10)

# ╔═╡ 04b69c26-73f1-4a6b-9440-4e7ffb564489
begin
	gr()
	T1= y1p_data.*1.
	t1= xp_data./60.
	model = loess(t1, T1, span=0.7)
	T1 = predict(model, t1)
	dTdt = [(T1[i]-T1[i-step])/(t1[i]-t1[i-step]) for i = (step+1):length(T1)]
end

# ╔═╡ f227cdd0-a248-4f86-9f6a-0106fb0fe38e
begin
	k = 2
	Snat = interpolate(t1, T1, BSplineOrder(k), Natural())  # with natural BCs
end

# ╔═╡ f5ebbe81-55cf-48b6-83e1-8635a6d53569
begin
	plot(xq_data,
			y1q_data, label="10.01 lpm **, 601 kW/m2")
	plot!(xr_data,
			y1r_data, label="13.15 lpm **, 601 kW/m2")
	plot!(xs_data,
			y1s_data, label="8.03 lpm **, 601 kW/m2")
	plot!(xn_data,
			y1n_data, label="13.50 lpm **, 308 kW/m2")
plot!(xl_data,
			y1l_data, label= "9.73 lpm **, 308 kW/m2", title= "Gas Steady State Temp.",  ylabel= "Temperature (°C)", xlabel=" Time (s)") 
	plot!(xe_data,
			y1e_data, label= "10.41 lpm **, 127 kW/m2")
	plot!(xh_data,
			y1h_data, label="7.63 lpm **, 127 kW/m2")

end

# ╔═╡ 8f562278-71b9-433c-94a9-4052a95cf84b
begin 
	plot(t1, T1, label="T")
	plot!(twinx(), t1[(step+1):end], dTdt, color=:red, label="dTdt")
	
end

# ╔═╡ ff9b0104-799f-4a91-9e3f-8aef84b46b1a
begin
	der_S=Derivative(1)*Snat
	t11=range(0, t1[end], step=1)
	plot( t1, Snat.(t1), color=:blue, label="spline")
	plot!(twinx(), t1, der_S.(t1), color=:green, label="deriv")

end

# ╔═╡ 0965bfd4-fa01-4d4a-9a64-78d1e82a3f31
md"""
# New Experimental Set Up Campaign
"""

# ╔═╡ 07a79684-e035-4617-8289-38ccb6d9654d
begin 
#Exp 67 - T3, T8
Z = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0067_231125_161757.xlsx")["Sheet 1 - Data_FPT0067_231125_1"]["A3:C3932"]
xz_data = Z[:,1]
y1z_data = Z[:,2]
y2z_data = Z[:,3]
end
#Exp 67 - T9, T10, T11, T12
begin
Z1 = XLSX.readxlsx("/Users/aishamelhim/Documents/Github/Aysha/SolarSimulator/EXCEL/Data_FPT0067_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0067_T9"]["A3:E3932"]
y3z_data = Z1[:,2]
y4z_data = Z1[:,3]
y5z_data = Z1[:,4]
y6z_data = Z1[:,5]
end

begin 
#Exp 68 - T3, T8
A1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0068_231126_115725.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:C5365"]
xa1_data = A1[:,1]
y1a1_data = A1[:,2]
y2a1_data = A1[:,3]
end
#Exp 68 - T9, T10, T11, T12
begin
A11 = XLSX.readxlsx("/Users/aishamelhim/Documents/Github/Aysha/SolarSimulator/EXCEL/Data_FPT0068_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:E5365"]
y3a_data = A11[:,2]
y4a_data = A11[:,3]
y5a_data = A11[:,4]
y6a_data = A11[:,5]
end

begin 
#Exp 69 - T3, T8
B1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0069_231126_140153.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:C5366"]
xb1_data = B1[:,1]
y1b1_data = B1[:,2]
y2b1_data = B1[:,3]
end
#Exp 69 - T9, T10, T11, T12
begin
B11 = XLSX.readxlsx("/Users/aishamelhim/Documents/Github/Aysha/SolarSimulator/EXCEL/Data_FPT0069_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:E5366"]
y3b_data = B11[:,2]
y4b_data = B11[:,3]
y5b_data = B11[:,4]
y6b_data = B11[:,5]
end

begin 
#Exp 70 - T3, T8
C1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0070_231127_090339.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:C6705"]
xc1_data = C1[:,1]
y1c1_data = C1[:,2]
y2c1_data = C1[:,3]
end
#Exp 70 - T9, T10, T11, T12
begin
C11 = XLSX.readxlsx("/Users/aishamelhim/Documents/Github/Aysha/SolarSimulator/EXCEL/Data_FPT0070_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:E6705"]
y3c_data = C11[:,2]
y4c_data = C11[:,3]
y5c_data = C11[:,4]
y6c_data = C11[:,5]
end
	
begin 
#Exp 71 - T3, T8
D1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0071_231128_102707.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:C7087"]
xd1_data = D1[:,1]
y1d1_data = D1[:,2]
y2d1_data = D1[:,3]
end
#Exp 71 - T9, T10, T11, T12
begin
D11 = XLSX.readxlsx("/Users/aishamelhim/Documents/Github/Aysha/SolarSimulator/EXCEL/Data_FPT0071_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:E7087"]
y3d_data = D11[:,2]
y4d_data = D11[:,3]
y5d_data = D11[:,4]
y6d_data = D11[:,5]
end

begin 
#Exp 72 - T3, T8
E1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0072_231129_104140.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:C3217"]
xe1_data = E1[:,1]
y1e1_data = E1[:,2]
y2e1_data = E1[:,3]
end
#Exp 72 - T9, T10, T11, T12
begin
E11 = XLSX.readxlsx("/Users/aishamelhim/Documents/Github/Aysha/SolarSimulator/EXCEL/Data_FPT0072_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:E3217"]
y3e_data = E11[:,2]
y4e_data = E11[:,3]
y5e_data = E11[:,4]
y6e_data = E11[:,5]
end

begin 
#Exp 73 - T3, T8
F1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0073_231129_132744.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:C4575"]
xf1_data = F1[:,1]
y1f1_data = F1[:,2]
y2f1_data = F1[:,3]
end
#Exp 73 - T9, T10, T11, T12
begin
F11 = XLSX.readxlsx("/Users/aishamelhim/Documents/Github/Aysha/SolarSimulator/EXCEL/Data_FPT0073_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:E4575"]
y3f_data = F11[:,2]
y4f_data = F11[:,3]
y5f_data = F11[:,4]
y6f_data = F11[:,5]
end

begin 
#Exp 74 - T3, T8
G1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0074_231130_123228.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:C6018"]
xg1_data = G1[:,1]
y1g1_data = G1[:,2]
y2g1_data = G1[:,3]
end


begin 
#Exp 75
H1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0075_231201_162138.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:C6353"]
xh1_data = H1[:,1]
y1h1_data = H1[:,2]
y2h1_data = H1[:,3]
end

# ╔═╡ 2b01c6f3-bc6d-4a75-981b-98dbcf6796c8
begin 
#Exp 76
I1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0076_231203_120521.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:C7141"]
xi1_data = I1[:,1]
y1i1_data = I1[:,2]
y2i1_data = I1[:,3]
end

# ╔═╡ c38189a6-5f72-4b8d-a047-7a96ded9cd10
begin 
#Exp 77
J1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0077_231203_161315.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:C3044"]
xj1_data = J1[:,1]
y1j1_data = J1[:,2]
y2j1_data = J1[:,3]
end

# ╔═╡ de7bcfe1-cf89-48ad-91fd-ebfcb406c112
begin 
#Exp 78
K1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0078_231204_132252.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:C5384"]
xk1_data = K1[:,1]
y1k1_data = K1[:,2]
y2k1_data = K1[:,3]
end

# ╔═╡ 888f2e04-ad6c-4955-bc92-9647df89a47f
begin 
#Exp 79
L1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0079_231204_172244.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:C5233"]
xl1_data = L1[:,1]
y1l1_data = L1[:,2]
y2l1_data = L1[:,3]
end

# ╔═╡ a859522f-f07c-42ac-a8fd-07b66798f411
begin 
#Exp 80
M1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0080_231205_095122.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:C5814"]
xm1_data = M1[:,1]
y1m1_data = M1[:,2]
y2m1_data = M1[:,3]
end

# ╔═╡ 42635db4-9771-4b53-80ba-84e866324e0e
begin 
#Exp 81
N1 = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0081_231205_135354.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:C5989"]
xn1_data = N1[:,1]
y1n1_data = N1[:,2]
y2n1_data = N1[:,3]
end

# ╔═╡ edb4063d-92a6-447e-a34a-5989e3c2961f
md"""
# Plots for New Experimental Set Up Campaign
"""

#T3 at 456 kW/m^2
begin
	plot(xz_data,
			y1z_data, label="15.27 lpm")
	plot!(xa1_data,
			y1a1_data, label="12.50 lpm")
	plot!(xb1_data,
			y1b1_data, label="10.50 lpm")
	plot!(xc1_data,
			y1c1_data, label="9.10 lpm")
plot!(xd1_data,
			y1d1_data, label= "7.12 lpm", title= "Gas Temperature Profile (TI-3) at 456 kW/m\u00B2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end

#T8 at 456 kW/m^2
begin
	plot(xz_data,
			y2z_data, label="15.27 lpm")
	plot!(xa1_data,
			y2a1_data, label="12.50 lpm")
	plot!(xb1_data,
			y2b1_data, label="10.50 lpm")
	plot!(xc1_data,
			y2c1_data, label="9.10 lpm")
plot!(xd1_data,
			y2d1_data, label= "7.12 lpm", title= "Solid Temperature Profile (TI-8) at 456 kW/m\u00B2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end
#T9 at 456 kW/m^2
begin
	plot(xz_data,
			y3z_data, label="15.27 lpm")
	plot!(xa1_data,
			y3a_data, label="12.50 lpm")
	plot!(xb1_data,
			y3b_data, label="10.50 lpm")
	plot!(xc1_data,
			y3c_data, label="9.10 lpm")
plot!(xd1_data,
			y3d_data, label= "7.12 lpm", title= "Solid Temperature Profile (TI-9) at 456 kW/m\u00B2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end
#T10 at 456 kW/m^2
begin
	plot(xz_data,
			y4z_data, label="15.27 lpm")
	plot!(xa1_data,
			y4a_data, label="12.50 lpm")
	plot!(xb1_data,
			y4b_data, label="10.50 lpm")
	plot!(xc1_data,
			y4c_data, label="9.10 lpm")
plot!(xd1_data,
			y4d_data, label= "7.12 lpm", title= "Solid Temperature Profile (TI-10) at 456 kW/m\u00B2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end
#T11 at 456 kW/m^2
begin
	plot(xz_data,
			y5z_data, label="15.27 lpm")
	plot!(xa1_data,
			y5a_data, label="12.50 lpm")
	plot!(xb1_data,
			y5b_data, label="10.50 lpm")
	plot!(xc1_data,
			y5c_data, label="9.10 lpm")
plot!(xd1_data,
			y5d_data, label= "7.12 lpm", title= "Solid Temperature Profile (TI-11) at 456 kW/m\u00B2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end

#T12 at 456 kW/m^2
begin
	plot(xz_data,
			y6z_data, label="15.27 lpm")
	plot!(xa1_data,
			y6a_data, label="12.50 lpm")
	plot!(xb1_data,
			y6b_data, label="10.50 lpm")
	plot!(xc1_data,
			y6c_data, label="9.10 lpm")
plot!(xd1_data,
			y6d_data, label= "7.12 lpm", title= "Solid Temperature Profile (TI-12) at 456 kW/m\u00B2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end

# ╔═╡ 6fe8e52b-4e2a-402f-b0ce-5771d98e7014
begin
	plot(xe1_data,
			y1e1_data, label="18.34 lpm")
	plot!(xf1_data,
			y1f1_data, label="13.16 lpm")
	plot!(xg1_data,
			y1g1_data, label="9.03 lpm")
	plot!(xh1_data,
			y1h1_data, label="6.95 lpm")
plot!(xi1_data,
			y1i1_data, label= "4.53 lpm", title= "Gas Steady State Temp. at 304 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end

# ╔═╡ 139b9668-47df-4bcb-ba2e-3e6ef66b251d
begin
	plot(xe1_data,
			y2e1_data, label="18.34 lpm")
	plot!(xf1_data,
			y2f1_data, label="13.16 lpm")
	plot!(xg1_data,
			y2g1_data, label="9.03 lpm")
	plot!(xh1_data,
			y2h1_data, label="6.95 lpm")
plot!(xi1_data,
			y2i1_data, label= "4.53 lpm", title= "Solid Steady State Temp. (T8) at 304 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)",legend=:bottomright) 

end

# ╔═╡ de622490-5ceb-40c7-a72f-9560542b8058
begin
	plot(xj1_data,
			y1j1_data, label="13.85 lpm")
	plot!(xk1_data,
			y1k1_data, label="10.02 lpm")
	plot!(xl1_data,
			y1l1_data, label="8.04 lpm")
	plot!(xm1_data,
			y1m1_data, label="6.62 lpm")
plot!(xn1_data,
			y1n1_data, label= "4.53 lpm", title= "Gas Steady State Temp. at 256 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end

# ╔═╡ e625719f-ccbc-4e57-947c-d321640e969f
begin
	plot(xj1_data,
			y2j1_data, label="13.85 lpm")
	plot!(xk1_data,
			y2k1_data, label="10.02 lpm")
	plot!(xl1_data,
			y2l1_data, label="8.04 lpm")
	plot!(xm1_data,
			y2m1_data, label="6.62 lpm")
plot!(xn1_data,
			y2n1_data, label= "4.53 lpm", title= "Solid Steady State Temp. at 256 kW/m2",  ylabel= "Temperature (°C)", xlabel=" Time (s)", legend=:bottomright) 

end
