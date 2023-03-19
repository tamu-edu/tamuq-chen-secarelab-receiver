### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 1e060c7e-886f-4208-b68e-09f2b08195fa
using PlutoUI

# ╔═╡ ad798898-61ea-4188-9c5f-d258e10cc47e
using Plots

# ╔═╡ a226d405-d910-4112-8b6f-a153420a38ca
using DelimitedFiles, LinearAlgebra, Statistics

# ╔═╡ 838cbfe0-f0ad-11eb-0146-a989eeaafda6
md"""
# Plotting COMSOL Sweep Results

##### Base Case (Uniform Reflectivity)

"""

# ╔═╡ 63b2aa9e-bdbf-4670-bc0a-ff86013a20e7
plotly()
# gr()

# ╔═╡ aa57714a-6ba3-48ec-afcd-31b44845a4e8
readdlm("Parametric Sweep Levels.txt")

# ╔═╡ 6317f5cc-bbf9-4c22-a813-c500a538d380
# Parameter Levels
begin
	R = [0.5, 1, 2, 4, 8]
	E = [0.1, 0.3, 0.5, 0.7, 0.9]
	LD = [25]
end

# ╔═╡ cdfd8cb1-7602-45ad-abc4-8d9ed53eee19
t = 0.4 	#channel wall half-thickness [mm]

# ╔═╡ 1f6244a6-09a4-4e7c-9468-8aa6ff376070
# Additional Parameters of Interest
begin
	# Re = readdlm("P:\\Chemical Engineering\\SECAReLab\\@Projects\\PJ.ABD.SolarAerosolSynthesis\\Comsol\\Cavity\\Sweep Results - Sep 15\\Reynolds No.txt")[5:end,1:3]
	# Re = readdlm("Reynolds No.txt")[6:end,1:3]
	# Re = Re[:,3]				#Reynolds number, not saved - used from diff. sweep
	# Re = [Re[i] for i = 1:length(R):length(R)*length(E)]
	ϕ = (R.^2) ./ (R .+ t).^2 	#Porosity

end

# ╔═╡ 5470fdd0-21cb-495c-9c69-ed9ba537e9bf


# ╔═╡ 1f1828f3-2f72-4da4-9058-7729e6316bfb
md"""

###### Checking if energy balance closes
"""

# ╔═╡ 3186cee0-226b-48d8-bad1-aba11756d9d2


# ╔═╡ 5c0de861-b973-4bf1-9b76-c2af3ce1a21b
md"""
###### Rest of processing data
"""

# ╔═╡ c037aea2-371d-45d1-b6e3-fb51aad8cb7c
 readdlm("Global Powers 1.txt")[5,:]

# ╔═╡ ce5f13f3-ec44-4cf1-8cba-fd9c331c7ca0
Eff = readdlm("Global Powers 1.txt"; comments=true, comment_char='%')

# ╔═╡ f82863f3-d2dc-46a8-9cb3-0b1a129a0508
P_ap_lin = readdlm("Power on Apperture.txt"; comments=true, comment_char='%')[:,end]

# ╔═╡ 1124a3c4-2b53-448a-9d63-278da7a2d834
begin
	n_param = 2
	
	Q_abs_lin = Eff[:,n_param+1]	# Heat absorbed by gas [W]
	BHS_lin = Eff[:,n_param+2]		# Boundary heat source [W]
	T_in_lin = Eff[:,n_param+3]		# Gas inflow temp. [K]
	T_out_lin = Eff[:,n_param+4]	# Gas outflow temp. [K]
	Q_rad_loss_lin = Eff[:,n_param+5]		# Radiative loss [W]
	Q_rad_loss_int_lin = Eff[:,n_param+6]	# Radiative loss from interior surfaces [W]
	Q_rad_loss_front_lin = Eff[:,n_param+7]	# Radiative loss from front surface [W]
	
	Q_abs_out_lin = Eff[:,n_param+8]	# Total energy flux through outlet [W]
	Q_abs_in_lin = Eff[:,n_param+9]	# Total energy flux through inlet [W]
	BHS_int_lin = Eff[:,n_param+10]	# Primary ray power absorbed by interior surfaces [W]
	BHS_front_lin = Eff[:,n_param+11]	# Primary ray power absorbed by front surface  [W]

end

# ╔═╡ da408ddb-5b01-4376-ad34-c19c3177910d
begin
	T_out = zeros(length(R),length(E))
	for i = 1:length(R), j = 1:length(E)
		lin_index = (i-1)*length(E) + j
		T_out[i,j] = T_out_lin[lin_index]
	end
end

# ╔═╡ 7b557345-bafa-44fb-b515-cd651a66ef7e
function interpolate_lin(U, x, x_int)
	#U = known 2D property array
	#x = levels at which 'U' is known (y is assumed the same for both U and U_int)
	#x_int = levels at which 'U' is desired
	Y,X = size(U)
	U_int = zeros(Y, length(x_int))
	
	for r = 1:Y, e = 1:length(x_int)
		slope = (U[r,e] - U[r,e+1])/(x[e] - x[e+1])
		
		U_int[r,e] = slope*(x_int[e] - x[e]) + U[r,e]
	end
	
	return U_int
end

# ╔═╡ 115b38d9-1c87-488d-8905-fd04ddc71682
# # Interpolating all relevant properties (for wall distribution cases)
# begin
# 	Eff_Ap_int = interpolate_lin(Eff_Ap, E, E_base)
# 	T_out_int = interpolate_lin(T_out, E, E_base)
# 	vol_eff_int = interpolate_lin(vol_eff[:,ld,:], E, E_base)
# 	Q_abs_tot_int = interpolate_lin(Q_abs_tot, E, E_base)
# 	BHS_tot_int = interpolate_lin(BHS_tot, E, E_base)
# 	Q_loss_tot_int = interpolate_lin(Q_rad_loss_tot, E, E_base)
# 	Q_loss_tot_front_int = interpolate_lin(Q_rad_loss_tot_front, E, E_base)
# 	Q_loss_tot_int_int = interpolate_lin(Q_rad_loss_tot_int, ε, E_base)
# 	T_s_in_int = interpolate_lin(T_s_in, E, E_base)
# 	T_s_avg_int = interpolate_lin(T_s_avg, E, E_base)
# 	interpolated_results = [
# 		Eff_Ap_int; T_out_int; vol_eff_int; Q_abs_tot_int; BHS_tot_int; Q_loss_tot_int_int; Q_loss_tot_front_int; Q_loss_tot_int; T_s_in_int; T_s_avg_int
# 	]
# end

# ╔═╡ 3f30e71e-c0ac-4c9a-9ee7-e3707400118c
# open("interpolated_results_base.txt", "w") do io
#            writedlm(io, interpolated_results)
#        end

# ╔═╡ 5ba04d8b-0c65-46dc-b080-51413388dc30
log2.(R)

# ╔═╡ 5bf03252-b019-431a-9106-7a1df1f56806
begin
	# adding an arrow for direction to maximize efficiency
	x1 = collect(range(-0.1, 3, step = 0.1));
	m1 = (0.9 - 0.3)/(-0.1 - 2.5); a = pi/1	;
	c1 = 0.9;
	y1 = m1*(x1 .- -0.1) .+ c1;
end

# ╔═╡ 42efe36f-025c-44ad-8890-a32a86f7b571
begin
	contourf(
		log2.(R), E, transpose(T_out),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(T_out),
		# xlabel = "Radius (mm)",
		# xscale=:log3,
		# xticks = R,
		# xlim = (2e-1, 2e3),
		#---------------------------------------------
		# Re, E, transpose(T_out),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(T_out),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Exit Gas Temperature (K)",
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		titlefont = 16,
		# colorbar_title = "η",
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 14,
	 	# guidefont = 16,
		# titlefont = 18,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (440,335),	#for word doc.
			
		)
	
	# adding lines for selected design
	hline!(
		[0.8],
		line = :dash,
		linewidth = 2,
		color = :black,
		legend=:false
	)

	vline!(
		[0.],
		line = :dash,
		linewidth = 2,
		color = :black,
		legend=:false
	)
end

# ╔═╡ 69c19e51-1271-443e-8a36-30268e7afdf4
md"""
#### Observations:
- For a given channel size/porosity, increasing emissivity improves both exit gas temperature and thermal efficiency
 - For a given emissivity, there is an optimal value for porosity
"""

# ╔═╡ a87b183a-aebb-4398-920d-0f4e4e8a51a9


# ╔═╡ 1a9f5a82-5ce7-470d-a70e-67c766f932fd
md"""
## Auxillary Results

Intermediate fluxes and temperature profiles
"""

# ╔═╡ 44cec87f-537d-426d-8089-934210aee5b3
function line_avg(u, x, m, n)
	u_dx = zeros(length(u))
	for i = m:n
		if i == m
			u_dx[i] = u[1]*(x[2]-x[1])/2
		elseif i == n
			u_dx[i] = u[n]*(x[n]-x[n-1])/2
		else
			u_dx[i] = u[i]*(x[i+1]-x[i-1])/2
		end
	end
		
	u_avg = sum(u_dx)/(x[n] - x[m])
	return u_avg
end

# ╔═╡ 03f25a88-efb5-40bb-9c13-4ca2b8b97c14
σ = 5.670374419e-8 #W.m-2.K-4 (Stefan-Boltzmann Const)

# ╔═╡ 2e280ca0-1971-42f4-9f92-c04d9c464847
h_nat = 5 #W/m^2.K

# ╔═╡ ebe6df24-c623-4ffd-9857-40110fb0336d
T_amb = 318 	#K

# ╔═╡ 0593c86b-fb0f-40c2-a38f-9c7330c7a484
begin
	D_tot =  140 #[mm] "Inscribed diameter of square SolAir-200 reciever module"
	A_tot = D_tot^2 #"Area of square receiver module"
	n_channels(R_channel, t_channel) = (D_tot/(R_channel*2+t_channel*2))^2 
	q_ap =  650 #[kW/m^2] "Flux density on apperture"
	P_to_m = 700 # [kJ/kg] "Power on apperture to mass flowrate ratio"
	m_tot = q_ap*A_tot/P_to_m # "Mass flowrate on module"
	m_channel(R_channel,t_channel) =  m_tot/n_channels(R_channel,t_channel) # "Mass flowrate per channel"
end

# ╔═╡ 61c6bcf4-cf2c-400c-9922-32c6eef9dec0
P_ap_tot = [P_ap_lin[r]*n_channels(R[r],t) for r = 1:length(R), e = 1:length(E)]

# ╔═╡ b53b7e81-b817-42d2-8170-74affd2f5e12
P_ap_lin*n_channels(R[1],t)

# ╔═╡ 722d42d3-5c99-4e2f-826b-892006fac5e3
begin
	Eff_Ap = zeros(length(R), length(E))
	
	Q_abs = zeros(length(R), length(E))
	Q_abs_spec = zeros(length(R),length(E))
	Q_abs_tot = zeros(length(R),length(E))

	for r = 1:length(R), e = 1:length(E)
		lin_index = (r-1)*length(E) + e
		
		Q_abs[r,e] = Q_abs_lin[lin_index]		
		#--------------------------------------------------------------------------
		m_gas = m_channel(R[r]*1e-3,t)
		Q_abs_spec[r,e] = Q_abs_lin[lin_index]/m_gas #[W/kg]
		#--------------------------------------------------------------------------
		Q_abs_tot[r,e] = Q_abs_lin[lin_index]*n_channels(R[r],t) #[W]vper mod
		
		#--------------------------------------------------------------------------

		Eff_Ap[r,e] = Q_abs_lin[lin_index]./P_ap_lin[r]
	end
end

# ╔═╡ 31ee1595-be66-4d8a-9be1-96b2514dcc0c
begin
	contourf(
		log2.(R), E, transpose(Eff_Ap),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(Eff_Ap),
		# xlabel = "Radius (mm)",
		# xaxis = :log2, 
		# xlim = [-1, 8],
		#---------------------------------------------
		# Re, E, transpose(Eff_Ap),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(Eff_Ap),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Thermal Efficiency, η",# - using Power on Apperture",
		#--------------------------------------------
		# Size formatting for ppt
		framestyle = :box,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		titlefont = 16,
		# colorbar_title = "η",
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 14,
	 #    guidefont = 16,
		# titlefont = 18,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (440,335),	#for word doc.
			
		)

	# adding lines for selected design
# 	hline!(
# 		[0.8],
# 		line = :dash,
# 		linewidth = 2,
# 		color = :black,
# 		legend=:false
# 		)

# 	vline!(
# 		[0.],
# 		line = :dash,
# 		linewidth = 2,
# 		color = :black,
# 		legend=:false
# 		)
	
	# adding an arrow for direction to maximize efficiency
	# x1 = collect(range(-0.1, 3, step = 0.1));
	# m1 = (0.9 - 0.3)/(-0.1 - 2.5); a = pi/1	;
	# c1 = 0.9;
	# y1 = m1*(x1 .- -0.1) .+ c1;
	
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :black,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :black,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :black,
	# 	linewidth = 2,
	# 	label="")
end

# ╔═╡ 2c0b3dea-5df3-450d-b5a9-8a89930bc417
plot(
	log2.(R), 
	[Eff_Ap[:,e] for e = 1:length(E)],
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Eff_Ap),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Eff_Ap),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Eff_Ap),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Thermal Efficiency, η",# - using Boundary Heat Source",
	label = permutedims(["ε = $(E[e])" for e = 1:length(E)]),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 5*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 11,
    guidefont = 12,
	titlefont = 14,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
		
	)

# ╔═╡ f226cd78-4380-4340-a3b9-9d7bb1940abb
md"""
## Additional plots

For heat transfer analysis
"""

# ╔═╡ e09d1bcb-66f5-4e1a-bc34-1f1dc936335a
md"""
#### Plots that use absolute power values (per channel)

"""

# ╔═╡ 1bc2fdd5-36cf-4cb8-a388-9ac949ceabf2
contourf(
	log2.(R), E, transpose(Q_abs),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_abs),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_abs),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_abs),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Heat Absorbed by the Gas (W)",
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 12*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
		
	)

# ╔═╡ b78bdbdd-ba48-4e44-b19b-19caffe2b307


# ╔═╡ e5300a51-d793-4fe2-88b4-8a117924f91b
md"""
#### Plots that use specific power values (per unit mass or area per area)

"""

# ╔═╡ dec3866e-881e-42a0-9236-e65fd6fd8599
contourf(
	log2.(R), E, transpose(Q_abs_spec),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_abs_spec),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_abs_spec),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_abs_spec),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Heat Absorbed by the Gas,  <i>Q<sub>abs,g</sub>/m<sub>g</sub></i> <br> per channel (W/kg)",
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 12*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
		
	)

# ╔═╡ 2235d2e4-a48c-4102-a947-b40c0b1bdeb2


# ╔═╡ 2a730b37-d2d3-4b32-9c47-24fd5b122d95
md"""
#### Plots that use absolute power values (per module)

"""

# ╔═╡ 76f0b4fd-dc94-499d-8d41-c58cc2d3ff2d
contourf(
	log2.(R), E, transpose(Q_abs_tot),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_abs_tot),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_abs_tot),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_abs_tot),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Heat Absorbed by the Gas,  <i>Q<sub>abs,g</sub> </i> <br> per Module (W)",
	#--------------------------------------------
	# Size formatting for ppt
	top_margin = 10*Plots.mm,
	left_margin = 3*Plots.mm,
	right_margin = 3*Plots.mm,
	tickfontsize = 12,
	guidefont = 14,
	legendfont = 11,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (389,300),	#for ppt
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	# top_margin = 10*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 10,
 	# guidefont = 11,
	# titlefont = 12,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (325,250),	#for word doc.
		
	)

# ╔═╡ 736b6056-df0f-4f5b-8b8a-736bee68370f
cat([1], collect(range(5, 80, step=5));dims=1)

# ╔═╡ 538ae313-d93f-4c38-b049-491d21e15ecc
md"""
Plotting powers as fractions of the power on aperture is not as clear as absolute powers 
"""

# ╔═╡ db524b4c-336b-4a78-983e-4e1a72798d88


# ╔═╡ 622c7a82-f772-4705-b72f-aac21141114b


# ╔═╡ d0513898-51d8-4871-b36a-35d5fcaa0629


# ╔═╡ a2c93673-67ae-4ece-918c-add42bf0255e


# ╔═╡ c944cbc9-a9b6-4538-836d-eafea086189f
md"""
##### Plotting Temperature Profiles and Volumetric Effect
"""

# ╔═╡ dbcf829c-06a9-4914-a846-c0dd7ebdabbb
Temp_profile = "Cutline Temp. Profiles"

# ╔═╡ 693c84e3-f6c2-48b4-98d3-cc1e9e8ddc60
T_gas_dict = [
			readdlm("$(Temp_profile)//T_gas (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,2:end]
	
	for r = 1:length(R), ld = 1:length(LD)
			]

# ╔═╡ 54f952ce-8f4b-4a0f-bf1e-117fc48d3996
x_gas_dict = [
	readdlm("$(Temp_profile)//T_gas (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,1]

	for r = 1:length(R), ld = 1:length(LD)
			]

# ╔═╡ ff09ca2a-54b9-4112-abcf-92b744a86f9b
T_solid_dict = [
	readdlm("$(Temp_profile)//T_solid (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,2:end]
			
	for r = 1:length(R), ld = 1:length(LD)
			]

# ╔═╡ aaa5252e-0714-476d-8267-f13d8a1f9e0d
begin
	Q_rad_loss = zeros(length(R),length(E))
	Q_rad_loss_spec = zeros(length(R),length(E))
	
	Q_cnat_loss_tot = zeros(length(R),length(E))

	Q_rad_loss_tot = zeros(length(R),length(E))
	Q_rad_loss_tot_front = zeros(length(R),length(E))
	Q_rad_loss_tot_int = zeros(length(R),length(E))

	A_surface_tot = zeros(length(R))
	A_surface_tot_front = zeros(length(R))
	A_surface_tot_int = zeros(length(R))
	
	for r = 1:length(R), e = 1:length(E), ld = 1:length(LD)
		lin_index = (r-1)*length(E) + e

		A_walls = 2*4*R[r]*(LD[ld]*2*R[r])	#[mm^2]
		A_front = 4*(R[r]+t)^2 - 4*R[r]^2 	#[mm^2]
			
		T_s = T_solid_dict[r,ld][1,e]
		
		#--------------------------------------------------------------------------	
		Q_rad_loss[r,e] = Q_rad_loss_lin[lin_index]	# [W] per channel	
		Q_rad_loss_front = Q_rad_loss_front_lin[lin_index] 	#[W]
		Q_rad_loss_int = Q_rad_loss_int_lin[lin_index]  	#[W]
		#--------------------------------------------------------------------------
		Q_rad_loss_spec[r,e] = Q_rad_loss_lin[lin_index]/((A_walls + A_front)*1e-6) 	#[W/m^2]	
		#--------------------------------------------------------------------------	
		Q_rad_loss_tot[r,e] = Q_rad_loss[r,e]*n_channels(R[r], t) #[W]
		Q_rad_loss_tot_front[r,e] = Q_rad_loss_front*n_channels(R[r], t) #[W]
		Q_rad_loss_tot_int[r,e] = Q_rad_loss_int*n_channels(R[r], t) #[W]
				
		A_surface_tot_front[r] = A_front*n_channels(R[r], t) #[m^2]
		A_surface_tot_int[r] = A_walls*n_channels(R[r], t) #[m^2]
		A_surface_tot[r] =  A_surface_tot_front[r] + A_surface_tot_int[r] #[m^2]

		Q_cnat_loss_tot[r,e] = h_nat*(T_s[1] - T_amb)*A_front*1e-6*n_channels(R[r],t) #[W]
		
	end
end

# ╔═╡ 69f256c3-c603-4797-9571-73993d519580
contourf(
	log2.(R), E, transpose(Q_rad_loss),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_rad_loss),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_rad_loss),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_rad_loss),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	# title = "Average Boundary Heat Source (W/m<sup>2</sup>)",
	title = "Radiative Heat Loss per Channel (W)",
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 12*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
		
	)

# ╔═╡ 2caf1b17-51cf-4d96-93d3-948104c92fa5
contourf(
	log2.(R), E, transpose(Q_rad_loss_spec),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_rad_loss_spec),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_rad_loss_spec),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_rad_loss_spec),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Radiative Heat Loss per Channel (W/m<sup>2</sup>)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 12*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
	
	)

# ╔═╡ c917cb94-36a6-49f1-bf2c-7681f6172310
begin
	contourf(
		log2.(R), E, transpose(Q_rad_loss_tot),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(Q_rad_loss_tot),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(Q_rad_loss_tot),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(Q_rad_loss_tot),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Radiative Heat Loss per <br> Module (W)",
		# c = cgrad(:thermal, rev = true),
		c = cgrad(:thermal),
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 12*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		legendfont = 11,
		titlefont = 16,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 10*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 10,
	 #    guidefont = 11,
		# titlefont = 12,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (325,250),	#for word doc.
			
		)
	hline!(
		[0.8],
		line = :dash,
		color = :white,
		linewidth = 2,
		legend = false
		)
	
	vline!(
		[0.],
		line = :dash,
		color = :white,
		linewidth = 2,
		legend = false
		)
	
end

# ╔═╡ b87236e7-2e21-4d78-8f39-31292cae6218
[Q_rad_loss_tot[r,e]/(P_ap_lin[r]*n_channels(R[r],t)) for r = 1:length(R),  e = 1:length(E)]

# ╔═╡ 7e095764-2bb6-44ed-afc7-38f5cdff4d92
contour(
	log2.(R), E,
	transpose([Q_rad_loss_tot[r,e]/(P_ap_lin[r]*n_channels(R[r],t)) for r = 1:length(R),  e = 1:length(E)]),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_rad_loss_tot),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_rad_loss_tot),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_rad_loss_tot),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Q<sub>rad, loss</sub> / Power on Aperture",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 5*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 10,
    guidefont = 11,
	titlefont = 12,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (325,250),	#for word doc.
		
	)

# ╔═╡ 1672b8ce-a4a5-44c5-aa7a-54af60a29451
contourf(
	log2.(R), E, transpose(Q_rad_loss_tot_front),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_rad_loss_tot_front),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_rad_loss_tot_front),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_rad_loss_tot_front),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Radiative Heat Loss per Module<br> from Front Surface (W)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 10*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 10,
    guidefont = 11,
	titlefont = 12,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (325,250),	#for word doc.
		
	)

# ╔═╡ 549fbc4f-b840-4a49-b0dc-2174d9d8468d
contourf(
	log2.(R), E, transpose(Q_rad_loss_tot_int),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_rad_loss_tot_int),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_rad_loss_tot_int),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_rad_loss_tot_int),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Radiative Heat Loss per Module<br> from Interior Surfaces (W)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 10*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 10,
    guidefont = 11,
	titlefont = 12,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (325,250),	#for word doc.
		
	)

# ╔═╡ 7a6b554f-d315-4184-93b0-879b6e987e41
begin
	contour(
		log2.(R), E, transpose(Q_rad_loss_tot_int./Q_rad_loss_tot_front),
		# log2.(R), ε, transpose(BHS_int_tot./BHS_tot),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, ε, transpose(BHS_tot),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, ε, transpose(BHS_tot),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, ε, transpose(BHS_tot),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Q<sub>rad, loss (int)</sub> / Q<sub>rad, loss (front)</sub>",
		# title = "Total Boundary Heat Source (W)",
		# c = cgrad(:thermal, rev = true),
		c = cgrad(:thermal),
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 10*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		legendfont = 11,
		titlefont = 16,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 10,
	 #    guidefont = 11,
		# titlefont = 12,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (325,250),	#for word doc.
			
		)
	
	hline!(
		[0.8],
		line = :dash,
		color = :black,
		linewidth = 2,
		legend = false
		)
	
	vline!(
		[0.],
		line = :dash,
		color = :black,
		linewidth = 2,
		legend = false
		)
end

# ╔═╡ f005f378-0f08-4e0a-902d-a49ba445643f
contourf(
	log2.(R), E, transpose(Q_cnat_loss_tot),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_rad_loss_tot),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_rad_loss_tot),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_rad_loss_tot),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Natural Convection Loss <br>per Module (W)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 10*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 10,
    guidefont = 11,
	titlefont = 12,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (325,250),	#for word doc.
		
	)

# ╔═╡ b3509e00-6bfe-4ea2-9552-90d6aacb19a1
begin
	plot(
		log2.(R),
		[
			A_surface_tot,
			A_surface_tot_int,
			A_surface_tot_front
		]./1e6,
		label = ["Total Irradiated Surface Area" "Interior Surface Area" "Frontal Surface Area"],
		color = permutedims(palette(:tab10)[1:3]),
		linewidth = 2,
		xlabel = "Channel Radius, <i>R<sub>ch</sub><i>",
		ylabel = "Surface Area, <i>A</i> (mm<sup>2</sup>)",
		title = "Net Irradiated Surface Areas per Module",
		legend = :outerbottomright,
		top_margin = 12*Plots.mm,
		left_margin = 5*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (624,250)

	)
end

# ╔═╡ 2c8cc361-2eb0-445f-95d7-7d65295a3087
begin
	T_s_in = zeros(length(R),length(E))
	for r = 1:length(R), ld = 1:length(LD), e = 1:length(E)
		T_s = T_solid_dict[r,ld][:,e]
		T_s_in[r,e] = T_s[1]
	end
end

# ╔═╡ 475c4317-427d-4558-923b-982662107ea2
contourf(
	log2.(R), E, transpose(T_s_in),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(T_s_in),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(T_s_in),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(T_s_in),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Front Solid Surface Temp. (K)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 5*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
	
	)

# ╔═╡ b3c91cbe-baa8-421a-8ff1-b2b0bca215fa
x_solid_dict = [
	readdlm("$(Temp_profile)//T_solid (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,1]

	for r = 1:length(R), ld = 1:length(LD)
			]

# ╔═╡ a71d51c0-efdb-45bb-bb1b-821b1ee4bb81
begin
	T_s_avg = zeros(length(R),length(E))
	for r = 1:length(R), ld = 1:length(LD), e = 1:length(E)
		x = x_solid_dict[r,ld]
		A_int = (x*1e-3)*(4*R[r]*2*1e-3) 	#Interior surface area per channel
		T_s = T_solid_dict[r,ld][:,e]
		m = 1	# Start length considered for integration
		n = length(x)
	
		T_s_avg[r,e] = line_avg(T_s, x, m, n)
	end
end

# ╔═╡ 1dab20f7-dd45-42c1-9565-156b9d4f5009
contourf(
	log2.(R), E, transpose(T_s_avg),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(T_s_avg),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(T_s_avg),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(T_s_avg),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Average Solid Surface Temp. (K)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 12*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
	
	)

# ╔═╡ be8c4bdf-501a-4f2c-8a0d-d7e034040b24
colors = permutedims(palette(:tab10)[1:length(E)])
# colors = permutedims(palette(:curl, length(E)))

# ╔═╡ da253a52-1954-4890-842b-b9444917deeb
@bind r1 PlutoUI.Slider(1:length(R), show_value=true)

# ╔═╡ 92bd2aae-62c7-4f84-9427-16e09a7bac14
begin
	# r1 = 2
	ld1 = 1
	plot(
			x_gas_dict[r1,ld1],

			[T_gas_dict[r1,ld1][:,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :dash,
		linewidth = 2,
		label = false
		)

	plot!(
			x_solid_dict[r1,ld1],

			[T_solid_dict[r1,ld1][:,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :solid,
		linewidth = 2,
		label = permutedims(["Emissivity = $(E[e])" for e = 1:length(E)]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Axial Temperature Profiles (R<sub>ch</sub> = $(R[r1])mm)",
		legend = :bottomright,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 10,
		guidefont = 14,
		titlefont = 16,
		size = (624,477)

		)

end

# ╔═╡ fc079169-d9ba-4f8d-a57d-613e8f2c5d49
md"""
#### Extracting and Plotting The Volumetric Effect/Efficiency
"""

# ╔═╡ bd7dd4ab-ee29-4780-92a6-af77c440f6ba
vol_eff = 	[
	T_solid_dict[r,ld][1,e1]\T_gas_dict[r,ld][end,e1]	
	for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)
			]

# ╔═╡ e9065f52-64bf-45d2-b5da-e75a4b408607
@bind ld PlutoUI.Slider(1:length(LD), show_value=true)

# ╔═╡ fa856550-9cd1-496b-a2b5-2cf79c22dd33
begin
	contourf(
		log2.(R), E, transpose(vol_eff[:,ld,:]),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, T<sub>g,out</sub> / T<sub>s,in</sub>",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		framestyle = :box,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		titlefont = 16,
		# colorbar_title = "η",
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 14,
	 #    guidefont = 16,
		# titlefont = 18,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (440,335),	#for word doc.
			
		)

	# adding lines for selected design
	# hline!(
	# 	[0.8],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend=:false
	# )

	# vline!(
	# 	[0.],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend=:false
	# )
end

# ╔═╡ 3c5eaa90-ec96-4757-9356-956a98840aef
begin
	BHS = zeros(length(R),length(E))
	BHS_spec = zeros(length(R),length(E))
	BHS_tot = zeros(length(R),length(E))
		
	BHS_int_tot = zeros(length(R),length(E))
	BHS_front_tot = zeros(length(R),length(E))
		
	for r = 1:length(R), e = 1:length(E)
		lin_index = (r-1)*length(E) + e
		BHS[r,e] = BHS_lin[lin_index]	# [W]
		A_walls = 2*4*R[r]*(LD[ld]*2*R[r])	#[mm^2]
		A_front = 4*(R[r]+t)^2 - 4*R[r]^2 	#[mm^2]
		BHS_spec[r,e] = BHS_lin[lin_index]/((A_walls + A_front)*1e-6) 	#[W/m^2]
		BHS_tot[r,e] = BHS_lin[lin_index] * n_channels(R[r],t) 	#[W]	
		
		BHS_int_tot[r,e] = BHS_int_lin[lin_index] * n_channels(R[r],t) 	#[W]
		BHS_front_tot[r,e] = BHS_front_lin[lin_index] * n_channels(R[r],t) 	#[W]
	end
end

# ╔═╡ 7c3fc9a5-3f6b-428e-a32e-9c44c92fdb08
maximum(BHS_tot .- Q_abs_tot .- Q_rad_loss_tot_front .- Q_rad_loss_tot_int)

# ╔═╡ ba014dd3-7b46-444e-8337-ba096d164127
minimum(BHS_tot .- Q_abs_tot .- Q_rad_loss_tot_front .- Q_rad_loss_tot_int)

# ╔═╡ e74945d2-06f7-4d36-bb27-37e5e10095ba
# EB = P_ap_tot .- (Q_abs_tot .+ Q_rad_loss_tot)
EB = BHS_tot .- (Q_abs_tot .+ Q_rad_loss_tot)

# ╔═╡ 464815e3-6f96-4b1d-8a22-a88afcd7aafb
contourf(R,E,transpose(EB))

# ╔═╡ f57af54b-3a7f-4337-982a-8c1b21554694
interpolated_results = [
		Eff_Ap; T_out; vol_eff[:,1,:]; Q_abs_tot; BHS_tot; Q_rad_loss_tot_int; Q_rad_loss_tot_front; Q_rad_loss_tot; T_s_in; T_s_avg
	]

# ╔═╡ 0cdb6335-31ce-4ee2-9e81-b3d025a11179
contourf(
	log2.(R), E, transpose(BHS),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(BHS),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(BHS),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(BHS),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	# title = "Total Boundary Heat Source (W/m<sup>2</sup>)",
	title = "Total Boundary Heat Source (W)",
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 12*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
		
	)

# ╔═╡ 8cb6d0d1-9301-419e-a242-817b0b718ae6
contourf(
	log2.(R), E, transpose(BHS_spec),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(BHS_spec),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(BHS_spec),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(BHS_spec),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Total Boundary Heat Source (W/m<sup>2</sup>)",
	# title = "Total Boundary Heat Source (W)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 12*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 12,
    guidefont = 14,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (624/1.5,477/1.5),	#for word doc.
		
	)

# ╔═╡ ee777a68-05da-4696-aa3e-d4fcddd3367a
begin
	contourf(
		log2.(R), E, transpose(BHS_tot),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(BHS_tot),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(BHS_tot),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(BHS_tot),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Total Boundary Heat Source <br> per Module (W)",
		# title = "Total Boundary Heat Source (W)",
		# c = cgrad(:thermal, rev = true),
		c = cgrad(:thermal),
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 12*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		legendfont = 11,
		titlefont = 16,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 10*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 10,
	 #    guidefont = 11,
		# titlefont = 12,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (325,250),	#for word doc.
			
		)

	hline!(
		[0.8],
		line = :dash,
		color = :black,
		linewidth = 2,
		legend = false
		)
	
	vline!(
		[0.],
		line = :dash,
		color = :black,
		linewidth = 2,
		legend = false
		)
end

# ╔═╡ 90b3c5be-feeb-4e94-b96f-49622ee88ec6
contourf(
	log2.(R), E, transpose(BHS_int_tot),
	# log2.(R), ε, transpose(BHS_int_tot./BHS_tot),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, ε, transpose(BHS_tot),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, ε, transpose(BHS_tot),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, ε, transpose(BHS_tot),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Interior Boundary Heat Source <br> per Module (W)",
	# title = "Total Boundary Heat Source (W)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	top_margin = 10*Plots.mm,
	left_margin = 3*Plots.mm,
	right_margin = 3*Plots.mm,
	tickfontsize = 12,
	guidefont = 14,
	legendfont = 11,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (389,300),	#for ppt
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	# top_margin = 10*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 10,
 #    guidefont = 11,
	# titlefont = 12,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (325,250),	#for word doc.
		
	)

# ╔═╡ 0a38ee74-9d7a-451f-a171-9aabdbc2939c
contourf(
	log2.(R), E, transpose(BHS_front_tot),
	# log2.(R), ε, transpose(BHS_front_tot./BHS_tot),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, ε, transpose(BHS_tot),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, ε, transpose(BHS_tot),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, ε, transpose(BHS_tot),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "Front Boundary Heat Source <br> per Module (W)",
	# title = "Total Boundary Heat Source (W)",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	top_margin = 10*Plots.mm,
	left_margin = 3*Plots.mm,
	right_margin = 3*Plots.mm,
	tickfontsize = 12,
	guidefont = 14,
	legendfont = 11,
	titlefont = 16,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (389,300),	#for ppt
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	# top_margin = 10*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 10,
 #    guidefont = 11,
	# titlefont = 12,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (325,250),	#for word doc.
		
	)

# ╔═╡ 238739ff-d8c0-4089-a87b-a435e2c08325
BHS_int_tot./BHS_front_tot

# ╔═╡ 416b7566-578e-4606-8fee-6305e3f5a11d
begin
	contour(
		log2.(R), E, transpose(BHS_int_tot./BHS_front_tot),
		# log2.(R), ε, transpose(BHS_int_tot./BHS_tot),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, ε, transpose(BHS_tot),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, ε, transpose(BHS_tot),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, ε, transpose(BHS_tot),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "BHS<sub>int</sub> / BHS<sub>front</sub>",
		# title = "Total Boundary Heat Source (W)",
		# c = cgrad(:thermal, rev = true),
		c = cgrad(:thermal),
		levels = range(1, 80, step=4),
		# levels = cat([1], collect(range(5, 80, step=5));dims=1),
		#--------------------------------------------
		# Size formatting for ppt
		framestyle = :box,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		legendfont = 11,
		titlefont = 16,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 10,
	 #    guidefont = 11,
		# titlefont = 12,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (325,250),	#for word doc.	
		)

	# hline!(
	# 	[0.8],
	# 	line = :dash,
	# 	color = :black,
	# 	linewidth = 2,
	# 	legend = false
	# 	)
	
	# vline!(
	# 	[0.],
	# 	line = :dash,
	# 	color = :black,
	# 	linewidth = 2,
	# 	legend = false
	# 	)
end

# ╔═╡ 410622b9-8213-4424-a09f-ff5eff7487c9
BHS_int_tot./BHS_front_tot

# ╔═╡ 31a18798-98f7-4d0b-935d-be8610eae64e
contour(
	log2.(R), E,
	transpose([BHS_tot[r,e]/(P_ap_lin[r]*n_channels(R[r],t)) for r = 1:length(R),  e = 1:length(E)]),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(Q_rad_loss_tot),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(Q_rad_loss_tot),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(Q_rad_loss_tot),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Emissivity",
	title = "BHS / Power on Aperture",
	# c = cgrad(:thermal, rev = true),
	c = cgrad(:thermal),
	#--------------------------------------------
	# Size formatting for ppt
	# top_margin = 12*Plots.mm,
	# left_margin = 3*Plots.mm,
	# tickfontsize = 12,
	# guidefont = 14,
	# titlefont = 16,
	# fontfamily = "ComputerModern",
	# colorbar_font = "ComputerModern",
	# contour_labels = true,
	# size = (624,477), #for ppt
	#--------------------------------------------
	# Size formatting for word
	top_margin = 5*Plots.mm,
	left_margin = 3*Plots.mm,
	tickfontsize = 10,
    guidefont = 11,
	titlefont = 12,
	fontfamily = "ComputerModern",
	colorbar_font = "ComputerModern",
	contour_labels = true,
	size = (325,250),	#for word doc.
		
	)

# ╔═╡ ed899759-bc40-4d2a-9b2f-323eeb4870a0
md"""
##### (1) Based on inlet and outlet temperatures
"""

# ╔═╡ e26b7863-94ed-44d1-a0e3-a30b2546cfd7
# Based on definiton of LMTD
begin
	∆T1 = [T_solid_dict[r,ld][1,e1] - T_gas_dict[r,ld][1,e1] for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)]
	
	∆T2 = [T_solid_dict[r,ld][end,e1] - T_gas_dict[r,ld][end,e1] for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)]
	
	LMTD = 	[
		(∆T1[r,ld,e1] - ∆T2[r,ld,e1])/log(∆T1[r,ld,e1] / ∆T2[r,ld,e1])
		
	for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)
			]

	# LMTD = LMTD_f(T_solid_dict, T_gas_dict)
	
	LMTD_A = 	[
		LMTD[r,ld,e1] * A_surface_tot_int[r]*1e-6
		
	for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)
			]
end

# ╔═╡ 00f75f99-b015-4c20-89a9-dd47b903cd60
function LMTD_f(T_solid_dict, T_gas_dict)
	
	∆T1 = [T_solid_dict[r,ld][1,e1] - T_gas_dict[r,ld][1,e1] for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)]
	
	∆T2 = [T_solid_dict[r,ld][end,e1] - T_gas_dict[r,ld][end,e1] for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)]
	
	LMTD = 	[
		(∆T1[r,ld,e1] - ∆T2[r,ld,e1])/log(∆T1[r,ld,e1] / ∆T2[r,ld,e1])
		
	for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)
			]
	
	return LMTD
end

# ╔═╡ f29302bb-ac90-4f6a-bbfb-ea22a8dae0db


# ╔═╡ 939a8c76-7fd1-43ca-ac90-8dd6fe4847df
md"""
##### (2) Based on power ratios
"""

# ╔═╡ 43d07b3b-19d4-410b-a02b-b21bd162d16c
# Based on ratio of powers
begin
	vol_eff_p1 = Q_abs_tot./Q_rad_loss_tot
	vol_eff_p2 = Q_rad_loss_tot./BHS_tot
	vol_eff_p3 = Q_abs_tot./BHS_tot
	vol_eff_p4 = [BHS_tot[r,e]/(P_ap_lin[r] * n_channels(R[r],t)) for r = 1:length(R), e = 1:length(E)]

end

# ╔═╡ 69be736e-08ac-4ee4-b291-f18c88cf6caf


# ╔═╡ bc815fec-4251-47b2-99e7-2f4407d5bb87
md"""
##### (3) Based on statistical evaluations of ∆T
"""

# ╔═╡ 4a3698e5-bb3a-4d79-95df-41226407873e
function interpolate_lin_1D(U, z, z_int)
	#U = known 1D property array
	#z = levels at which 'U' is known 
	#z_int = levels at which 'U' is desired
	
	Z = length(U)
	U_int = zeros(length(z_int))
	
	for i = 1:length(z_int)
		
		if z_int[i] == z[1]
			U_int[i] = U[1]
		elseif z_int[i] == z[end]
			U_int[i] = U[end]
		else
			r_ind = findfirst(z_int[i] .< z) - 1
			slope = (U[r_ind] - U[r_ind+1])/(z[r_ind] - z[r_ind+1])
			U_int[i] = slope*(z_int[i] - z[r_ind]) + U[r_ind]
		end
	end
	
	return U_int
end

# ╔═╡ 20d1daa7-4015-403c-95b2-716b31137dc2
# Extracting temperature profiles at the same 'z' points
begin
	T_s_int = Array{Any, 2}(undef, (length(R),length(LD)))
	T_g_int = Array{Any, 2}(undef, (length(R),length(LD)))
	Z_int = Array{Any, 2}(undef, (length(R),length(LD)))
	
	for r = 1:length(R), ld = 1:length(LD),  e = 1:length(E)
		L_channel = 2*R[r]*LD[ld] #mm
		z_int = range(0, L_channel, step = 1)
		Z_int[r,ld] = z_int
		
		T_s_int[r,ld] = [interpolate_lin_1D(T_solid_dict[r,ld][:,e], x_solid_dict[r,ld], z_int) for e = 1:length(E)]

		# println("r = ", r, "  ld = ", ld)
		T_g_int[r,ld] = [interpolate_lin_1D(T_gas_dict[r,ld][:,e], x_gas_dict[r,ld], z_int) for e = 1:length(E)]

		
	end
end
		

# ╔═╡ 91e63de9-3227-41b9-8a7a-d44cfbd24309
begin
	# r1 = 2
	# ld1 = 1
	plot(
			x_gas_dict[r1,ld1],

			[T_gas_dict[r1,ld1][:,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :dash,
		linewidth = 2,
		label = false
		)

	plot!(
			x_solid_dict[r1,ld1],

			[T_solid_dict[r1,ld1][:,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :solid,
		linewidth = 2,
		label = permutedims(["Emissivity = $(E[e])" for e = 1:length(E)]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Axial Temperature Profiles (R<sub>ch</sub> = $(R[r1])mm)",
		legend = :bottomright,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 10,
		guidefont = 14,
		titlefont = 16,
		size = (624,477)

		)

	scatter!(
			Z_int[r1,ld1],

			[T_g_int[r1,ld1][e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		linewidth = 2,
		label = false
		)
	
	scatter!(
			Z_int[r1,ld1],

			[T_s_int[r1,ld1][e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :dash,
		linewidth = 2,
		label = false
		)
end

# ╔═╡ 60362194-1571-44f9-89c1-8ae09f35e9e6


# ╔═╡ aa8713e4-f665-4177-adae-9e9560aa4c81
# Based on statistical evaluations of ∆T
begin
	RMSE = Array{Any, 2}(undef, (length(R),length(E)))
	R_corr = Array{Any, 2}(undef, (length(R),length(E)))
	NMB = Array{Any, 2}(undef, (length(R),length(E)))
	NMSD = Array{Any, 2}(undef, (length(R),length(E)))

	SIG_S = Array{Any, 2}(undef, (length(R),length(E)))
	SIG_G = Array{Any, 2}(undef, (length(R),length(E)))
	T_S_BAR = Array{Any, 2}(undef, (length(R),length(E)))
	T_G_BAR = Array{Any, 2}(undef, (length(R),length(E)))

	
	for r = 1:length(R), ld = 1:length(LD), e = 1:length(E)
		T_s = T_s_int[r,ld][e]
		T_g = T_g_int[r,ld][e]
		Z = Z_int[r,ld]
		
		N = length(Z)
		
		T_s_bar = line_avg(T_s, Z, 1, N)
		T_g_bar = line_avg(T_g, Z, 1, N)
		
		sig_s = sqrt(1/N * sum((T_s .- T_s_bar).^2))
		sig_g = sqrt(1/N * sum((T_g .- T_g_bar).^2))
		
		T_S_BAR[r,e] = T_s_bar
		T_G_BAR[r,e] = T_g_bar
		SIG_S[r,e] = sig_s
		SIG_G[r,e] = sig_g

		
		RMSE[r,e] = sqrt(1/N * sum((T_s .- T_g).^2))
		R_corr[r,e] = (sum((T_s .- T_s_bar).*(T_g .- T_g_bar))			
						)/(
				sqrt(sum((T_s .- T_s_bar).^2)) * sqrt(sum((T_g .- T_g_bar).^2))
						)
		NMB[r,e] = (T_s_bar - T_g_bar)/(T_g_bar)
		NMSD[r,e] = (sig_s - sig_g)/(sig_g)
					
	end
end

# ╔═╡ 78512ef2-878e-4980-afef-a464e45cf03c


# ╔═╡ fbb59f1a-3ca6-4a8a-a5af-6914a9bbe5ee
md"""
#### Plotting all vol. Effect Definitions
"""

# ╔═╡ da058b09-797e-42c2-8614-61df883a2ea0
md"""
##### (1) Based on inlet and outlet temperatures
"""

# ╔═╡ 7b080495-6d5c-455d-91ec-e243476d7001
begin
	p1 = 
		contourf(
		log2.(R), E, transpose(vol_eff[:,ld,:]),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, T<sub>g,out</sub>/T<sub>s,in</sub>",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
	# adding an arrow for direction to maximize efficiency
	plot!(
		x1,y1,
		linewidth = 2,
		line = (:dash, :arrow),
		color = :black,
		# label="Direction to max(η)",
		label="",
		)
	plot!(
		[-0.1, 0.046], 
		[0.9, 0.832],
		color = :black,
		linewidth = 2,
		label="")
	plot!(
		[-0.1, 0.20], 
		[0.9, 0.875],
		color = :black,
		linewidth = 2,
		label="")
	
	# plot!(
	# 	[-1, 3], 
	# 	[0.1, 0.9],
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="Direction to max(E<sub>vol</sub>)",
	# 	legend = (0.4,-0.3),
	# 	)
	# plot!(
	# 	[-1, -0.85], 
	# 	[0.1, 0.15],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-1, -0.8], 
	# 	[0.1, 0.12],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
end

# ╔═╡ 957c0e3f-b8be-4058-a781-d77242249312
begin
	p2 = 
		contourf(
		log2.(R), E, transpose(LMTD[:,ld,:]),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, LMTD (K)",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
		# adding an arrow for direction to maximize efficiency
		plot!(
			x1,y1,
			linewidth = 2,
			line = (:dash, :arrow),
			color = :black,
			label = false,
			)
		plot!(
			[-0.1, 0.046], 
			[0.9, 0.832],
			color = :black,
			linewidth = 2,
			label="")
		plot!(
			[-0.1, 0.20], 
			[0.9, 0.875],
			color = :black,
			linewidth = 2,
			label="")
end

# ╔═╡ 428978e3-8551-4423-840e-75fbe67da471


# ╔═╡ b2c0a4e7-565b-41ac-8f52-08640820ee39
md"""
##### (2) Based on power ratios
"""

# ╔═╡ c95fe5ff-5a5c-4164-bf1d-c53e9ad99b14
begin
	contourf(
		log2.(R), E, transpose(vol_eff_p1),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, Q<sub>abs,g</sub>/Q<sub>rad,loss</sub>",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
		# adding an arrow for direction to maximize efficiency
		plot!(
			x1,y1,
			linewidth = 2,
			line = (:dash, :arrow),
			color = :white,
			label = false,
			)
		plot!(
			[-0.1, 0.046], 
			[0.9, 0.832],
			color = :white,
			linewidth = 2,
			label="")
		plot!(
			[-0.1, 0.20], 
			[0.9, 0.875],
			color = :white,
			linewidth = 2,
			label="")
	
		# adding an arrow for direction to maximize volumetric effect defs.
		plot!(
			[-1, 3], 
			[0.1, 0.9],
			line = (:dash, :arrow),
			color = :white,
			linewidth = 2,
			label="",
			)
		plot!(
			[-1, -0.85], 
			[0.1, 0.15],
			color = :white,
			linewidth = 2,
			label="")
		plot!(
			[-1, -0.8], 
			[0.1, 0.12],
			color = :white,
			linewidth = 2,
			label="")
end

# ╔═╡ 09d41f66-871d-47cf-a852-b0bb794870c4
begin
	contourf(
		log2.(R), E, transpose(vol_eff_p2),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, Q<sub>rad,loss</sub>/Q<sub>BHS</sub>",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
		# adding an arrow for direction to maximize efficiency
		plot!(
			x1,y1,
			linewidth = 2,
			line = (:dash, :arrow),
			color = :white,
			label = false,
			)
		plot!(
			[-0.1, 0.046], 
			[0.9, 0.832],
			color = :white,
			linewidth = 2,
			label="")
		plot!(
			[-0.1, 0.20], 
			[0.9, 0.875],
			color = :white,
			linewidth = 2,
			label="")
end

# ╔═╡ bfccd459-c94e-4502-bad6-b4d2dd1694b3
begin
	contourf(
		log2.(R), E, transpose(vol_eff_p3),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, Q<sub>abs,g</sub>/Q<sub>BHS</sub>",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
	# adding an arrow for direction to maximize efficiency
	plot!(
		x1,y1,
		linewidth = 2,
		line = (:dash, :arrow),
		color = :white,
		label = false,
		)
	plot!(
		[-0.1, 0.046], 
		[0.9, 0.832],
		color = :white,
		linewidth = 2,
		label="")
	plot!(
		[-0.1, 0.20], 
		[0.9, 0.875],
		color = :white,
		linewidth = 2,
		label="")
end

# ╔═╡ d9268572-df0b-40ba-b3ef-7e908562bc34
begin
	contourf(
		log2.(R), E, transpose(vol_eff_p4),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, Q<sub>BHS</sub>/P<sub>ap</sub>",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
	# adding an arrow for direction to maximize efficiency
	plot!(
		x1,y1,
		linewidth = 2,
		line = (:dash, :arrow),
		color = :white,
		label = false,
		)
	plot!(
		[-0.1, 0.046], 
		[0.9, 0.832],
		color = :white,
		linewidth = 2,
		label="")
	plot!(
		[-0.1, 0.20], 
		[0.9, 0.875],
		color = :white,
		linewidth = 2,
		label="")
end

# ╔═╡ 8ec3b3a2-83cd-469e-b34e-2f2c20ca4eb0


# ╔═╡ 46c1346c-1cd9-41d1-8216-13b882b1b8bf
md"""
##### (3) Based on statistical evaluations of ∆T driving force
"""

# ╔═╡ ff999ca8-b613-43a4-9718-8545bf61f184
begin
	contourf(
		log2.(R), E, transpose(RMSE),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, RMSE",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		size = (440,335),	#for word doc.
			
		)
	
	# adding an arrow for direction to maximize efficiency
	plot!(
		x1,y1,
		linewidth = 2,
		line = (:dash, :arrow),
		color = :white,
		label = false,
		)
	plot!(
		[-0.1, 0.046], 
		[0.9, 0.832],
		color = :white,
		linewidth = 2,
		label="")
	plot!(
		[-0.1, 0.20], 
		[0.9, 0.875],
		color = :white,
		linewidth = 2,
		label="")

	# adding an arrow for direction to maximize volumetric effect defs.
	x2 = range(-1, 3, length = 20)
	# y2 = 0.1*log10.(x2 .+ 1) .+ 0.2
	y2 = -10 * exp.(-1.75 * (x2 .+ 2.6)) .+ 0.7
	# plot!(
	# 	[-1, 3], 
	# 	[0.1, 0.9],
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="",
	# 	)	#linear direction
	plot!(
		x2, 
		y2,
		line = (:dash, :arrow),
		color = :black,
		linewidth = 2,
		label="",
	)  		#curved direction
	plot!(
		[-1, -0.85], 
		[0.1, 0.15],
		color = :white,
		linewidth = 2,
		label="")
	plot!(
		[-1, -0.8], 
		[0.1, 0.12],
		color = :white,
		linewidth = 2,
		label="")
end

# ╔═╡ e3c2b3b1-6faf-482d-b880-583dbd177912
begin
	contourf(
		log2.(R), E, -transpose(R_corr),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, R",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
	
		# adding an arrow for direction to maximize efficiency
		plot!(
			x1,y1,
			linewidth = 2,
			line = (:dash, :arrow),
			color = :white,
			label = false,
			)
		plot!(
			[-0.1, 0.046], 
			[0.9, 0.832],
			color = :white,
			linewidth = 2,
			label="")
		plot!(
			[-0.1, 0.20], 
			[0.9, 0.875],
			color = :white,
			linewidth = 2,
			label="")
end

# ╔═╡ 354093b7-0205-4906-8c1f-8a6cdc4212b6
begin
	contourf(
		log2.(R), E, transpose(NMB),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, NMB",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
	
	# adding an arrow for direction to maximize efficiency
	plot!(
		x1,y1,
		linewidth = 2,
		line = (:dash, :arrow),
		color = :white,
		label = false,
		)
	plot!(
		[-0.1, 0.046], 
		[0.9, 0.832],
		color = :white,
		linewidth = 2,
		label="")
	plot!(
		[-0.1, 0.20], 
		[0.9, 0.875],
		color = :white,
		linewidth = 2,
		label="")

	# adding an arrow for direction to maximize volumetric effect defs.
	# plot!(
	# 	x2, 
	# 	y2,
	# 	line = (:dash, :arrow),
	# 	color = :black,
	# 	linewidth = 2,
	# 	label="",
	# )  		#curved direction
	plot!(
		[-1, -0.85], 
		[0.1, 0.15],
		color = :black,
		linewidth = 2,
		label="")
	plot!(
		[-1, -0.8], 
		[0.1, 0.12],
		color = :black,
		linewidth = 2,
		label="")
end

# ╔═╡ e5db2608-f53d-4fa0-9e72-510b954eee8e
begin
	contourf(
		log2.(R), E, transpose(NMSD),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, NMSD",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
		# adding an arrow for direction to maximize efficiency
		plot!(
			x1,y1,
			linewidth = 2,
			line = (:dash, :arrow),
			color = :white,
			label = false,
			)
		plot!(
			[-0.1, 0.046], 
			[0.9, 0.832],
			color = :white,
			linewidth = 2,
			label="")
		plot!(
			[-0.1, 0.20], 
			[0.9, 0.875],
			color = :white,
			linewidth = 2,
			label="")
	
		# adding an arrow for direction to maximize volumetric effect defs.
		plot!(
			[-1, 3], 
			[0.1, 0.9],
			line = (:dash, :arrow),
			color = :white,
			linewidth = 2,
			label="",
			)
		plot!(
			[-1, -0.85], 
			[0.1, 0.15],
			color = :white,
			linewidth = 2,
			label="")
		plot!(
			[-1, -0.8], 
			[0.1, 0.12],
			color = :white,
			linewidth = 2,
			label="")
end

# ╔═╡ ea8f5356-81b7-4bd1-896a-686aa331c45e
md"""
Changing the basis of normalization (observed values) for $NMSD$ and $NMB$ does not significantly change the trends. Only the values change.

"""

# ╔═╡ 9908d8b6-61bf-4f73-a0c3-7a4f7a3a9016


# ╔═╡ cde431b4-2a17-4507-9124-cf3f8b953abe
md"""
#### Assessing Bases for Vol. Effect Justification
"""

# ╔═╡ 40c4ff27-670a-4c80-b588-7919f48bc5bd
md"""
##### (1) Based on minimizing local solid temperature gradients
"""

# ╔═╡ a83af44b-ad18-48c4-a600-deb5f5f42f29
function line_grad(u, x)
	n = length(u)
	du_dx = zeros(n)
	for i = 1:n
		if i == 1
			du_dx[i] = (u[2]-u[1]) / (x[2]-x[1])
		elseif i == n
			du_dx[i] = (u[n]-u[n-1]) / (x[n]-x[n-1])
		else
			du_dx[i] = (u[i+1]-u[i-1]) / (x[i+1]-x[i-1])
		end
	end
		
	# u_avg = sum(u_dx)/(x[n] - x[m])
	return du_dx
end

# ╔═╡ 9631161e-0e1f-41f5-9c9f-92756fa5d736
# Based on minimizing local solid temperature gradients
begin
	
	IOTS = Array{Any, 2}(undef, (length(R),length(E)))
	IOTG = Array{Any, 2}(undef, (length(R),length(E)))

	MAXTS = Array{Any, 2}(undef, (length(R),length(E)))
	MAXTG = Array{Any, 2}(undef, (length(R),length(E)))	
	
	DT_S_DZ = Array{Any, 2}(undef, (length(R),length(E)))
	DT_G_DZ = Array{Any, 2}(undef, (length(R),length(E)))

	
	for r = 1:length(R), ld = 1:length(LD), e = 1:length(E)
		T_s = T_s_int[r,ld][e]
		T_g = T_g_int[r,ld][e]
		Z = Z_int[r,ld]
		
		N = length(Z)
		
		dT_s_dz = line_grad(T_s, Z)
		dT_g_dz = line_grad(T_g, Z)
		
		# display(plot(Z, dT_s_dz, title = "R=$(R[r])mm, ε=$(E[e])"))
		
		DT_S_DZ[r,e] = line_avg(dT_s_dz, Z, 1, N)
		DT_G_DZ[r,e] = line_avg(dT_g_dz, Z, 1, N)
		
		IOTS[r,e] = maximum(T_s) - minimum(T_s)
		IOTG[r,e] = maximum(T_g) - minimum(T_g)
					
	end
end

# ╔═╡ ce1da142-ba57-489f-8762-fb0be86acc30
md"""
##### (2) Based on maximizing solid temperature uniformity
"""

# ╔═╡ 6c6a8bbe-d011-4a77-80b7-c211abd29b2f
md"""
Uniformity can be expressed using the standard deviation in the axial solid temperature

"""

# ╔═╡ bf1196cd-e33c-4aaa-8ed3-294afe7be6e5
md"""
##### (3) Based on maximizing exergy efficiency
"""

# ╔═╡ 994fef1c-56eb-4410-8009-caf9a47c111d


# ╔═╡ 1b43324f-e9ca-473e-b093-7bf24c93f426
begin
	contourf(
		log2.(R), E, transpose(SIG_S./maximum(SIG_S)),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Solid Temp. Non-Uniformity, σ<sub>s</sub>",# <br>(L/D = $(LD[ld]))",
		
		color = cgrad(:thermal, rev = true),
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
	
			# adding an arrow for direction to maximize efficiency
			plot!(
				x1,y1,
				linewidth = 2,
				line = (:dash, :arrow),
				color = :white,
				label = false,
				)
			plot!(
				[-0.1, 0.046], 
				[0.9, 0.832],
				color = :white,
				linewidth = 2,
				label="")
			plot!(
				[-0.1, 0.20], 
				[0.9, 0.875],
				color = :white,
				linewidth = 2,
				label="")
end

# ╔═╡ 8b20cb06-7688-4d69-9d04-9fea0986ffbd
begin
	contourf(
		log2.(R), E, transpose(DT_S_DZ),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Avg. Gradient, <i>dT<sub>s</sub> / dz</i> (K/mm)",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
			# adding an arrow for direction to maximize efficiency
			plot!(
				x1,y1,
				linewidth = 2,
				line = (:dash, :arrow),
				color = :white,
				label = false,
				)
			plot!(
				[-0.1, 0.046], 
				[0.9, 0.832],
				color = :white,
				linewidth = 2,
				label="")
			plot!(
				[-0.1, 0.20], 
				[0.9, 0.875],
				color = :white,
				linewidth = 2,
				label="")
end

# ╔═╡ a9737cdc-0f2b-4b0e-92a0-b0fb463563ca


# ╔═╡ b6a330b6-9d9c-4f63-83ad-516873f99057
md"""
##### (4) Based on maximizing driving force uniformity
"""

# ╔═╡ 722da3ba-286d-4d0b-bf6e-598e623afe28
σ_DF = std.(
	[T_s_int[r,1][e] .- T_g_int[r,1][e]  for r = 1:length(R), e = 1:length(E)]
)#./SIG_G

# ╔═╡ 0b56fb23-ad3d-4d97-93d4-f7433beca6da
begin
	contourf(
		# log2.(R), E, transpose(σ_DF./maximum(σ_DF)),
		log2.(R), E, transpose(σ_DF),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(vol_eff[:,1,:]),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Emissivity",
		title = "Volumetric Effect, σ<sub>DF</sub>",# <br>(L/D = $(LD[ld]))",
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
		)
			# adding an arrow for direction to maximize efficiency
			plot!(
				x1,y1,
				linewidth = 2,
				line = (:dash, :arrow),
				color = :white,
				label = false,
				)
			plot!(
				[-0.1, 0.046], 
				[0.9, 0.832],
				color = :white,
				linewidth = 2,
				label="")
			plot!(
				[-0.1, 0.20], 
				[0.9, 0.875],
				color = :white,
				linewidth = 2,
				label="")
end

# ╔═╡ b73d7ae7-75ff-44da-bf34-21b3ab3b36e9


# ╔═╡ 41277ca6-db9f-4481-98e4-fb44674014e5
md"""
#### Plotting temperature profiles along line that maximizes η
"""

# ╔═╡ 9e40ab4c-64b7-4ddb-8bdb-a2796bfec9ed
E_int = reverse([0.185, 0.415, 0.646, 0.877]); R_int = reverse([8, 4, 2, 1])

# ╔═╡ e853f4ec-7cf4-48e1-bd34-ef57ae1b4772
cat([vec(T_s_int[2+1,ld][e]) for e =1:length(E)]; dims = 2)

# ╔═╡ 157f8136-29c4-44b4-ad9e-e5aa1c291cb2
begin
	T_s_inter = Array{Any, 2}(undef, (length(R_int),length(LD)))
	T_g_inter = Array{Any, 2}(undef, (length(R_int),length(LD)))
	
	for r = 1:length(R_int), ld = 1:length(LD)
		T_s = cat([vec(T_s_int[r+1,ld][e]) for e =1:length(E)]; dims = 2)
		T_s = transpose([T_s[e][z] for e =1:length(E), z = 1:length(Z_int[r+1])])
		T_s_inter[r,ld] = interpolate_lin(T_s, E, E_int[r])
		
		T_g = cat([vec(T_g_int[r+1,ld][e]) for e =1:length(E)]; dims = 2)
		T_g = transpose([T_g[e][z] for e =1:length(E), z = 1:length(Z_int[r+1])])
		T_g_inter[r,ld] = interpolate_lin(T_g, E, E_int[r])
	end
end

# ╔═╡ 4a39856f-87ae-472e-becb-4d2e4a911d61


# ╔═╡ 6aa20591-3d86-441b-a480-4cc54dcb79dd
md"""
#### Plotting temperature profiles along line that maximizes volumetric effect definitions
"""

# ╔═╡ 9eb8ea54-25de-4690-b1ab-f73435049442
m_vol = (0.9 - 0.1)/(3 - -1)

# ╔═╡ d9543116-23a2-495c-833e-eab8e1b378e1
R_vol = [1, 2, 4, 8]

# ╔═╡ 94430cd4-ec29-4b33-9f86-e3d8b36b6729
# E_vol = m_vol*(log2.(R_vol) .- 3) .+ 0.9
E_vol = E[1:end-1]

# ╔═╡ 2e210e9a-e7d2-4032-bdac-da50ab1599f2
begin
	T_s_vol = Array{Any, 2}(undef, (length(R_vol),length(LD)))
	T_g_vol = Array{Any, 2}(undef, (length(R_vol),length(LD)))
	
	for r = 1:length(R_vol), ld = 1:length(LD)
		T_s = cat([vec(T_s_int[r+1,ld][e]) for e =1:length(E)]; dims = 2)
		T_s = transpose([T_s[e][z] for e =1:length(E), z = 1:length(Z_int[r+1])])
		T_s_vol[r,ld] = interpolate_lin(T_s, E, E_vol[r])
		
		T_g = cat([vec(T_g_int[r+1ld][e]) for e =1:length(E)]; dims = 2)
		T_g = transpose([T_g[e][z] for e =1:length(E), z = 1:length(Z_int[r+1])])
		T_g_vol[r,ld] = interpolate_lin(T_g, E, E_vol[r])
	end
end

# ╔═╡ af7bb812-33ab-4b15-9ac6-baea5b3b81b2
colors_max_eta = reverse(permutedims(palette(:reds, 4)[:]))

# ╔═╡ ad09b9cb-2fe1-47f6-abcc-1ac66aaae37e
begin
	# r1 = 2
	# ld1 = 1
	plot(
			[Z_int[r+1,ld1] for r = 1:length(R_int)],

			[T_g_inter[r,ld1] for r = 1:length(R_int)],
		color = colors_max_eta,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :dash,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_int[r]) mm" for r = 1:length(R_int)]),

		# label = false
		)

	plot!(
			[Z_int[r+1,ld1] for r = 1:length(R_int)],

			[T_s_inter[r,ld1] for r = 1:length(R_int)],
		color = colors_max_eta,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :solid,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_int[r]) mm" for r = 1:length(R_int)]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Axial Temperature Profiles for <i>max(η)</i>",
		legend = :bottomright,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 10,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (624,477)
		)

end

# ╔═╡ d01fa8da-857d-4436-96ad-4eddc6e11c23
colors_max_vol = reverse(permutedims(palette(:blues, 4)[:]))

# ╔═╡ 21af0363-49a6-4b71-b303-4da665444174
begin	
	plot(
			[Z_int[r+1,ld1] for r = 1:length(R_int)],

			[T_g_vol[r,ld1] for r = 1:length(R_int)],
		color = colors_max_vol,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :dash,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_int[r]) mm" for r = 1:length(R_int)]),

		# label = false
		)

	plot!(
			[Z_int[r+1,ld1] for r = 1:length(R_int)],

			[T_s_vol[r,ld1] for r = 1:length(R_int)],
		color = colors_max_vol,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :solid,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_int[r]) mm" for r = 1:length(R_int)]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Axial Temperature Profiles for <i>max(E<sub>vol</sub>)</i>",
		legend = :bottomright,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 10,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (624,477)
		)

end

# ╔═╡ cbba62c9-4b99-4724-a283-d3f92c162aa1


# ╔═╡ f2ebf95d-1dc7-4d43-b569-3c4619019bcc
md"""
#### Combining temperature profile plots
"""

# ╔═╡ 36f7c1be-1ba3-4b5e-bae7-2f6ad5dd6c46
md"""
Darker color values ≡ Greater Value for $η$ and $E_{vol}$
"""

# ╔═╡ 91fde27d-92b3-4dda-b4ac-9fe5aa49cac1
begin
	plot(
			[Z_int[(r+1),ld1] for r = 1:length(R_int)],

			[T_g_inter[r,ld1] for r = 1:length(R_int)],
		color = colors_max_eta,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :dash,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_int[r]) mm" for r = 1:length(R_int)]),

		# label = false
		)

	plot!(
			[Z_int[(r+1),ld1] for r = 1:length(R_int)],

			[T_s_inter[r,ld1] for r = 1:length(R_int)],
		color = colors_max_eta,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :solid,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_int[r]) mm" for r = 1:length(R_int)]),
		
		)
	
	plot!(
			[Z_int[(r+1),ld1] for r = 1:length(R_int)],

			[T_g_vol[r,ld1] for r = 1:length(R_int)],
		color = colors_max_vol,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :dash,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_vol[r]) mm" for r = 1:length(R_int)]),

		# label = false
		)

	plot!(
			[Z_int[(r+1),ld1] for r = 1:length(R_int)],

			[T_s_vol[r,ld1] for r = 1:length(R_int)],
		color = colors_max_vol,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :solid,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_vol[r]) mm" for r = 1:length(R_int)]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Comparison of Temperature Profiles for <i>max(E<sub>vol</sub>)</i> vs. <i>max(η)</i>",
		legend = (1.0,0.9),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 10,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (624,477)
		)

end

# ╔═╡ 36f94bae-3200-413a-b1fc-fc22e1413d8c
md"""

The intersection point up till which increasing volumetric effect increases efficiency is seen somewhere around $R_{ch} = 2 [mm]$
"""

# ╔═╡ 9dec2b52-79ea-4588-ac89-237d3fa5b302


# ╔═╡ 16ebda6d-7bb5-426d-9e0d-a792c6202eb8
md"""
### Plotting _η_ vs. Volumetric Effect Definitions along the lines the Maximize _η_
"""

# ╔═╡ afec87c3-dbd2-4671-bfd2-3442a39cf4e4
E_vol2 = [T_g_inter[r,ld][end] / T_s_inter[r,ld][1] for r = 1:length(R_int), ld = 1:length(LD)]

# ╔═╡ 573f9fb9-c5e3-4657-9aca-fa2530402ba9
begin
	dT1 = [T_s_inter[r,ld][1] - T_g_inter[r,ld][1] for r = 1:length(R_int), ld = 1:length(LD)]
	
	dT2 = [T_s_inter[r,ld][end] - T_g_inter[r,ld][end] for r = 1:length(R_int), ld = 1:length(LD)]
	
	LMTD2 = (dT1 .- dT2)./log.(dT1 ./ dT2)
	""
end

# ╔═╡ 4a415610-63c7-45f3-b805-cf7fc0f4ae41
σ_DF2 = std.(
	[T_s_inter[r,1] .- T_g_inter[r,1]  for r = 1:length(R_int)]
)#./SIG_G

# ╔═╡ f0fc001e-7c67-4bea-9e0c-44949f023e9c
function interpolate_lin2(U, x, x_int)
	#U = known 2D property array
	#x = levels at which 'U' is known (y is assumed the same for both U and U_int)
	#x_int = levels at which 'U' is desired
	Y,X = size(U)
	U_int = zeros(Y, length(x_int))
	
	for r = 1:Y, e = 1:length(x_int)
		upper = findfirst(x .> x_int[e])
		
	 	if(x[upper - 1] .< x_int[e]) 
			lower = (upper - 1) 
		else 
			lower = (upper - 2)
		end
		slope = (U[r,lower] - U[r,upper])/(x[lower] - x[upper])
		
		U_int[r,e] = slope*(x_int[e] - x[upper]) + U[r,upper]
	end
	
	return U_int
end

# ╔═╡ c92a01ef-909b-4948-8e8a-74f459f435c0
begin
	# Q_abs_inter = interpolate_lin(Q_abs_tot[2:end,:], E, E_int)
	# Q_rad_loss_inter = interpolate_lin(Q_rad_loss_tot[2:end,:], E, E_int)
	# BHS_inter = interpolate_lin(BHS_tot[2:end,:], E, E_int)
	
	Q_abs_inter = Array{Any, 1}(undef, (length(R_int)))
	Q_rad_loss_inter = Array{Any, 1}(undef, (length(R_int)))
	BHS_inter =  Array{Any, 1}(undef, (length(R_int)))
	
	for r = 1:length(R_int)
		Q_abs = reshape(Q_abs_tot[r+1,:], (1,length(E)))
		Q_abs_inter[r] = interpolate_lin2(Q_abs, E, E_int[r])

		Q_rad = reshape(Q_rad_loss_tot[r+1,:], (1,length(E)))
		Q_rad_loss_inter[r] = interpolate_lin2(Q_rad, E, E_int[r])

		Q_BHS = reshape(BHS_tot[r+1,:], (1,length(E)))
		BHS_inter[r] = interpolate_lin2(Q_BHS, E, E_int[r])		
		
	end

	Q_abs_inter = [Q_abs_inter[r][1] for r = 1:length(R_int)]
	Q_rad_loss_inter = [Q_rad_loss_inter[r][1] for r = 1:length(R_int)]
	BHS_inter = [BHS_inter[r][1] for r = 1:length(R_int)]

end

# ╔═╡ 9c6ad57b-48d6-40b5-82db-fb64e1e2c858
begin
	Eff_inter = Array{Any, 1}(undef, (length(R_int)))
	
	for r = 1:length(R_int)
		Eff = reshape(Eff_Ap[r+1,:], (1,length(E)))
		Eff_inter[r] = interpolate_lin2(Eff, E, E_int[r])
	end

	Eff_inter = [Eff_inter[r][1] for r = 1:length(R_int)]

end

# ╔═╡ 0f5ef458-a444-4c86-b3fb-86b7b82b878c
begin
	surface(
		log2.(R), E, transpose(Eff_Ap),
		
		xlabel = "log<sub>2</sub>(R<sub>ch</sub> - mm)",
		ylabel = "Emissivity",
		zlabel = "Thermal Efficinecy",

		gridalpha = 0.8,
		
		zaxis = 0:0.2:1.0,
		#------------------------------------------------
		# tickfontsize = 10,
		# guidefont = 12,
		# # titlefont = 18,
		# fontfamily = "ComputerModern",
		# colorbar = :false,
		# size = 2 .* (320,250)
		#------------------------------------------------
		tickfontsize = 9,
	    guidefont = 11,
		# titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar = :false,
		size = 1.4 .* (389,300),	#for ppt
		)
	scatter3d!(log2.(R_int), E_int, Eff_inter, color = :black, label = "")
	# plot3d!(log2.(R_int), E_int, Eff_inter, color = :black, label = "")
end

# ╔═╡ 96e1d173-f24d-4682-bfe1-5871f372fec0
begin
	vol_eff_p12 = Q_abs_inter./Q_rad_loss_inter
	vol_eff_p22 = Q_rad_loss_inter./BHS_inter
	vol_eff_p32 = Q_abs_inter./BHS_inter
end

# ╔═╡ 83e96667-2543-42b4-826a-ddcf2d8ba252
# Based on statistical evaluations of ∆T
begin
	RMSE2 = Array{Any, 1}(undef, (length(R_int)))
	R_corr2 = Array{Any, 1}(undef, (length(R_int)))
	NMB2 = Array{Any, 1}(undef, (length(R_int)))
	NMSD2 = Array{Any, 1}(undef, (length(R_int)))

	SIG_S2 = Array{Any, 1}(undef, (length(R_int)))
	SIG_G2 = Array{Any, 1}(undef, (length(R_int)))
	T_S_BAR2 = Array{Any, 1}(undef, (length(R_int)))
	T_G_BAR2 = Array{Any, 1}(undef, (length(R_int)))

	
	for r = 1:length(R_int), ld = 1:length(LD)
		T_s = T_s_inter[r,ld]
		T_g = T_g_inter[r,ld]
		Z = Z_int[r+1,ld]
		
		N = length(Z)
		
		T_s_bar = line_avg(T_s, Z, 1, N)
		T_g_bar = line_avg(T_g, Z, 1, N)
		
		sig_s = sqrt(1/N * sum((T_s .- T_s_bar).^2))
		sig_g = sqrt(1/N * sum((T_g .- T_g_bar).^2))
		
		T_S_BAR2[r] = T_s_bar
		T_G_BAR2[r] = T_g_bar
		SIG_S2[r] = sig_s
		SIG_G2[r] = sig_g

		
		RMSE2[r] = sqrt(1/N * sum((T_s .- T_g).^2))
		R_corr2[r] = (sum((T_s .- T_s_bar).*(T_g .- T_g_bar))			
						)/(
				sqrt(sum((T_s .- T_s_bar).^2)) * sqrt(sum((T_g .- T_g_bar).^2))
						)
		NMB2[r] = (T_s_bar - T_g_bar)/(T_g_bar)
		NMSD2[r] = (sig_s - sig_g)/(sig_g)
					
	end
end

# ╔═╡ d1626a4b-4c56-4085-831f-01a4efc1e7dc
labels = ["dTₛ/dz", "σₛ" , "T<sub>g,out</sub>/T<sub>s,in</sub>" , "LMTD" ,"R_corr", "RMSE", "NMB"  ,"NMSD", "Q<sub>abs,g</sub>/Q<sub>rad,loss</sub>", "Q<sub>rad,loss</sub>/Q<sub>BHS</sub>", "Q<sub>abs,g</sub>/Q<sub>BHS</sub>", "σ_∆T" ] 


# ╔═╡ bd58561d-03c6-47cd-87ee-c80ba7b54601
begin
	plot(
		log2.(R_int),
		Eff_inter,

		# ylims = [0.55, 0.9],
		# xlim = (0.8,1.0),
				
		color = permutedims(palette(:tab10)[1:4]),

		ylabel = "Thermal Efficiency, η",
		xlabel = "log2(R)",	
		title = "Efficiency Corss-Section for Selected Designs",
		#--------------------------------------------
		markers = :true,		# markers = marker_value[1],
		# markeralpha = marker_value[2],
		# markersize = marker_value[3],
		xgrid = :none,
		linewidth = 2,
		# legend = :right,
		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 5*Plots.mm,
		legendfontsize = 12,
		tickfontsize = 14,
		guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		size = (446,340),	#for word doc.
		dpi = 400,

		)
	
	hline!(
		[maximum(Eff_inter)],		color = :black,
		linewidth = 0.5,
		line = :dash,
		label = "max(η)",
		)
end

# ╔═╡ 83520631-205e-4a3e-bff0-6bc33ef25d01
begin
	plot(
		Eff_inter,
		[
			# BHS_inter, 
			# Q_abs_inter, 
			# Q_rad_loss_inter,

			BHS_inter./BHS_inter.*100, 
			Q_abs_inter./BHS_inter.*100, 
			Q_rad_loss_inter./BHS_inter.*100, 
			
		],

		label = ["Q<sub>BHS</sub>" "Q<sub>abs,g</sub>" "Q<sub>rad,loss</sub>"],
		
		color = permutedims(palette(:tab10)[1:3]),

		xlims = [0.55, 0.9],

		xlabel = "Thermal Efficiency, η",
		ylabel = "Powers(W)",	
		#--------------------------------------------
		markers = :false,
		# markers = marker_value[1],
		# markeralpha = marker_value[2],
		# markersize = marker_value[3],
		xgrid = :none,
		linewidth = 2,
		# legend = :right,
		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 5*Plots.mm,
		legendfontsize = 12,
		tickfontsize = 14,
		guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		size = (446,340),	#for word doc.
		dpi = 400,

		)
	
	vline!(
		[maximum(Eff_inter)],
		color = :black,
		linewidth = 0.5,
		line = :dash,
		label = "max(η)",
		)
end

# ╔═╡ d9f767c2-d1b4-4998-a113-1524d3762ce8
begin
	plot(
		[
			vol_eff_p12 ./ maximum(vol_eff_p12), 
			vol_eff_p22 ./ maximum(vol_eff_p22), 
			vol_eff_p32 ./ maximum(vol_eff_p32)
		],
		Eff_inter,

		label = permutedims(labels[end-2:end]),
		
		color = permutedims(palette(:tab10)[1:3]),

		ylims = [0.55, 0.9],

		ylabel = "Thermal Efficiency, η",
		xlabel = "Volumetric Effect (norm)",	
		title = "Based on Power Ratios",
		#--------------------------------------------
		markers = :false,
		# markers = marker_value[1],
		# markeralpha = marker_value[2],
		# markersize = marker_value[3],
		xgrid = :none,
		linewidth = 2,
		# legend = :right,
		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 5*Plots.mm,
		legendfontsize = 12,
		tickfontsize = 14,
		guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		size = (446,340),	#for word doc.
		dpi = 400,

		)
	
	hline!(
		[maximum(Eff_inter)],
		color = :black,
		linewidth = 0.5,
		line = :dash,
		label = "max(η)",
		)
end

# ╔═╡ 6c76cc1a-ded6-4224-acee-f140e9fae43c
begin
	plot(
		[
			E_vol2,# ./ maximum(E_vol2), 
			LMTD2 ./ maximum(LMTD2), 
		],
		Eff_inter,

		label = permutedims(labels[3:4]),
		
		color = permutedims(palette(:tab10)[1:2]),

		ylims = [0.55, 0.9],
		
		ylabel = "Thermal Efficiency, η",
		xlabel = "Volumetric Effect (norm)",	
		title = "Based on Inlet & Outlet Temp.",
		#--------------------------------------------
		markers = :false,
		# markers = marker_value[1],
		# markeralpha = marker_value[2],
		# markersize = marker_value[3],
		xgrid = :none,
		linewidth = 2,
		# legend = :right,
		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 5*Plots.mm,
		legendfontsize = 12,
		tickfontsize = 14,
		guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		size = (446,340),	#for word doc.
		dpi = 400,

		)
	
	hline!(
		[maximum(Eff_inter)],		color = :black,
		linewidth = 0.5,
		line = :dash,
		label = "max(η)",
		)
end

# ╔═╡ e37aa4e5-0574-4e20-842c-c46df54e3ca0
begin
	plot(
		[
			E_vol2,# ./ maximum(E_vol2), 
		],
		Eff_inter,

		label = permutedims(labels[3:4]),
		
		color = permutedims(palette(:tab10)[1:2]),

		ylims = [0.55, 0.9],
		
		ylabel = "Thermal Efficiency, η",
		xlabel = "Volumetric Effect",	
		# title = "Based on Inlet & Outlet Temp.",
		# title = "(A)", titlelocation = :left,
		#--------------------------------------------
		xgrid = :none,
		linewidth = 2,
		# legend = :outerright,
		legend = (1.0, 0.7),
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 10*Plots.mm,
		left_margin = 5*Plots.mm,
		bottom_margin = 1*Plots.mm,
		right_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		legendfont = 10,
		titlefont = 16,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt
		dpi = 400,
		)
	
	hline!(
		[maximum(Eff_inter)],		color = :black,
		linewidth = 0.5,
		line = :dash,
		label = "max(η)",
		)
end

# ╔═╡ beab92a2-f8b8-440d-b2ba-3c88a55a3441
begin
	plot(
		[
			# R_corr2,
			# RMSE2 ./ maximum(RMSE2), 
			# # RMSE2,
			# NMB2, 
			# NMSD2,
			
			# R_corr2 ./ maximum(R_corr2), 
			RMSE2 ./ maximum(RMSE2), 
			NMB2 ./ maximum(NMB2), 
			NMSD2 ./ maximum(NMSD2), 
		],
		Eff_inter,

		# label = permutedims(labels[5:8]),
		label = permutedims(labels[6:8]),
		
		color = permutedims(palette(:tab10)[1:4]),

		ylims = [0.55, 0.9],
		xlims = [-1,1],
		
		ylabel = "Thermal Efficiency, η",
		xlabel = "Volumetric Effect (norm)",	
		title = "Based on Statistical Evaluations of ∆T",
		#--------------------------------------------
		markers = :false,
		# markers = marker_value[1],
		# markeralpha = marker_value[2],
		# markersize = marker_value[3],
		xgrid = :none,
		linewidth = 2,
		# legend = :right,
		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 5*Plots.mm,
		legendfontsize = 12,
		tickfontsize = 14,
		guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		size = (446,340),	#for word doc.
		dpi = 400,

		)
	
	hline!(
		[maximum(Eff_inter)],		color = :black,
		linewidth = 0.5,
		line = :dash,
		label = "max(η)",
		)
	
	vline!(
		[0],
		color = :black,
		linewidth = 0.5,
		line = :dot,
		label = "Target",
		)
end

# ╔═╡ 838076f5-04ab-4eca-8763-98a27b8637de
begin
	plot(
		[
			# σ_DF2, 
			σ_DF2 ./ maximum(σ_DF2), 

		],
		Eff_inter,

		ylims = [0.55, 0.9],
		xlim = (0.8,1.0),
		
		label = "σ<sub>∆T</sub>",
		
		color = permutedims(palette(:tab10)[1:4]),

		ylabel = "Thermal Efficiency, η",
		xlabel = "Volumetric Effect (norm)",	
		title = "Based on ∆T Non-Uniformity",
		#--------------------------------------------
		markers = :false,
		# markers = marker_value[1],
		# markeralpha = marker_value[2],
		# markersize = marker_value[3],
		xgrid = :none,
		linewidth = 2,
		# legend = :right,
		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 5*Plots.mm,
		legendfontsize = 12,
		tickfontsize = 14,
		guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		size = (446,340),	#for word doc.
		dpi = 400,

		)
	
	hline!(
		[maximum(Eff_inter)],		color = :black,
		linewidth = 0.5,
		line = :dash,
		label = "max(η)",
		)
end

# ╔═╡ 67e10d4c-94fa-4586-a138-fe3548d42468


# ╔═╡ ffedd213-09d3-414f-aa7a-619d80b6d748
md"""
### Plotting the design for the maximum of each definition over efficiency
"""

# ╔═╡ 4235aef9-9789-48c7-a548-fdceeebbd3ac
begin
	# target = Array{CartesianIndex{2}, 1}(undef, 13)

	# target[1] = findall(abs.(DT_S_DZ) .== minimum(abs.(DT_S_DZ)))[1]
	# target[2] = findall(SIG_S .== minimum(SIG_S))[1]
	# target[3] = findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1]
	# target[4] = findall(LMTD[:,1,:] .== maximum(LMTD[:,1,:]))[1]
	# target[5] = findall(abs.(R_corr) .== maximum(abs.(R_corr)))[1]
	# target[6] = findall(RMSE .== minimum(RMSE))[1]
	# target[7] = findall(NMB .== minimum(NMB))[1]
	# target[8] = findall(abs.(NMSD) .== minimum(abs.(NMSD)))[1]
	# target[9] = findall(vol_eff_p1 .== maximum(vol_eff_p1))[1]
	# target[10] = findall(vol_eff_p2 .== minimum(vol_eff_p2))[1]
	# target[11] = findall(vol_eff_p3 .== maximum(vol_eff_p3))[1]
	# target[12] = findall(vol_eff_p4 .== maximum(vol_eff_p4))[1]
	# target[13] = findall(σ_DF .== minimum(σ_DF))[1]

	# target = Array{CartesianIndex{2}, 1}(undef, 12)

	# target[1] = findall(abs.(DT_S_DZ) .== minimum(abs.(DT_S_DZ)))[1]
	# target[2] = findall(SIG_S .== minimum(SIG_S))[1]
	# target[3] = findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1]
	# target[4] = findall(LMTD[:,1,:] .== maximum(LMTD[:,1,:]))[1]
	# target[5] = findall(abs.(R_corr) .== maximum(abs.(R_corr)))[1]
	# target[6] = findall(RMSE .== minimum(RMSE))[1]
	# target[7] = findall(NMB .== minimum(NMB))[1]
	# target[8] = findall(abs.(NMSD) .== minimum(abs.(NMSD)))[1]
	# target[9] = findall(vol_eff_p1 .== maximum(vol_eff_p1))[1]
	# target[10] = findall(vol_eff_p2 .== minimum(vol_eff_p2))[1]
	# target[11] = findall(vol_eff_p3 .== maximum(vol_eff_p3))[1]
	# target[12] = findall(σ_DF .== minimum(σ_DF))[1]

	
	target = Array{CartesianIndex{2}, 1}(undef, 10)
	
	target[1] = findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1]
	target[2] = findall(LMTD[:,1,:] .== maximum(LMTD[:,1,:]))[1]
	target[3] = findall(RMSE .== minimum(RMSE))[1] #closest value to 0
	target[4] = findall(NMB .== minimum(NMB))[1] #closest value to 0
	target[5] = findall(abs.(NMSD) .== minimum(abs.(NMSD)))[1] #closest value to 0
	target[6] = findall(vol_eff_p1 .== maximum(vol_eff_p1))[1]
	target[7] = findall(vol_eff_p2 .== minimum(vol_eff_p2))[1]
	target[8] = findall(vol_eff_p3 .== maximum(vol_eff_p3))[1]
	target[9] = findall(vol_eff_p4 .== maximum(vol_eff_p4))[1]
	target[10] = findall(σ_DF .== minimum(σ_DF))[1]

end

# ╔═╡ c092187a-a148-4092-bae2-450bd9e1d3d4
labels3 = cat(labels[3:4], labels[6:11], ["Q<sub>BHS,g</sub>/P<sub>ap</sub>"], "σ<sub>∆T</sub>"; dims = 1)

# ╔═╡ bca491f5-9446-4fcc-af0e-ca5835c14286
labels3[8] = "η<sub>conv</sub>"; labels3[9] = "η<sub>opt</sub>"

# ╔═╡ 85097c68-d655-4ecd-ba48-aadeb4f3368d
labels2 = cat(copy(labels[1:11]), ["Q<sub>BHS</sub>/P<sub>ap</sub>"], [labels[12]]; dims = 1)

# ╔═╡ dee5a117-26d1-4371-a765-be6d0eab8037
labels2[13] = "σ<sub>∆T</sub>";  labels2[5] = "R<sub>corr</sub>";

# ╔═╡ 0707634c-6a9e-465e-96be-1598a9336da5
# marker_shapes = [:hex  :hex  :square  :square  :diamond  :diamond  :diamond  :diamond  :circle  :circle  :circle  :circle  :utriangle]

# marker_shapes = [:hex  :hex  :square  :square  :diamond  :diamond  :diamond  :diamond  :circle  :circle  :circle  :utriangle]

marker_shapes = [:circle  :circle  :square  :square  :square  :diamond  :diamond  :diamond  :diamond  :utriangle]

# ╔═╡ 8fe1a290-24fa-416e-949c-11448b9a1e71
begin
	contourf(
		log2.(R), E, transpose(Eff_Ap),
		# log2.(R), E_base, transpose(Eff_Ap_int),
		#---------------------------------------------
		xlabel = "log<sub>2</sub>(Radius - mm)",
		ylabel = "Emissivity",
		
		# title = "Thermal Efficiency, η",
		# title = "(A)", titlelocation = :left,

		colorbar_title = " η<sub> thermal<sub>",
		colorbar_titlefontsize = 16,
	    colorbar_title_location = :right,           # also :left or :right
	    colorbar_titlefontfamily = :match,
		# color = :grays,
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		framestyle = :box,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		titlefont = 16,
		# colorbar_title = "η",
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt	
		)
	
	scatter!(
		[log2(R[target[1][1]])],
		[E[target[1][2]]],

		# yaxis = :false,
		# xaxis = :false,
		
		# label = permutedims(labels2),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = "max($(labels3[1]))",
		
		color = :gray,
		# color = colors2,
		markershape = marker_shapes,
		markerstrokewidth = 1,
		markerstrokecolor  = :white,
		markersize = 8,

		legendfontsize = 13,
		legend = (1.3,0.9)
		)

	scatter!(
		[log2(R[2])],
		[E[5]],

		marker = :xcross,
		color = :black,
		markersize = 5,
		label = "max(η<sub>th</sub>)",
		legend = false,
		)
end

# ╔═╡ 2dd60715-9a7d-4892-aa45-c695b143a6da
begin
	contour(
		log2.(R), E, transpose(BHS_int_tot./BHS_front_tot),
		levels = range(1,80, step=4),
		# log2.(R), E_base, transpose(Eff_Ap_int),
		#---------------------------------------------
		xlabel = "log<sub>2</sub>(Radius - mm)",
		ylabel = "Emissivity",
		
		# title = "Thermal Efficiency, η",
		# title = "(A)", titlelocation = :left,

		colorbar_title = "BHS<sub>int</sub>/BHS<sub>front</sub>",
		colorbar_titlefontsize = 16,
	    colorbar_title_location = :right,           # also :left or :right
	    colorbar_titlefontfamily = :match,
		# color = :grays,
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		framestyle = :box,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		titlefont = 16,
		# colorbar_title = "η",
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (389,300),	#for ppt
		# size = (624,477), #for ppt	
		)
	
	scatter!(
		[log2(R[target[1][1]])],
		[E[target[1][2]]],

		# yaxis = :false,
		# xaxis = :false,
		
		# label = permutedims(labels2),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = "max($(labels3[1]))",
		
		color = :gray,
		# color = colors2,
		markershape = marker_shapes,
		markerstrokewidth = 1,
		markerstrokecolor  = :white,
		markersize = 8,

		legendfontsize = 13,
		legend = (1.3,0.9)
		)

	scatter!(
		[log2(R[2])],
		[E[5]],

		marker = :xcross,
		color = :black,
		markersize = 5,
		label = "max(η<sub>th</sub>)",
		legend = false,
		)
end

# ╔═╡ 7b11fbbd-f3a9-42b0-aea7-a5c5f63b95e8
# colors2 = permutedims(palette([:red, :blue], 12)[:])
# colors2 = permutedims(palette(:solar, 13)[:])
# colors2 = permutedims(palette(:solar, 12)[:])
colors2 = permutedims(palette(:solar, 10)[:])

# ╔═╡ 4893e9f2-0482-4f51-8581-f015135e09e7
begin
	contourf(
		log2.(R), E, transpose(Eff_Ap),
		# log2.(R), E_base, transpose(Eff_Ap_int),
		#---------------------------------------------
		xlabel = "log<sub>2</sub>(Radius - mm)",
		ylabel = "Emissivity",
		
		# title = "Thermal Efficiency, η",
		# title = "(A)", titlelocation = :left,

		colorbar_title = " η<sub> thermal<sub>",
		colorbar_titlefontsize = 16,
	    colorbar_title_location = :right,           # also :left or :right
	    colorbar_titlefontfamily = :match,
		color = :grays,
		#--------------------------------------------
		# Size formatting for ppt
		# top_margin = 12*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		framestyle = :box,
		top_margin = 15*Plots.mm,
		left_margin = 7*Plots.mm,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		# size = (440,335),	#for word doc.
			
		)
	
	scatter!(
		[[log2(R[target[i][1]])] for i = 1:length(target)],
		[[E[target[i][2]]] for i = 1:length(target)],

		# yaxis = :false,
		# xaxis = :false,
		
		# label = permutedims(labels2),		
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = permutedims(labels3),		
		
		color = colors2,
		markershape = marker_shapes,
		markerstrokewidth = 1,
		markerstrokecolor  = :white,
		markersize = 8,

		legendfontsize = 13,
		legend = (1.3,0.9)
		)

	scatter!(
		[log2(mean(R[1:2]))],
		[E[5]],

		marker = :xcross,
		color = :black,
		markersize = 5,
		label = "η<sub>max</sub>",
		)
end

# ╔═╡ 9e43acb4-e0a4-4658-a4f7-a4154d60ea89
begin
	bar(
		# [[labels2[i]] for i = 1:length(target)],
		# [[labels2[i]] for i = cat(1:11, [13]; dims=1)],
		[[labels3[i]] for i = 1:length(labels3)],
		[[Eff_Ap[target[i]]] for i = 1:length(target)],

		ylabel = "Thermal Efficiency, η",
		ylim = [0.2, 1.0],
		xgrid = :none,
		
		xrotation = 90,
		bottom_margin = -120*Plots.mm,
	
		color = colors2,
		legend = false,

		# title = "(B)", titlelocation = :left, 
		topmargin = 15*Plots.mm,
		#--------------------------------------------
		# Size formatting for word
		tickfontsize = 14,
		guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440, 400),	#for word doc.
				
	)
	
	hline!(
		[maximum(Eff_Ap)],
	
		label = "η<sub>max</sub>",
		
		color = :black,
		line = :dash,
		lineswidth = 2,
		
		# annotations = (14.5, maximum(Eff_Ap)*1.05, Plots.text("<i> η<sub>max</sub> </i>", :left, "ComputerModern")),	
		annotations = (11.5, maximum(Eff_Ap)*1.06, Plots.text("<i> η<sub>max</sub> </i>", :left, "ComputerModern")),	
	)
end

# ╔═╡ 4f1e689a-f775-4de9-83ba-9515fc1db364


# ╔═╡ 05c33ad2-e2ae-4916-a1ec-6c885a424c1e
begin
	plot(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_g_int[target[i][1],1][target[i][2]] for i = 1:length(target)],
		
		color = colors2,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = permutedims(labels3),

		)
	

	plot!(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_s_int[target[i][1],1][target[i][2]]  for i = 1:length(target)],
		
		color = colors2,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = permutedims(labels3),

		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Axial Temperature Profiles for <i>max(E<sub>vol</sub>)</i> <br> (Base Case)",
		legend = (1.0, 1.0),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 10,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (624,477)
		)

	plot!(
		Z_int[1,1],

		[T_s_int[1,1][5], T_g_int[1,1][5]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		
		)
	

end

# ╔═╡ 222030b4-7074-4d42-aa09-089e5f54b644
findall(Eff_Ap .== maximum(Eff_Ap))

# ╔═╡ 050b289a-c636-4944-b411-ce272deb1e2c
md"""
selecting relevant definitions
"""

# ╔═╡ cd78642a-8ab5-41b6-a5ac-ee07a4ac9f45
begin
	plot(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_g_int[target[i][1],1][target[i][2]] for i = 1:length(target)],
		
		color = colors2,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = permutedims(["" for i = 1:length(target)]),
		)
	

	plot!(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_s_int[target[i][1],1][target[i][2]]  for i = 1:length(target)],
		
		color = colors2,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),

		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Axial Temperature Profiles for <i>max(E<sub>vol</sub>)</i> <br> (Base Case)",
		legend = (1.0, 1.0),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (460,240).*1.5,
		)

	plot!(
		Z_int[1,1],

		[T_s_int[1,1][5], T_g_int[1,1][5]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		
		)
	

end

# ╔═╡ 32fe31d3-d242-4604-ae7c-7e99009eadd5
md"""
Separating grouped definitions
"""

# ╔═╡ 7d783db2-a669-4a09-b7d6-b21a26b94326
# colors3 = permutedims(palette(:Spectral_4)[:])
colors3 = permutedims(palette(:Dark2_4)[:])

# ╔═╡ f4f3dbbb-f202-4c8d-bd02-99a3d6d932b2
begin
	plot(
		[Z_int[target[i][1],1] for i = 1:2],

		[T_g_int[target[i][1],1][target[i][2]] for i = 1:2],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 1:2]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i = 1:2],

		[T_s_int[target[i][1],1][target[i][2]]  for i = 1:2],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3[1:2]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		title = "(A)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (320,250).*1.5,
		)

	plot!(
		Z_int[1,1],

		[T_s_int[1,1][5], T_g_int[1,1][5]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		
		)
	

end

# ╔═╡ 3adba9e9-fa6b-4e8a-891c-82bbd9f48f45
begin
	plot(
		[Z_int[target[i][1],1] for i = 6:9],

		[T_g_int[target[i][1],1][target[i][2]] for i = 6:9],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 6:9]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i = 6:9],

		[T_s_int[target[i][1],1][target[i][2]]  for i = 6:9],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3[6:9]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		title = "(B)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (320,250).*1.5,
		)

	plot!(
		Z_int[1,1],

		[T_s_int[1,1][5], T_g_int[1,1][5]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		
		)
	

end

# ╔═╡ 16484697-c3b3-4509-8640-7e8121df9d23
begin
	plot(
		[Z_int[target[i][1],1] for i = 3:5],

		[T_g_int[target[i][1],1][target[i][2]] for i = 3:5],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 3:5]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i = 3:5],

		[T_s_int[target[i][1],1][target[i][2]]  for i = 3:5],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3[3:5]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		title = "(C)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (320,250).*1.5,
		)

	plot!(
		Z_int[1,1],

		[T_s_int[1,1][5], T_g_int[1,1][5]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		
		)
	

end

# ╔═╡ 58073901-4f8d-4b37-8a0a-c6eeb32d3461
begin
	plot(
		[Z_int[target[i][1],1] for i in [10]],

		[T_g_int[target[i][1],1][target[i][2]] for i in [10]],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i in [10]]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i in [10]],

		[T_s_int[target[i][1],1][target[i][2]]  for i in [10]],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = labels3[10],
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		title = "(D)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (320,250).*1.5,
		)

	plot!(
		Z_int[1,1],

		[T_s_int[1,1][5], T_g_int[1,1][5]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		
		)
	

end

# ╔═╡ aed47d09-bd50-4ce5-8771-f700874efddc
begin
	plot(
		Z_int[target[1][1],1]/maximum(Z_int[target[1][1],1]),

		T_g_int[target[1][1],1][target[1][2]],
		
		color = :gray,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims([""]),		
		)

	plot!(
		Z_int[target[1][1],1]/maximum(Z_int[target[1][1],1]),

		T_s_int[target[1][1],1][target[1][2]],
		
		color = :gray,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = "max($(labels3[1]))",
		
		xlabel = "Normalized Receiver Length",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		# title = "(A)", titlelocation = :left, topmargin = 8*Plots.mm,

		yaxis = [300, 600, 900],
		grid = :none,
		framestyle = :box,
		
		legend = (0.58, 0.66),
		# legend = :right,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
		tickfontsize = 12,
		legendfont = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (320,250).*1.2,
		)

	plot!(
		Z_int[1,1]/maximum(Z_int[1,1]),

		[T_s_int[1,1][5], T_g_int[1,1][5]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		
		)
	

end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Plots = "~1.22.3"
PlutoUI = "~0.7.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c2178cfbc0a5a552e16d097fae508f2024de61a3"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.59.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "ef49a187604f865f4708c90e3f431890724e9012"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.59.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "2537ed3c0ed5e03896927187f5f2ee6a4ab342db"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.14"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "cfbd033def161db9494f86c5d18fbf874e09e514"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.3"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "d1fb76655a95bf6ea4348d7197b22e889a4375f4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.14"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─838cbfe0-f0ad-11eb-0146-a989eeaafda6
# ╠═1e060c7e-886f-4208-b68e-09f2b08195fa
# ╠═ad798898-61ea-4188-9c5f-d258e10cc47e
# ╠═63b2aa9e-bdbf-4670-bc0a-ff86013a20e7
# ╠═a226d405-d910-4112-8b6f-a153420a38ca
# ╠═aa57714a-6ba3-48ec-afcd-31b44845a4e8
# ╠═6317f5cc-bbf9-4c22-a813-c500a538d380
# ╠═cdfd8cb1-7602-45ad-abc4-8d9ed53eee19
# ╠═1f6244a6-09a4-4e7c-9468-8aa6ff376070
# ╟─5470fdd0-21cb-495c-9c69-ed9ba537e9bf
# ╟─1f1828f3-2f72-4da4-9058-7729e6316bfb
# ╠═7c3fc9a5-3f6b-428e-a32e-9c44c92fdb08
# ╠═ba014dd3-7b46-444e-8337-ba096d164127
# ╟─61c6bcf4-cf2c-400c-9922-32c6eef9dec0
# ╠═e74945d2-06f7-4d36-bb27-37e5e10095ba
# ╠═464815e3-6f96-4b1d-8a22-a88afcd7aafb
# ╟─3186cee0-226b-48d8-bad1-aba11756d9d2
# ╟─5c0de861-b973-4bf1-9b76-c2af3ce1a21b
# ╠═c037aea2-371d-45d1-b6e3-fb51aad8cb7c
# ╠═ce5f13f3-ec44-4cf1-8cba-fd9c331c7ca0
# ╠═f82863f3-d2dc-46a8-9cb3-0b1a129a0508
# ╠═b53b7e81-b817-42d2-8170-74affd2f5e12
# ╠═1124a3c4-2b53-448a-9d63-278da7a2d834
# ╠═722d42d3-5c99-4e2f-826b-892006fac5e3
# ╠═da408ddb-5b01-4376-ad34-c19c3177910d
# ╠═7b557345-bafa-44fb-b515-cd651a66ef7e
# ╟─115b38d9-1c87-488d-8905-fd04ddc71682
# ╠═f57af54b-3a7f-4337-982a-8c1b21554694
# ╠═3f30e71e-c0ac-4c9a-9ee7-e3707400118c
# ╠═5ba04d8b-0c65-46dc-b080-51413388dc30
# ╟─31ee1595-be66-4d8a-9be1-96b2514dcc0c
# ╠═5bf03252-b019-431a-9106-7a1df1f56806
# ╟─2c0b3dea-5df3-450d-b5a9-8a89930bc417
# ╟─42efe36f-025c-44ad-8890-a32a86f7b571
# ╟─fa856550-9cd1-496b-a2b5-2cf79c22dd33
# ╟─69c19e51-1271-443e-8a36-30268e7afdf4
# ╟─a87b183a-aebb-4398-920d-0f4e4e8a51a9
# ╟─1a9f5a82-5ce7-470d-a70e-67c766f932fd
# ╠═3c5eaa90-ec96-4757-9356-956a98840aef
# ╟─44cec87f-537d-426d-8089-934210aee5b3
# ╟─03f25a88-efb5-40bb-9c13-4ca2b8b97c14
# ╟─2e280ca0-1971-42f4-9f92-c04d9c464847
# ╟─ebe6df24-c623-4ffd-9857-40110fb0336d
# ╠═aaa5252e-0714-476d-8267-f13d8a1f9e0d
# ╠═0593c86b-fb0f-40c2-a38f-9c7330c7a484
# ╟─f226cd78-4380-4340-a3b9-9d7bb1940abb
# ╟─e09d1bcb-66f5-4e1a-bc34-1f1dc936335a
# ╟─0cdb6335-31ce-4ee2-9e81-b3d025a11179
# ╟─1bc2fdd5-36cf-4cb8-a388-9ac949ceabf2
# ╟─69f256c3-c603-4797-9571-73993d519580
# ╟─b78bdbdd-ba48-4e44-b19b-19caffe2b307
# ╟─e5300a51-d793-4fe2-88b4-8a117924f91b
# ╟─dec3866e-881e-42a0-9236-e65fd6fd8599
# ╟─8cb6d0d1-9301-419e-a242-817b0b718ae6
# ╟─2caf1b17-51cf-4d96-93d3-948104c92fa5
# ╟─2235d2e4-a48c-4102-a947-b40c0b1bdeb2
# ╟─2a730b37-d2d3-4b32-9c47-24fd5b122d95
# ╟─76f0b4fd-dc94-499d-8d41-c58cc2d3ff2d
# ╟─ee777a68-05da-4696-aa3e-d4fcddd3367a
# ╟─90b3c5be-feeb-4e94-b96f-49622ee88ec6
# ╟─0a38ee74-9d7a-451f-a171-9aabdbc2939c
# ╠═238739ff-d8c0-4089-a87b-a435e2c08325
# ╠═736b6056-df0f-4f5b-8b8a-736bee68370f
# ╠═416b7566-578e-4606-8fee-6305e3f5a11d
# ╠═410622b9-8213-4424-a09f-ff5eff7487c9
# ╟─c917cb94-36a6-49f1-bf2c-7681f6172310
# ╠═b87236e7-2e21-4d78-8f39-31292cae6218
# ╟─7e095764-2bb6-44ed-afc7-38f5cdff4d92
# ╟─31a18798-98f7-4d0b-935d-be8610eae64e
# ╟─538ae313-d93f-4c38-b049-491d21e15ecc
# ╟─1672b8ce-a4a5-44c5-aa7a-54af60a29451
# ╟─549fbc4f-b840-4a49-b0dc-2174d9d8468d
# ╟─7a6b554f-d315-4184-93b0-879b6e987e41
# ╟─f005f378-0f08-4e0a-902d-a49ba445643f
# ╟─db524b4c-336b-4a78-983e-4e1a72798d88
# ╟─b3509e00-6bfe-4ea2-9552-90d6aacb19a1
# ╟─622c7a82-f772-4705-b72f-aac21141114b
# ╟─d0513898-51d8-4871-b36a-35d5fcaa0629
# ╠═a71d51c0-efdb-45bb-bb1b-821b1ee4bb81
# ╠═2c8cc361-2eb0-445f-95d7-7d65295a3087
# ╟─1dab20f7-dd45-42c1-9565-156b9d4f5009
# ╟─475c4317-427d-4558-923b-982662107ea2
# ╟─a2c93673-67ae-4ece-918c-add42bf0255e
# ╟─c944cbc9-a9b6-4538-836d-eafea086189f
# ╠═dbcf829c-06a9-4914-a846-c0dd7ebdabbb
# ╠═693c84e3-f6c2-48b4-98d3-cc1e9e8ddc60
# ╠═54f952ce-8f4b-4a0f-bf1e-117fc48d3996
# ╠═ff09ca2a-54b9-4112-abcf-92b744a86f9b
# ╠═b3c91cbe-baa8-421a-8ff1-b2b0bca215fa
# ╠═be8c4bdf-501a-4f2c-8a0d-d7e034040b24
# ╠═da253a52-1954-4890-842b-b9444917deeb
# ╟─92bd2aae-62c7-4f84-9427-16e09a7bac14
# ╟─fc079169-d9ba-4f8d-a57d-613e8f2c5d49
# ╠═bd7dd4ab-ee29-4780-92a6-af77c440f6ba
# ╠═e9065f52-64bf-45d2-b5da-e75a4b408607
# ╟─ed899759-bc40-4d2a-9b2f-323eeb4870a0
# ╠═e26b7863-94ed-44d1-a0e3-a30b2546cfd7
# ╠═00f75f99-b015-4c20-89a9-dd47b903cd60
# ╟─f29302bb-ac90-4f6a-bbfb-ea22a8dae0db
# ╟─939a8c76-7fd1-43ca-ac90-8dd6fe4847df
# ╠═43d07b3b-19d4-410b-a02b-b21bd162d16c
# ╟─69be736e-08ac-4ee4-b291-f18c88cf6caf
# ╟─bc815fec-4251-47b2-99e7-2f4407d5bb87
# ╟─4a3698e5-bb3a-4d79-95df-41226407873e
# ╠═20d1daa7-4015-403c-95b2-716b31137dc2
# ╟─91e63de9-3227-41b9-8a7a-d44cfbd24309
# ╟─60362194-1571-44f9-89c1-8ae09f35e9e6
# ╠═aa8713e4-f665-4177-adae-9e9560aa4c81
# ╟─78512ef2-878e-4980-afef-a464e45cf03c
# ╟─fbb59f1a-3ca6-4a8a-a5af-6914a9bbe5ee
# ╟─da058b09-797e-42c2-8614-61df883a2ea0
# ╟─7b080495-6d5c-455d-91ec-e243476d7001
# ╟─957c0e3f-b8be-4058-a781-d77242249312
# ╟─428978e3-8551-4423-840e-75fbe67da471
# ╟─b2c0a4e7-565b-41ac-8f52-08640820ee39
# ╟─c95fe5ff-5a5c-4164-bf1d-c53e9ad99b14
# ╟─09d41f66-871d-47cf-a852-b0bb794870c4
# ╟─bfccd459-c94e-4502-bad6-b4d2dd1694b3
# ╟─d9268572-df0b-40ba-b3ef-7e908562bc34
# ╟─8ec3b3a2-83cd-469e-b34e-2f2c20ca4eb0
# ╟─46c1346c-1cd9-41d1-8216-13b882b1b8bf
# ╟─ff999ca8-b613-43a4-9718-8545bf61f184
# ╟─e3c2b3b1-6faf-482d-b880-583dbd177912
# ╟─354093b7-0205-4906-8c1f-8a6cdc4212b6
# ╟─e5db2608-f53d-4fa0-9e72-510b954eee8e
# ╟─ea8f5356-81b7-4bd1-896a-686aa331c45e
# ╟─9908d8b6-61bf-4f73-a0c3-7a4f7a3a9016
# ╟─cde431b4-2a17-4507-9124-cf3f8b953abe
# ╟─40c4ff27-670a-4c80-b588-7919f48bc5bd
# ╟─a83af44b-ad18-48c4-a600-deb5f5f42f29
# ╟─9631161e-0e1f-41f5-9c9f-92756fa5d736
# ╟─ce1da142-ba57-489f-8762-fb0be86acc30
# ╟─6c6a8bbe-d011-4a77-80b7-c211abd29b2f
# ╟─bf1196cd-e33c-4aaa-8ed3-294afe7be6e5
# ╠═994fef1c-56eb-4410-8009-caf9a47c111d
# ╟─1b43324f-e9ca-473e-b093-7bf24c93f426
# ╟─8b20cb06-7688-4d69-9d04-9fea0986ffbd
# ╟─a9737cdc-0f2b-4b0e-92a0-b0fb463563ca
# ╟─b6a330b6-9d9c-4f63-83ad-516873f99057
# ╠═722da3ba-286d-4d0b-bf6e-598e623afe28
# ╟─0b56fb23-ad3d-4d97-93d4-f7433beca6da
# ╟─b73d7ae7-75ff-44da-bf34-21b3ab3b36e9
# ╟─41277ca6-db9f-4481-98e4-fb44674014e5
# ╠═9e40ab4c-64b7-4ddb-8bdb-a2796bfec9ed
# ╠═e853f4ec-7cf4-48e1-bd34-ef57ae1b4772
# ╠═157f8136-29c4-44b4-ad9e-e5aa1c291cb2
# ╟─ad09b9cb-2fe1-47f6-abcc-1ac66aaae37e
# ╟─4a39856f-87ae-472e-becb-4d2e4a911d61
# ╟─6aa20591-3d86-441b-a480-4cc54dcb79dd
# ╟─9eb8ea54-25de-4690-b1ab-f73435049442
# ╠═d9543116-23a2-495c-833e-eab8e1b378e1
# ╠═94430cd4-ec29-4b33-9f86-e3d8b36b6729
# ╟─2e210e9a-e7d2-4032-bdac-da50ab1599f2
# ╟─21af0363-49a6-4b71-b303-4da665444174
# ╠═af7bb812-33ab-4b15-9ac6-baea5b3b81b2
# ╠═d01fa8da-857d-4436-96ad-4eddc6e11c23
# ╟─cbba62c9-4b99-4724-a283-d3f92c162aa1
# ╟─f2ebf95d-1dc7-4d43-b569-3c4619019bcc
# ╟─36f7c1be-1ba3-4b5e-bae7-2f6ad5dd6c46
# ╟─91fde27d-92b3-4dda-b4ac-9fe5aa49cac1
# ╟─36f94bae-3200-413a-b1fc-fc22e1413d8c
# ╟─9dec2b52-79ea-4588-ac89-237d3fa5b302
# ╟─16ebda6d-7bb5-426d-9e0d-a792c6202eb8
# ╠═afec87c3-dbd2-4671-bfd2-3442a39cf4e4
# ╠═573f9fb9-c5e3-4657-9aca-fa2530402ba9
# ╠═c92a01ef-909b-4948-8e8a-74f459f435c0
# ╠═4a415610-63c7-45f3-b805-cf7fc0f4ae41
# ╠═9c6ad57b-48d6-40b5-82db-fb64e1e2c858
# ╟─0f5ef458-a444-4c86-b3fb-86b7b82b878c
# ╟─f0fc001e-7c67-4bea-9e0c-44949f023e9c
# ╠═96e1d173-f24d-4682-bfe1-5871f372fec0
# ╟─83e96667-2543-42b4-826a-ddcf2d8ba252
# ╠═d1626a4b-4c56-4085-831f-01a4efc1e7dc
# ╟─bd58561d-03c6-47cd-87ee-c80ba7b54601
# ╟─83520631-205e-4a3e-bff0-6bc33ef25d01
# ╟─d9f767c2-d1b4-4998-a113-1524d3762ce8
# ╟─6c76cc1a-ded6-4224-acee-f140e9fae43c
# ╟─e37aa4e5-0574-4e20-842c-c46df54e3ca0
# ╟─beab92a2-f8b8-440d-b2ba-3c88a55a3441
# ╟─838076f5-04ab-4eca-8763-98a27b8637de
# ╟─67e10d4c-94fa-4586-a138-fe3548d42468
# ╟─ffedd213-09d3-414f-aa7a-619d80b6d748
# ╠═4235aef9-9789-48c7-a548-fdceeebbd3ac
# ╠═c092187a-a148-4092-bae2-450bd9e1d3d4
# ╠═bca491f5-9446-4fcc-af0e-ca5835c14286
# ╠═85097c68-d655-4ecd-ba48-aadeb4f3368d
# ╠═dee5a117-26d1-4371-a765-be6d0eab8037
# ╠═0707634c-6a9e-465e-96be-1598a9336da5
# ╠═4893e9f2-0482-4f51-8581-f015135e09e7
# ╠═8fe1a290-24fa-416e-949c-11448b9a1e71
# ╟─2dd60715-9a7d-4892-aa45-c695b143a6da
# ╠═7b11fbbd-f3a9-42b0-aea7-a5c5f63b95e8
# ╟─9e43acb4-e0a4-4658-a4f7-a4154d60ea89
# ╟─4f1e689a-f775-4de9-83ba-9515fc1db364
# ╟─05c33ad2-e2ae-4916-a1ec-6c885a424c1e
# ╠═222030b4-7074-4d42-aa09-089e5f54b644
# ╟─050b289a-c636-4944-b411-ce272deb1e2c
# ╟─cd78642a-8ab5-41b6-a5ac-ee07a4ac9f45
# ╟─32fe31d3-d242-4604-ae7c-7e99009eadd5
# ╠═7d783db2-a669-4a09-b7d6-b21a26b94326
# ╟─f4f3dbbb-f202-4c8d-bd02-99a3d6d932b2
# ╟─3adba9e9-fa6b-4e8a-891c-82bbd9f48f45
# ╟─16484697-c3b3-4509-8640-7e8121df9d23
# ╟─58073901-4f8d-4b37-8a0a-c6eeb32d3461
# ╟─aed47d09-bd50-4ce5-8771-f700874efddc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
