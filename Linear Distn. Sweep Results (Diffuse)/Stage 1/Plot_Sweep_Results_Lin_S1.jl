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
using Plots; plotly()

# ╔═╡ a226d405-d910-4112-8b6f-a153420a38ca
using DelimitedFiles, LinearAlgebra, Statistics

# ╔═╡ 838cbfe0-f0ad-11eb-0146-a989eeaafda6
md"""
# Plotting COMSOL Sweep Results

##### Stepped Distribution Case - Stage 1: Effect of $R$ and $L_e$
"""

# ╔═╡ aa57714a-6ba3-48ec-afcd-31b44845a4e8
readdlm("Parametric Sweep Levels.txt")

# ╔═╡ 6317f5cc-bbf9-4c22-a813-c500a538d380
# Parameter Levels
begin
	R = [0.5, 1, 2, 4, 8]
	E = [1e-5, 0.1, 0.3, 0.5, 0.7, 0.9]
	LD = [25]
	M = [1000]	#m^-1
end

# ╔═╡ cdfd8cb1-7602-45ad-abc4-8d9ed53eee19
t = 0.4 	#channel wall half-thickness [mm]

# ╔═╡ 1f6244a6-09a4-4e7c-9468-8aa6ff376070
# Additional Parameters of Interest
begin
	# Re = readdlm("Reynolds No.txt")[6:end,1:3]
	# Re = Re[:,3]				#Reynolds number
	
	ϕ = (R.^2) ./ (R .+ t).^2 	#Porosity
	
	const emis_low = 0.1
	const emis_high = 0.9

	A_front = 4*t * (2*R .+ t)	#[mm^2] per channel
	A_int = 16*R.^2 * LD[1] 		#[mm^2] per channel

	ε = [
		(emis_high*(A_front[r] + (1 - E[e])*A_int[r]) + emis_low*E[e]*A_int[r]) / (A_front[r] + A_int[r])
		for r = 1:length(R), e = 1:length(E)
	]

	ε = vec(ε[end,:])
end

# ╔═╡ b6491ea6-b497-4c2b-b3af-4cfbf85aa007
readdlm("Global Powers.txt")[5,:]

# ╔═╡ c62cb59a-864e-4a16-996f-5423fd08569b
Eff = readdlm("Global Powers.txt"; comments=true, comment_char='%')

# ╔═╡ 1124a3c4-2b53-448a-9d63-278da7a2d834
P_ap_lin = readdlm("Power on Apperture.txt"; comments=true, comment_char='%')[:,end]

# ╔═╡ 75527708-4ca3-4277-a4fe-a8de6aa888a9
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
	BHS_int_lin = Eff[:,n_param+11]	# Primary ray power absorbed by interior surfaces [W]
	BHS_front_lin = Eff[:,n_param+10]	# Primary ray power absorbed by front surface  [W]

end

# ╔═╡ da408ddb-5b01-4376-ad34-c19c3177910d
begin
	T_out = zeros(length(R),length(E))
	for i = 1:length(R), j = 1:length(E)
		lin_index = (i-1)*length(E) + j
		T_out[i,j] = T_out_lin[lin_index]
	end
end

# ╔═╡ 42efe36f-025c-44ad-8890-a32a86f7b571
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
	ylabel = "Refelctive Region,  L<sub>e</sub>/L<sub>ch</sub>",
	title = "Exit Gas Temperature (K)",
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

# ╔═╡ 8c024801-84be-4fff-aaee-3a55cc84b396
surface(
	log2.(R), E, transpose(T_out),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, E, transpose(T_out),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, E, transpose(T_out),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, E, transpose(T_out),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Refelctive Region,  L<sub>e</sub>/L<sub>ch</sub>",
	title = "Exit Gas Temperature (K)",
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

# ╔═╡ d33a36b0-1e13-4476-9e8b-bd1b99cd5061
md"""
Replotting based on average emissivity

"""

# ╔═╡ 743091d1-3f24-4899-b694-2b5e3d7e9ed8
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

# ╔═╡ b1b2f925-8f6f-4b23-bae0-c2c2550c5b53
# # Interpolating all relevant properties
# begin
# 	Eff_Ap_int = interpolate_lin(Eff_Ap, ε, E_base)
# 	T_out_int = interpolate_lin(T_out, ε, E_base)
# 	vol_eff_int = interpolate_lin(vol_eff[:,ld,:], ε, E_base)
# 	Q_abs_tot_int = interpolate_lin(Q_abs_tot, ε, E_base)
# 	BHS_tot_int = interpolate_lin(BHS_tot, ε, E_base)
# 	Q_loss_tot_int = interpolate_lin(Q_rad_loss_tot, ε, E_base)
# 	Q_loss_tot_front_int = interpolate_lin(Q_rad_loss_tot_front, ε, E_base)
# 	Q_loss_tot_int_int = interpolate_lin(Q_rad_loss_tot_int, ε, E_base)
# 	T_s_in_int = interpolate_lin(T_s_in, ε, E_base)
# 	T_s_avg_int = interpolate_lin(T_s_avg, ε, E_base)
# 	interpolated_results = [
# 		Eff_Ap_int; T_out_int; vol_eff_int; Q_abs_tot_int; BHS_tot_int; Q_loss_tot_int_int; Q_loss_tot_front_int; Q_loss_tot_int; T_s_in_int; T_s_avg_int
# 	]
# end

# ╔═╡ 335c4f83-d00f-4b07-ab28-37271fc5308c
md"""
Value of $L_e/L_{ch}$ for length optimization $= p_x$
"""

# ╔═╡ 142b3096-9165-4470-bdc1-43ffb28b8eef
begin
	contourf(
		log2.(R), ε, transpose(T_out),
		# log2.(R), E_base, transpose(T_out_int),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(T_out),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, E, transpose(T_out),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(T_out),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		yaxis = [0.2, 0.4, 0.6, 0.8],
		ylabel = "Average Emissivity",
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
		# # Size formatting for word
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
	hline!(
		[0.8],
		line = :dash,
		linewidth = 2,
		color = :black,
		legend=:false
	)

	vline!(
		[-01],
		line = :dash,
		linewidth = 2,
		color = :black,
		legend=:false
	)
end

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

# ╔═╡ ebe6df24-c623-4ffd-9857-40110fb0336d
T_amb = 318 	#K

# ╔═╡ 0c49cfb9-2787-4bc7-a25a-5869e33e6457
h_nat = 10  #[W/m^2.K]

# ╔═╡ 5ccfecd6-c794-4ca8-ab1e-d21ed06e69dd
F_amb_min = 0.01 #Lowest view factor to be considered

# ╔═╡ 20b8d53d-d6ca-4b69-b409-775e9722e746
function emissivity(refl_reg, m, z)
	emis = zeros(length(z))
	L_channel = z[end]
	Len = refl_reg * L_channel
	
	Lex_hypo = ((emis_high-emis_low)/m)+Len  	
	Lex = (if (Lex_hypo<L_channel) Lex_hypo else L_channel end)
		
	i_start = findfirst(isapprox.(z, Len; atol=0.5))
	i_end = findfirst(isapprox.(z, Lex; atol=0.5))
		
	emis[1:i_start] .= emis_low
	emis[i_start+1:i_end] = m * (z[i_start+1:i_end] .- Len) .+ emis_low
	emis[i_end+1:end] .= emis_high
		
	return emis
end

# ╔═╡ d9eab463-1477-4703-bc73-c3d8fcc9314c
begin
	plot(
		range(0,50,length=100),
		[emissivity(E[e],M[1],range(0,50,length=100)) for e = 1:length(E)
			],
		xlabel = "Channel Depth, z",
		ylabel = "Emissivity, <i>ε(z)</i>",
		label = permutedims(["L<sub>e</sub>/L<sub>ch</sub> = $(E[e])" for e = 1:length(E)]),
		xlims = (-5,50),
		ylim = (0,1),
		color = permutedims(palette(:tab10)[1:length(R)]),
		linewidth = 2,
		legend = :right,
		top_margin = 12*Plots.mm,
		left_margin = 3*Plots.mm,
		right_margin = 5*Plots.mm,
		tickfontsize = 12,
		guidefont = 14,
		titlefont = 16,
		fontfamily = "ComputerModern",
		size = (624/1.2,477/1.2)
	)
	
	plot!(
		range(-5,50,length=100),
		[ones(100)*emis_low,
		ones(100)*emis_high],
		color = :black,
		line = :dash,
		linealpha = 0.3,
		linewidth = 2,
		label = "",
		
		)
end

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

# ╔═╡ 6536390f-6cee-44c2-9a88-c23df8205992
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
		Q_abs_tot[r,e] = Q_abs_lin[lin_index]*n_channels(R[r],t) #[W] per mod
		
		#--------------------------------------------------------------------------

		Eff_Ap[r,e] = Q_abs_lin[lin_index]./P_ap_lin[r]
	end
end

# ╔═╡ 31ee1595-be66-4d8a-9be1-96b2514dcc0c
contourf(
	log2.(R), E, transpose(Eff_Ap),
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
	ylabel = "Refelctive Region,  L<sub>e</sub>/L<sub>ch</sub>",
	title = "Thermal Efficiency, η",# - using Power on Apperture",
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

# ╔═╡ 6e2d8d5b-627e-416b-8e84-79a3af3b384a
surface(
	log2.(R), E, transpose(Eff_Ap),
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
	ylabel = "Refelctive Region,  L<sub>e</sub>/L<sub>ch</sub>",
	title = "Thermal Efficiency, η",# - using Boundary Heat Source",
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

# ╔═╡ faec19c6-a584-4421-991b-3cb01abbbc8c
begin
	contourf(
		log2.(R), ε, transpose(Eff_Ap),
		# log2.(R), E_base, transpose(Eff_Ap_int),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, E, transpose(Eff_Ap),
		# xlabel = "Radius (mm)",
		# xscale=:log3,
		# xticks = R,
		# xlim = (2e-1, 2e3),
		#---------------------------------------------
		# Re, E, transpose(Eff_Ap),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, E, transpose(Eff_Ap),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Average Emissivity",
		title = "Thermal Efficiency, η",
		yaxis = range(0.2, 0.8, step=0.2),
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
		# # Size formatting for word
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
	# 	[-1.],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend=:false
	# )
end

# ╔═╡ b9ae10ff-f59d-4be3-b7da-99fb678310cb
ε

# ╔═╡ f226cd78-4380-4340-a3b9-9d7bb1940abb
md"""
## Additional plots

For heat transfer analysis
"""

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

# ╔═╡ 2a730b37-d2d3-4b32-9c47-24fd5b122d95
md"""
#### Plots that use absolute power values (per module)

"""

# ╔═╡ 76f0b4fd-dc94-499d-8d41-c58cc2d3ff2d
contourf(
	log2.(R), ε, transpose(Q_abs_tot),
	xlabel = "log<sub>2</sub>(Radius - mm)",
	#---------------------------------------------
	# R, ε, transpose(Q_abs_tot),
	# xlabel = "Radius (mm)",
	#---------------------------------------------
	# Re, ε, transpose(Q_abs_tot),
	# xlabel = "Reynolds Number",
	#---------------------------------------------
	# ϕ, ε, transpose(Q_abs_tot),
	# xlabel = "Porosity, ϕ",
	#---------------------------------------------
	ylabel = "Averag Emissivity",
	title = "Heat Absorbed by the Gas,  <i>Q<sub>abs,g</sub> </i> <br> per Module (J/kg)",
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

# ╔═╡ b0c6d7b6-13eb-4dc9-9b31-bc3e0a85da57


# ╔═╡ 29909be9-8863-4e09-8219-528adbcdb0d3
E

# ╔═╡ c944cbc9-a9b6-4538-836d-eafea086189f
md"""
##### Plotting Temperature Profiles and Volumetric Effect
"""

# ╔═╡ d9902d65-9c44-4d64-8f55-71a6bb4e555a
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

# ╔═╡ 95d1da60-d807-4ebc-8633-105fa948ddf4
Q_rad_loss_tot

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
	title = "Radiative Heat Loss per Channel using<br> Elements with <i>F<sub>amb</sub> > $(round(F_amb_min; digits=2)) </i> (W/m<sup>2</sup>)",
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

# ╔═╡ c917cb94-36a6-49f1-bf2c-7681f6172310
begin
	contourf(
		log2.(R), ε, transpose(Q_rad_loss_tot),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, ε, transpose(Q_rad_loss_tot),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, ε, transpose(Q_rad_loss_tot),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, ε, transpose(Q_rad_loss_tot),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Average Emissivity",
		yaxis = [0.2, 0.4, 0.6, 0.8],
		# title = "Radiative Heat Loss per Module<br>w/ Elements with <i>F<sub>amb</sub> > $(round(F_amb_min; digits=2)) </i> (W)",
		title = "Radiative Heat Loss per <br> Module (W)",
		# c = cgrad(:thermal, rev = true),
		c = cgrad(:thermal),
		levels = range(0, 2500, step = 100),
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
	 	# guidefont = 11,
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
		legend = false,
		)
	
	vline!(
		[-1.],
		line = :dash,
		color = :white,
		linewidth = 2,
		legend = false,
		)
end

# ╔═╡ 1672b8ce-a4a5-44c5-aa7a-54af60a29451
contour(
	log2.(R), ε, transpose(Q_rad_loss_tot_front),
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
	ylabel = "Average Emissivity",
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
begin
	contour(
		log2.(R), ε, transpose(Q_rad_loss_tot_int),
		xlabel = "log<sub>2</sub>(Radius - mm)",
		#---------------------------------------------
		# R, ε, transpose(Q_rad_loss_tot_int),
		# xlabel = "Radius (mm)",
		#---------------------------------------------
		# Re, ε, transpose(Q_rad_loss_tot_int),
		# xlabel = "Reynolds Number",
		#---------------------------------------------
		# ϕ, ε, transpose(Q_rad_loss_tot_int),
		# xlabel = "Porosity, ϕ",
		#---------------------------------------------
		ylabel = "Average Emissivity",
		title = "Radiative Heat Loss per Module<br> from Interior Surfaces (W)",
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
	 	# guidefont = 11,
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
		legend = false,
		)
	
	vline!(
		[-1.],
		line = :dash,
		color = :black,
		linewidth = 2,
		legend = false,
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
	log2.(R), ε, transpose(T_s_in),
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
	ylabel = "Average Emissivity",
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
	log2.(R), ε, transpose(T_s_avg),
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
	ylabel = "Average Emissivity",
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

# ╔═╡ be8c4bdf-501a-4f2c-8a0d-d7e034040b24
colors = permutedims(palette(:tab10)[1:length(E)])

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
		# label = false,
		label = permutedims(["L<sub>e</sub>/L<sub>ch</sub> = $(E[e])" for e = 1:length(E)]),

		)

	plot!(
			x_solid_dict[r1,ld1],

			[T_solid_dict[r1,ld1][:,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :solid,
		linewidth = 2,
		label = permutedims(["L<sub>e</sub>/L<sub>ch</sub> = $(E[e])" for e = 1:length(E)]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Axial Temp. Profiles (R<sub>ch</sub> = $(R[r1])mm)",
		legend = :outerbottomright,
		# legend = false,
		
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 5*Plots.mm,
		left_margin = 6*Plots.mm,
		legendfont = 10,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		# size = (440,335),	#for word doc.
		size = (440*1.4,335),	#for word doc.
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 5*Plots.mm,
		# left_margin = 5*Plots.mm,
		# tickfontsize = 10,
		# guidefont = 11,
		# titlefont = 12,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (325,250),	#for word doc.
			
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
	ylabel = "Refelctive Region,  L<sub>e</sub>/L<sub>ch</sub>",
	title = "Volumetric Effect, T<sub>g,out</sub> / T<sub>s,in</sub> <br>(L/D = $(LD[ld]))",
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

# ╔═╡ 92f3e394-6119-4f2c-96d9-ee5057bea2f1
begin
	contourf(
		log2.(R), ε, transpose(vol_eff[:,ld,:]),
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
		ylabel = "Average Emissivity",
		yaxis = [0.2, 0.4, 0.6, 0.8],
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
	# 	[-1.],
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

# ╔═╡ ee777a68-05da-4696-aa3e-d4fcddd3367a
begin
	contourf(
		log2.(R), ε, transpose(BHS_tot),
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
		ylabel = "Average Emissivity",
		yaxis = [0.2, 0.4, 0.6, 0.8],
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
		legend = false,
		)
	
	vline!(
		[-1.],
		line = :dash,
		color = :black,
		linewidth = 2,
		legend = false,
		)
end

# ╔═╡ 5efa2d9f-d099-4c78-90f1-252f8f576a1a
contourf(
	log2.(R), ε, transpose(BHS_int_tot),
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
	ylabel = "Average Emissivity",
	title = "Interior Boundary Heat Source <br> per Module (W)",
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

# ╔═╡ 0d83b674-5971-4b0d-b237-56d3c08fe0dd
contourf(
	log2.(R), ε, transpose(BHS_front_tot),
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
	ylabel = "Average Emissivity",
	title = "Front Boundary Heat Source <br> per Module (W)",
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

# ╔═╡ 004a727f-1214-4b9b-a90b-482270633700
begin
	contour(
		log2.(R), ε, transpose(BHS_int_tot./BHS_front_tot),
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
		ylabel = "Average Emissivity",
		# title = "BHS<sub>int</sub> / BHS<sub>front</sub> per Module",
		title = "BHS<sub>int</sub> / BHS<sub>front</sub>",
		yaxis = range(0.2, 0.8, step=0.2),
		# levels = range(1,20,step=2),
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
		# tickfontsize = 10,
	 #    guidefont = 11,
		# titlefont = 12,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (325,250),	#for word doc.
			
		)
	# plot!(
	# 	log2.(range(R[1], R[end], length=10)),
	# 	mean(ε[1:3])*ones(10),
	# 	linewidth = 2,
	# 	line = :dash,
	# 	color = :black,
	# 	legend = false,
	
	# 	)
	
	# plot!(
	# 	log2.(mean(R[1]))*ones(10),
	# 	range(ε[1], ε[end], length=10),
	# 	linewidth = 2,
	# 	line = :dash,
	# 	color = :black,
	# 	legend = false,
	
	# 	)
end

# ╔═╡ 38e0e42f-db6a-4269-8aa8-0693b2650e17
BHS_int_tot./BHS_front_tot

# ╔═╡ 0c4c7b1b-f5bd-46c1-acd0-173ca3c8ab02
[BHS_tot[i,3] - BHS_tot[i,end] for i =1:length(R)]

# ╔═╡ 7466e9ae-0a09-4c7f-8067-baed695e80ec
md"""
#### Extracting and Plotting The BHS Profiles
"""

# ╔═╡ 700df4c1-d7c3-4698-8f53-b3a91f993a44
BHS_profiles = "Cutline BHS Profiles"

# ╔═╡ d274af72-1018-4061-b9de-2723b0fba417
begin
	x_bhs_dcit = [
		readdlm("$(BHS_profiles)//BHS_comb (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,1]
		for r = 1:length(R), ld = 1:length(LD)
				]
	
	bhs_z_dict = [
				readdlm("$(BHS_profiles)//BHS_comb (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,2:end]
		for r = 1:length(R), ld = 1:length(LD)
				]
end

# ╔═╡ f0f553de-2103-4f00-8eda-23fccba63545
begin
	# r1 = 2
	# ld1 = 1
	plot(
		x_bhs_dcit[r1,ld1],
		[bhs_z_dict[r1,ld1][:,e1] for e1 = 1:length(E)]./1e3,
		color = colors,
		line = :solid,
		linewidth = 2,
		label = permutedims(["L<sub>e</sub>/L<sub>ch</sub> = $(E[e])" for e = 1:length(E)]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "BHS Density (kW/m<sup>2</sup>)",
		title = "Axial BHS Profiles (R<sub>ch</sub> = $(R[r1])mm)",
		legend = :right,

		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 5*Plots.mm,
		left_margin = 6*Plots.mm,
		legendfont = 12,
		tickfontsize = 14,
	    guidefont = 16,
		titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
		#--------------------------------------------
		# Size formatting for word
		# top_margin = 5*Plots.mm,
		# left_margin = 5*Plots.mm,
		# tickfontsize = 10,
		# guidefont = 11,
		# titlefont = 12,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (325,250),	#for word doc.
			
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
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
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
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

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
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

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
# ╠═a226d405-d910-4112-8b6f-a153420a38ca
# ╠═aa57714a-6ba3-48ec-afcd-31b44845a4e8
# ╠═6317f5cc-bbf9-4c22-a813-c500a538d380
# ╠═cdfd8cb1-7602-45ad-abc4-8d9ed53eee19
# ╠═1f6244a6-09a4-4e7c-9468-8aa6ff376070
# ╠═b6491ea6-b497-4c2b-b3af-4cfbf85aa007
# ╠═c62cb59a-864e-4a16-996f-5423fd08569b
# ╠═1124a3c4-2b53-448a-9d63-278da7a2d834
# ╠═75527708-4ca3-4277-a4fe-a8de6aa888a9
# ╠═6536390f-6cee-44c2-9a88-c23df8205992
# ╠═da408ddb-5b01-4376-ad34-c19c3177910d
# ╟─31ee1595-be66-4d8a-9be1-96b2514dcc0c
# ╟─6e2d8d5b-627e-416b-8e84-79a3af3b384a
# ╟─42efe36f-025c-44ad-8890-a32a86f7b571
# ╟─8c024801-84be-4fff-aaee-3a55cc84b396
# ╟─fa856550-9cd1-496b-a2b5-2cf79c22dd33
# ╟─d33a36b0-1e13-4476-9e8b-bd1b99cd5061
# ╠═743091d1-3f24-4899-b694-2b5e3d7e9ed8
# ╠═b1b2f925-8f6f-4b23-bae0-c2c2550c5b53
# ╟─335c4f83-d00f-4b07-ab28-37271fc5308c
# ╟─faec19c6-a584-4421-991b-3cb01abbbc8c
# ╟─142b3096-9165-4470-bdc1-43ffb28b8eef
# ╟─92f3e394-6119-4f2c-96d9-ee5057bea2f1
# ╠═3c5eaa90-ec96-4757-9356-956a98840aef
# ╟─44cec87f-537d-426d-8089-934210aee5b3
# ╟─03f25a88-efb5-40bb-9c13-4ca2b8b97c14
# ╟─ebe6df24-c623-4ffd-9857-40110fb0336d
# ╟─0c49cfb9-2787-4bc7-a25a-5869e33e6457
# ╠═5ccfecd6-c794-4ca8-ab1e-d21ed06e69dd
# ╟─20b8d53d-d6ca-4b69-b409-775e9722e746
# ╟─d9eab463-1477-4703-bc73-c3d8fcc9314c
# ╠═aaa5252e-0714-476d-8267-f13d8a1f9e0d
# ╠═0593c86b-fb0f-40c2-a38f-9c7330c7a484
# ╠═95d1da60-d807-4ebc-8633-105fa948ddf4
# ╠═b9ae10ff-f59d-4be3-b7da-99fb678310cb
# ╟─f226cd78-4380-4340-a3b9-9d7bb1940abb
# ╟─e5300a51-d793-4fe2-88b4-8a117924f91b
# ╟─dec3866e-881e-42a0-9236-e65fd6fd8599
# ╟─8cb6d0d1-9301-419e-a242-817b0b718ae6
# ╟─2caf1b17-51cf-4d96-93d3-948104c92fa5
# ╟─2a730b37-d2d3-4b32-9c47-24fd5b122d95
# ╟─76f0b4fd-dc94-499d-8d41-c58cc2d3ff2d
# ╟─ee777a68-05da-4696-aa3e-d4fcddd3367a
# ╟─5efa2d9f-d099-4c78-90f1-252f8f576a1a
# ╟─0d83b674-5971-4b0d-b237-56d3c08fe0dd
# ╟─004a727f-1214-4b9b-a90b-482270633700
# ╠═38e0e42f-db6a-4269-8aa8-0693b2650e17
# ╟─b0c6d7b6-13eb-4dc9-9b31-bc3e0a85da57
# ╠═0c4c7b1b-f5bd-46c1-acd0-173ca3c8ab02
# ╠═29909be9-8863-4e09-8219-528adbcdb0d3
# ╟─c917cb94-36a6-49f1-bf2c-7681f6172310
# ╟─1672b8ce-a4a5-44c5-aa7a-54af60a29451
# ╟─549fbc4f-b840-4a49-b0dc-2174d9d8468d
# ╟─a71d51c0-efdb-45bb-bb1b-821b1ee4bb81
# ╟─2c8cc361-2eb0-445f-95d7-7d65295a3087
# ╟─1dab20f7-dd45-42c1-9565-156b9d4f5009
# ╟─475c4317-427d-4558-923b-982662107ea2
# ╟─c944cbc9-a9b6-4538-836d-eafea086189f
# ╠═d9902d65-9c44-4d64-8f55-71a6bb4e555a
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
# ╟─7466e9ae-0a09-4c7f-8067-baed695e80ec
# ╠═700df4c1-d7c3-4698-8f53-b3a91f993a44
# ╠═d274af72-1018-4061-b9de-2723b0fba417
# ╟─f0f553de-2103-4f00-8eda-23fccba63545
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
