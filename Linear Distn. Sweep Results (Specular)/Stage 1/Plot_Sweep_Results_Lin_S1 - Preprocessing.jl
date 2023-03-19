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

##### _Fully_ Specular Stepped Distribution Case - Stage 1: Effect of $R$ and $L_e$
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

# ╔═╡ 7e296429-ca31-449c-a1c7-514001586851
L_channel = R.*2*LD[1]

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

# ╔═╡ 1d0c4527-1a3b-4767-b9ed-06df5054e5ae
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
	BHS_front_lin = Eff[:,n_param+10]	# Primary ray power absorbed by front surface  [W]
	BHS_int_lin = Eff[:,n_param+11]	# Primary ray power absorbed by interior surfaces [W]

end

# ╔═╡ 819592b5-8152-4db4-8ced-7472c0381e38
plot(Q_abs_lin - (-(Q_abs_out_lin) .+ (Q_abs_in_lin)))

# ╔═╡ da408ddb-5b01-4376-ad34-c19c3177910d
begin
	T_out = zeros(length(R),length(E))
	for i = 1:length(R), j = 1:length(E)
		lin_index = (i-1)*length(E) + j
		T_out[i,j] = T_out_lin[lin_index]
	end
end

# ╔═╡ 80f21ac4-6478-40f1-aad1-7f169a20dfcf
E[4]

# ╔═╡ a7d27b67-7b6d-49ef-b267-cb6942bd2d02
R[4]

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

# ╔═╡ 382543f1-c574-47c2-93d0-e3643288c9a6
mean(T_out[4,4:5])

# ╔═╡ 1c39de78-d278-4b89-924e-113693e6bfb8
md"""
Optimal design parameters drawn with dashed lines: R = $(round(2^1.7;digits=1))mm and ε = 0.45
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
		ylabel = "Average Emissivity",
		title = "Exit Gas Temperature (K)",
		yaxis = range(0.2, 0.8, step=0.2),
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

	hline!(
		[0.45],
		line = :dash,
		linewidth = 2,
		color = :black,
		legend = false,
		)

	vline!(
		[1.7],
		line = :dash,
		linewidth = 2,
		color = :black,
		legend = false,
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
	q_ap_lin = [q_ap*1e3/n_channels(R[r],t)*(A_tot*1e-6) for r = 1:length(R)] #[W]
	
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
		
		# Eff_Ap[r,e] = Q_abs_lin[lin_index]./q_ap_lin[r]
		# Eff_Ap[r,e] = Q_abs_out_lin[lin_index]./P_ap_lin[r]
		# Eff_Ap[r,e] = (Q_abs_out_lin[lin_index] - Q_abs_in_lin[lin_index])./P_ap_lin[r]
	end
end

# ╔═╡ 0799d6a5-84fc-431c-83b1-249ce9805e12
Eff_Ap

# ╔═╡ da2add34-b3f8-4cab-915a-bfe7e378246a
maximum(Eff_Ap)

# ╔═╡ b79e969b-f842-4eb9-a4a6-b6585d33dcef
Eff_Ap[4,4]

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

# ╔═╡ 929c0b8e-50b0-4143-85d9-30b0810b1d01
Eff_Ap[4,5]

# ╔═╡ 2240b754-2cd1-4056-897d-e2f99fa37e86
mean(Eff_Ap[4,4:5])

# ╔═╡ d9c232e5-bd1f-47bc-80c8-15f540aeadc7
surface(
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
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 14,
	 #    guidefont = 16,
		# titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		contour_labels = true,
		size = (440,335),	#for word doc.
			
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

	# hline!(
	# 	[0.45],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend = false,
	# 	)

	# vline!(
	# 	[1.7],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend = false,
	# 	)
	
	# plot!(
	# 	log2.(range(R[1], R[end], length=10)),
	# 	mean(ε[4:5])*ones(10),
	# 	linewidth = 2,
	# 	line = :dash,
	# 	color = :black,
	# 	legend = false,
	
	# 	)
	
	# plot!(
	# 	log2.(mean(R[3:4]))*ones(10),
	# 	range(ε[1], ε[end], length=10),
	# 	linewidth = 2,
	# 	line = :dash,
	# 	color = :black,
	# 	legend = false,
	
	# 	)
end

# ╔═╡ b9ae10ff-f59d-4be3-b7da-99fb678310cb
ε

# ╔═╡ f226cd78-4380-4340-a3b9-9d7bb1940abb
md"""
## Additional plots

For heat transfer analysis
"""

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
	title = "Heat Absorbed by the Gas,  <i>Q<sub>abs,g</sub> </i> <br> per Module (W)",
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

# ╔═╡ c944cbc9-a9b6-4538-836d-eafea086189f
md"""
##### Plotting Temperature Profiles and Volumetric Effect
"""

# ╔═╡ d9902d65-9c44-4d64-8f55-71a6bb4e555a
Temp_profile = "Cutline Temp. Profiles (unprocessed)"

# ╔═╡ 54f952ce-8f4b-4a0f-bf1e-117fc48d3996
begin
	x_gas = [
		readdlm("$(Temp_profile)//T_gas (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,1]
		for r = 1:length(R), ld = 1:length(LD)
				]
	
	T_gas = [
				readdlm("$(Temp_profile)//T_gas (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,2:end]
		for r = 1:length(R), ld = 1:length(LD)
				]
	#---------------------------------------------------------------------
	end_ind_gas = cat([
		[
			findall(x_gas[r,1] .== L_channel[r]) for r = 1:length(R)
		][r] 
		for r =1:length(R)]; dims=2)
	#---------------------------------------------------------------------
	x_gas_dict = [
		x_gas[r,1][(if (e == 1) 1 else (1+end_ind_gas[r][e-1]) end):end_ind_gas[r][e]]
		for r = 1:length(R), e = 1:length(E)]
	
	T_gas_dict = [
		T_gas[r,1][(if (e == 1) 1 else (1+end_ind_gas[r][e-1]) end):end_ind_gas[r][e]]
		for r = 1:length(R), e = 1:length(E)]
end

# ╔═╡ b3c91cbe-baa8-421a-8ff1-b2b0bca215fa
begin
	x_solid = 	[
		readdlm("$(Temp_profile)//T_solid (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,1]
		for r = 1:length(R), ld = 1:length(LD)
				]
	
	T_solid = [
		readdlm("$(Temp_profile)//T_solid (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,2:end]
		for r = 1:length(R), ld = 1:length(LD)
				]
	#---------------------------------------------------------------------
	end_ind_sol = cat( 	[
		[
			findall(x_solid[r,1] .== L_channel[r]) for r = 1:length(R)
		][r] for r =1:length(R)
						]; dims=2)
	#---------------------------------------------------------------------
	x_solid_dict = 	[
		x_solid[r,1][(if (e == 1) 1 else (1+end_ind_sol[r][e-1]) end):end_ind_sol[r][e]] for r = 1:length(R), e = 1:length(E)
					]
	
	T_solid_dict = 	[
		T_solid[r,1][(if (e == 1) 1 else (1+end_ind_sol[r][e-1]) end):end_ind_sol[r][e]] for r = 1:length(R), e = 1:length(E)]
	
end

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
			
		T_s = T_solid_dict[r,e][1]
		
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

# ╔═╡ c917cb94-36a6-49f1-bf2c-7681f6172310
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
	# title = "Radiative Heat Loss per Module<br>w/ Elements with <i>F<sub>amb</sub> > $(round(F_amb_min; digits=2)) </i> (W)",
	title = "Radiative Heat Loss per <br> Module (W)",
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

# ╔═╡ 1672b8ce-a4a5-44c5-aa7a-54af60a29451
contourf(
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
contourf(
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

# ╔═╡ a71d51c0-efdb-45bb-bb1b-821b1ee4bb81
begin
	T_s_avg = zeros(length(R),length(E))
	for r = 1:length(R), ld = 1:length(LD), e = 1:length(E)
		x = x_solid_dict[r,e]
		local A_int = (x*1e-3)*(4*R[r]*2*1e-3) 	#Interior surface area per channel
		T_s = T_solid_dict[r,e]
		m = 1	# Start length considered for integration
		n = length(x)
	
		T_s_avg[r,e] = line_avg(T_s, x, m, n)
	end
	
	md""" 
	Contains $T_{s, avg}$ calculations
	"""
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

# ╔═╡ 2c8cc361-2eb0-445f-95d7-7d65295a3087
begin
	T_s_in = zeros(length(R),length(E))
	for r = 1:length(R), ld = 1:length(LD), e = 1:length(E)
		T_s = T_solid_dict[r,e]
		T_s_in[r,e] = T_s[1]
	end
		
	md""" 
	Contains $T_{s, front}$ calculations
	"""
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

# ╔═╡ be8c4bdf-501a-4f2c-8a0d-d7e034040b24
colors = permutedims(palette(:tab10)[1:length(E)])

# ╔═╡ da253a52-1954-4890-842b-b9444917deeb
@bind r1 PlutoUI.Slider(1:length(R), show_value=true)

# ╔═╡ 92bd2aae-62c7-4f84-9427-16e09a7bac14
begin
	# r1 = 2
	ld1 = 1
	plot(
			[x_gas_dict[r1,e1] for e1 = 1:length(E)],

			[T_gas_dict[r1,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :dash,
		linewidth = 2,
		# label = false
		label = permutedims(["L<sub>e</sub>/L<sub>ch</sub> = $(E[e])" for e = 1:length(E)]),

		)

	plot!(
			[x_solid_dict[r1,e1] for e1 = 1:length(E)],

			[T_solid_dict[r1,e1] for e1 = 1:length(E)],
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

# ╔═╡ 7466e9ae-0a09-4c7f-8067-baed695e80ec
md"""
#### Extracting and Plotting The BHS Profiles
"""

# ╔═╡ 05d998a1-18ae-47bc-ad52-b931c987ef43
BHS_profiles = "Cutline BHS Profiles"

# ╔═╡ c45b508a-f7e0-4a55-a16f-1a3ebcb88df6
begin
	x_bhs = [
		readdlm("$(BHS_profiles)//BHS_comb (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,1]
		for r = 1:length(R), ld = 1:length(LD)
				]
	
	bhs_z = [
				readdlm("$(BHS_profiles)//BHS_comb (R=$(if (R[r] < 1) R[r] else Int(R[r]) end)mm, LD=$(LD[ld])).txt"; comments=true, comment_char='%')[:,2:end]
		for r = 1:length(R), ld = 1:length(LD)
				]
	#---------------------------------------------------------------------
	end_ind_bhs = cat([
		[
			findall(x_bhs[r,1] .== L_channel[r]) for r = 1:length(R)
		][r] 
		for r =1:length(R)]; dims=2)
	#---------------------------------------------------------------------
	x_bhs_dcit = [
		x_bhs[r,1][(if (e == 1) 1 else (1+end_ind_bhs[r][e-1]) end):end_ind_bhs[r][e]]
		for r = 1:length(R), e = 1:length(E)]
	
	bhs_z_dict = [
		bhs_z[r,1][(if (e == 1) 1 else (1+end_ind_bhs[r][e-1]) end):end_ind_bhs[r][e]]
		for r = 1:length(R), e = 1:length(E)]

	#---------------------------------------------------------------------
	# Removing NaNs
	for r = 1:length(R), e = 1:length(E)
		bhs = bhs_z_dict[r,e]
		local x_bhs = x_bhs_dcit[r,e]
		
		non_nans = findall(.!(isnan.(bhs)))
		
		if any(isnan.(bhs))
			bhs_z_dict[r,e] = bhs[non_nans]
			x_bhs_dcit[r,e] = x_bhs[non_nans]
		else
			continue
		end
	end

		
	md""" 
	Processing cutline $Q_{BHS}(z)$ data
	"""
end

# ╔═╡ fdffa90d-b353-4837-9a24-c1a69afd5191
begin
	# r1 = 2
	# ld1 = 1
	plot(
			[x_bhs_dcit[r1,e1] for e1 = 1:length(E)],

			[bhs_z_dict[r1,e1] for e1 = 1:length(E)]./1e3,
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
		# size = (440*1.4,335),	#for word doc.
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

	# hline!(
	# 	# [x_bhs_dcit[r1,e1] for e1 = 1:length(E)],
	# 	transpose([mean(bhs_z_dict[r1,e1]) for e1 = 1:length(E)]./1e3),
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = colors,
	# 	# label = "", 
	# 	label = permutedims(["L<sub>e</sub>/L<sub>ch</sub> = $(E[e])" for e = 1:length(E)]),
	# 	legend = (1.0,1.0)
	# 	)

end

# ╔═╡ fc079169-d9ba-4f8d-a57d-613e8f2c5d49
md"""
#### Extracting and Plotting The Volumetric Effect/Efficiency
"""

# ╔═╡ e9065f52-64bf-45d2-b5da-e75a4b408607
@bind ld PlutoUI.Slider(1:length(LD), show_value=true)

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

# ╔═╡ ee777a68-05da-4696-aa3e-d4fcddd3367a
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
	title = "Total Boundary Heat Source <br> per Module (W)",
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

# ╔═╡ bd03bcfb-bd4a-43ce-bb38-8cbbb402fe55
isapprox.((BHS_int_tot .+ BHS_front_tot), BHS_tot; atol = 1e-6)

# ╔═╡ 22ead363-fe5a-4b5a-8718-9a0c000d0c4d
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

# ╔═╡ d996de9d-406f-4cb2-930b-b33909befc74
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

# ╔═╡ 46f66376-05a1-4047-b572-adee3b37d519
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
		# levels = range(1,20,step=2),
		yaxis = range(0.2, 0.8, step=0.2),
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
		# tickfontsize = 14,
		# guidefont = 16,
		# titlefont = 18,
		# fontfamily = "ComputerModern",
		# colorbar_font = "ComputerModern",
		# contour_labels = true,
		# size = (440,335),	#for word doc.
			
		)
	
	# hline!(
	# 	[0.45],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend = false,
	# 	)

	# vline!(
	# 	[1.7],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend = false,
	# 	)
	
end

# ╔═╡ 7b9f1ec8-d09e-477d-ba5d-f2cc9c46503e
mean((BHS_int_tot./BHS_front_tot)[3:4,5])

# ╔═╡ 0c4c7b1b-f5bd-46c1-acd0-173ca3c8ab02
[BHS_tot[i,3] - BHS_tot[i,end] for i =1:length(R)]

# ╔═╡ 4a8e9b11-9aea-4dd7-b03f-b362860148aa
BHS_int_tot

# ╔═╡ a943d8e2-7dfe-4e92-80ea-bb0a8caa15c1
md"""
##### (1) Based on inlet and outlet temperatures
"""

# ╔═╡ bd7dd4ab-ee29-4780-92a6-af77c440f6ba
begin
	vol_eff = 	[
		T_solid_dict[r,e1][1]\T_gas_dict[r,e1][end]	
		for r = 1:length(R), ld = 1:length(LD),  e1 = 1:length(E)
	];
	
	md"""
	
	Calculating vol. effect definitions based on Luque (2017) definition
	
	"""
end

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

# ╔═╡ dc08bd04-8b58-4cdc-97b8-45d268739074
mean(vol_eff[4,:,4:5])

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
		title = "Volumetric Effect, T<sub>g,out</sub> / T<sub>s,in</sub>",# <br>(L/D = $(LD[ld]))",
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

	# hline!(
	# 	[0.45],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend = false,
	# 	)

	# vline!(
	# 	[1.7],
	# 	line = :dash,
	# 	linewidth = 2,
	# 	color = :black,
	# 	legend = false,
	# 	)
end

# ╔═╡ ec8a5e12-9eaf-4f42-989b-064d27169fdc
# Based on definiton of LMTD
begin
	∆T1 = [T_solid_dict[r,e1][1] - T_gas_dict[r,e1][1] for r = 1:length(R),  e1 = 1:length(E)]
	
	∆T2 = [T_solid_dict[r,e1][end] - T_gas_dict[r,e1][end] for r = 1:length(R), e1 = 1:length(E)]
	
	LMTD = 	[
		(∆T1[r,e1] - ∆T2[r,e1])/log(∆T1[r,e1] / ∆T2[r,e1])
		
	for r = 1:length(R),  e1 = 1:length(E)
			]

	# LMTD = LMTD_f(T_solid_dict, T_gas_dict)
	
	LMTD_A = 	[
		LMTD[r,e1] * A_surface_tot_int[r]*1e-6
		
	for r = 1:length(R),  e1 = 1:length(E)
			]

	md"""
	Calculation of $LMTD$

	"""
end

# ╔═╡ bea1991d-b1cd-4de2-90b8-2030d147d13a
function LMTD_f(T_solid_dict, T_gas_dict)
	
	∆T1 = [T_solid_dict[r,e1][1] - T_gas_dict[r,e1][1] for r = 1:length(R),  e1 = 1:length(E)]
	
	∆T2 = [T_solid_dict[r,e1][end] - T_gas_dict[r,e1][end] for r = 1:length(R), e1 = 1:length(E)]
	
	LMTD = 	[
		(∆T1[r,e1] - ∆T2[r,e1])/log(∆T1[r,e1] / ∆T2[r,e1])
		
	for r = 1:length(R),  e1 = 1:length(E)
			]

	
	return LMTD
end

# ╔═╡ ccbfa443-7eea-4986-a3c2-a00c26d779d8


# ╔═╡ 67ad79e5-dd7c-4361-8bb0-3078bee29a67
md"""
##### (2) Based on power ratios
"""

# ╔═╡ b82ca70a-b701-42b2-a774-25c187811e85
# Based on ratio of powers
begin
	vol_eff_p1 = Q_abs_tot./Q_rad_loss_tot
	vol_eff_p2 = Q_rad_loss_tot./BHS_tot
	vol_eff_p3 = Q_abs_tot./BHS_tot
	vol_eff_p4 = [BHS_tot[r,e]/(P_ap_lin[r] * n_channels(R[r],t)) for r = 1:length(R), e = 1:length(E)]
	md"""
	Calculating vol. effect definitions based on power ratios
	"""
end

# ╔═╡ 0e3968c3-da01-497d-832f-cb517797b8ce


# ╔═╡ b68d4847-221f-460c-b8a9-c8688bce14ac
md"""
##### (3) Based on statistical evaluations of ∆T
"""

# ╔═╡ 7409dde6-64ba-4209-94ab-fd7e3a97feb7
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

# ╔═╡ d066f9e0-768c-45ba-b029-154773517c34
# Extracting temperature profiles at the same 'z' points
begin
	T_s_int = Array{Any, 2}(undef, (length(R),length(E)))
	T_g_int = Array{Any, 2}(undef, (length(R),length(E)))
	Z_int = Array{Any, 2}(undef, (length(R),length(E)))
	
	for r = 1:length(R),  e = 1:length(E)
		L_channel = 2*R[r]*LD[ld] #mm
		z_int = range(0, L_channel, step = 1)
		Z_int[r,e] = z_int
		
		T_s_int[r,e] = interpolate_lin_1D(T_solid_dict[r,e], x_solid_dict[r,e], z_int)

		# println("r = ", r, "  ld = ", ld)
		T_g_int[r,e] = interpolate_lin_1D(T_gas_dict[r,e], x_gas_dict[r,e], z_int)

		
	end

			
	md""" 
	Extracting all temperature profiles at the same '$z$' points/locations
	"""
end		

# ╔═╡ ad913c29-00e1-46ef-b079-66aa0d37c1b7
# Saving Processed Temperature Profiles
processed_folder = "Cutline Temp. Profiles (processed)"

# ╔═╡ 3dc007bd-200a-4419-9122-86dd41e471f8
# for  r = 1:length(R), e = 1:length(E)
# 	open("$(processed_folder)\\T_gas (R=$(R[r])mm, LD=25, E=$(E[e])).txt", "w") do io
# 		   writedlm(io, [Z_int[r,e] T_g_int[r,e]])
# 			end
# 	open("$(processed_folder)\\T_solid (R=$(R[r])mm, LD=25, E=$(E[e])).txt", "w") do io
# 		   writedlm(io, [Z_int[r,e] T_s_int[r,e]])
# 			end
# end

md"""

 Saving processed temperature profiles
"""

# ╔═╡ 66c0abdf-9820-4478-9a9b-c1f9927fe7d2
begin
	# r1 = 2
	# ld1 = 1
	plot(
			[x_gas_dict[r1,e1] for e1 = 1:length(E)],

			[T_gas_dict[r1,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		line = :dash,
		linewidth = 2,
		# label = false
		label = permutedims(["L<sub>e</sub>/L<sub>ch</sub> = $(E[e])" for e = 1:length(E)]),

		)

	plot!(
			[x_solid_dict[r1,e1] for e1 = 1:length(E)],

			[T_solid_dict[r1,e1] for e1 = 1:length(E)],
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

	scatter!(
			[Z_int[r1,e1] for e1 = 1:length(E)],

			[T_g_int[r1,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		# linewidth = 2,
		label = false
		)
	
	scatter!(
			[Z_int[r1,e1] for e1 = 1:length(E)],

			[T_s_int[r1,e1] for e1 = 1:length(E)],
		color = colors,
		# color = permutedims(palette(:tab10)[1:length(E)]),
		# linewidth = 2,
		label = false
		)

end

# ╔═╡ c809e1e3-2edb-4640-bd5a-8e7c668efd0f
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

	
	for r = 1:length(R), e = 1:length(E)
		T_s = T_s_int[r,e]
		T_g = T_g_int[r,e]
		Z = Z_int[r,e]
		
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
		# R_corr[r,e] = (sum((T_s .- T_s_bar).*(T_g .- T_g_bar))			
		# 				)/(
		# 		sqrt(sum((T_s .- T_s_bar).^2)) * sqrt(sum((T_g .- T_g_bar).^2))
		# 				)		
		R_corr[r,e] = cor(T_s, T_g)
		
		NMB[r,e] = (T_s_bar - T_g_bar)/(T_g_bar)
		NMSD[r,e] = (sig_s - sig_g)/(sig_g)
					
	end

	md"""
	Calculating vol. effect definitions based on statistical evaluations of $∆T$
	"""
end

# ╔═╡ dc0f8c5c-b709-4aa1-a96c-9494e3f8ce58


# ╔═╡ 1d75576d-06aa-48dc-8ed7-d2e992e4beaf
md"""
#### Plotting all vol. Effect Definitions
"""

# ╔═╡ a9130f68-3908-4822-adfa-2bcce7d1c96d
md"""
##### (1) Based on inlet and outlet temperatures
"""

# ╔═╡ 1daeca94-9774-4c6f-ac51-1d104c071d81
findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1][1]

# ╔═╡ 2efd401c-ab68-4161-bda5-461754ea69b0
ε[findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1][2]]


# ╔═╡ 02f6c434-b6b4-457f-bc9d-6f682566aa5b
log2.(R[findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1][1]]) 


# ╔═╡ f4435862-38c7-40d9-bd4e-8958673ba89e
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

	# adding a point for the design that's closest to the definition's target
	scatter!(
		[log2.(R[findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1][1]])] ,
		[ε[findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1][2]]],
		label = ""

		)
	
	# adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :black,
	# 	# label="Direction to max(η)",
	# 	label="",
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

# ╔═╡ 26d01d52-129d-4ee2-97f6-3c6f4a0a4e4c
begin
	p2 = 
		contourf(
		log2.(R), ε, transpose(LMTD),
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

# ╔═╡ 3b15f409-c60d-4d82-aaa0-b79d3c6668b0
md"""
##### (2) Based on power ratios
"""

# ╔═╡ 7f2bd49f-803b-4a6b-9b27-6fb704c0692f
begin
	contourf(
		log2.(R), ε, transpose(vol_eff_p1),
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
		# plot!(
		# 	x1,y1,
		# 	linewidth = 2,
		# 	line = (:dash, :arrow),
		# 	color = :white,
		# 	label = false,
		# 	)
		# plot!(
		# 	[-0.1, 0.046], 
		# 	[0.9, 0.832],
		# 	color = :white,
		# 	linewidth = 2,
		# 	label="")
		# plot!(
		# 	[-0.1, 0.20], 
		# 	[0.9, 0.875],
		# 	color = :white,
		# 	linewidth = 2,
		# 	label="")
	
		# adding an arrow for direction to maximize volumetric effect defs.
		# plot!(
		# 	[-1, 3], 
		# 	[0.1, 0.9],
		# 	line = (:dash, :arrow),
		# 	color = :white,
		# 	linewidth = 2,
		# 	label="",
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

# ╔═╡ 891cfd4a-2470-423d-b3d4-d3101f6bf470
begin
	contourf(
		log2.(R), ε, transpose(vol_eff_p2),
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
		# # adding an arrow for direction to maximize efficiency
		# plot!(
		# 	x1,y1,
		# 	linewidth = 2,
		# 	line = (:dash, :arrow),
		# 	color = :white,
		# 	label = false,
		# 	)
		# plot!(
		# 	[-0.1, 0.046], 
		# 	[0.9, 0.832],
		# 	color = :white,
		# 	linewidth = 2,
		# 	label="")
		# plot!(
		# 	[-0.1, 0.20], 
		# 	[0.9, 0.875],
		# 	color = :white,
		# 	linewidth = 2,
		# 	label="")
end

# ╔═╡ 735b39ad-8c45-497b-b758-13f2663f7adc
begin
	contourf(
		log2.(R), ε, transpose(vol_eff_p3),
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
	# # adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
end

# ╔═╡ d7220c55-6d7a-40ea-97b7-51c8b6882d66
begin
	contourf(
		log2.(R), ε, transpose(vol_eff_p4),
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
	# # adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
end

# ╔═╡ 33cc5659-032d-49d9-b002-638330a82213
md"""
##### (3) Based on statistical evaluations of ∆T driving force
"""

# ╔═╡ 15fc0f33-fe14-4264-b23e-10703148c4f0
begin
	contourf(
		log2.(R), ε, transpose(RMSE),
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
	
	# # adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")

	# # adding an arrow for direction to maximize volumetric effect defs.
	# x2 = range(-1, 3, length = 20)
	# # y2 = 0.1*log10.(x2 .+ 1) .+ 0.2
	# y2 = -10 * exp.(-1.75 * (x2 .+ 2.6)) .+ 0.7
	# # plot!(
	# # 	[-1, 3], 
	# # 	[0.1, 0.9],
	# # 	line = (:dash, :arrow),
	# # 	color = :white,
	# # 	linewidth = 2,
	# # 	label="",
	# # 	)	#linear direction
	# plot!(
	# 	x2, 
	# 	y2,
	# 	line = (:dash, :arrow),
	# 	color = :black,
	# 	linewidth = 2,
	# 	label="",
	# )  		#curved direction
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

# ╔═╡ b2feba5d-d69a-4913-a2a7-717b5648403c
begin
	contourf(
		log2.(R), ε, -transpose(R_corr),
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
	
		# # adding an arrow for direction to maximize efficiency
		# plot!(
		# 	x1,y1,
		# 	linewidth = 2,
		# 	line = (:dash, :arrow),
		# 	color = :white,
		# 	label = false,
		# 	)
		# plot!(
		# 	[-0.1, 0.046], 
		# 	[0.9, 0.832],
		# 	color = :white,
		# 	linewidth = 2,
		# 	label="")
		# plot!(
		# 	[-0.1, 0.20], 
		# 	[0.9, 0.875],
		# 	color = :white,
		# 	linewidth = 2,
		# 	label="")
end

# ╔═╡ 5b4f24d1-7857-4b38-b059-7ef4af18ff1c
begin
	contourf(
		log2.(R), ε, transpose(NMB),
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
	
# 	# adding an arrow for direction to maximize efficiency
# 	plot!(
# 		x1,y1,
# 		linewidth = 2,
# 		line = (:dash, :arrow),
# 		color = :white,
# 		label = false,
# 		)
# 	plot!(
# 		[-0.1, 0.046], 
# 		[0.9, 0.832],
# 		color = :white,
# 		linewidth = 2,
# 		label="")
# 	plot!(
# 		[-0.1, 0.20], 
# 		[0.9, 0.875],
# 		color = :white,
# 		linewidth = 2,
# 		label="")

# 	# adding an arrow for direction to maximize volumetric effect defs.
# 	# plot!(
# 	# 	x2, 
# 	# 	y2,
# 	# 	line = (:dash, :arrow),
# 	# 	color = :black,
# 	# 	linewidth = 2,
# 	# 	label="",
# 	# )  		#curved direction
# 	plot!(
# 		[-1, -0.85], 
# 		[0.1, 0.15],
# 		color = :black,
# 		linewidth = 2,
# 		label="")
# 	plot!(
# 		[-1, -0.8], 
# 		[0.1, 0.12],
# 		color = :black,
# 		linewidth = 2,
# 		label="")
end

# ╔═╡ e3244257-4baf-47d3-8511-efb7b0aafca8
begin
	contourf(
		log2.(R), ε, transpose(NMSD),
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
	# # adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")

	# # adding an arrow for direction to maximize volumetric effect defs.
	# plot!(
	# 	[-1, 3], 
	# 	[0.1, 0.9],
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="",
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

# ╔═╡ 17f04945-b670-495c-bbd4-e2780da9cfc2


# ╔═╡ 3c620e52-15ef-46e2-8679-ce250782ccb6
md"""
#### Assessing Bases for Vol. Effect Justification
"""

# ╔═╡ b57b927e-e2a6-4496-97ba-7e34765c89bb
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

# ╔═╡ 45003d50-62c3-47a4-a638-51dfc75d41aa
# Based on minimizing local solid temperature gradients
begin
	
	IOTS = Array{Any, 2}(undef, (length(R),length(E)))
	IOTG = Array{Any, 2}(undef, (length(R),length(E)))

	MAXTS = Array{Any, 2}(undef, (length(R),length(E)))
	MAXTG = Array{Any, 2}(undef, (length(R),length(E)))	
	
	DT_S_DZ = Array{Any, 2}(undef, (length(R),length(E)))
	DT_G_DZ = Array{Any, 2}(undef, (length(R),length(E)))

	
	for r = 1:length(R),  e = 1:length(E)
		T_s = T_s_int[r,e]
		T_g = T_g_int[r,e]
		Z = Z_int[r,e]
		
		N = length(Z)
		
		dT_s_dz = line_grad(T_s, Z)
		dT_g_dz = line_grad(T_g, Z)
		
		# display(plot(Z, dT_s_dz, title = "R=$(R[r])mm, ε=$(E[e])"))
		
		DT_S_DZ[r,e] = line_avg(dT_s_dz, Z, 1, N)
		DT_G_DZ[r,e] = line_avg(dT_g_dz, Z, 1, N)
		
		IOTS[r,e] = maximum(T_s) - minimum(T_s)
		IOTG[r,e] = maximum(T_g) - minimum(T_g)
					
	end

	md"""

	Calculation of avg. local solid temperature gradients
	"""
end

# ╔═╡ 8da0dadc-bc72-454b-8388-1f20fa102884


# ╔═╡ 06d9b17f-4cc3-4010-b292-1a65a997d5e1
begin
	contourf(
		log2.(R), ε, transpose(SIG_S./maximum(SIG_S)),
		# log2.(R), ε, transpose(SIG_S),
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
	
	# # adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
end

# ╔═╡ 7e883bd7-4c50-47a5-8310-6f9735bfa0cc
begin
	contourf(
		log2.(R), ε, transpose(DT_S_DZ),
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
	# # adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
end

# ╔═╡ 10672baf-dea8-4812-a178-ec2771b8dd1c
md"""
##### (4) Based on maximizing driving force uniformity
"""

# ╔═╡ adb9e601-973e-4223-8c92-cf09b4c5d638
σ_DF = std.(
	[T_s_int[r,e] .- T_g_int[r,e]  for r = 1:length(R), e = 1:length(E)]
)#./SIG_G

# ╔═╡ 99307a96-a115-4335-8631-1c37b6faae87
begin
	contourf(
		# log2.(R), E, transpose(σ_DF./maximum(σ_DF)),
		log2.(R), ε, transpose(σ_DF),
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
	# # adding an arrow for direction to maximize efficiency
	# plot!(
	# 	x1,y1,
	# 	linewidth = 2,
	# 	line = (:dash, :arrow),
	# 	color = :white,
	# 	label = false,
	# 	)
	# plot!(
	# 	[-0.1, 0.046], 
	# 	[0.9, 0.832],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
	# plot!(
	# 	[-0.1, 0.20], 
	# 	[0.9, 0.875],
	# 	color = :white,
	# 	linewidth = 2,
	# 	label="")
end

# ╔═╡ 9e771134-bb10-40e5-aa61-93afdea594ee
md"""
#### Plotting temperature profiles along line that maximizes η
"""

# ╔═╡ 017acb44-bf48-4ae8-8f13-7e4a29a76a22
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

# ╔═╡ 726a6a6a-ac12-4f86-a894-65aee6429f29
colors_max_eta = (permutedims(palette(:reds, 4)[:]))

# ╔═╡ 662ec8de-38f5-4835-b8f0-db6b67a68fe6
colors_max_vol = reverse(permutedims(palette(:blues, 4)[:]))

# ╔═╡ b58a6822-a678-45c8-a0a6-57356d0c8d0e
md"""
### Plotting _η_ vs. Volumetric Effect Definitions along a line the Maximizes _η_
"""

# ╔═╡ c6681122-ddab-45cb-b44c-4b9972a2f689
E_int = 0.5001 * ones(length(R)-1); R_int = R[1:end-1]

# ╔═╡ e8ced466-7773-41f8-bed6-20e3057134e4
begin
	T_s_inter = Array{Any, 1}(undef, (length(R_int)))
	T_g_inter = Array{Any, 1}(undef, (length(R_int)))
	
	for r = 1:length(R_int)
		Z = Z_int[r,1]
		
		T_s = cat([vec(T_s_int[r,e]) for e = 1:length(E)]; dims = 2)
		T_s = transpose([T_s[e][z] for e = 1:length(E), z = 1:length(Z)])
		T_s_inter[r] = interpolate_lin2(T_s, E, E_int[r])
		
		T_g = cat([vec(T_g_int[r,e]) for e = 1:length(E)]; dims = 2)
		T_g = transpose([T_g[e][z] for e = 1:length(E), z = 1:length(Z)])
		T_g_inter[r] = interpolate_lin2(T_g, E, E_int[r])
	end

	md"""
	Interpolating for temperature profiles of selected designs
	"""
end

# ╔═╡ 5623dd49-2018-4e83-8557-20e043c4d319
begin
	# r1 = 2
	# ld1 = 1
	plot(
			[Z_int[r,1] for r = 1:length(R_int)],

			[T_g_inter[r] for r = 1:length(R_int)],
		color = colors_max_eta,
		# color = permutedims(palette(:tab10)[1:length(E_int)]),
		line = :dash,
		linewidth = 2,
		label = permutedims(["R<sub>ch</sub> = $(R_int[r]) mm" for r = 1:length(R_int)]),

		# label = false
		)

	plot!(
			[Z_int[r,1] for r = 1:length(R_int)],

			[T_s_inter[r] for r = 1:length(R_int)],
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

# ╔═╡ 6658af7c-02f9-4655-b30e-5c948cda3317
E_vol2 = [T_g_inter[r][end] / T_s_inter[r][1] for r = 1:length(R_int)]

# ╔═╡ 7b5ed0d5-d2be-4709-905f-b2fec6589d89
begin
	dT1 = [T_s_inter[r][1] - T_g_inter[r][1] for r = 1:length(R_int)]
	
	dT2 = [T_s_inter[r][end] - T_g_inter[r][end] for r = 1:length(R_int)]
	
	LMTD2 = (dT1 .- dT2)./log.(dT1 ./ dT2)
	
	md"""
	LMTD for selected designs
	"""
end

# ╔═╡ 1361caf4-eb56-4745-b39d-8b0a495927dd
begin
	# Q_abs_inter = interpolate_lin(Q_abs_tot[2:end,:], E, E_int)
	# Q_rad_loss_inter = interpolate_lin(Q_rad_loss_tot[2:end,:], E, E_int)
	# BHS_inter = interpolate_lin(BHS_tot[2:end,:], E, E_int)
	
	Q_abs_inter = Array{Any, 1}(undef, (length(R_int)))
	Q_rad_loss_inter = Array{Any, 1}(undef, (length(R_int)))
	BHS_inter =  Array{Any, 1}(undef, (length(R_int)))
	
	for r = 1:length(R_int)
		Q_abs = reshape(Q_abs_tot[r,:], (1,length(E)))
		Q_abs_inter[r] = interpolate_lin2(Q_abs, E, E_int[r])

		Q_rad = reshape(Q_rad_loss_tot[r,:], (1,length(E)))
		Q_rad_loss_inter[r] = interpolate_lin2(Q_rad, E, E_int[r])

		Q_BHS = reshape(BHS_tot[r,:], (1,length(E)))
		BHS_inter[r] = interpolate_lin2(Q_BHS, E, E_int[r])		
		
	end

	Q_abs_inter = [Q_abs_inter[r][1] for r = 1:length(R_int)]
	Q_rad_loss_inter = [Q_rad_loss_inter[r][1] for r = 1:length(R_int)]
	BHS_inter = [BHS_inter[r][1] for r = 1:length(R_int)]

	md"""
	Calculation of interpolated powers for selected designs
	"""
end

# ╔═╡ ca0bc822-741a-476f-b0e1-e1a97edadbfe
σ_DF2 = std.(
	[T_s_inter[r] .- T_g_inter[r]  for r = 1:length(R_int)]
)#./SIG_G

# ╔═╡ 98ecc826-152d-4cd5-af90-9e6c26ad9fea
begin
	Eff_inter = Array{Any, 1}(undef, (length(R_int)))
	
	for r = 1:length(R_int)
		Eff = reshape(Eff_Ap[r,:], (1,length(E)))
		Eff_inter[r] = interpolate_lin2(Eff, E, E_int[r])
	end

	Eff_inter = [Eff_inter[r][1] for r = 1:length(R_int)]

	md"""
	Calculation of $η$ for selected designs


	"""
end

# ╔═╡ 6fbe124b-2a1d-4f61-9256-d24c8b7ebca3
begin
	surface(
		log2.(R), ε, transpose(Eff_Ap),
		
		xlabel = "log<sub>2</sub>(R<sub>ch</sub> - mm)",
		ylabel = "Avg. Emissivity",
		zlabel = "Thermal Efficinecy",

		zaxis = 0:0.2:1.0,
		tickfontsize = 10,
	    guidefont = 12,
		# titlefont = 18,
		fontfamily = "ComputerModern",
		colorbar = :false,
		size = 2 .* (320,250)
		)
	scatter3d!(log2.(R_int), E_int, Eff_inter, color = :black, label = "")
	# plot3d!(log2.(R_int), E_int, Eff_inter, color = :black, label = "")
end

# ╔═╡ d547c6b0-09ae-4cc0-8d8e-a21aadf0773b
begin
	vol_eff_p12 = Q_abs_inter./Q_rad_loss_inter
	vol_eff_p22 = Q_rad_loss_inter./BHS_inter
	vol_eff_p32 = Q_abs_inter./BHS_inter

	md"""
	Calculating vol. effect definitions based on power ratios for selected designs

	"""
end

# ╔═╡ 56caee93-7377-4827-927b-5eba2899ccf6
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

	
	for r = 1:length(R_int) 
		T_s = T_s_inter[r]
		T_g = T_g_inter[r]
		Z = Z_int[r,ld]
		
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

	md"""
	Calculating statistical vol. effect definitions for selected designs

	"""
end

# ╔═╡ e2da42c7-7ac2-47e9-b13f-2b14649dbea1


# ╔═╡ a86195ea-2bc5-4784-94e3-8cc342adc505
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

# ╔═╡ c8bb1d22-0718-4652-aab7-725c0bfa08d6
labels = ["dTₛ/dz", "σₛ" , "T<sub>g,out</sub>/T<sub>s,in</sub>" , "LMTD" ,"R_corr", "RMSE", "NMB"  ,"NMSD", "Q<sub>abs,g</sub>/Q<sub>rad,loss</sub>", "Q<sub>rad,loss</sub>/Q<sub>BHS</sub>", "Q<sub>abs,g</sub>/Q<sub>BHS</sub>", "σ_∆T" ] 


# ╔═╡ 0c673fda-0a19-44a4-abd0-da8a502ee223
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

		xaxis = range(0.90, 0.96, step=0.02),
		xlim = [0.90, 0.96],

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

# ╔═╡ ed983369-71fd-4553-8eda-0e782addd1b2
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

		yaxis = range(0.90, 0.96, step=0.02),
		ylim = [0.90, 0.96],
		
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

# ╔═╡ 0d6b1cca-051b-4f3f-a82d-3bb833c3d0bf
begin
	plot(
		[
			E_vol2,# ./ maximum(E_vol2), 
			LMTD2 ./ maximum(LMTD2), 
		],
		Eff_inter,

		label = permutedims(labels[3:4]),
		
		color = permutedims(palette(:tab10)[1:2]),

		yaxis = range(0.90, 0.96, step=0.02),
		ylim = [0.90, 0.96],	
		
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

# ╔═╡ 18cc2de4-73b9-4800-89cd-ca4f8a02fb89
begin
	plot(
		[
			# R_corr2,
			RMSE2 ./ maximum(RMSE2), 
			# RMSE2,
			NMB2, 
			NMSD2,
			
			# R_corr2 ./ maximum(R_corr2), 
			# RMSE2 ./ maximum(RMSE2), 
			# NMB2 ./ maximum(NMB2), 
			# NMSD2 ./ maximum(NMSD2), 
		],
		Eff_inter,

		# label = permutedims(labels[5:8]),
		label = permutedims(labels[6:8]),
		
		color = permutedims(palette(:tab10)[1:4]),


		yaxis = range(0.90, 0.95, step=0.02),
		ylim = [0.90, 0.96],
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

# ╔═╡ 42c6d9f2-a2d4-4bf4-a95b-bb439de63266
begin
	plot(
		[
			# σ_DF2, 
			σ_DF2 ./ maximum(σ_DF2), 

		],
		Eff_inter,


		yaxis = range(0.90, 0.95, step=0.02),
		ylim = [0.90, 0.96],

		xaxis = range(0.85, 1.0, step=0.05),
		xlim = (0.85,1.0),
		
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

# ╔═╡ b4ed0cf6-d0c5-4024-a042-78724c111689
md"""

All definitions except the three _selected_ in the previous stage follow  the trend of $η$.

However, $T_{g,out}/T_{s,in}$ follows the trend along other $η(R,ε)$ cutlines. 

This might be an indicator that the design selection method is not suitable for comparison.
"""

# ╔═╡ 7c0de0d2-2878-407b-90d7-a7daaa4ba2a4


# ╔═╡ 4f23bf74-a90f-40df-a621-32e9bea7b64c
md"""
An alternative to selecting a set of designs (a cutline) is to overlay the maximum/closest target for each definition on top of the efficiency map
"""

# ╔═╡ c3d18dd2-aaaf-4aca-a068-8461c063df59
md"""
### Plotting the design for the maximum of each definition over efficiency
"""

# ╔═╡ f8a35c85-f238-439d-ba0e-04605a36e1d2
begin
	# target = Array{CartesianIndex{2}, 1}(undef, 13)

	# target[1] = findall(abs.(DT_S_DZ) .== minimum(abs.(DT_S_DZ)))[1]
	# target[2] = findall(SIG_S .== minimum(SIG_S))[1]
	# target[3] = findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1]
	# target[4] = findall(LMTD .== maximum(LMTD))[1]
	# target[5] = findall(abs.(R_corr) .== maximum(abs.(R_corr)))[1] #closest value to 1
	# target[6] = findall(RMSE .== minimum(RMSE))[1] #closest value to 0
	# target[7] = findall(NMB .== minimum(NMB))[1] #closest value to 0
	# target[8] = findall(abs.(NMSD) .== minimum(abs.(NMSD)))[1] #closest value to 0
	# target[9] = findall(vol_eff_p1 .== maximum(vol_eff_p1))[1]
	# target[10] = findall(vol_eff_p2 .== minimum(vol_eff_p2))[1]
	# target[11] = findall(vol_eff_p3 .== maximum(vol_eff_p3))[1]
	# target[12] = findall(vol_eff_p4 .== maximum(vol_eff_p4))[1]
	# target[13] = findall(σ_DF .== minimum(σ_DF))[1]

	
	# target = Array{CartesianIndex{2}, 1}(undef, 12)

	# target[1] = findall(abs.(DT_S_DZ) .== minimum(abs.(DT_S_DZ)))[1]
	# target[2] = findall(SIG_S .== minimum(SIG_S))[1]
	# target[3] = findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1]
	# target[4] = findall(LMTD .== maximum(LMTD))[1]
	# target[5] = findall(abs.(R_corr) .== maximum(abs.(R_corr)))[1] #closest value to 1
	# target[6] = findall(RMSE .== minimum(RMSE))[1] #closest value to 0
	# target[7] = findall(NMB .== minimum(NMB))[1] #closest value to 0
	# target[8] = findall(abs.(NMSD) .== minimum(abs.(NMSD)))[1] #closest value to 0
	# target[9] = findall(vol_eff_p1 .== maximum(vol_eff_p1))[1]
	# target[10] = findall(vol_eff_p2 .== minimum(vol_eff_p2))[1]
	# target[11] = findall(vol_eff_p3 .== maximum(vol_eff_p3))[1]
	# target[12] = findall(σ_DF .== minimum(σ_DF))[1]

	
	target = Array{CartesianIndex{2}, 1}(undef, 10)
	
	target[1] = findall(vol_eff[:,ld,:] .== maximum(vol_eff[:,ld,:]))[1]
	target[2] = findall(LMTD .== maximum(LMTD))[1]
	target[3] = findall(RMSE .== minimum(RMSE))[1] #closest value to 0
	target[4] = findall(NMB .== minimum(NMB))[1] #closest value to 0
	target[5] = findall(abs.(NMSD) .== minimum(abs.(NMSD)))[1] #closest value to 0
	target[6] = findall(vol_eff_p1 .== maximum(vol_eff_p1))[1]
	target[7] = findall(vol_eff_p2 .== minimum(vol_eff_p2))[1]
	target[8] = findall(vol_eff_p3 .== maximum(vol_eff_p3))[1]
	target[9] = findall(vol_eff_p4 .== maximum(vol_eff_p4))[1]
	target[10] = findall(σ_DF .== minimum(σ_DF))[1]
end

# ╔═╡ 160f1526-f6e5-42bd-b7eb-c33e1fdd0dc5
labels3 = cat(labels[3:4], labels[6:11], ["Q<sub>BHS,g</sub>/P<sub>ap</sub>"], "σ<sub>∆T</sub>"; dims = 1)

# ╔═╡ 54d9e648-aac6-4ce0-a37a-ca68e79e7e01
labels3[8] = "η<sub>conv</sub>"; labels3[9] = "η<sub>opt</sub>"

# ╔═╡ 5d501574-063b-401c-bdf1-60dee8224e09
labels2 = cat(copy(labels[1:11]), ["Q<sub>BHS</sub>/P<sub>ap</sub>"], [labels[12]]; dims = 1)

# ╔═╡ 7284d794-de2d-4ec1-aff9-f87aface8c98
labels2[13] = "σ<sub>∆T</sub>";  labels2[5] = "R<sub>corr</sub>";

# ╔═╡ ad0762c2-a0ed-433c-b35a-f038f2f6fcda
# marker_shapes = [:hex  :hex  :square  :square  :diamond  :diamond  :diamond  :diamond  :circle  :circle  :circle  :circle  :utriangle]

# marker_shapes = [:hex  :hex  :square  :square  :diamond  :diamond  :diamond  :diamond  :circle  :circle  :circle  :utriangle]

marker_shapes = [:circle  :circle  :square  :square  :square  :diamond  :diamond  :diamond  :diamond  :utriangle]

# ╔═╡ 6c4de3e9-0803-4e56-bd53-64e74cef62de
[[ε[target[i][2]]] for i = 1:length(target)]

# ╔═╡ c55d4d1d-c24e-4126-a286-2836332f7a16
begin
	path3d(
		R, 
		[ε[i]*ones(length(R)) for i = 1:length(E)],
		[Eff_Ap[:, i] for i = 1:length(E)],
		#---------------------------------------------
		label = permutedims(["ε = $(round(ε[i]; digits=3))" for i = 1:length(ε)]),

		ylabel = "ε",
		xlabel = "log<sub>2</sub>(R - mm)",
		zlabel = "η",

		# size = (624,477), #for ppt
		#--------------------------------------------
		# Size formatting for word
		linewidth = 2,
		top_margin = 5*Plots.mm,
		left_margin = 3*Plots.mm,
	
		# legend = (1.0, 0.5),
		fontfamily = "ComputerModern",
		colorbar_font = "ComputerModern",
		size = (440,335),	#for word doc.
			
		)

end

# ╔═╡ 4abf0f3c-e989-440f-8e3a-e44d5c6ee05a
vol_eff[:,1,:]

# ╔═╡ 8d4d1129-7273-4ec7-8e2b-22619dd7e022
begin
	contourf(
		log2.(R), ε, transpose(Eff_Ap),
		# log2.(R), E_base, transpose(Eff_Ap_int),
		#---------------------------------------------
		# xaxis = :log2,
		xlabel = "log<sub>2</sub>(Radius - mm)",
		ylabel = "Average Emissivity",
		
		# title = "Thermal Efficiency, η",
		# title = "(A)", titlelocation = :left,

		colorbar_title = " η<sub> thermal<sub>",
		colorbar_titlefontsize = 16,
	    colorbar_title_location = :right,           # also :left or :right
	    colorbar_titlefontfamily = :match,
		# color = :grays,
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
		)
	
	scatter!(
		[log2(R[target[1][1]])],
		[ε[target[1][2]]],

		# yaxis = :false,
		# xaxis = :false,
		
		# label = permutedims(labels2),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = "max($(labels3[1]))",
		
		color = :gray,
		markershape = marker_shapes,
		markerstrokewidth = 1,
		markerstrokecolor  = :white,
		markersize = 8,

		legendfontsize = 13,
		legend = (1.3,0.9)
		)

	scatter!(
		[log2(R_int[end])],
		[ε[4]],

		marker = :xcross,
		color = :black,
		markersize = 5,
		label = "max(η<sub>th</sub>)",
		legend = false,
		)
end

# ╔═╡ b2166f6a-0725-4604-89ae-b6f03593117c
begin
	contour(
		log2.(R), ε, transpose(BHS_int_tot./BHS_front_tot),
		# log2.(R), E_base, transpose(Eff_Ap_int),
		#---------------------------------------------
		# xaxis = :log2,
		xlabel = "log<sub>2</sub>(Radius - mm)",
		ylabel = "Average Emissivity",
		
		# title = "Thermal Efficiency, η",
		# title = "(A)", titlelocation = :left,

		colorbar_title = " BHS<sub>int</sub>/BHS<sub>front</sub>",
		colorbar_titlefontsize = 16,
	    colorbar_title_location = :right,           # also :left or :right
	    colorbar_titlefontfamily = :match,
		# color = :grays,
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
		)
	
	scatter!(
		[log2(R[target[1][1]])],
		[ε[target[1][2]]],

		# yaxis = :false,
		# xaxis = :false,
		
		# label = permutedims(labels2),
		# label = permutedims(cat(labels2[1:11],labels2[13]; dims=1)),
		label = "max($(labels3[1]))",
		
		color = :gray,
		markershape = marker_shapes,
		markerstrokewidth = 1,
		markerstrokecolor  = :white,
		markersize = 8,

		legendfontsize = 13,
		legend = (1.3,0.9)
		)

	scatter!(
		[log2(R_int[end])],
		[ε[4]],

		marker = :xcross,
		color = :black,
		markersize = 5,
		label = "max(η<sub>th</sub>)",
		legend = false,
		)
end

# ╔═╡ ca02df8c-48cb-4c9a-a307-cd0e702ec4bc
# colors2 = permutedims(palette([:red, :blue], 12)[:])
# colors2 = permutedims(palette(:solar, 13)[:])
# colors2 = permutedims(palette(:solar, 12)[:])
colors2 = permutedims(palette(:solar, 11)[2:end])

# ╔═╡ 5e083ba9-2ee3-449f-8045-d620e071366c
begin
	contourf(
		log2.(R), ε, transpose(Eff_Ap),
		# log2.(R), E_base, transpose(Eff_Ap_int),
		#---------------------------------------------
		# xaxis = :log2,
		xlabel = "log<sub>2</sub>(Radius - mm)",
		ylabel = "Average Emissivity",
		
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
		[[ε[target[i][2]]] for i = 1:length(target)],

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
		[log2(R_int[end])],
		[ε[4]],

		marker = :xcross,
		color = :black,
		markersize = 5,
		label = "η<sub>max</sub>",
		)
end

# ╔═╡ d9acabeb-a742-41b3-a14d-777b915a321e
begin
	sel_range = cat(1:11, [13]; dims = 1)
	non_sel = setdiff(collect(1:length(target)), sel_range)
	
	bar(
		# [[labels2[i]] for i in sel_range],
		# [[Eff_Ap[target[i]]] for i in sel_range],
		[[labels3[i]] for i = 1:length(labels3)],
		[[Eff_Ap[target[i]]] for i = 1:length(target)],

		ylabel = "Thermal Efficiency, η",
		ylim = [0.5, 1.0],
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
		
		# annotations = (14.5, maximum(Eff_Ap)*1.035, Plots.text("<i> η<sub>max</sub> </i>", :left, "ComputerModern")),	
		annotations = (11.5, maximum(Eff_Ap)*1.035, Plots.text("<i> η<sub>max</sub> </i>", :left, "ComputerModern")),	
	)

	# bar!(
	# 	[[labels2[i]] for i in non_sel],
	# 	[[Eff_Ap[target[i]]] for i in non_sel],

	# 	ylabel = "Thermal Efficiency, η",
	# 	ylim = [0.5, 1.0],
	# 	xgrid = :none,
		
	# 	xrotation = 90,
	# 	bottom_margin = -120*Plots.mm,
	
	# 	color = reverse(colors2),
	# 	legend = false,
		
	# )
	
end

# ╔═╡ 7224e04f-4073-441d-9418-833d65d2409a
md"""
##### Observations:
- Statistical definitions (except $R_{corr}$) aim for very different extremes in terms of desgin, but all of them end up acheiving the same efficiency (→ _see temp. profiles_)

- It might make more sense to put $Q_{BHS}/P_{ap}$ and $Q_{rad,loss}/P_{ap}$ inplace of $Q_{abs}/Q_{rad,loss}$ and $Q_{rad,loss}/Q_{BHS}$ 

  - Apparently they're all different ( >_<)

- 
"""

# ╔═╡ f6d66ed3-b6da-4d0e-9ac4-8e155261c25e
begin
	plot(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_g_int[target[i]] for i = 1:length(target)],
		
		color = colors2,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3),
		)

	plot!(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_s_int[target[i]] for i = 1:length(target)],
		
		color = colors2,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		title = "Axial Temperature Profiles for <i>max(E<sub>vol</sub>)</i> <br> (Specular/Volumetric Case)",
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
		Z_int[4,1],

		[T_s_int[4,4], T_g_int[4,4]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		
		)
	

end

# ╔═╡ b689cb78-d2b8-436c-801f-a970c8d78d76
findall(Eff_Ap .== maximum(Eff_Ap))

# ╔═╡ ed834410-8ee5-47ad-910c-c98a68db6b14
md"""
selecting relevant definitions
"""

# ╔═╡ f22e1991-f04e-4008-a65c-a7f626b33795
begin
	plot(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_g_int[target[i]] for i = 1:length(target)],
		
		color = colors2,
		line = :dash,
		linewidth = 3,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 1:length(target)]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i = 1:length(target)],

		[T_s_int[target[i]] for i = 1:length(target)],
		
		color = colors2,
		line = :solid,
		linewidth = 3,
		# label = permutedims(labels2),
		label = permutedims(labels3),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Axial Temperature Profiles for <i>max(E<sub>vol</sub>)</i> <br> (Specular/Volumetric Case)",
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
		Z_int[4,1],

		[T_s_int[4,4], T_g_int[4,4]],

		color = :black,
		line = [:solid :dash],
		linewidth = 1.5,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		)
	

end

# ╔═╡ 1ebc7b96-b275-4957-8661-7c54159a5339
# colors3 = permutedims(palette(:Spectral_4)[:])
colors3 = permutedims(palette(:Dark2_7)[:])

# ╔═╡ 5c282dc4-22de-41d9-a867-1b65333923e8
begin
	plot(
		[Z_int[target[i][1],1] for i = 1:2],

		[T_g_int[target[i]] for i = 1:2],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 1:2]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i = 1:2],

		[T_s_int[target[i]] for i = 1:2],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3[1:2]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		# title = "(A)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		#--------------------------------------------
		# Size formatting for manuscript	
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# legendfont = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# size = (320,250).*1.5,
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 10*Plots.mm,
		left_margin = 3*Plots.mm,
		bottom_margin = 1*Plots.mm,
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
		)

	plot!(
		Z_int[4,1],

		[T_s_int[4,4], T_g_int[4,4]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		)
	

end

# ╔═╡ 60c3e5ab-4f74-4f6e-b408-787858a40bbc
begin
	plot(
		[Z_int[target[i][1],1] for i = 6:9],

		[T_g_int[target[i]] for i = 6:9],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 6:9]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i = 6:9],

		[T_s_int[target[i]] for i = 6:9],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3[6:9]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		# title = "(B)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		#--------------------------------------------
		# Size formatting for manuscript	
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# legendfont = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# size = (320,250).*1.5,
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 10*Plots.mm,
		left_margin = 5*Plots.mm,
		bottom_margin = 1*Plots.mm,
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
		)

	plot!(
		Z_int[4,1],

		[T_s_int[4,4], T_g_int[4,4]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		)
	

end

# ╔═╡ 704aaec2-9fa2-447b-8853-818245ae6f28
begin
	plot(
		[Z_int[target[i][1],1] for i = 3:5],

		[T_g_int[target[i]] for i = 3:5],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 3:5]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i = 3:5],

		[T_s_int[target[i]] for i = 3:5],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3[3:5]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		# title = "(C)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		#--------------------------------------------
		# Size formatting for manuscript	
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# legendfont = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# size = (320,250).*1.5,
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 10*Plots.mm,
		left_margin = 3*Plots.mm,
		bottom_margin = 1*Plots.mm,
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
		)

	plot!(
		Z_int[4,1],

		[T_s_int[4,4], T_g_int[4,4]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		)
	

end

# ╔═╡ 560ddbcb-3733-4d40-b983-a5e4f4fe80c3
begin
	plot(
		[Z_int[target[i][1],1] for i in [10]],

		[T_g_int[target[i]] for i in [10]],
		
		color = colors3,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 1:2]),		
		)

	plot!(
		[Z_int[target[i][1],1] for i in [10]],

		[T_s_int[target[i]] for i in [10]],
		
		color = colors3,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(labels3[10:end]),
		
		xlabel = "Channel Length (L<sub>ch</sub> - mm)",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		# title = "(D)", titlelocation = :left, topmargin = 8*Plots.mm,

		legend = (1.0, 0.7),
		#--------------------------------------------
		# Size formatting for manuscript	
		# top_margin = 5*Plots.mm,
		# left_margin = 3*Plots.mm,
		# tickfontsize = 12,
		# legendfont = 12,
		# guidefont = 14,
		# titlefont = 16,
		# fontfamily = "ComputerModern",
		# size = (320,250).*1.5,
		#--------------------------------------------
		# Size formatting for ppt
		top_margin = 10*Plots.mm,
		left_margin = 3*Plots.mm,
		bottom_margin = 1*Plots.mm,
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
		# size = (320,250).*1.5,
		)

	plot!(
		Z_int[4,1],

		[T_s_int[4,4], T_g_int[4,4]],

		color = :black,
		line = [:solid :dash],
		linewidth = 2,
		# label = ["<i>max(η)</i>"  "<i>max(η)</i>"],
		label = ["<i>max(η)</i>"  ""],
		)
	

end

# ╔═╡ 335e921c-13d7-4402-909a-e2ec85e98d9e
labels[end]= "σ_∆T"

# ╔═╡ 3079d456-c3aa-4c2a-9c7b-1ada1e64a54c
begin
	plot(
		Z_int[target[1][1],1]/maximum(Z_int[target[1][1],1]),

		T_g_int[target[1]],
		
		color = :gray,
		line = :dash,
		linewidth = 2,
		# label = permutedims(labels2),
		label = permutedims(["" for i = 1:2]),		
		)

	plot!(
		Z_int[target[1][1],1]/maximum(Z_int[target[1][1],1]),
		
		T_s_int[target[1]],
		
		color = :gray,
		line = :solid,
		linewidth = 2,
		# label = permutedims(labels2),
		label = "max($(labels3[1]))",
		
		xlabel = "Normalized Receiver Length",
		ylabel = "Temperature (K)",
		# title = "Based on Power ratios",
		# title = "(A)", titlelocation = :left, topmargin = 8*Plots.mm,

		grid = :none,
		framestyle = :box,
		
		legend = (0.58, 0.45),
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
		Z_int[4,1]/maximum(Z_int[4,1]),

		[T_s_int[4,4], T_g_int[4,4]],

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

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

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

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

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
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

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
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

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

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

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
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

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
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "d189c6d2004f63fd3c91748c458b09f26de0efaa"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.61.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

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
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

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
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "4f00cc36fede3c04b8acf9b2e2763decfdcecfa6"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.13"

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

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "db0eee9b3bb2b38ab2d94349a3b0272d0a68e21f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.8"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

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
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

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
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

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
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

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
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "e7523dd03eb3aaac09f743c23c1a553a8c834416"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.7"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

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
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

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
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

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
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

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
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

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
# ╠═7e296429-ca31-449c-a1c7-514001586851
# ╠═cdfd8cb1-7602-45ad-abc4-8d9ed53eee19
# ╠═1f6244a6-09a4-4e7c-9468-8aa6ff376070
# ╠═1d0c4527-1a3b-4767-b9ed-06df5054e5ae
# ╠═c62cb59a-864e-4a16-996f-5423fd08569b
# ╠═1124a3c4-2b53-448a-9d63-278da7a2d834
# ╠═75527708-4ca3-4277-a4fe-a8de6aa888a9
# ╠═819592b5-8152-4db4-8ced-7472c0381e38
# ╠═6536390f-6cee-44c2-9a88-c23df8205992
# ╠═0799d6a5-84fc-431c-83b1-249ce9805e12
# ╠═da408ddb-5b01-4376-ad34-c19c3177910d
# ╠═da2add34-b3f8-4cab-915a-bfe7e378246a
# ╠═b79e969b-f842-4eb9-a4a6-b6585d33dcef
# ╠═80f21ac4-6478-40f1-aad1-7f169a20dfcf
# ╠═a7d27b67-7b6d-49ef-b267-cb6942bd2d02
# ╟─31ee1595-be66-4d8a-9be1-96b2514dcc0c
# ╟─6e2d8d5b-627e-416b-8e84-79a3af3b384a
# ╟─42efe36f-025c-44ad-8890-a32a86f7b571
# ╟─8c024801-84be-4fff-aaee-3a55cc84b396
# ╟─fa856550-9cd1-496b-a2b5-2cf79c22dd33
# ╟─d33a36b0-1e13-4476-9e8b-bd1b99cd5061
# ╠═743091d1-3f24-4899-b694-2b5e3d7e9ed8
# ╟─b1b2f925-8f6f-4b23-bae0-c2c2550c5b53
# ╠═929c0b8e-50b0-4143-85d9-30b0810b1d01
# ╠═2240b754-2cd1-4056-897d-e2f99fa37e86
# ╠═382543f1-c574-47c2-93d0-e3643288c9a6
# ╠═dc08bd04-8b58-4cdc-97b8-45d268739074
# ╟─d9c232e5-bd1f-47bc-80c8-15f540aeadc7
# ╟─1c39de78-d278-4b89-924e-113693e6bfb8
# ╠═faec19c6-a584-4421-991b-3cb01abbbc8c
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
# ╟─2a730b37-d2d3-4b32-9c47-24fd5b122d95
# ╟─76f0b4fd-dc94-499d-8d41-c58cc2d3ff2d
# ╟─ee777a68-05da-4696-aa3e-d4fcddd3367a
# ╠═bd03bcfb-bd4a-43ce-bb38-8cbbb402fe55
# ╟─22ead363-fe5a-4b5a-8718-9a0c000d0c4d
# ╟─d996de9d-406f-4cb2-930b-b33909befc74
# ╠═46f66376-05a1-4047-b572-adee3b37d519
# ╠═7b9f1ec8-d09e-477d-ba5d-f2cc9c46503e
# ╠═0c4c7b1b-f5bd-46c1-acd0-173ca3c8ab02
# ╟─c917cb94-36a6-49f1-bf2c-7681f6172310
# ╟─1672b8ce-a4a5-44c5-aa7a-54af60a29451
# ╟─549fbc4f-b840-4a49-b0dc-2174d9d8468d
# ╟─a71d51c0-efdb-45bb-bb1b-821b1ee4bb81
# ╟─2c8cc361-2eb0-445f-95d7-7d65295a3087
# ╟─1dab20f7-dd45-42c1-9565-156b9d4f5009
# ╟─475c4317-427d-4558-923b-982662107ea2
# ╠═4a8e9b11-9aea-4dd7-b03f-b362860148aa
# ╟─c944cbc9-a9b6-4538-836d-eafea086189f
# ╠═d9902d65-9c44-4d64-8f55-71a6bb4e555a
# ╟─54f952ce-8f4b-4a0f-bf1e-117fc48d3996
# ╟─b3c91cbe-baa8-421a-8ff1-b2b0bca215fa
# ╠═be8c4bdf-501a-4f2c-8a0d-d7e034040b24
# ╠═da253a52-1954-4890-842b-b9444917deeb
# ╟─92bd2aae-62c7-4f84-9427-16e09a7bac14
# ╟─7466e9ae-0a09-4c7f-8067-baed695e80ec
# ╠═05d998a1-18ae-47bc-ad52-b931c987ef43
# ╟─c45b508a-f7e0-4a55-a16f-1a3ebcb88df6
# ╟─fdffa90d-b353-4837-9a24-c1a69afd5191
# ╟─fc079169-d9ba-4f8d-a57d-613e8f2c5d49
# ╠═e9065f52-64bf-45d2-b5da-e75a4b408607
# ╟─a943d8e2-7dfe-4e92-80ea-bb0a8caa15c1
# ╟─bd7dd4ab-ee29-4780-92a6-af77c440f6ba
# ╟─ec8a5e12-9eaf-4f42-989b-064d27169fdc
# ╟─bea1991d-b1cd-4de2-90b8-2030d147d13a
# ╟─ccbfa443-7eea-4986-a3c2-a00c26d779d8
# ╟─67ad79e5-dd7c-4361-8bb0-3078bee29a67
# ╠═b82ca70a-b701-42b2-a774-25c187811e85
# ╟─0e3968c3-da01-497d-832f-cb517797b8ce
# ╟─b68d4847-221f-460c-b8a9-c8688bce14ac
# ╟─7409dde6-64ba-4209-94ab-fd7e3a97feb7
# ╟─d066f9e0-768c-45ba-b029-154773517c34
# ╟─ad913c29-00e1-46ef-b079-66aa0d37c1b7
# ╟─3dc007bd-200a-4419-9122-86dd41e471f8
# ╟─66c0abdf-9820-4478-9a9b-c1f9927fe7d2
# ╟─c809e1e3-2edb-4640-bd5a-8e7c668efd0f
# ╟─dc0f8c5c-b709-4aa1-a96c-9494e3f8ce58
# ╟─1d75576d-06aa-48dc-8ed7-d2e992e4beaf
# ╟─a9130f68-3908-4822-adfa-2bcce7d1c96d
# ╠═1daeca94-9774-4c6f-ac51-1d104c071d81
# ╠═2efd401c-ab68-4161-bda5-461754ea69b0
# ╠═02f6c434-b6b4-457f-bc9d-6f682566aa5b
# ╟─f4435862-38c7-40d9-bd4e-8958673ba89e
# ╟─26d01d52-129d-4ee2-97f6-3c6f4a0a4e4c
# ╟─3b15f409-c60d-4d82-aaa0-b79d3c6668b0
# ╟─7f2bd49f-803b-4a6b-9b27-6fb704c0692f
# ╟─891cfd4a-2470-423d-b3d4-d3101f6bf470
# ╟─735b39ad-8c45-497b-b758-13f2663f7adc
# ╟─d7220c55-6d7a-40ea-97b7-51c8b6882d66
# ╟─33cc5659-032d-49d9-b002-638330a82213
# ╟─15fc0f33-fe14-4264-b23e-10703148c4f0
# ╟─b2feba5d-d69a-4913-a2a7-717b5648403c
# ╟─5b4f24d1-7857-4b38-b059-7ef4af18ff1c
# ╟─e3244257-4baf-47d3-8511-efb7b0aafca8
# ╠═17f04945-b670-495c-bbd4-e2780da9cfc2
# ╟─3c620e52-15ef-46e2-8679-ce250782ccb6
# ╟─b57b927e-e2a6-4496-97ba-7e34765c89bb
# ╟─45003d50-62c3-47a4-a638-51dfc75d41aa
# ╟─8da0dadc-bc72-454b-8388-1f20fa102884
# ╟─06d9b17f-4cc3-4010-b292-1a65a997d5e1
# ╟─7e883bd7-4c50-47a5-8310-6f9735bfa0cc
# ╟─10672baf-dea8-4812-a178-ec2771b8dd1c
# ╠═adb9e601-973e-4223-8c92-cf09b4c5d638
# ╟─99307a96-a115-4335-8631-1c37b6faae87
# ╟─9e771134-bb10-40e5-aa61-93afdea594ee
# ╟─017acb44-bf48-4ae8-8f13-7e4a29a76a22
# ╟─e8ced466-7773-41f8-bed6-20e3057134e4
# ╠═726a6a6a-ac12-4f86-a894-65aee6429f29
# ╠═662ec8de-38f5-4835-b8f0-db6b67a68fe6
# ╟─5623dd49-2018-4e83-8557-20e043c4d319
# ╟─b58a6822-a678-45c8-a0a6-57356d0c8d0e
# ╠═c6681122-ddab-45cb-b44c-4b9972a2f689
# ╠═6658af7c-02f9-4655-b30e-5c948cda3317
# ╟─7b5ed0d5-d2be-4709-905f-b2fec6589d89
# ╟─1361caf4-eb56-4745-b39d-8b0a495927dd
# ╠═ca0bc822-741a-476f-b0e1-e1a97edadbfe
# ╟─98ecc826-152d-4cd5-af90-9e6c26ad9fea
# ╟─6fbe124b-2a1d-4f61-9256-d24c8b7ebca3
# ╟─d547c6b0-09ae-4cc0-8d8e-a21aadf0773b
# ╟─56caee93-7377-4827-927b-5eba2899ccf6
# ╟─e2da42c7-7ac2-47e9-b13f-2b14649dbea1
# ╟─a86195ea-2bc5-4784-94e3-8cc342adc505
# ╠═c8bb1d22-0718-4652-aab7-725c0bfa08d6
# ╟─0c673fda-0a19-44a4-abd0-da8a502ee223
# ╟─ed983369-71fd-4553-8eda-0e782addd1b2
# ╟─0d6b1cca-051b-4f3f-a82d-3bb833c3d0bf
# ╟─18cc2de4-73b9-4800-89cd-ca4f8a02fb89
# ╟─42c6d9f2-a2d4-4bf4-a95b-bb439de63266
# ╟─b4ed0cf6-d0c5-4024-a042-78724c111689
# ╟─7c0de0d2-2878-407b-90d7-a7daaa4ba2a4
# ╟─4f23bf74-a90f-40df-a621-32e9bea7b64c
# ╟─c3d18dd2-aaaf-4aca-a068-8461c063df59
# ╠═f8a35c85-f238-439d-ba0e-04605a36e1d2
# ╠═160f1526-f6e5-42bd-b7eb-c33e1fdd0dc5
# ╠═54d9e648-aac6-4ce0-a37a-ca68e79e7e01
# ╟─5d501574-063b-401c-bdf1-60dee8224e09
# ╠═7284d794-de2d-4ec1-aff9-f87aface8c98
# ╠═ad0762c2-a0ed-433c-b35a-f038f2f6fcda
# ╠═6c4de3e9-0803-4e56-bd53-64e74cef62de
# ╟─c55d4d1d-c24e-4126-a286-2836332f7a16
# ╠═4abf0f3c-e989-440f-8e3a-e44d5c6ee05a
# ╟─5e083ba9-2ee3-449f-8045-d620e071366c
# ╟─8d4d1129-7273-4ec7-8e2b-22619dd7e022
# ╟─b2166f6a-0725-4604-89ae-b6f03593117c
# ╠═ca02df8c-48cb-4c9a-a307-cd0e702ec4bc
# ╟─d9acabeb-a742-41b3-a14d-777b915a321e
# ╟─7224e04f-4073-441d-9418-833d65d2409a
# ╟─f6d66ed3-b6da-4d0e-9ac4-8e155261c25e
# ╠═b689cb78-d2b8-436c-801f-a970c8d78d76
# ╟─ed834410-8ee5-47ad-910c-c98a68db6b14
# ╟─f22e1991-f04e-4008-a65c-a7f626b33795
# ╠═1ebc7b96-b275-4957-8661-7c54159a5339
# ╟─5c282dc4-22de-41d9-a867-1b65333923e8
# ╟─60c3e5ab-4f74-4f6e-b408-787858a40bbc
# ╟─704aaec2-9fa2-447b-8853-818245ae6f28
# ╟─560ddbcb-3733-4d40-b983-a5e4f4fe80c3
# ╠═335e921c-13d7-4402-909a-e2ec85e98d9e
# ╟─3079d456-c3aa-4c2a-9c7b-1ada1e64a54c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
