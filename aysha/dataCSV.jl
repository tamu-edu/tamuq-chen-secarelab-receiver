### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 316456f8-bb2a-11ed-0c9e-ab106d878adc
begin
	import Pkg
	Pkg.activate(".")
	
	using CSV
	using DataFrames
	using Plots, ColorSchemes
	plotly()

end

# ╔═╡ 9f193079-ca1e-4dd1-b077-8b208b31c70d
begin
	
	file="Data_FPT0026_230302_141014.csv"
	E26=CSV.read(file, DataFrame, delim=";", header=["Time","T3","T8"], skipto=2, ntasks = 1)
	dropmissing(E26)
	file="Data_FPT0027_230302_131154.csv"
	E27=CSV.read(file, DataFrame, delim=";", header=["Time","T3","T8"], skipto=2, ntasks = 1)
	dropmissing(E27)
	file="Data_FPT0028_230302_111716.csv"
	E28=CSV.read(file, DataFrame, delim=";", header=["Time","T3","T8"], skipto=2, ntasks = 1)
	dropmissing(E28)
end

# ╔═╡ 320b21e4-9b9f-429d-ba62-1496af802af8
begin
		T3=[E26.T3, E27.T3, E28.T3],
		T8=[E26.T8, E27.T8, E28.T8],
		time=[E26.Time, E27.Time, E28.Time]
end

# ╔═╡ 5b9b6cd4-92a8-4090-b808-aaa8d4eea773
plot(T3, T8, time, labels=string.(transpose(HR_A.rate)),lw=2, legend = :outertopright, xlabel = "Time [s]", ylabel = "Temperature [C]", color_palette= colors_A)

# ╔═╡ Cell order:
# ╠═316456f8-bb2a-11ed-0c9e-ab106d878adc
# ╠═9f193079-ca1e-4dd1-b077-8b208b31c70d
# ╠═320b21e4-9b9f-429d-ba62-1496af802af8
# ╠═5b9b6cd4-92a8-4090-b808-aaa8d4eea773
