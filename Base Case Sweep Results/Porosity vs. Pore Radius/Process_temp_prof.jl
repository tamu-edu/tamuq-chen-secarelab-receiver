### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 9ac747d0-a46b-11ec-0f2a-6fcbe88bf23d
using DelimitedFiles, Statistics, LinearAlgebra

# ╔═╡ ccfedbeb-2448-438f-97d5-04306145a2a1
begin
	R = [0.5,  1.0,  2.0,  4.0,  8.0]
	LD = 25
	ϕ = [0.1,  0.3,  0.5,  0.7,  0.9]
end

# ╔═╡ 2d54c9e0-8e71-4d9c-a237-d1a471184bb0
folder = "Cutline Temp Profiles (unprocessed)"

# ╔═╡ 863d8415-456c-4a20-b86e-e30c57a9f856
begin
	T_gas = Array{Vector, 2}(undef, length(R),length(ϕ))
	x_gas = Array{Vector, 2}(undef, length(R),length(ϕ))
	
	T_sol = Array{Vector, 2}(undef, length(R),length(ϕ))
	x_sol = Array{Vector, 2}(undef, length(R),length(ϕ))
	
	for r = 1:length(R)
		A_g = readdlm("$(folder)\\T_gas (R=$(R[r])mm, LD=25).txt", comments= true, comment_char = '%')
		x_g = A_g[:,1]; T_g = A_g[:,2]
		
		A_s = readdlm("$(folder)\\T_solid (R=$(R[r])mm, LD=25).txt", comments= true, comment_char = '%')
		x_s = A_s[:,1]; T_s = A_s[:,2]
	
		gas_starts = findall(x_g .== 0.0)
		sol_starts = findall(x_s .== 0.0)
		for phi = 1:length(ϕ)
			if phi != length(ϕ)
				x_gas[r,phi] = x_g[gas_starts[phi]:gas_starts[phi+1]-1]
				T_gas[r,phi] = T_g[gas_starts[phi]:gas_starts[phi+1]-1]
				
				x_sol[r,phi] = x_s[sol_starts[phi]:sol_starts[phi+1]-1]
				T_sol[r,phi] = T_s[sol_starts[phi]:sol_starts[phi+1]-1]
			else
				x_gas[r,phi] = x_g[gas_starts[phi]:end]
				T_gas[r,phi] = T_g[gas_starts[phi]:end]
				
				x_sol[r,phi] = x_s[sol_starts[phi]:end]
				T_sol[r,phi] = T_s[sol_starts[phi]:end]
			end
		end
	end

end

# ╔═╡ fef03fc7-d94c-41c9-a590-7c4993e06fa5
function interpolate_lin_1D(U, z, z_int)
	#U = known 1D property array
	#z = levels at which 'U' is known 
	#z_int = levels at which 'U' is desired
	
	Z = length(U)
	U_int = zeros(length(z_int))
	
	for i = 1:length(z_int)
		
		if z_int[i] == z[1]
			U_int[i] = U[1]
		elseif z_int[i] ≈ z[end]
			U_int[i] = U[end]
		else
			r_ind = findfirst(z_int[i] .< z) - 1
			slope = (U[r_ind] - U[r_ind+1])/(z[r_ind] - z[r_ind+1])
			U_int[i] = slope*(z_int[i] - z[r_ind]) + U[r_ind]
		end
	end
	
	return U_int
end

# ╔═╡ 9175a4ba-ab05-415d-ac36-f595f8654816
# Extracting temperature profiles at the same 'z' points
begin
	T_s_int = Array{Any, 2}(undef, (length(R),length(ϕ)))
	T_g_int = Array{Any, 2}(undef, (length(R),length(ϕ)))
	Z_int = Array{Any, 1}(undef, length(R))
	
	for r = 1:length(R), e = 1:length(ϕ)
		L_channel = 2*R[r]*LD #mm
		z_int = range(0, L_channel, step = 1)
		Z_int[r] = z_int
		
		T_s_int[r,e] = interpolate_lin_1D(T_sol[r,e], x_sol[r,e], z_int)

		println("r = ", r, "  phi = ", e)
		println(length(z_int))
		
		println("---------------------- in loop")
		T_g_int[r,e] = interpolate_lin_1D(T_gas[r,e], x_gas[r,e], z_int)
	end
	
	#-----------------------------------------------------------------
	# Assembling vactors of the same length into arrays
	T_g = Array{Array,1}(undef, length(R))
	T_s = Array{Array,1}(undef, length(R))
	
	for r = 1:length(R)
		T_g[r] = cat(T_g_int[r,1], T_g_int[r,2], T_g_int[r,3], T_g_int[r,4], T_g_int[r,5]; dims = 2)

		T_s[r] = cat(T_s_int[r,1], T_s_int[r,2], T_s_int[r,3], T_s_int[r,4], T_s_int[r,5]; dims = 2)
	end
end

# ╔═╡ 57d39072-96b1-426c-bd77-d5b02a31d443
begin
	# for  r = 1:length(R), phi = 1:length(ϕ)
	# 	open("Cutline Temp Profiles\\T_gas (R=$(R[r])mm, phi=$(ϕ[phi])).txt", "w") do io
	# 		   writedlm(io, [x_gas[r,phi] T_gas[r,phi]])
	# 			end
	# 	open("Cutline Temp Profiles\\T_solid (R=$(R[r])mm, phi=$(ϕ[phi])).txt", "w") do io
	# 		   writedlm(io, [x_sol[r,phi] T_sol[r,phi]])
	# 			end
	# end

	md"""
	Writing Separarted Raw Temperature Profiles 
	"""
end

# ╔═╡ 5bee07b0-c958-4296-83e6-240786f41863
begin
	col_labels = cat(["%z[mm]"], ["ϕ=$(ϕ[e])" for e = 1:length(ϕ)]; dims = 1)
	col_labels = permutedims(col_labels)
	
	# for  r = 1:length(R)
	# 	open("Cutline Temp Profiles (processed)\\T_gas (R=$(R[r])mm, LD=25).txt", "w") do io
	# 		   writedlm(io, [col_labels; Z_int[r] T_g[r]])
	# 			end
		
	# 	open("Cutline Temp Profiles (processed)\\T_solid (R=$(R[r])mm, LD=25).txt", "w") do io
	# 		   writedlm(io, [col_labels; Z_int[r] T_s[r]])
	# 			end
	# end

	md"""
	Writing Assembled Temperature Profiles at 1mm increments

	"""
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╠═9ac747d0-a46b-11ec-0f2a-6fcbe88bf23d
# ╠═ccfedbeb-2448-438f-97d5-04306145a2a1
# ╠═2d54c9e0-8e71-4d9c-a237-d1a471184bb0
# ╠═863d8415-456c-4a20-b86e-e30c57a9f856
# ╟─fef03fc7-d94c-41c9-a590-7c4993e06fa5
# ╠═9175a4ba-ab05-415d-ac36-f595f8654816
# ╟─57d39072-96b1-426c-bd77-d5b02a31d443
# ╟─5bee07b0-c958-4296-83e6-240786f41863
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
