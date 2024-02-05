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

# ╔═╡ e5ea9ee0-6954-11eb-2942-e76574442366
begin
	import Pkg
	Pkg.activate(".")
	
	using CSV
	using DataFrames
	using Plots, ColorSchemes
	using Loess
	using Optim
	using LsqFit
	using Statistics, StatsBase
	using DifferentialEquations
	using PlutoUI
	using GLM
	plotly()

	#using Flux, DiffEqFlux, Sundials, DiffEqCallbacks, LinearAlgebra
end

# ╔═╡ 51322bfd-c9b6-4a18-ac4d-7d9add4a004e
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 0fb798e8-de08-424a-88bf-1f0727faa828
md"# Get the data"

# ╔═╡ b733f4fe-0b25-4248-8f6f-d872f3fd868f
begin
	
	#files with different HRate
	
	file="DataTGA3/DataTGA3_28_15min_a.csv"
	HR15=CSV.read(file, DataFrame, delim=";", header=["time","TF","TS", "TG"], skipto=10, ntasks = 1)
	dropmissing(HR15)
	file="DataTGA3/DataTGA3_29_12min_a.csv"
	HR12=CSV.read(file, DataFrame, delim=";", header=["time","TF","TS", "TG"], skipto=10, ntasks = 1)
	dropmissing(HR12)
	file="DataTGA3/DataTGA3_30_02min_a.csv"
	HR02=CSV.read(file, DataFrame, delim=";", header=["time","TF","TS", "TG"], skipto=10, ntasks = 1)
	dropmissing(HR02)
	file="DataTGA3/DataTGA3_44_06min_a.csv"
	HR06=CSV.read(file, DataFrame, delim=";", header=["time","TF","TS", "TG"], skipto=10, ntasks = 1)
	dropmissing(HR06)
	file="DataTGA3/DataTGA3_45_04min_b.csv"
	HR04=CSV.read(file, DataFrame, delim=";", header=["time","TF","TS", "TG"], skipto=10, ntasks = 1)
	dropmissing(HR04)

end

# ╔═╡ 2d594622-4bd9-4a34-b938-7830aff31dd3
md"# Common Functions"

# ╔═╡ fe678aa2-6954-11eb-08d3-bdd150d83c49
begin

	function smooth(hrate, HR_0)
	
		# #HR.mass=ones(hrate)
		HR = deepcopy(HR_0)
		HR.massdiff=zeros(hrate)
		
		HR.a=Array{Any}(undef, hrate)
		HR.dadT=Array{Any}(undef, hrate)
			
		for i=1:hrate
			T = HR.T[i]
			TG = HR.TG[i]
			mTGA = HR.mTGA[i]
			
			limlow=findfirst(collect(T.>110)) #lower limit temperature
			limhigh=findfirst(collect(T.>550)) #upper limit temperature
			T=T[limlow:limhigh] .+ 273.
			TG=TG[limlow:limhigh]		
			a=(TG .- TG[1]) ./ mTGA
			a ./=minimum(a) #re normalize to account for small differences at reactio ending masses
			mT=convert(Array{Float64},T)
			mod=loess(mT, a, span=0.05, degree=1)
			a=Loess.predict(mod, mT)
			
			dadT= diff(a)./diff(T)
			append!(dadT, 0)
			df = DataFrame(T = T, TG = TG, a = a, dadT = dadT)
			
			filter!(row -> all(x -> !(x isa Number && isnan(x)), row), df)

			HR.T[i] = df.T
			HR.TG[i] = df.TG
			HR.a[i] = df.a
			HR.dadT[i] = df.dadT
		end
		return HR
	end
end

# ╔═╡ 09d02877-479d-480b-b503-19de160f8ac0
function plot_dadT(HR, SelPeaks, colors)
		t= plot(HR.T, HR.dadT, labels=string.(transpose(HR.rate)),lw=2, legend=:outertopright, xlabel="T [K]", ylabel="dadT", color_palette= colors)
		for i in 1:length(HR.rate)
			vline!(SelPeaks[i], label=string("Peak:", HR.rate[i]), color_palette=colors)
		end
		return	t
end

# ╔═╡ cb3d56c7-b08d-4de5-98eb-fc01bdcae591
md"# Block A - two peaks"

# ╔═╡ 47e9064b-0554-4928-ba89-372e0c1f4d02
begin
	hrate_A=2
	HR_0A= DataFrame(rate=[2.0, 4.0],
		T=[HR02.TS, HR04.TS], #select TS or TF
		TG=[HR02.TG, HR04.TG],
		time=[HR02.time, HR04.time])
	HR_0A.mTGA=[17.6, 17.11] # mass of sample as per TGA
end

# ╔═╡ 48293b1c-ce1f-4010-9d85-0bef313df73a
HR_A = smooth(hrate_A, HR_0A)


# ╔═╡ 0537a765-4048-4280-a40b-3ee64004f862
colors_A = colorschemes[:tol_muted][1:hrate_A]

# ╔═╡ 9e71b200-6955-11eb-1cc7-7dd941326057
plot(HR_A.T, HR_A.a, labels=string.(transpose(HR_A.rate)),lw=2, legend = :outertopright, xlabel = "Temperature [K]", ylabel = "alpha reaction extent [-]", color_palette= colors_A)

# ╔═╡ a6e6c3d0-6955-11eb-0231-9fb731da375d
begin
	pn1=2 #number of peaks
	wi_peak=[1., 1.]
	SelPeaks=[
		[597, 681],
		[646, 715],
		]

		plot_dadT(HR_A, SelPeaks, colors_A)

end

# ╔═╡ 2bd08850-6957-11eb-2578-ad33186981d0
begin
	function FS(xdata, p) #Fraser-Suzuki 
		#h: Amplitude, w: halfwidth, s: assymetry, xc: position
		h, xc, w, s =p
		@. f(x)= h * exp(- log(2)/s^2 * (log(max(1+2*s*(x-xc),0))/w)^2)
		return f(xdata)
	end
	function FSOrigin(xdata, p) #Fraser-Suzuki as in Origin, symmetric
		yo, A, sig, xc =p
		#xc=position
		@. f(x)= yo + A*exp(-(x-xc)^2/(2*sig*(1+A)))/(sig*sqrt(2*π))
		return f(xdata)
	end
	function FSn(xdata, p) #Fraser-Suzuki for multiple peaks
		ydata=zeros(length(xdata))
		pn= Int(length(p)/4)
		for i=1:pn
			ptemp =p[1+4*(i-1) : 4*i]
			ydata+=FS(xdata, ptemp) .* wi_peak[i]
		end
		return ydata
	end
	function SBred(αi, mi, ni, ci) #Sestak-Berggren reduced with c
		@. ci*(1-αi)^ni *αi^mi
	end
	function SBred(αi, mi, ni) #Sestak-Berggren reduced without c
		@. (1-αi)^ni *αi^mi
	end
	function SB(αi, mi, ni, pi) #Sestak-Berggren
		@. (1-αi)^ni * αi^mi *(- log(1-αi))^pi
	end
	function kT(T, A, E) #Arrehnius
		R=8.314
		A*exp(-E/R/(T))
	end
	md"Define Theory Functions"
end

# ╔═╡ 39570e90-6957-11eb-28d4-af76c02c6de1
#Function to fit the FS peaks
begin
	function fitFS2(SelPeaks, DF, pn)
		hrate=size(SelPeaks,1)
		thT=10.
		th=Inf
		asym=0.1
		hw=50
		FSpar=fill(zeros(4,pn), hrate)
		covar=zeros(hrate)
		for i=1:hrate
		#i=5
			SelPeak=SelPeaks[i]
			pp2=vcat([[0.17, SelPeak[ipk], 30., 0.001...] for ipk=1:pn]...)
			lb2=vcat([[0, SelPeak[ipk]-thT, -hw, -asym...] for ipk=1:pn]...)
			ub2=vcat([[th, SelPeak[ipk]+thT, hw, asym...] for ipk=1:pn]...)
			#pp2=[0.17, 30., 0.001, SelPeak[1],
			#	0.14, 30., 0.001, SelPeak[2]]
			#lb2=[pp2[1]*(1-th), -hw, -asym, pp2[4]*(1-thT),
			#	pp2[5]*(1-th), -hw, -asym, pp2[8]*(1-thT)]
			#ub2=[pp2[1]*(1+th), hw, asym, pp2[4]*(1+thT),
			#	pp2[5]*(1+th), hw, asym, pp2[8]*(1+thT)]	


			fitFSn=curve_fit(FSn, DF.T[i], DF.dadT[i], pp2, lower=lb2, upper=ub2)
			sigman = stderror(fitFSn)
			confidence_intern = confidence_interval(fitFSn, 0.05)
			FSpar[i]=transpose(reshape(fitFSn.param,4,pn))
			covar[i] = StatsBase.rmsd(DF.dadT[i], FSn(DF.T[i], fitFSn.param);normalize=true)
		end
		return covar, FSpar
	end
end

# ╔═╡ 43805fa0-69d1-11eb-1ec1-95a700d9e47f
#fit the SB model to the selected peak
begin
	function fitFa(p)
		mm, nn, pp = p
		LHS=log.(max.(tot_dadt, 0.)) - log.(max.(SB(tot_a, mm, nn, pp), 0.))		
		return StatsBase.cor(1 ./ tot_TF, LHS)
	end
	function LHS_T_corr(DFin)
		DF=copy(DFin)
		hrt=copy(DF.rate)
		#prepare the data
		for ihr=1:length(DF.T)
			TF=DF.T[ihr]
			dadT=DF.dadT[ihr]
			a_int=DF.a[ihr]
			a_int=cumsum(dadT[1:end-1] .*diff(TF)) #reconstruct a from diff
			push!(a_int, a_int[end])
			#c1=1/maximum(a_int)
			#a_int*=c1 
			dt=fit(UnitRangeTransform, a_int) #normalize a [0,1]
			a_int=StatsBase.transform(dt, a_int)
			dadT=diff(a_int) .*diff(TF) #reconstruct dadT from normalized
			push!(dadT, dadT[end])
			#modlhs= Loess.loess(TS, dadT, span=.1, degree=1)
        	#dadT= Loess.predict(modlhs, float(TS)) #smooth date to correct for negative rates
			idx=((a_int .>0.1*maximum(a_int)) .* (a_int .<0.9*maximum(a_int)))
			dadt=abs.(dadT[idx]*hrt[ihr]) #ABS is not verrified
			a_int=a_int[idx]
			TF=TF[idx]
			append!(tot_TF, TF)
			append!(tot_a, a_int)
			append!(tot_dadt, dadt)
			append!(lens, length(tot_TF))
		end
		upper=[1., 2., 5.] #upper limits of m, n, p
		lower=[0, 0., 0.] #lower limits of m, n, p
		init=[.5, .5, .5]
		res_lin= Optim.optimize(fitFa,lower, upper, init)
		m1, n1, p1 = res_lin.minimizer
		return m1, n1, p1, res_lin.minimum 
	end
	
	function KineticModel(HR)

		#DF=copy(HR)
		m1, n1, p1, res = LHS_T_corr(HR)
		tot_data = DataFrame(X=1 ./tot_TF, Y=log.(max.(tot_dadt, 0.)) - log.(max.(SB(tot_a, m1, n1, p1), 0.)))
		modln= GLM.lm(@formula(Y ~ X), tot_data)
		intercept, slope =GLM.coef(modln)
		A=exp(intercept)
		E=slope*8.314/1000
		return  m1, n1, p1, res, A, E, GLM.r2(modln)	
	end
	
	global tot_TF=[]
	global tot_a=[]
	global tot_dadt=[]
	global lens=[]
	

	md"Kinetic's fitting Function"

end

# ╔═╡ 05bce23a-c877-4dca-8570-b348a51eebfc
	function plot_deconv(HR, fullPar, FSPar, pn, hrt)
		pl_c=plot(HR.T[hrt],HR.dadT[hrt],lw=2, label="exp", legend=:outertopright, xlabel="T [K]", ylabel="dadT")
		plot!(HR.T[hrt], FSn(HR.T[hrt], fullPar),lw=3,ls=:dash, labels="num")
		for i=1:pn
			plot!(HR.T[hrt], FS(HR.T[hrt], FSPar[hrt][i,:]),label="peak " * string(i))
		end
		pl_c
	end

# ╔═╡ c3976722-a898-4261-bfdb-534042b4f4d4
md"## Deconvolution"

# ╔═╡ 48a38938-156d-44a1-9d03-0ac88bf1560b
covar, FSPar_A = fitFS2(SelPeaks, HR_A, pn1)

# ╔═╡ cf2287b0-6957-11eb-1d0e-4befdd8df2b3
md"Heating Rate: $(@bind hrt Slider(1:hrate_A, show_value=true))"

# ╔═╡ d6217440-6957-11eb-0573-6f4ddc7e5042
FSPar_A[hrt]

# ╔═╡ e8c61920-6957-11eb-37ab-4f54568b6b41
begin
	fullPar=reshape(transpose(FSPar_A[hrt]), pn1*4)

	plot_deconv(HR_A, fullPar, FSPar_A, pn1, hrt)

end

# ╔═╡ dc7b9552-6957-11eb-2fc7-55f9eab720f5
StatsBase.cor(HR_A.dadT[hrt],FSn(HR_A.T[hrt], fullPar))

# ╔═╡ 09472c70-6958-11eb-0f78-bbdffc69d8e0
md"Select Peak: $(@bind pk Slider(1:pn1, show_value=true))"

# ╔═╡ 1065d600-6958-11eb-2d59-8d98b52e0788
begin
	n=plot()
	for i=1:hrate_A
		n=plot!(HR_A.T[i], FS(HR_A.T[i], FSPar_A[i][pk,:]), label=HR_A.rate[i],xlabel="T", ylabel="dadT (original)")
	end
	n
end

# ╔═╡ 14098f8c-25e3-4263-b5c3-e341e45b9418
md"## Peak 1"

# ╔═╡ b816f330-6aec-11eb-3402-01980c895d22
begin
	#for peak 1
	empty!(tot_TF)
	empty!(tot_a)
	empty!(tot_dadt)
	empty!(lens)
	
	DF=copy(HR_A)
	DF.dadT=[FS(HR_A.T[hrt], FSPar_A[hrt][1,:]) for hrt in 1:hrate_A]
	m1, n1, p1, res1, A1, E1, lnres1 = KineticModel(DF)
end

# ╔═╡ 0db7cf46-66af-40b7-a2dd-6989acd9d4f5
DF.dadT

# ╔═╡ fae4e250-6aae-11eb-0763-1b5a759b100a
begin
	
	plot(xlabel="1/T", ylabel="LHS")
	yall=log.(tot_dadt) - log.(SB(tot_a, m1, n1, p1))
	for i=1:hrate_A
		i==1 ? start=1 : start=lens[i-1]+1
		x=1 ./tot_TF[start:lens[i]]
		y=yall[start:lens[i]]
		scatter!(x, y, label=string(HR_A.rate[i])) 
	end
	
	plot!(1 ./tot_TF, (log(A1) .+(E1/8.314*1000 ./tot_TF)),lw=4)
end

# ╔═╡ c65d41a0-07e9-4f0c-8a79-ddf0498e237f
md"## Peak 2"

# ╔═╡ 7f734920-6aec-11eb-09f2-49a95f3ad10d
begin
	#for peak 2
	empty!(tot_TF)
	empty!(tot_a)
	empty!(tot_dadt)
	empty!(lens)
	DF2=copy(HR_A)
	DF2.dadT=[FS(HR_A.T[hrt], FSPar_A[hrt][2,:]) for hrt in 1:hrate_A]
	m2, n2, p2, res2, A2, E2, lnres2 = KineticModel(DF2)
end

# ╔═╡ fc6505f0-6aeb-11eb-3482-2ba2e564c4de
begin
	
	plot(xlabel="1/T", ylabel="LHS")
	#scatter!(1 ./tot_TF, log.(tot_dadt) - log.(SB(tot_a, m2, n2, p2)))
	yall2=log.(tot_dadt) - log.(SB(tot_a, m2, n2, p2))
	
	for i=1:hrate_A
		i==1 ? start=1 : start=lens[i-1]+1
		x=1 ./tot_TF[start:lens[i]]
		y=yall2[start:lens[i]]
		scatter!(x, y, label=string(HR_A.rate[i])) 
	end
	plot!(1 ./tot_TF, (log(A2) .+(E2/8.314*1000 ./tot_TF)),lw=4)
end

# ╔═╡ 2c5c7157-a544-4274-9e8a-a7e133903f09
md"# Block B - one peak"

# ╔═╡ c8bd7c14-9279-4d84-b343-b812a01d8017
begin
	#combine files to a DataFrame
	hrate_B=3
	HR_0B= DataFrame(rate=[6.0, 12., 15.0],
		T=[HR06.TS, HR12.TS, HR15.TS], #select TS or TF
		TG=[HR06.TG, HR12.TG, HR15.TG],
		time=[HR06.time, HR12.time, HR15.time])
	HR_0B.mTGA=[17.3, 18.6, 15.4] # mass of sample as per TGA
end

# ╔═╡ 64b42735-0bda-4473-a35d-a74e6a49dafe
begin
	HR_B = smooth(hrate_B, HR_0B)
	colors_B = colorschemes[:tol_muted][1:hrate_B]
	plot(HR_B.T, HR_B.a, labels=string.(transpose(HR_B.rate)),lw=2, legend = :outertopright, xlabel = "Temperature [K]", ylabel = "alpha reaction extent [-]", color_palette= colors_B)
end


# ╔═╡ 49dc7661-8202-4620-8250-9495f4bef954
begin
	pnB=1 #number of peaks
	wi_peak_B=[1.]
	SelPeaks_B=[
		[734],
		[747],
		[751]
		]

	plot_dadT(HR_B, SelPeaks_B, colors_B)
end

# ╔═╡ 60a31d5e-3dbf-4475-bc0b-b49a2748818a
md"## Deconvolution"

# ╔═╡ 9e0eed6c-90f1-443f-86d4-7d8921a7a655
covar_B, FSPar_B = fitFS2(SelPeaks_B, HR_B, pnB)

# ╔═╡ 1e3fff9f-a2ad-4e57-b020-c0b515639a82
md"Heating Rate: $(@bind hrt_B Slider(1:hrate_B, show_value=true))"

# ╔═╡ a3bfa0c4-e4fe-46a2-984c-52f3d61be98b
begin
	fullParB=reshape(transpose(FSPar_B[hrt_B]), pnB*4)
	plot_deconv(HR_B, fullParB, FSPar_B, pnB, hrt_B)
end

# ╔═╡ e2fa7da1-47fb-4b7a-bbea-17ed12902e8a
md"## Peak 1"

# ╔═╡ 22b403b6-dfbf-4069-850d-52758dc3239f
begin
	empty!(tot_TF)
	empty!(tot_a)
	empty!(tot_dadt)
	empty!(lens)
	DF_B=deepcopy(HR_B)
	DF_B.dadT=[FS(HR_B.T[hrt], FSPar_B[hrt_B][1,:]) for hrt in 1:hrate_B]
	Bm1, Bn1, Bp1, Bres1, BA1, BE1, Blnres1 = KineticModel(DF_B)
end

# ╔═╡ ac2d5ee7-604e-4e43-9a4f-0224c136c296
md" SB model parameters $Bm1, $Bn1, $Bp1
with a residual $Bres1"

# ╔═╡ 61f9c4af-f22d-4546-bfef-3f00eb812b5b
md" pre exponential factor $BA1 and activation energy $BE1 with a residual $Blnres1"

# ╔═╡ b05c11fd-4a0f-4679-9efb-2fd287d5dbff
begin
	plot(xlabel="1/T", ylabel="LHS")
	#scatter!(1 ./tot_TF, log.(tot_dadt) - log.(SB(tot_a, m2, n2, p2)))
	yallB=log.(tot_dadt) - log.(SB(tot_a, Bm1, Bn1, Bp1))
	
	for i=1:hrate_B
		i==1 ? start=1 : start=lens[i-1]+1
		x=1 ./tot_TF[start:lens[i]]
		y=yallB[start:lens[i]]
		scatter!(x, y, label=string(HR_B.rate[i])) 
	end
	plot!(1 ./tot_TF, (log(BA1) .+(BE1/8.314*1000 ./tot_TF)),lw=4, label="overall")
end

# ╔═╡ ab378cf0-6aed-11eb-14d3-e520239acfbf
md"# Verify with the differential model"

# ╔═╡ Cell order:
# ╠═e5ea9ee0-6954-11eb-2942-e76574442366
# ╠═51322bfd-c9b6-4a18-ac4d-7d9add4a004e
# ╟─0fb798e8-de08-424a-88bf-1f0727faa828
# ╠═b733f4fe-0b25-4248-8f6f-d872f3fd868f
# ╟─2d594622-4bd9-4a34-b938-7830aff31dd3
# ╟─fe678aa2-6954-11eb-08d3-bdd150d83c49
# ╠═2bd08850-6957-11eb-2578-ad33186981d0
# ╟─39570e90-6957-11eb-28d4-af76c02c6de1
# ╠═43805fa0-69d1-11eb-1ec1-95a700d9e47f
# ╟─09d02877-479d-480b-b503-19de160f8ac0
# ╠═05bce23a-c877-4dca-8570-b348a51eebfc
# ╟─cb3d56c7-b08d-4de5-98eb-fc01bdcae591
# ╠═47e9064b-0554-4928-ba89-372e0c1f4d02
# ╟─48293b1c-ce1f-4010-9d85-0bef313df73a
# ╟─0537a765-4048-4280-a40b-3ee64004f862
# ╠═9e71b200-6955-11eb-1cc7-7dd941326057
# ╠═a6e6c3d0-6955-11eb-0231-9fb731da375d
# ╠═c3976722-a898-4261-bfdb-534042b4f4d4
# ╠═48a38938-156d-44a1-9d03-0ac88bf1560b
# ╠═cf2287b0-6957-11eb-1d0e-4befdd8df2b3
# ╠═d6217440-6957-11eb-0573-6f4ddc7e5042
# ╠═e8c61920-6957-11eb-37ab-4f54568b6b41
# ╠═dc7b9552-6957-11eb-2fc7-55f9eab720f5
# ╠═09472c70-6958-11eb-0f78-bbdffc69d8e0
# ╟─1065d600-6958-11eb-2d59-8d98b52e0788
# ╠═14098f8c-25e3-4263-b5c3-e341e45b9418
# ╠═b816f330-6aec-11eb-3402-01980c895d22
# ╠═0db7cf46-66af-40b7-a2dd-6989acd9d4f5
# ╠═fae4e250-6aae-11eb-0763-1b5a759b100a
# ╠═c65d41a0-07e9-4f0c-8a79-ddf0498e237f
# ╠═7f734920-6aec-11eb-09f2-49a95f3ad10d
# ╠═fc6505f0-6aeb-11eb-3482-2ba2e564c4de
# ╠═2c5c7157-a544-4274-9e8a-a7e133903f09
# ╠═c8bd7c14-9279-4d84-b343-b812a01d8017
# ╠═64b42735-0bda-4473-a35d-a74e6a49dafe
# ╠═49dc7661-8202-4620-8250-9495f4bef954
# ╠═60a31d5e-3dbf-4475-bc0b-b49a2748818a
# ╠═9e0eed6c-90f1-443f-86d4-7d8921a7a655
# ╠═1e3fff9f-a2ad-4e57-b020-c0b515639a82
# ╠═a3bfa0c4-e4fe-46a2-984c-52f3d61be98b
# ╠═e2fa7da1-47fb-4b7a-bbea-17ed12902e8a
# ╠═22b403b6-dfbf-4069-850d-52758dc3239f
# ╠═ac2d5ee7-604e-4e43-9a4f-0224c136c296
# ╠═61f9c4af-f22d-4546-bfef-3f00eb812b5b
# ╠═b05c11fd-4a0f-4679-9efb-2fd287d5dbff
# ╠═ab378cf0-6aed-11eb-14d3-e520239acfbf
