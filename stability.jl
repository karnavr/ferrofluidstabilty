### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ ce56ba82-79e8-11ee-1b30-272520e58195
begin
	using Plots
	using LinearAlgebra
	using SpecialFunctions
	using LaTeXStrings
	using Trapz
	using DelimitedFiles
	using PlutoUI
	using ProgressLogging
	TableOfContents()
end

# ╔═╡ 6009f94c-6eb6-400f-8f73-f0433e82e42d
md"## Stability of a Ferrofluid Jet"

# ╔═╡ 8eeaba0e-8d13-4dea-9534-b45c53e9847f
md"First define constants common to all problems, and a few helper functions."

# ╔═╡ 1e018b6a-6e6f-470d-8e41-deaa52570d7f
function β(n, k, b, S0)

	beta1 = besseli.(1, k*b) .* besselk.(n, k.*S0)
	beta2 = (-1)^n .* besselk.(1, k*b) .* besseli.(n, k.*S0)

	return beta1 + beta2
end

# ╔═╡ 48249892-d7de-4048-88bb-dcd93e81da62
function c0(k, b, B)
	# wave speed for small amplitude waves, depending on the wave-number k
	
	c0  = sqrt.((1 ./ k).*((-β(1,k,b,1) ./ β(0,k,b,1)) .* (k.^2 .- 1 .+ B)))

	return c0
	
end

# ╔═╡ d60a4f84-af0b-48ff-95ff-5ddaef40c034
begin

	Bond = 1.5
	b = 0.1
	ϵ = 1 - Bond/2

	c1 = c0(1, b, Bond) # wave speed at k = 1
	
end

# ╔═╡ 0c13895c-6bd9-4377-929a-285d2b59843c
function fourierToProfile(coeffs, domain)
	
	N = length(coeffs)	# number of fourier coeffs (a0, a1, ...)
	
	profile = 0 		# initialize profile
	
	for i = 0:(N-1)
		profile = profile .+ coeffs[i+1] .* cos.(i .* domain) 
	end

	return profile
end

# ╔═╡ 042be35a-81e4-45ca-b1cf-2023be8092bb
md"### Equilibrium Jet

First consider the stability of a the ferrofluid jet at equilibrium $S(z) = 1$."

# ╔═╡ cc94cd11-da4f-4e31-800c-3053d7bfb2fd
md"##### Analytically

The eigenvalues of the spectral problem for the equilibrium case are given by:

$$\lambda_{\mu + m} = ic(\mu + m) \pm ic_0 (\mu + m)$$ 

where "

# ╔═╡ e3a9db1c-3387-4a39-9c97-c5a44c7f576a
md"m = $(@bind m PlutoUI.Slider(1:5, show_value = true, default=1))"

# ╔═╡ 5fe20173-ce5d-4cbd-860d-f36833c1fdeb
begin
	μ = collect(range(0.001,1.0,500))

	λ1 = im .* c1 .* (μ .+ m) .+ im .* c0(μ .+ m, b, Bond) .* (μ .+ m)
	λ2 = im .* c1 .* (μ .+ m) .- im .* c0(μ .+ m, b, Bond) .* (μ .+ m)
	
end

# ╔═╡ a0401fc6-9ddd-4bdc-98b9-a1543bc010fa
begin
	scatter(real(λ1),imag(λ1), label=L"\lambda_+")
	scatter!(real(λ2),imag(λ2), label=L"\lambda_-")

	xlabel!(L"Re\{\lambda\}")
	ylabel!(L"Im\{\lambda\}")

	xlims!(-0.5,0.5)
	ylims!(-5,7)

	# scatter(λ1)
end

# ╔═╡ 7a13978f-9f3d-4206-8262-3a8929dbe957
md"##### Numerically

We can also solve the problem numerically and see if we get the same results"

# ╔═╡ 9196cc07-f410-4843-8a4c-5c716f36fa4b
md"### Generalized Eigenvalue Problem"

# ╔═╡ 1b148894-7680-4f9e-b489-eec0d89db4a5
md"Import numerical solution data and convert from fourier coeffs to profiles:"

# ╔═╡ bb883c2d-8049-414f-adfe-919bacc1b6a9
begin
	# import results as array
	msolutions = readdlm("test & misc/matlab_solutions.csv", ',', Float64)

	mcoeffs = msolutions[:,2:end]
	mspeeds = msolutions[:,1]

	nsols = length(msolutions[:,1])

	# convert profiles + extract speeds
	mdomain = collect(range(-pi,pi,100))
	mprofiles = zeros(nsols,length(mdomain))

	for i = 1:nsols
		mprofiles[i,:] = fourierToProfile(mcoeffs[i,:], mdomain)
	end

	# reflect profiles 
	mprofiles = [mprofiles[:,Int(end/2):end] mprofiles[:,1:Int(end/2)-1]]; nothing

end

# ╔═╡ 25c765a0-dda0-468d-ae7b-a8d7019fce7c
function fourierSeries(coefficients::Vector{Float64}, domain, L::Number)
	
    N = length(coefficients) - 1
	
    S = zeros(length(domain))  		# profile S
    Sz = zeros(length(domain))  	# first derivative Sz
    Szz = zeros(length(domain))  	# second derivative Szz

    for (i, x) in enumerate(domain)
		
        # Calculate the series and its derivatives at each point x
        S[i] = coefficients[1] 
		
        for n in 1:N

			k = n*π/L

            S[i] += coefficients[n + 1] * cos(k * x)
            Sz[i] -= k * coefficients[n + 1] * sin(k * x)
            Szz[i] -= k^2 * coefficients[n + 1] * cos(k * x)
        end
		
    end

    return S, Sz, Szz
end

# ╔═╡ ef4bcf38-abc2-416f-a377-547930c69b62
@bind mprofileindex PlutoUI.Slider(1:nsols, default=1)

# ╔═╡ 6d9b1e2d-032f-4f78-ad59-84c408de490b
begin
	# plot profiles 
	profile_plot = plot(mdomain, mprofiles[mprofileindex,:], legend=false, title = "a1 = $(round(mcoeffs[mprofileindex,2], digits=3))", lw=2)
	ylims!(0.45,1.3)
	xlabel!(L"z"); ylabel!(L"S")

	# plot coeffs 
	first_coeff = 0
	coeff_plot = scatter(abs.(mcoeffs[mprofileindex,first_coeff+1:end]), legend=false, title="k = $(mprofileindex)", xticks = :all, yaxis=:log)
	xlabel!("a$(first_coeff) to a$(length(mcoeffs[1,:])-1)")

	# plot branch
	branch_plot = scatter(mspeeds[1:mprofileindex], mcoeffs[1:mprofileindex,2], legend = false, markersize=4)
	xlabel!(L"c"); ylabel!(L"a_1")
	xlims!(0.73,0.82); ylims!(0.04,0.34)
	
	plot(profile_plot, branch_plot, size=(700,350))
end

# ╔═╡ 30364597-1f57-4e30-ade6-31ad6f11949d
md"Create ideal data for now:"

# ╔═╡ 458b05b7-5655-415d-a9bf-f64d157ac891
begin
	z = collect(range(-π,+π,100))
	
	# newprofiles = zeros(nsols,length(mdomain))
	S0, S0z, S0zz = fourierSeries(mcoeffs[40,:], z, π)

	# S0 = mprofiles[66,:]
	# S0z = gradient(S0, z)
	# S0zz = derivative(S0z, z)

	# S0 = - 0.001 .* cos.(z) .+ 1
	# S0 = - 0.00 .* cos.(z) .+ 1
	# profile1 = mprofiles[1,:]
	# S0z = zeros(length(z))
	# S0zz = zeros(length(z))

	S0sq = 1 .+ S0z.^2

	scatter(z,S0)
	# scatter!(z,profile1)
	# scatter!(z,S0z)

	c = mspeeds[40]
	N = 4
	# q0z = c0 + c0

	κ = - (S0zz./(S0sq.^(3/2))) .+ (1 ./ (S0.*S0sq.^(1/2)))
	q0z = c .+ (1 ./ S0) .* sqrt.(S0sq .* ( (c^2 + 2 .* ϵ .- 2 .* κ).*(S0.^2) .+ Bond))
end

# ╔═╡ 70f948f9-e526-47c3-98c9-d745493a5e18
begin
	# c = c1
	# c = mspeeds[20]
	# N = 4
	# # q0z = c0 + c0

	# κ = - (S0zz./(S0sq.^(3/2))) .+ (1 ./ (S0.*S0sq.^(1/2)))
	# q0z = c .+ (1 ./ S0) .* sqrt.(S0sq .* ( (c^2 + 2 .* ϵ .- 2 .* κ).*(S0.^2) .+ Bond))
end

# ╔═╡ 50917e2b-48b0-41cb-95cd-a2d7b3a8bf7b
# ╠═╡ disabled = true
#=╠═╡
begin

	λ = zeros(ComplexF64, length(μ))
	# λ = zeros(ComplexF64, 0)
	
for i = 1:length(μ)

	# create matrices
	A = Ag(N, z, S0z, q0z, c)
	B = Bg(N, z)
	E = Eg(N, z, S0, S0z, S0zz, q0z, c, Bond, μ[i])
	F = Fg(N, z, S0, S0z, q0z, c, μ[i])

	C = Cg(N, z, S0, b, μ[i])
	D = zeros(2N+1,2N+1)
	G = Gg(N, z, S0, b, c, μ[i])
	H = Hg(N, z, S0, b, μ[i])

	lhs = [A B; C D]
    rhs = [E F; G H]
	
	# solve problem 
	solutions = eigen(lhs, rhs)
	
	# save solution
	λ[i] = (solutions.values)[5]

	# append!(λ, solutions.values)
end
end
  ╠═╡ =#

# ╔═╡ 2ed210ad-a828-40ed-bbc0-dffb2e808b17
λ = zeros(ComplexF64, 4*N+2, length(μ))

# ╔═╡ 43165603-47c4-4107-8c7a-6bd4f056939c
begin
	scatter(vec(λ))
	# xlims!(0.3,0.7)
	ylims!(0.325,0.35)
end

# ╔═╡ 7d535fe6-f885-4ab6-8660-392887cb8e4b
begin
	maxrealλ = zeros(length(μ))
	for i = 1:length(μ)
		maxrealλ[i] = maximum(real(λ[:,i]))
	end
	scatter(μ,maxrealλ, label = "Re{λ}", markersize = 3)
	# scatter!(μ,imag(λ), label = "Im{λ}")
	# xlims!(0.5,0.6)

	xlabel!(L"\mu")
	# ylims!(-0.1,0.1)
	# xlims!(-0.1,0.1)
end

# ╔═╡ f9244069-3f08-4446-8857-4d824ae13834
md"### Full solver function"

# ╔═╡ 615b3f79-3a72-49ba-bb1d-90a52311c3cd
function stabilityPlots(λ, Nmu)

	μ = collect(range(0.001,1.0,Nmu))

	# plot λ on complex plane
		complexPlot = scatter(vec(λ), markersize = 1, legend = false)
		# xlims!(0.3,0.7)
		# ylims!(0.325,0.35)

		# plot max real λ vs μ
		maxrealλ = zeros(Nmu)
		for i = 1:Nmu
			maxrealλ[i] = maximum(real(λ[:,i]))
		end

		muPlot = scatter(μ,maxrealλ, label = "Re{λ}", markersize = 1)
		xlabel!(L"\mu")

		# combine into one plot
		plot(complexPlot, muPlot, size=(700,350))
		title!("Nmu = $(Nmu)")
		

	return plot(complexPlot, muPlot, size=(700,350))
end

# ╔═╡ d2eeccbe-1f59-4e5e-9e64-342a38b8f477
md"##### Matrix definitions"

# ╔═╡ 191f060d-4302-4b18-82f8-6e20224ee201
md"###### Local 

In this case, $j = - m$ so we only iterate over the $m$ index and the matrices are left-right flipped diagonal."

# ╔═╡ 9bc57aaf-5cba-4569-975c-a504730b8008
function Ag(N, z, S0z, q0z, c)

	# initialize matrix
	A = zeros(2*N+1, 2*N+1) .+ 0.0im
	
	# constants
	α = 1 ./ (1 .+ S0z.^2)
	f = S0z .* (q0z .- c) .* α

	# populate matrix
	for mm = 1:(2*N + 1)

		# working with both code and matrix indicies
		m = mm
		j = m
		jj = j

		# define integrand specifc to matrix
		term = f
		
		# populate matrix entries
		A[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* j .* z))

	end

	return A 
	
end

# ╔═╡ 6f8174f6-ee58-4c76-81fd-d6c522491339
function Atest(N, z, S0z, q0z, c)

	# initialize matrix
	A = zeros(2*N+1, 2*N+1) .+ 0.0im
	
	# constants
	α = 1 ./ (1 .+ S0z.^2)
	f = S0z .* (q0z .- c) .* α

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

		for jj = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = jj - 1
	
			# define integrand specifc to matrix
			term = f
			
			# populate matrix entries
			A[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* j .* z))

		end

	end

	return A 
	
end

# ╔═╡ 4ee55b92-512b-431e-9973-d0c283aa13d2
Ag(N, z, S0z, q0z, c1)

# ╔═╡ fe9a9fca-f098-441e-b51f-8c27be65785d
Atest(N, z, S0z, q0z, c1)

# ╔═╡ 8d7afa78-57ae-4fb9-8899-4254f22a11f8
function Bg(N, z)

	# initialize matrix
	B = zeros(2*N+1, 2*N+1) .+ 0.0im

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

		# working with both code and matrix indicies
		m = mm - 1
		
		# j = -m
		# jj = -j + 2*N - 2*(mm-1)

		j = m
		jj = j + 1

		# define integrand specifc to matrix
		term = -1
		
		# populate matrix entries
		B[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* j .* z))

	end
	

	return B
end

# ╔═╡ 2f4b4efd-a892-4797-a9f1-e5f94222af33
function Btest(N, z)

	# initialize matrix
	B = zeros(2*N+1, 2*N+1) .+ 0.0im

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

		for jj = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = jj - 1
	
			# define integrand specifc to matrix
			term = -1
			
			# populate matrix entries
			B[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m-j) .* z))

		end

	end
	
	return B
end

# ╔═╡ 1eabaff3-9776-4280-a4dc-5686441544f4
Bg(N, z)

# ╔═╡ 17fa2078-d17e-43b6-89e9-cc4396582676
real(Btest(N,z))

# ╔═╡ c2ec1e3a-7ed5-43de-867a-bfa2f25de5af
function Eg(N, z, S0, S0z, S0zz, q0z, c, B, μ)

	# initialize matrix
	E = zeros(2*N+1, 2*N+1) .+ 0.0im
	
	# constants
	α = 1 ./ (1 .+ S0z.^2)
	f = S0z .* (q0z .- c) .* α
	γ = c .- q0z .+ f.*S0z

	# loop over index to populate matrix
	for mm = 1:(2*N + 1)

		# working with both code and matrix indicies
		m = mm
		
		j = -m
		jj = -j + 2*N - 2*(mm-1)

		# define integrand specifc to matrix
		term = ( γ.*f .+ 3 .*S0zz.*S0z.*(α.^(5/2)) .- ((α.^(3/2)) ./ S0) .* S0z .- (α.^(3/2)) .* im .* (μ .+ m) ) .* (im.*(μ .+ m)) .- (α.^(1/2)) ./ S0.^2 .+ B ./ S0.^3
		
		# populate matrix entries
		E[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* j .* z))

	end

	return E
end

# ╔═╡ a61c0de8-5ec1-4f58-bad5-488931f53b60
function Etest(N, z, S0, S0z, S0zz, q0z, c, B, μ)

	# initialize matrix
	E = zeros(2*N+1, 2*N+1) .+ 0.0im
	
	# constants
	α = 1 ./ (1 .+ S0z.^2)
	f = S0z .* (q0z .- c) .* α
	γ = c .- q0z .+ f.*S0z

	for mm = 1:(2*N+1)

		for jj = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = jj - 1
	
			# define integrand specifc to matrix
			term = ( γ.*f .+ 3 .*S0zz.*S0z.*(α.^(5/2)) .- ((α.^(3/2)) ./ S0) .* S0z .- (α.^(3/2)) .* im .* (μ .+ m) ) .* (im.*(μ .+ m)) .- (α.^(1/2)) ./ S0.^2 .+ B ./ S0.^3
			
			# populate matrix entries
			E[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m-j) .* z))

		end

	end

	return E
end

# ╔═╡ ce187832-a000-4d42-bdfd-da238e644f59
Eg(N, z, S0, S0z, S0zz, q0z, c1, Bond, 0.1)

# ╔═╡ 187edc0f-e9a5-406d-8a33-a326a6ccfa99
real(Etest(N, z, S0, S0z, S0zz, q0z, c1, Bond, 0.1))

# ╔═╡ cf7a32d7-59ca-46dc-9206-dfcd0d030e37
(0.1 + 1).^2 - 1 + 1.5

# ╔═╡ 0b757a49-b1a3-4d6d-b482-fe7adce2c499
function Fg(N, z, S0, S0z, q0z, c, μ)

	# initialize matrix
	F = zeros(2*N+1, 2*N+1) .+ 0.0im
	
	# constants
	α = 1 ./ (1 .+ S0z.^2)
	f = S0z .* (q0z .- c) .* α
	γ = c .- q0z .+ f.*S0z

	# loop over index to populate matrix
	for mm = 1:(2*N + 1)

		# working with both code and matrix indicies
		m = mm
		
		j = -m
		jj = -j + 2*N - 2*(mm-1)

		# define integrand specifc to matrix
		term = - γ .* im .* (μ .+ m)
		
		# populate matrix entries
		F[jj,mm] =  1/(2*π) .* round(trapz(z,term .* exp.(- im .* j .* z)))

	end
	

	return F
end

# ╔═╡ 9338bc36-9f62-45a9-8805-bd7964c3e090
function Ftest(N, z, S0, S0z, q0z, c, μ)

	# initialize matrix
	F = zeros(2*N+1, 2*N+1) .+ 0.0im
	
	# constants
	α = 1 ./ (1 .+ S0z.^2)
	f = S0z .* (q0z .- c) .* α
	γ = c .- q0z .+ f.*S0z

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

		for jj = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = jj - 1
	
			# define integrand specifc to matrix
			term = - γ .* im .* (μ .+ m)
			# term = 5.0 * im
			
			# populate matrix entries
			F[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m-j) .* z))

		end

	end
	
	return F
end

# ╔═╡ 5a36d3de-8341-4500-bff2-b32a88939377
term = - (-c1) .* im .* (0 .+ 1)

# ╔═╡ 84d5dba8-ac87-49eb-b7ca-0b186c8ba203
trapz(z,term .* exp.(- im .* (0-1) .* z))

# ╔═╡ 73828383-6b1f-4f8b-ab22-b2c5e1f581a0
Fg(N, z, S0, S0z, q0z, c1, μ)

# ╔═╡ 3cf16498-0770-4308-9962-eac685a43599
imag(Ftest(N, z, S0, S0z, q0z, c1, 0.1))

# ╔═╡ 86128a17-7089-4f55-9be3-c7da169b6ff6
imag(- im .* (-c1) .* (0.1 + 1))

# ╔═╡ b676b7ee-9d47-41dd-a80d-60fa8556a38e
md"###### Non-local

In this case, $j$ is a ''free'' index, so we iterate over both $j$ and $m$."

# ╔═╡ 01b0c32a-a664-4d1c-9519-17d05a8dddb8
md"Always seems to be off by a factor of m + μ..."

# ╔═╡ 9c6745e6-0757-43f4-89da-ae6b04b1803c
function Cg(N, z, S0, b, μ)

	# initialize matrix
	C = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N+1)

		for jj = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = jj - 1

			k = m .+ μ

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = - k .* S0 .* β0
			
			# populate matrix entries
			C[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m - j).* z))
		
		end
	end

	return C
	
end

# ╔═╡ ef4b6992-1053-4b64-a659-0f45407ce8ef
function Ctest(N, z, S0, b, μ)

	# initialize matrix
	C = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N+1)

		for jj = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = jj - 1

			k = m .+ μ

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = - k .* S0 .* β0
			
			# populate matrix entries
			C[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m - j).* z)) ./ k
		
		end
	end

	return C
	
end

# ╔═╡ 36d65f3b-ef45-4960-b4ae-5e30b94e7500
function Ctest2(N, z, S0, b, μ)

	# initialize matrix
	C = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N+1)

		for jj = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = jj - 1

			k = j .+ μ

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = - (j .+ μ) .* S0 .* β0
			
			# populate matrix entries
			C[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m - j).* z))
		
		end
	end

	return C
	
end

# ╔═╡ 39d542a2-0ecb-4ed0-b0cf-4eca6bd67094
real(Cg(N, z, S0, b, 0.1))

# ╔═╡ 20a266d9-7eaa-448c-bef0-209dc5771959
-β(0, (3 + 0.1), b, S0) .* (3 + 0.1)

# ╔═╡ a45b4050-302b-4842-a13e-3934041d75c4
real(Ctest(N, z, S0, b, 0.1))

# ╔═╡ dbea8be9-f915-4221-a71b-c98e4ca40a92
-β(0, (1 + 0.1), b, S0)

# ╔═╡ c75a4af4-ec70-41ae-848b-0d4464047936
function Gg(N, z, S0, b, c, μ)

	# initialize matrix
	G = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N + 1)

		# converting from code to matrix index
		m = mm 

		for jj = 1:(2*N + 1)

			j = jj 

			k = μ .+ m

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = S0.*S0z.*c .* k.^2 .* β1 .- S0z.*c.*k.*β0 .+ im.*S0.*q0z.* k.^2 .* β0 .- im.*S0.*c.* k.^2 .* β0
			
			# populate matrix entries
			G[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m - j).* z))
		
		end
	end

	return G
	
end

# ╔═╡ a44953e2-352b-4524-9191-ac0f0d3d9a43
function Gtest(N, z, S0, b, c, μ)

	# initialize matrix
	G = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N + 1)

		# converting from code to matrix index
		m = mm 

		for jj = 1:(2*N + 1)

			j = jj 

			k = μ .+ m

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = S0.*S0z.*c .* k.^2 .* β1 .- S0z.*c.*k.*β0 .+ im.*S0.*q0z.* k.^2 .* β0 .- im.*S0.*c.* k.^2 .* β0
			
			# populate matrix entries
			G[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m - j).* z)) ./ k
		
		end
	end

	return G
	
end

# ╔═╡ bdecea25-78ce-43d1-9ad9-17785889f722
function Gtest2(N, z, S0, b, c, μ)

	# initialize matrix
	G = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N + 1)

		# converting from code to matrix index
		m = mm 

		for jj = 1:(2*N + 1)

			j = jj 

			k = μ .+ j

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = S0.*S0z.*c .* k.^2 .* β1 .- S0z.*c.*k.*β0 .+ im.*S0.*q0z.* k.^2 .* β0 .- im.*S0.*c.* k .* (μ .+ m) .* β0
			
			# populate matrix entries
			G[jj,mm] =  1/(2*π) .* trapz(z,term .* exp.(- im .* (m - j).* z))
		
		end
	end

	return G
	
end

# ╔═╡ 8c00452c-2585-4edc-9ed6-b5befaa7991d
imag(Gg(N, z, S0, b, c1, 0.1))

# ╔═╡ 8a542361-959e-4c6d-b15d-2c5700b2ef3c
imag((im .* (0.1 + 2) .* β(0, (0.1 + 2), b, S0)) .* (c1)) .* (0.1 + 2) 

# ╔═╡ 04077a45-a8a3-4d4b-a558-174c059d413d
imag(Gtest(N, z, S0, b, c1, 0.1))

# ╔═╡ 9bb89559-2490-4f7a-8256-e02239ccff29
imag((im .* (0.1 + 3) .* β(0, (0.1 + 3), b, S0)) .* (c1))

# ╔═╡ b586cc56-cef0-4a0a-b31d-b7a9b37ecffa
function Hg(N, z, S0, b, μ)

	# initialize matrix
	H = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N + 1)

		# converting from code to matrix index
		m = mm

		for jj = 1:(2*N + 1)

			j = jj 

			k = μ .+ m

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = S0 .* k.^2 .* β1
			
			# populate matrix entries
			H[jj,mm] =  1/(2*π) .* round(trapz(z,term .* exp.(- im .* (m - j).* z)))
		
		end
	end

	return H
	
end

# ╔═╡ 727e16e0-989a-4168-b705-7512e52fbf83
function Htest(N, z, S0, b, μ)

	# initialize matrix
	H = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N + 1)

		# converting from code to matrix index
		m = mm

		for jj = 1:(2*N + 1)

			j = jj 

			k = μ .+ m

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = S0 .* k.^2 .* β1
			
			# populate matrix entries
			H[jj,mm] =  1/(2*π) .* round(trapz(z,term .* exp.(- im .* (m - j).* z))) ./ k
		
		end
	end

	return H
	
end

# ╔═╡ e7ef8aa5-8668-4489-97e8-fb46c80a501c
function Htest2(N, z, S0, b, μ)

	# initialize matrix
	H = zeros(2*N+1, 2*N+1) .+ 0.0im

	for mm = 1:(2*N + 1)

		# converting from code to matrix index
		m = mm

		for jj = 1:(2*N + 1)

			j = jj 

			k = μ .+ j

			β0 = β(0, k, b, S0)
			β1 = β(1, k, b, S0)

			# define integrand specifc to matrix
			term = S0 .* k .* (μ .+ m) .* β1
			
			# populate matrix entries
			H[jj,mm] =  1/(2*π) .* round(trapz(z,term .* exp.(- im .* (m - j).* z)))
		
		end
	end

	return H
	
end

# ╔═╡ ad51ce37-7a2d-4741-b69b-f85683ec4b28
begin

	# λ = zeros(ComplexF64, length(μ))
	# λ = zeros(ComplexF64, 0)

	# λ = zeros(4*N+2, length(μ))
	# λ = zeros(ComplexF64, 4*N+2, length(μ))
	
@progress for i = 1:length(μ)

	# create matrices
	A = Atest(N, z, S0z, q0z, c)
	B = Btest(N, z)
	E = Etest(N, z, S0, S0z, S0zz, q0z, c, Bond, μ[i])
	F = Ftest(N, z, S0, S0z, q0z, c, μ[i])

	C = Ctest2(N, z, S0, b, μ[i])
	D = zeros(2N+1,2N+1)
	G = Gtest2(N, z, S0, b, c, μ[i])
	H = Htest2(N, z, S0, b, μ[i])

	lhs = [A B; C D]
    rhs = [E F; G H]
	
	# solve problem 
	solutions = eigen(rhs, lhs)
	
	# save solution
	λ[:,i] = (solutions.values)

	# append!(λ, solutions.values)
end
end

# ╔═╡ 5f6b4556-54f2-42e8-841f-684c00a19ee9
function solveGenEig(solution, Nmodes, Nmu, plotting = true)

	# un-pack solution 
	coeffs = solution[2:end]
	c = solution[1]
	N = Nmodes

	# create domain and convert to real space
	z = collect(range(-π,+π,100))
	S0, S0z, S0zz = fourierSeries(coeffs, z, π)

	# commonly used constants
	S0sq = 1 .+ S0z.^2
	κ = - (S0zz./(S0sq.^(3/2))) .+ (1 ./ (S0.*S0sq.^(1/2)))
	q0z = c .+ (1 ./ S0) .* sqrt.(S0sq .* ( (c^2 + 2 .* ϵ .- 2 .* κ).*(S0.^2) .+ Bond))

	# set up stability stuff
	μ = collect(range(0.001,1.0,Nmu))
	λ = zeros(ComplexF64, 4*N+2, length(μ))

	# create matrices that stay constant (don't depend on μ)
	A = Atest(N, z, S0z, q0z, c)
	B = Btest(N, z)
	D = zeros(2N+1,2N+1)
	
	
	Threads.@threads for i = 1:length(μ)

		# create matrices
		
		E = Etest(N, z, S0, S0z, S0zz, q0z, c, Bond, μ[i])
		F = Ftest(N, z, S0, S0z, q0z, c, μ[i])
	
		C = Ctest2(N, z, S0, b, μ[i])
		G = Gtest2(N, z, S0, b, c, μ[i])
		H = Htest2(N, z, S0, b, μ[i])
	
		lhs = [A B; C D]
	    rhs = [E F; G H]
		
		# solve problem 
		solutions = eigen(rhs, lhs)
		
		# save solution
		λ[:,i] = (solutions.values)
	
	end

	if plotting == true

		# plot λ on complex plane
		complexPlot = scatter(vec(λ), markersize = 1, legend = false)
		# xlims!(0.3,0.7)
		# ylims!(0.325,0.35)

		# plot max real λ vs μ
		maxrealλ = zeros(length(μ))
		for i = 1:length(μ)
			maxrealλ[i] = maximum(real(λ[:,i]))
		end

		muPlot = scatter(μ,maxrealλ, label = "Re{λ}", markersize = 1)
		xlabel!(L"\mu")

		# combine into one plot
		plot(complexPlot, muPlot, size=(700,350))
		
	end

	return λ
end

# ╔═╡ 001a6f78-7e07-46d4-9a58-d53f86fa0239
begin
	λ21000 = solveGenEig(msolutions[40,:], 2, 1000, true)
	stabilityPlots(λ21000, 1000)
end

# ╔═╡ f1c23bf3-15ba-4db7-be8e-49c5abcdf6a8
begin
	λ41000 = solveGenEig(msolutions[40,:], 4, 1000, true)
	stabilityPlots(λ41000, 1000)
end

# ╔═╡ cd603fe7-bcac-4609-91dc-f32ef2a64c83
begin
	λ45000 = solveGenEig(msolutions[40,:], 4, 5000, true)
	stabilityPlots(λ45000, 5000)
end

# ╔═╡ 69ad7b3f-45ad-4a3b-8381-f6ac0cd6daf5
λ410000 = solveGenEig(msolutions[40,:], 4, 10000, true)

# ╔═╡ ddb3bfb8-4c1f-4a6e-b72a-a96e99e75859
stabilityPlots(λ410000, 10000)

# ╔═╡ cebac023-d60e-464c-89bf-ed15f26c432f
λ420000 = solveGenEig(msolutions[40,:], 4, 20000, true)

# ╔═╡ 30b461ba-a69f-4de9-bf84-ecca8a6da742
begin 
	stabilityPlots(λ420000, 20000)
	ylims!(-1e-10, 1e-10)
	xlims!(0.04, 0.045)
end

# ╔═╡ 157fc67a-93cf-4bec-8a2d-8be2230dcb51
λ810000 = solveGenEig(msolutions[40,:], 8, 10000, true)

# ╔═╡ 134fd2c0-61ac-4ba3-9246-da571dfdb352
begin
	stabilityPlots(λ810000, 10000)
	ylims!(-1e-7, 1e-7)
end

# ╔═╡ d89e6f1b-697e-4ebb-90cb-2d8649b49a31
function solveGenEignothreads(solution, Nmodes, Nmu, plotting = true)

	# un-pack solution 
	coeffs = solution[2:end]
	c = solution[1]
	N = Nmodes

	# create domain and convert to real space
	z = collect(range(-π,+π,100))
	S0, S0z, S0zz = fourierSeries(coeffs, z, π)

	# commonly used constants
	S0sq = 1 .+ S0z.^2
	κ = - (S0zz./(S0sq.^(3/2))) .+ (1 ./ (S0.*S0sq.^(1/2)))
	q0z = c .+ (1 ./ S0) .* sqrt.(S0sq .* ( (c^2 + 2 .* ϵ .- 2 .* κ).*(S0.^2) .+ Bond))

	# set up stability stuff
	μ = collect(range(0.001,1.0,Nmu))
	λ = zeros(ComplexF64, 4*N+2, length(μ))

	# create matrices that stay constant (don't depend on μ)
	A = Atest(N, z, S0z, q0z, c)
	B = Btest(N, z)
	D = zeros(2N+1,2N+1)
	
	
	@progress for i = 1:length(μ)

		# create matrices
		
		E = Etest(N, z, S0, S0z, S0zz, q0z, c, Bond, μ[i])
		F = Ftest(N, z, S0, S0z, q0z, c, μ[i])
	
		C = Ctest2(N, z, S0, b, μ[i])
		G = Gtest2(N, z, S0, b, c, μ[i])
		H = Htest2(N, z, S0, b, μ[i])
	
		lhs = [A B; C D]
	    rhs = [E F; G H]
		
		# solve problem 
		solutions = eigen(rhs, lhs)
		
		# save solution
		λ[:,i] = (solutions.values)
	
	end

	if plotting == true

		# plot λ on complex plane
		complexPlot = scatter(vec(λ), markersize = 1, legend = false)
		# xlims!(0.3,0.7)
		# ylims!(0.325,0.35)

		# plot max real λ vs μ
		maxrealλ = zeros(length(μ))
		for i = 1:length(μ)
			maxrealλ[i] = maximum(real(λ[:,i]))
		end

		muPlot = scatter(μ,maxrealλ, label = "Re{λ}", markersize = 1)
		xlabel!(L"\mu")

		# combine into one plot
		plot(complexPlot, muPlot, size=(700,350))
		
		
	end

	return plot(complexPlot, muPlot, size=(700,350))
end

# ╔═╡ 7bc4c8dc-0134-47fd-bc6a-18fbe81aff82
real(Hg(N, z, S0, b, 0.1))

# ╔═╡ ac93e9d5-6097-4ecf-b35f-e60952f1284b
β(1, (0.1 + 2), b, S0) .* (0.1 + 2) .* (0.1 + 2)

# ╔═╡ 729abf56-2a75-48af-bb72-9168fcf9560a
real(Htest(N, z, S0, b, 0.1))

# ╔═╡ 20683dd5-d1f1-4969-865b-a12c201c3c76
β(1, (0.1 + 3), b, S0) .* π 

# ╔═╡ 6823e7bf-66f5-4a1f-b063-b19fadf071db
md"#### Equilibrium"

# ╔═╡ 3ff4a644-f633-4d30-b36d-8b56094bfb20
md"###### Local "

# ╔═╡ d52ecd4c-9755-491b-94f7-db96fae1eb05
function Be(N)

	# initialize matrix
	B = zeros(2*N+1, 2*N+1) .+ 0.0im

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = m
			
			# populate matrix entries
			B[mm,mm] =  -1

	end
	
	return B
end

# ╔═╡ f5b244e8-0a8e-41a1-8ead-b0f2922fb125
real(Be(N))

# ╔═╡ 87d3da6f-7537-48b1-9c9e-6ffd9e45f602
function Ee(N, B, μ)

	# initialize matrix
	E = zeros(2*N+1, 2*N+1) .+ 0.0im

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = m
			
			# populate matrix entries
			E[mm,mm] = (μ .+ m).^2 - 1 + B

	end
	
	return E
end

# ╔═╡ 3cf7fb88-db37-4d34-a13f-efde9000e8b7
real(Ee(N, 1.5, 0.1))

# ╔═╡ 570496eb-a486-41b4-a441-9c5efc095fab
function Fe(N, S0z, q0z, c, μ)

	# initialize matrix
	F = zeros(2*N+1, 2*N+1) .+ 0.0im

	# constants
	α = 1 ./ (1 .+ S0z.^2)
	f = S0z .* (q0z .- c) .* α
	γ = c .- q0z .+ f.*S0z

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = m
			
			# populate matrix entries
			F[mm,mm] = -im .* γ[1] .* (μ .+ m)

	end
	
	return F
end

# ╔═╡ 75f39b14-c053-471c-bea8-6657bd70656d
imag(Fe(N, S0z, q0z, c1, 0.1))

# ╔═╡ bc5cdcfe-cf3e-4e83-b423-712e3a4d13e6
md"###### Non-local "

# ╔═╡ 9904c284-ef99-4ee6-877c-9e06eb416f1f
function Ce(N, z, S0, b, μ)

	# initialize matrix
	C = zeros(2*N+1, 2*N+1) .+ 0.0im

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = m

			k = μ .+ m
			
			# populate matrix entries
			C[mm,mm] =  -β(0, k, b, S0)[1]

	end
	
	return C
	
end

# ╔═╡ c72c3408-a7c1-4427-8ebc-1572fe40d1c3
real(Ce(N, z, S0, b, 0.1))

# ╔═╡ 19004d6e-b4a5-4754-8564-8bde3e2c2532
function Ge(N, z, S0, b, c, μ)

	# initialize matrix
	G = zeros(2*N+1, 2*N+1) .+ 0.0im

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = m

			k = μ .+ m
			
			# populate matrix entries
			G[mm,mm] = ((im .* k .* β(0, k, b, S0)) .* (q0z .- c))[1]

	end

	return G
	
end

# ╔═╡ 41e552a4-0976-40d9-a1b1-e40f241dea4a
imag(Ge(N, z, S0, b, c1, 0.1))

# ╔═╡ 471c8398-1306-4514-ab4c-2e37c04cfe9b
function He(N, z, S0, b, μ)

	# initialize matrix
	H = zeros(2*N+1, 2*N+1) .+ 0.0im

	# loop over index to populate matrix
	for mm = 1:(2*N+1)

			# working with both code and matrix indicies
			m = mm - 1
			j = m

			k = μ .+ m
		
			β1 = β(1, k, b, S0)
			
			# populate matrix entries
			H[mm,mm] = (-β1 .* k)[1]

	end

	return H
	
end

# ╔═╡ eda75fa2-e051-44a5-9cb7-c041ffd322ef
begin

	λe = zeros(ComplexF64, length(μ))
	# λ = zeros(ComplexF64, 0)
	
@progress for i = 1:length(μ)

	# create matrices
	Aeq = zeros(2N+1,2N+1)
	Beq = Be(N)
	Eeq = Ee(N, Bond, μ[i])
	Feq = Fe(N, S0z, q0z, c, μ[i])

	Ceq = Ce(N, z, S0, b, μ[i])
	Deq = zeros(2N+1,2N+1)
	Geq = Ge(N, z, S0, b, c, μ[i])
	Heq = He(N, z, S0, b, μ[i])

	lhse = [Aeq Beq; Ceq Deq]
    rhse = [Eeq Feq; Geq Heq]
	
	# solve problem 
	solutionse = eigen(rhse, lhse)
	
	# save solution
	# λe[i] = (solutions.values)[5]

	append!(λe, solutionse.values)
end
end

# ╔═╡ 50ecf6a6-1c16-48c3-9597-d5967402b0c9
scatter(λe)

# ╔═╡ f0b14fd1-f089-4e5d-bf07-9e66bc97e2e9
begin
	scatter(μ,real(λe), label = "Re{λ}")
	scatter!(μ,imag(λe), label = "Im{λ}")

	xlabel!(L"\mu")
	# ylims!(-0.1,0.1)
	# xlims!(-0.1,0.1)
end

# ╔═╡ 329a421b-a0aa-488a-a4a3-1c3d8bb153e6
real(He(N, z, S0, b, 0.1))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Trapz = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"

[compat]
LaTeXStrings = "~1.3.1"
Plots = "~1.39.0"
PlutoUI = "~0.7.52"
ProgressLogging = "~0.1.4"
SpecialFunctions = "~2.3.1"
Trapz = "~2.0.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "eea4a68aea7773a7f6d9f8f333849d5dd9377431"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e0af648f0692ec1691b5d094b8724ba1346281cf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.18.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "2fba81a302a7be671aefe194f0525ef231104e7f"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.8"
weakdeps = ["InverseFunctions"]

    [deps.ChangesOfVariables.extensions]
    ChangesOfVariablesInverseFunctionsExt = "InverseFunctions"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "8cfa272e8bdedfa88b6aefbbca7c19f1befac519"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "9fb0b890adab1c0a4a475d4210d51f228bfc250d"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "7c29f0e8c575428bd84dc3c72ece5178caa67336"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.2+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "49cbf7c74fafaed4c529d47d48c8f7da6a19eb75"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.1"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Trapz]]
git-tree-sha1 = "79eb0ed763084a3e7de81fe1838379ac6a23b6a0"
uuid = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
version = "2.0.3"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "24b81b59bd35b3c42ab84fa589086e19be919916"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.11.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "47cf33e62e138b920039e8ff9f9841aafe1b733e"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.35.1+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "f7c281e9c61905521993a987d38b5ab1d4b53bef"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═ce56ba82-79e8-11ee-1b30-272520e58195
# ╟─6009f94c-6eb6-400f-8f73-f0433e82e42d
# ╟─8eeaba0e-8d13-4dea-9534-b45c53e9847f
# ╟─d60a4f84-af0b-48ff-95ff-5ddaef40c034
# ╟─48249892-d7de-4048-88bb-dcd93e81da62
# ╟─1e018b6a-6e6f-470d-8e41-deaa52570d7f
# ╟─0c13895c-6bd9-4377-929a-285d2b59843c
# ╟─042be35a-81e4-45ca-b1cf-2023be8092bb
# ╟─cc94cd11-da4f-4e31-800c-3053d7bfb2fd
# ╠═5fe20173-ce5d-4cbd-860d-f36833c1fdeb
# ╟─e3a9db1c-3387-4a39-9c97-c5a44c7f576a
# ╠═a0401fc6-9ddd-4bdc-98b9-a1543bc010fa
# ╟─7a13978f-9f3d-4206-8262-3a8929dbe957
# ╠═eda75fa2-e051-44a5-9cb7-c041ffd322ef
# ╠═50ecf6a6-1c16-48c3-9597-d5967402b0c9
# ╠═f0b14fd1-f089-4e5d-bf07-9e66bc97e2e9
# ╟─9196cc07-f410-4843-8a4c-5c716f36fa4b
# ╟─1b148894-7680-4f9e-b489-eec0d89db4a5
# ╠═bb883c2d-8049-414f-adfe-919bacc1b6a9
# ╟─25c765a0-dda0-468d-ae7b-a8d7019fce7c
# ╟─ef4bcf38-abc2-416f-a377-547930c69b62
# ╟─6d9b1e2d-032f-4f78-ad59-84c408de490b
# ╟─30364597-1f57-4e30-ade6-31ad6f11949d
# ╠═458b05b7-5655-415d-a9bf-f64d157ac891
# ╠═70f948f9-e526-47c3-98c9-d745493a5e18
# ╠═50917e2b-48b0-41cb-95cd-a2d7b3a8bf7b
# ╠═2ed210ad-a828-40ed-bbc0-dffb2e808b17
# ╠═ad51ce37-7a2d-4741-b69b-f85683ec4b28
# ╠═43165603-47c4-4107-8c7a-6bd4f056939c
# ╠═7d535fe6-f885-4ab6-8660-392887cb8e4b
# ╟─f9244069-3f08-4446-8857-4d824ae13834
# ╟─5f6b4556-54f2-42e8-841f-684c00a19ee9
# ╟─615b3f79-3a72-49ba-bb1d-90a52311c3cd
# ╟─d89e6f1b-697e-4ebb-90cb-2d8649b49a31
# ╠═001a6f78-7e07-46d4-9a58-d53f86fa0239
# ╠═f1c23bf3-15ba-4db7-be8e-49c5abcdf6a8
# ╠═cd603fe7-bcac-4609-91dc-f32ef2a64c83
# ╠═69ad7b3f-45ad-4a3b-8381-f6ac0cd6daf5
# ╠═ddb3bfb8-4c1f-4a6e-b72a-a96e99e75859
# ╠═cebac023-d60e-464c-89bf-ed15f26c432f
# ╠═30b461ba-a69f-4de9-bf84-ecca8a6da742
# ╠═157fc67a-93cf-4bec-8a2d-8be2230dcb51
# ╠═134fd2c0-61ac-4ba3-9246-da571dfdb352
# ╟─d2eeccbe-1f59-4e5e-9e64-342a38b8f477
# ╟─191f060d-4302-4b18-82f8-6e20224ee201
# ╟─9bc57aaf-5cba-4569-975c-a504730b8008
# ╟─6f8174f6-ee58-4c76-81fd-d6c522491339
# ╠═4ee55b92-512b-431e-9973-d0c283aa13d2
# ╠═fe9a9fca-f098-441e-b51f-8c27be65785d
# ╟─8d7afa78-57ae-4fb9-8899-4254f22a11f8
# ╟─2f4b4efd-a892-4797-a9f1-e5f94222af33
# ╠═1eabaff3-9776-4280-a4dc-5686441544f4
# ╠═17fa2078-d17e-43b6-89e9-cc4396582676
# ╟─c2ec1e3a-7ed5-43de-867a-bfa2f25de5af
# ╟─a61c0de8-5ec1-4f58-bad5-488931f53b60
# ╠═ce187832-a000-4d42-bdfd-da238e644f59
# ╠═187edc0f-e9a5-406d-8a33-a326a6ccfa99
# ╠═cf7a32d7-59ca-46dc-9206-dfcd0d030e37
# ╟─0b757a49-b1a3-4d6d-b482-fe7adce2c499
# ╟─9338bc36-9f62-45a9-8805-bd7964c3e090
# ╠═5a36d3de-8341-4500-bff2-b32a88939377
# ╠═84d5dba8-ac87-49eb-b7ca-0b186c8ba203
# ╠═73828383-6b1f-4f8b-ab22-b2c5e1f581a0
# ╠═3cf16498-0770-4308-9962-eac685a43599
# ╠═86128a17-7089-4f55-9be3-c7da169b6ff6
# ╟─b676b7ee-9d47-41dd-a80d-60fa8556a38e
# ╟─01b0c32a-a664-4d1c-9519-17d05a8dddb8
# ╟─9c6745e6-0757-43f4-89da-ae6b04b1803c
# ╟─ef4b6992-1053-4b64-a659-0f45407ce8ef
# ╟─36d65f3b-ef45-4960-b4ae-5e30b94e7500
# ╠═39d542a2-0ecb-4ed0-b0cf-4eca6bd67094
# ╠═20a266d9-7eaa-448c-bef0-209dc5771959
# ╠═a45b4050-302b-4842-a13e-3934041d75c4
# ╠═dbea8be9-f915-4221-a71b-c98e4ca40a92
# ╟─c75a4af4-ec70-41ae-848b-0d4464047936
# ╟─a44953e2-352b-4524-9191-ac0f0d3d9a43
# ╟─bdecea25-78ce-43d1-9ad9-17785889f722
# ╠═8c00452c-2585-4edc-9ed6-b5befaa7991d
# ╠═8a542361-959e-4c6d-b15d-2c5700b2ef3c
# ╠═04077a45-a8a3-4d4b-a558-174c059d413d
# ╠═9bb89559-2490-4f7a-8256-e02239ccff29
# ╟─b586cc56-cef0-4a0a-b31d-b7a9b37ecffa
# ╟─727e16e0-989a-4168-b705-7512e52fbf83
# ╟─e7ef8aa5-8668-4489-97e8-fb46c80a501c
# ╠═7bc4c8dc-0134-47fd-bc6a-18fbe81aff82
# ╠═ac93e9d5-6097-4ecf-b35f-e60952f1284b
# ╠═729abf56-2a75-48af-bb72-9168fcf9560a
# ╠═20683dd5-d1f1-4969-865b-a12c201c3c76
# ╟─6823e7bf-66f5-4a1f-b063-b19fadf071db
# ╟─3ff4a644-f633-4d30-b36d-8b56094bfb20
# ╟─d52ecd4c-9755-491b-94f7-db96fae1eb05
# ╠═f5b244e8-0a8e-41a1-8ead-b0f2922fb125
# ╟─87d3da6f-7537-48b1-9c9e-6ffd9e45f602
# ╠═3cf7fb88-db37-4d34-a13f-efde9000e8b7
# ╟─570496eb-a486-41b4-a441-9c5efc095fab
# ╠═75f39b14-c053-471c-bea8-6657bd70656d
# ╟─bc5cdcfe-cf3e-4e83-b423-712e3a4d13e6
# ╟─9904c284-ef99-4ee6-877c-9e06eb416f1f
# ╠═c72c3408-a7c1-4427-8ebc-1572fe40d1c3
# ╟─19004d6e-b4a5-4754-8564-8bde3e2c2532
# ╠═41e552a4-0976-40d9-a1b1-e40f241dea4a
# ╟─471c8398-1306-4514-ab4c-2e37c04cfe9b
# ╠═329a421b-a0aa-488a-a4a3-1c3d8bb153e6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
