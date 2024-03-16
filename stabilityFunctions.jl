using Plots
using LinearAlgebra
using SpecialFunctions
using LaTeXStrings
using Trapz
using DelimitedFiles

## Solvers 

function solveGenEig(solution, Nmodes, Nmu)
    
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
	A = Ag(N, z, S0z, q0z, c)
	B = Bg(N, z)
	D = zeros(2N+1,2N+1)
	
	
	Threads.@threads for i = 1:length(μ)

		# create matrices
		
		E = Eg(N, z, S0, S0z, S0zz, q0z, c, Bond, μ[i])
		F = Fg(N, z, S0, S0z, q0z, c, μ[i])
	
		C = Cg(N, z, S0, b, μ[i])
		G = Gg(N, z, S0, b, c, μ[i])
		H = Hg(N, z, S0, b, μ[i])
	
		lhs = [A B; C D]
	    rhs = [E F; G H]
		
		# solve problem 
		solutions = eigen(rhs, lhs)
		
		# save solution
		λ[:,i] = (solutions.values)
	
	end

	return λ
end

## Helper functions

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

## Matrices

# LOCAL 
function Ag(N, z, S0z, q0z, c)

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

function Bg(N, z)

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

function Eg(N, z, S0, S0z, S0zz, q0z, c, B, μ)

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

function Fg(N, z, S0, S0z, q0z, c, μ)

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

# NON-LOCAL

function Cg(N, z, S0, b, μ)

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

function Gg(N, z, S0, b, c, μ)

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

function Hg(N, z, S0, b, μ)

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