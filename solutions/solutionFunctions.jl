using Plots
using LaTeXStrings

using SpecialFunctions
using Trapz
using LinearAlgebra

using Base.Threads

using DelimitedFiles

## Constants struct

struct Constants
	N::Int64  					# number of modes for solution S(z)
	L::Number                   # half-domain length
	
	# domain definition
	dz::Float64 				# domain spacing
    z::Vector{Float64} 			# domain vector of values (2N + 2 ponts)

	# magnetic constants 
	B::Float64 					# Bond number 
	b::Float64 					# inner rod radius
    E::Float64                  # Bernoulli constant 
	
	
	function Constants(N::Int64, L::Number, B::Float64, b::Float64)
        dz = 2*L / (2*N+1)
        z = collect(-L:dz:L)
		
        E = 1 - B/2

        new(N, L, dz, z, B, b, E)
    end
end


## Main functions

# SOLVERS

function mySolver(f, initial_guess::Vector{Float64}; solver = :NewtonRaphson, tol::Float64 = 1e-8, max_iter::Int64 = 1000)

	x = initial_guess

	if solver == :Newton
		
	    for i in 1:max_iter
	        J = finite_diff_jacobian(f, x)
	        δx = -J \ f(x)  # Newton's update step
	        x += δx
	        if norm(δx) < tol  # Check for convergence
	            return x, i
	        end
	    end
	    error("Failed to converge after $max_iter iterations")

	elseif solver == :NewtonRaphson 
	    
	    # alpha = 1.0  # Initial step size
	    c = 1e-4  # Sufficient decrease parameter
	    rho = 0.5  # Step size reduction factor
		
	    for i in 1:max_iter
	        J = finite_diff_jacobian(f, x)
	        δx = -J \ f(x)  # Newton's update step
	        t = 1.0  # Initialize step size
	
	        # Backtracking line search
	        while norm(f(x + t * δx)) > norm(f(x)) + c * t * dot(f(x), δx)
	            t *= rho
	        end
	
	        x += t * δx  # Update with the found step size
	        if norm(δx) < tol  # Check for convergence
	            return x, i
	        end
	    end
	    error("Failed to converge after $max_iter iterations")

	elseif solver == :Broyden
		J = finite_diff_jacobian(f, x)
		for i in 1:max_iter
			δx = -J \ f(x)  # Broyden's update step
			x_new = x + δx
			if norm(δx) < tol  # Check for convergence
				return x, i
			end
			Δf = f(x_new) - f(x)
			J += (Δf - J * δx) * δx' / (δx' * δx)  # Update Jacobian approximation
			x = x_new
		end
		error("Failed to converge after $max_iter iterations")

	elseif solver == :LevenbergMarquardt
		λ = 1e-3  # Initial damping factor
		ν = 2.0   # Damping factor increment

		for i in 1:max_iter
			J = finite_diff_jacobian(f, x)
			A = J' * J + λ * I  # Levenberg-Marquardt update
			g = J' * f(x)
			δx = -A \ g

			x_new = x + δx
			if norm(f(x_new)) < norm(f(x))  # If the step improves the solution
				x = x_new
				λ /= ν
			else
				λ *= ν
			end

			if norm(δx) < tol  # Check for convergence
				return x, i
			end
		end
		error("Failed to converge after $max_iter iterations")

	elseif solver == :Secant
		δx = ones(length(x)) * 1e-4
		for i in 1:max_iter
			f_val = f(x)
			if norm(f_val) < tol
				return x
			end

			J_approx = finite_diff_jacobian(f, x)
			δx = -J_approx \ f_val
			x += δx

			if norm(δx) < tol  # Check for convergence
				return x, i
			end
		end
		error("Failed to converge after $max_iter iterations")

	else
		error("Enter which algorithm you want to use!")
	end
end

function finite_diff_jacobian(f, x)
    h = 1e-8  # Small perturbation (typically sqrt(machine epsilon))
    n = length(x)
    J = zeros(n, n)
    fx = f(x)

    Threads.@threads for i in 1:n
        x_perturbed = copy(x)
        x_perturbed[i] += h
        J[:, i] = (f(x_perturbed) - fx) / h
    end
    return J
end

# PROBLEM FUNCTIONS

function equations(unknowns::Vector{Float64}, constants::Constants, a₁::Float64, a₀::Float64)

	# Returns the $N + 2$ equations that we want to solve for:

	# problem constants 
	z = constants.z
	N = constants.N
	B = constants.B
	b = constants.b
	E = constants.E
	L = constants.L

	c = unknowns[1]
	coeffs = unknowns[2:N+2] # N + 1 coeffs

	a0 = coeffs[1]
	a1 = coeffs[2]

	S, Sz, Szz = fourierSeries(coeffs, z, L)

	integrands = zeros(N, length(z)) # N integrands for k = 1:N
	integrals = zeros(N) # N integrals (array gets condensed on the z-axis)
	eqs = zeros(N+2) # integrands + 2 extra equations I define later

	# define common factors in equations 
	Szsq = 1 .+ (Sz.^2);

	one_p = (Szsq).*((c.^2)./2 .- 1 ./ (S.*sqrt.(Szsq)) .+ Szz./(Szsq.^(3/2)) .+ B./(2 .* S.^2) .+ E);

	Threads.@threads for n = 1:N

		k = n*π/L 
		
	    one = k .* S .* sqrt.(Complex.(one_p))
	    two = besselk.(1, k * b) .* besseli.(1, k .* S) .- besseli.(1, k * b) .* besselk.(1, k .* S)
		
	    integrands[n, :] = real.(one .* two)
		
	    # Normalize the integrand before integration to prevent numerical issues
	    integrands[n, :] ./= maximum(abs.(integrands[n, :]))
	    integrals[n] = trapz(z, integrands[n, :] .* cos.(k .* z))
	end

	eqs[1:N] = real.(integrals)
	eqs[N+1] = abs(a0 - a₀)
	eqs[N+2] = abs(a1 - a₁)

	return eqs
	
end

function bifurcation(initial_guess, a1Vals, branchN, constants, tol = 1e-8, solver = :NewtonRaphson, max_iter = 1000)

	# compute the birfucation branch for brancN branch points and provided a₁ values, starting at the given intial guess

	
	# initialize solution array
	solutions = zeros(branchN, constants.N+2)

	# initialize convergence array
	iterations = zeros(branchN)
	
	for i = 1:branchN

		f(u::Vector{Float64}) = equations(u, constants, a1Vals[i], 1.0)

		# solve for the current branch point + capture
		solutions[i,:], iterations[i] = mySolver(f, initial_guess[i,:], tol = tol, solver = solver, max_iter = max_iter)

		# update intial guess 
		initial_guess[i+1,:] = solutions[i,:]

        # print progress for every 10% of branch points and;
		# zero the last 20% of coefficients on every (0.1*branchN)th iteration
        if i % Int(round(0.1*branchN)) == 0
            println("Branch point $i of $branchN, $(Int(iterations[i])) iterations.")
			initial_guess[i+1, end - Int(round(0.2*length(initial_guess[1,:]))):end] .= 0
        end
		
	end
	
	return solutions, iterations 
end



## Helper functions

function c0(k, constants::Constants)
	# linearized wave speed for small amplitude waves c(k)

	B = constants.B
	b = constants.b
	
	c0 = sqrt.((1 ./ k).*((-β(1,k,b,1) ./ β(0,k,b,1)) .* (k.^2 .- 1 .+ B)))

	return c0
	
end

function β(n, k, b, S0)

	beta1 = besseli.(1, k*b) .* besselk.(n, k.*S0)
	beta2 = (-1)^n .* besselk.(1, k*b) .* besseli.(n, k.*S0)

	return beta1 + beta2
end

function fourierSeries(coefficients::Vector{Float64}, domain, L::Number)
	
    N = length(coefficients) - 1
	
    S = zeros(length(domain))  		# profile S
    Sz = zeros(length(domain))  	# first derivative Sz
    Szz = zeros(length(domain))  	# second derivative Szz

    Threads.@threads for i in 1:length(domain)

		x = domain[i]
		
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

function plotting(solutions, index::Int, constants::Constants, shift_profiles = true)

	ii = index 

	# un-pack constants 
	z = constants.z
	N = constants.N
	L = constants.L

	branchN = length(solutions[:,1])

	# seperate coeffs and speeds
	coeffs = solutions[:,2:end]
	speeds = solutions[:,1];

	# create array for profiles
	profiles = zeros(branchN,length(z))
	
	# convert profiles
	Threads.@threads for i = 1:branchN
		profiles[i,:] .= fourierSeries(coeffs[i,:], z, L)[1]
	end

	# shift to the right profiles over by L
	if shift_profiles == true
		profiles = [profiles[:,Int(end/2)+1:end] profiles[:,1:Int(end/2)]]; nothing
	end

	# plot profiles 
	profile_plot = plot(z, profiles[ii,:], legend=false, title = "a1 = $(round(coeffs[ii,2], digits=3))", lw=2)
	xlabel!(L"z"); ylabel!(L"S")

	# plot coeffs 
	first_coeff = 0
	coeff_plot = scatter(abs.(coeffs[ii,first_coeff+1:end]), legend=false, title="k = $(round(ii*π/L))", xticks = :all, yaxis=:log)
	xlabel!("a$(first_coeff) to a$(length(coeffs[1,:])-1)")

	# plot branch
	branch_plot = scatter(speeds[1:ii], coeffs[1:ii,2], legend = false, markersize=4)
	xlabel!(L"c"); ylabel!(L"a_1")

	return profile_plot, branch_plot, coeff_plot
	
end