#= Program for Problem_Set_4
    Fengfan Xiang & Brice Gueyap    
    January 2022
    NGM with leisure in utility function
=#
using Plots
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
using Roots # Pkg.add("Roots")
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
using Optim # Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
    using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots # Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Parameters
    # Generate structure for parameters using Parameters module
    # We can default values for our parameters
    @with_kw struct Par
        # Model Parameters
        α::Float64 = 1/3        ; # Productivity function
        z::Float64 = 1          ; # Productivity
        σ::Float64 = 2          ; # utility function consumption
        η::Float64 = 1          ; # utility function labor
        δ::Float64 = 0.05          ; # depreciation rate
        β::Float64 = 0.98       ; # discount factor
        # VFI Parameters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-8  ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-8  ; # Tolerance for policy function iteration
    end

    # Allocate paramters to object p for future calling
    p = Par()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#= Part 4
    solve for \chi that l_ss = 0.4
=#
kl_ss_ratio = (1/(p.α*p.z) * (1/p.β-1+p.σ))^(1/(p.α-1));
# c_ss when l_ss = 0.4
l_ss = 0.4;
c_ss = (p.z*kl_ss_ratio^p.α - p.σ*kl_ss_ratio)*l_ss;
χ = c_ss^(-p.σ)*(1-p.α)*p.z*kl_ss_ratio^p.α * (1/(l_ss^p.η));

println(" ")
println("------------------------")
println("Result for Part 4")
println("χ that gives a steady state labor as 0.4 = ", χ)
println("------------------------")
println(" ")

k_ss = kl_ss_ratio*l_ss;

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#= Part 5
        Solve the planner's problem using: Plain VFI
=#

#-------------------------------------------------------------------------------
#=
    useful functions for Part 5
=#

# Grid for k_grid
function Make_K_Grid(n_k, kss)
    k_grid = range(1e-5,2*kss,length=n_k) # Equally spaced grid between 0 and 2*kss
    return k_grid
end

# utility function
function utility(k,kp,χ,p::Par)
    @unpack σ, η, α, δ, z = p
    # use bisection method solve for l
    f(l_temp) = (1-α)*(z*k^α*l_temp^(1-α) + (1-δ)*k - kp)^(-σ)*z*k^α*l_temp^(-α)-χ*l_temp^η
    l_temp = find_zero(f,(0,10),Bisection())
    if 0<= l_temp <=1
        l = l_temp
    else
        l = 1
    end
    # # Given l, then it is easy to get c, thus the utility
    c = z*k^α*l^(1-α) + (1-δ)*k - kp
    if c>0
        return c^(1-σ)/(1-σ) - χ* l^(1+η)/(1+η)
    else
        return -Inf
    end
end

# Euler Error
function Euler_Error(k,kp,kpp,l,p::Par)
    # Return percentage error in Euler Equation
    @unpack z, α, δ, β, σ = p
    LHS = (z*k^α*l^(1-α) + (1-δ)*k - kp)^(-σ)
    RHS =  β*(z*kp^α^l^(1-α) +(1-δ)*kp - kpp)^(-σ) * (α*z*kp^(α-1)*l^(1-α)+1-σ)
    return (RHS/LHS-1)*100
end

# Spline interpolation
function spline_linear_interpolation(x,x_grid,y_grid)
    N = length(x_grid)
    ind = findmax(sign.(x_grid .- x))[2] - 1
    A_x = (x_grid[ind+1] - x)/(x_grid[ind+1]-x_grid[ind])

    V_hat = A_x*y_grid[ind] + (1-A_x)*y_grid[ind+1]

    return V_hat
end

#-------------------------------------------------------------------------------
#=
    Plain VFI
=#
# VFI - Grid Search - Matrices
println(" "); println(" "); println(" ")
println("----------------------------------------------------")
println("----------------------------------------------------")
println("Value Function Iteration - Grid Search with Matrices")

# Define function for value update and policy function
function T_grid_mat(V_old, U_mat, p::Par)
    @unpack β = p
    n_k = length(V_old)
    V, G_kp = findmax(U_mat .+  β*repeat(V_old', n_k, 1) , dims = 2 )
    G_kp = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
    return V, G_kp
end

# Solve VFI with grid search and loops
function VFI_grid(T::Function, k_grid, p::Par)
    # VFI parameters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    n_k    = length(k_grid) ; # Number of grid nodes
    V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
    V_dist = 1              ; # Initialize distance
    for iter = 1:max_iter
        # Update value function
        V_new, G_kp = T(V_old)
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # Update old function
        V_old  = V_new
        if V_dist<= dist_tol
            println(" ")
            println("------------------------")
            println("VFI - Grid Search - n_k=$n_k")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            return V_new, G_kp
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in VFI - Grid Search - Solution not found")
end

# Solve VFI with grid search using Matrices
function Solve_VFI_mat(n_k,χ,k_ss,p::Par)
    # Get Grid
    k_grid = Make_K_Grid(n_k,k_ss)
    # Utility matrix
    U_mat = [utility(k_grid[i],k_grid[j],χ,p) for i in 1:n_k, j in 1:n_k]
    # Solve VFI
    V, G_kp = VFI_grid(x->T_grid_mat(x,U_mat,p),k_grid,p)
    return V, G_kp, k_grid, U_mat
end


#-------------------------------------------------------------------------------
#=
    Plain VFI + Interpolations & Optmization
=#

function VFI_interpolation(n_k,χ,k_ss,p::Par)
    @unpack β, z, α = p
    # VFI in discrete grids
    Sol = Solve_VFI_mat(n_k,χ,k_ss,p)
    V_discrete = Sol[1]
    kgrid = Sol[3]
    G_kp = zeros(n_k)
    V_new = zeros(n_k)
    V_p(kp) = spline_linear_interpolation(kp,kgrid,V_discrete)


    for i = 1:n_k
        # Define objective function
        F(kp) = - utility(kgrid[i],kp,χ,p) - β*V_p(kp)
        dF(x) = ForwardDiff.derivative(F,x)

        # Find the feasible range of kp
        k_min = kgrid[1] ; k_max = kgrid[n_k]
        #z*kgrid[i]^α


        # Check for corner solutions
        opt = optimize(F,k_min,k_max)

        if -dF(k_min) < 0
            kp = k_min
        elseif -dF(k_max) > 0
            kp = k_max
        else
            kp = opt.minimizer
        end
        # Solve minimization
        G_kp[i] = kp
        V_new[i] = -F(kp)
    end

    return G_kp, V_new
end


@time vfi_res = VFI_interpolation(50,χ,k_ss,p)



#-------------------------------------------------------------------------------
#=
    Report Euler Error Results
=#

function euler_error_report(n_k,χ,k_ss,p::Par)
    @unpack σ, η, α, δ, z = p
    V, G_kp, k_grid = Solve_VFI_mat(n_k,χ,k_ss,p)
    e_error = zeros(n_k)
    for i = 1:n_k
        k = k_grid[i]
        kp = k_grid[G_kp[i]]
        kpp = k_grid[G_kp[G_kp[i]]]
        f(l_temp) = (1-α)*(z*k^α*l_temp^(1-α) + (1-δ)*k - kp)^(-σ)*z*k^α*l_temp^(-α)-χ*l_temp^η
        l_temp = find_zero(f,(0,10),Bisection())
        if 0<= l_temp <=1
            l = l_temp
        else
            l = 1
        end
        e_error[i] = Euler_Error(k,kp,kpp,l,p)
    end
    # return
    return e_error
end

test_e_error = euler_error_report(50,χ,k_ss,p)
