#= Program for Problem_Set_2
    Fengfan Xiang
    January 2022
    NGM with leisure in utility function
=#
using Plots
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl

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
        δ::Float64 = 1          ; # depreciation rate
        β::Float64 = 0.98       ; # discount factor
        # VFI Parameters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-9  ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
    end

    # Allocate paramters to object p for future calling
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
χ = c_ss^(-p.σ)*(1-p.α)*p.z*kl_ss_ratio^p.α;

println(" ")
println("------------------------")
println("Result for Part 4")
println("χ that gives a steady state labor as 0.4 = ", χ)
println("------------------------")
println(" ")
