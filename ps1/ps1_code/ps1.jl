# This is the code for Problem set 1 part 7
α = 1/3;
# Assume that β = 0.96 for this question
β = 0.96;
z = 1.0;

# Policy function analytical solution
function ngm_policy_fn(k,α,β,z)
    return α*β*z*k^α;
end

# Solve for the steady state given α = 1/3 and z = 1
k_ss_original = (α*β*z)^(1/(1-α));

# (a) Captial decreases to 80% of its steady state value
k_eighty_percentage = 0.8 *  k_ss_original;

# Start with k_eighty_percentage
k_path_a_temp = [k_eighty_percentage]; # Path a is the path of capital for part a, temp means after the capital decreases
k_temp_in = k_eighty_percentage;
k_temp_out = 0;

while abs(k_ss_original - k_temp_out) > 1e-6
    k_temp_out = ngm_policy_fn(k_temp_in,α,β,z);
    k_path_a_temp = push!(k_path_a_temp,k_temp_out);
    k_temp_in = k_temp_out;
end

k_path_a = [k_ss_original;k_path_a_temp]; # Result for capital path in (a)

# (b) Now the productivity increases by 5% permanetly
z = z*1.05; 

# New steady state with higher productivity
k_ss_high_productivity = (α*β*z)^(1/(1-α));

# Start with original steady state
k_path_b = [k_ss_original]; # Path b is the path of capital for part b
k_temp_in = k_ss_original;
k_temp_out = 0;

while abs(k_ss_high_productivity - k_temp_out) > 1e-6
    k_temp_out = ngm_policy_fn(k_temp_in,α,β,z);
    k_path_b = push!(k_path_b,k_temp_out);
    k_temp_in = k_temp_out;
end

# Graphs
using Plots
# (a)
size_t_a = length(k_path_a);
t_a = 1:1:size_t_a;
gr()
plot(t_a,k_path_a,linewidth=2,marker=(:diamond,9),markercolor=RGB(0.1,0.1,0.1))
plot(t_a,k_path_a)

# (b）