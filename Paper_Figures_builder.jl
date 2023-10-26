include("Moment_Estimation_Functions.jl")
using Plots
using JSON
using GLMakie

function scenario_string_name(name, degree_limit_d)
    density_name = name*"_"*string(degree_limit_d)
    return density_name
end


#Occupation density pictures
#First we load the data that was computed in Paper_Data_Calculator.jl
occupation_data_filename = "occupation_densities.json"
densities_dict = JSON.parsefile(occupation_data_filename)
name = "Brownian"
degree_limit_d = 2
density_name = scenario_string_name(name, degree_limit_d)
Data = densities_dict[density_name]
X = Vector{Float64}(Data[1])
Y = Vector{Float64}(Data[2]) 
Z = reduce(hcat,Vector{Vector{Float64}}(Data[3])) 
#Create the contourPlot
f = Figure(resolution= (1024,1024))
ax = Axis(f[1, 1],aspect = 1)
hm = contourf!(ax, X, Y, Z, levels=20)
f
picture_filename = "occupation_"*density_name*".png"
save(picture_filename*".png",f)


#We evaluate the function at a large grid of points for the picture
g = Figure(resolution= (1024,1024))
ax = Axis3(g[1, 1])
surf = surface!(ax, X, Y, Z)
g

#contourf
f = Figure(resolution= (1024,1024))
ax = Axis(f[1, 1],aspect = 1)
hm = contourf!(ax, X, Y, Z, levels=30)
f


"""The figures show the behavior of three stochastic processes in the unit disc, specified by the infinitesimal generators below"""
function Brownian_infinitesimal_operator(p)
    #this function implements the differential operator L for Brownian motion, which is negative ONE HALF of the laplacian
    list = []
    vars = variables(p)
    for var in vars
        push!(list, differentiate(differentiate(p,var),var))
    end
    return((-1/2)*sum(list))
end


function Brownian_with_drift_infinitesimal_operator(p)
    # this function implements the differential operator L for Brownian motion with drift:
    # dX_1 = dW_1
    # dX_2 = 2dt + dW_2 
    vars = variables(p)
    list = [(-2)*differentiate(p,x[2])]
    for var in vars
        push!(list, (-1/2)*differentiate(differentiate(p,var),var))
    end
    return(sum(list))
end


function Sq_bessel_infinitesimal_operator(p)
    #this function implements the differential operator L, in this changes
    #for the two dimensional diffusion which is brownian in one component
    #and the squared bessel in the other
    # dX_1 = dW_1
    # dX_2 = dt + 2 \sqrt{X_2}dW_2 
    list = [-differentiate(p,x[2])]
    push!(list, (-1/2)*differentiate(differentiate(p,x[1]),x[1]))
    push!(list, -2*x[2]*differentiate(differentiate(p,x[2]),x[2]))
    return(sum(list))
end





"""
#Testing the polynomial inverse image map:
n = 2 #dimension two
@polyvar x[1:n]
space_vars_vector = x
differential_operator_L = Brownian_with_drift_infinitesimal_operator #this function is defined above
target_g = x[1]^4*x[2]^4 #This is the function whose boundary moments we wish to estimate
initial_location = [0.0,0.0]
target_functions_list = [target_g]
error_tolerance = 1e-8
maximum_allowed_degree_increases = 7
find_a_low_degree_inverse_image_under_differential_operator(target_g, differential_operator_L)

#lower_bounds, upper_bounds = compute_certified_exit_location_moments(space_vars_vector, target_functions_list, differential_operator_L, initial_location, error_tolerance, maximum_allowed_degree_increases)
"""




