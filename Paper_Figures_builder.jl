#The functions in this file construct the pictures of the article using data previously computed in Paper_Data_Calculator.jl
#which is where most of the computation takes place.
include("Moment_Estimation_Functions.jl")
using JSON
using GLMakie

function scenario_string_name(name, degree_limit_d)
    density_name = name*"_"*string(degree_limit_d)
    return density_name
end


#1. Occupation density estimations
#First we load the data that was computed via Paper_Data_Calculator.jl
occupation_data_filename = "occupation_densities.json"
densities_dict = JSON.parsefile(occupation_data_filename)
name = "Bessel" #
degree_limit_d = 2
density_name = scenario_string_name(name, degree_limit_d)
Data = densities_dict[density_name]
X = Vector{Float64}(Data[1])
Y = Vector{Float64}(Data[2]) 
Z = reduce(hcat,Vector{Vector{Float64}}(Data[3])) 
#The we create the contourPlot
f = Figure(resolution= (1024,1024))
ax = Axis(f[1, 1],aspect = 1)
hm = GLMakie.contourf!(ax, X, Y, Z, levels=50)
f
picture_filename = "occupation_"*density_name
save(picture_filename*".png",f)


#We evaluate the function at a large grid of points for the picture
g = Figure(resolution= (1024,1024))
ax = Axis3(g[1, 1])
surf = GLMakie.surface!(ax, X, Y, Z)
g

#contourf
f = Figure(resolution= (1024,1024))
ax = Axis(f[1, 1],aspect = 1)
hm = GLMakie.contourf!(ax, X, Y, Z, levels=30)
f


#2. Exit location density estimations:
#First we load the data that was computed via Paper_Data_Calculator.jl
using Plots
data_filename = "exit_densities.json"
densities_dict = JSON.parsefile(data_filename)
thetas = [t for t in range(-pi/2, 3*pi/2, length=300)]
points = [[cos(theta),sin(theta)] for theta in thetas]

#We specify the figures we wish to plot via a dictionary specified in the following function,
function chosen_pictures_dict()
    #Each picture is specified by its filename 
    #and a dictionary of parameters specifying the operator name and the approximation_degrees we wish to plot.
    pictures_dict=Dict(
        "Exit_location_density_1" => Dict("name"=>"Brownian", "desired_degrees"=>[6,8]), 
        "Exit_location_density_2" => Dict("name"=>"Brownian_drift", "desired_degrees"=> [6,8]),
        "Exit_location_density_3" => Dict("name"=>"Bessel", "desired_degrees"=> [6,8])
        )
    #This function also verifies that the desired data exists in our json file with data
    for (key,value) in pictures_dict
        @assert(value["name"] in keys(operators_dict))
        for k in value["desired_degrees"]
            @assert(scenario_string_name(value["name"], k) in keys(densities_dict), scenario_string_name(value["name"], k)*" not present in computed densities.")
        end
    end
    return pictures_dict
end
#We verify that we have data for the desired scenariosi
pc_dict = chosen_pictures_dict()
xlabels_array = [round(x,sigdigits=2) for x in range(-pi/2, 3*pi/2, length=10)]
using Plots
for (picture_filename, specification_dict) in pc_dict
    Plots.plot(dpi=400, xticks=xlabels_array)
    name = specification_dict["name"]
    for k in specification_dict["desired_degrees"]
        data_label = scenario_string_name(name, k) 
        y = densities_dict[data_label]
        y = Vector{Float32}(y)
        Plots.plot!(thetas, y, label = data_label, linewidth = 2)
    end
    savefig(picture_filename*".png")
end






