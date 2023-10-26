include("Moment_Estimation_Functions.jl")
using JSON
"""In this file we carry out the computations needed for the figures. The resulting data is saved in a .json file
which is then loaded by Paper_Figures_builder.jl to actually make the figures of the article"""

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

function scenario_string_name(name, degree_limit_d)
    density_name = name*"_"*string(degree_limit_d)
    return density_name
end

function compute_occupation_measure_data(scenario_name, differential_operator_L, degree_limit_d, occupation_data_filename = "occupation_densities.json")
    densities_dict = Dict()
    try 
        densities_dict = JSON.parsefile(occupation_data_filename)
    catch e
        #File does not exist, creating...
    end

    n = 2 #dimension two
    @polyvar x[1:n]
    space_vars_vector = x
    degree_limit_d = 2
    initial_condition = [0.0,0.5]
    required_error_bound = 1e-5
    #Run over the computation
    println("computing moments for "*name)
    moments_condition_number, maxgap, density_estimation_func = compute_certified_occupation_density_estimation_ball(degree_limit_d,space_vars_vector, differential_operator_L, initial_condition, required_error_bound)
    println(name)
    println("Moments condition number "*string(moments_condition_number))
    println("MaxGap "*string(maxgap))
    println("______________________________________")
    #
    NCD_func = density_estimation_func
    #Next we want to evaluate the density with respect to the Lebesgue measure so we need to factor out the equilibrium density
    function Lebesgue_density_CD(a,b)
        #This function computes the estimated Lebesgue density of the occupation measure
        density_estimate = NCD_func(x=>[b,a])
        if (a^2+b^2) < 0.95
            m = 1/(pi*sqrt(1-(a^2+b^2)))
        else
            m=0
        end
        return minimum([1,m*density_estimate])
    end    
    #The following function sets up the evaluation
    function density_wrt_Equilibrium(; n=49)
        #This function evalutes
        X = range(-1.0,1.0,n)
        Y = range(-1.0,1.0,n)
        #Interesting notation...
        a = Lebesgue_density_CD.(X',Y)    
        return (X,Y,a)
    end
    #We sample the data and save the result in our JSON
    X,Y,Z = density_wrt_Equilibrium(; n=400)
    density_name = scenario_string_name(name, degree_limit_d)
    #and write the result in our densities_dict (possibly overriding previous computations)
    densities_dict[density_name] = [X,Y,Z]        
    #and add this to the file
    open(occupation_data_filename, "w") do f
        JSON.print(f,densities_dict)    
    end
    return X, Y, Z
end




#First we compute the data for the evaluations of our density estimate for the occupation measure
#we will update older computations with same parameters with new data but also will keep old data 
#computed with different parameters

operators_dict=Dict("Brownian" => Brownian_infinitesimal_operator, "Brownian_drift"=> Brownian_with_drift_infinitesimal_operator, "Bessel"=>Sq_bessel_infinitesimal_operator)
occupation_data_filename = "occupation_densities.json"
name = "Brownian"
differential_operator_L = operators_dict[name]
degree_limit_d = 2
res = compute_occupation_measure_data(name,differential_operator_L,degree_limit_d)
#moments_condition_number, maxgap, density_estimation_func = compute_certified_occupation_density_estimation_ball(degree_limit_d,space_vars_vector, differential_operator_L, initial_condition, required_error_bound)


