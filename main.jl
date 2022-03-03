using Gurobi
using JuMP
using YAML
using Configurations


include("make_instance.jl")
include("load_instance.jl")

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "threads", 8) 

@option struct Configuration
    instance_file::String
    instance_name::String
end

function load_config(file_name)
    data = YAML.load_file(file_name, dicttype=Dict{String, Any})
    from_dict(Configuration, data)
end



struct ModelVariables 
    x
    y
    z
end

function build_model!(model, A, Γ, c, Φ, Π, Σ, κ, t, Λ)
    n, m = size(A)
    @variable(model, 0 <= x[1:n], Int)
    @variable(model, z[1:n], Bin)
    @variable(model, y[i = 1:n, j = 1:m], Bin)

    @constraint(model, c0, Γ' * z + (c .+ Φ)' * x <= t)
    @constraint(model, c1[i = 1:n],  x[i] <= Σ[i] * z[i])
    @constraint(model, c2, x .>= Π)
    @constraint(model, c3[i  = 1:n], κ * x[i] <=  Λ[i])
    @constraint(model, c4[i = 1:n], A[i, :]' * y[i, :] <= κ * x[i])

    @objective(model, Max, A[:]' * y[:])
    ModelVariables(x, y, z)
end

config = load_config("config.yml")

instance = load_istance(config.instance_file, config.instance_name)

A, Λ = make_instance(instance.clients, instance.stations)
vars = build_model!(model, A, instance.Γ, instance.c, instance.ϕ, instance.Π, instance.Σ, instance.k, instance.T, Λ)
optimize!(model)
if termination_status(model) == OPTIMAL
    println(value.(vars.x))
    println(value.(vars.z))
end
