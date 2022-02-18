using Cbc
using HDF5
using JuMP
using YAML
using Configurations

include("make_instance.jl")

model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 8) 

@option struct Configuration
    client_file::String
    client_resource::String
    station_file::String    
    station_resource::String
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


clients = h5read(config.client_file, config.client_resource)
stations = h5read(config.station_file, config.station_resource)


A, Λ = make_instance(clients, stations)
n, m = size(A)

ϕ = 1000rand(n) .+ 1500
ν = rand(1:6, n)
K = rand(50:100)

vars = build_model!(model, A, 1500ones(n), 15000, ϕ, zeros(n), ν, K, 25000, Λ)
optimize!(model)
if termination_status(model) == OPTIMAL
    println(value.(vars.x))
    # println(value.(vars.y))
    println(value.(vars.z))
end
