using Cbc
using HDF5
using JuMP

model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 8) 

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


n = 10
m = 10
ϕ = 1000rand(n) .+ 1500
δ = rand(0:1, n, m)
α = rand(n, m) .* δ
λ = 100rand(n)
ν = rand(1:6, n)
K = rand(50:100)

vars = build_model!(model, α, 1500ones(n), 15000, ϕ, zeros(n), ν, K, 25000, λ)
# optimize!(model)
if termination_status(model) == OPTIMAL
    println(value.(vars.x))
    println(value.(vars.y))
    println(value.(vars.z))
end
