using Cbc
using HDF5
using JuMP

model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 8) 

n = 100
m = 1000
ϕ = 1000rand(n) .+ 1500
δ = rand(0:1, n, m)
α = rand(n, m) .* δ
λ = 100rand(n)
ν = rand(1:6, n)
K = rand(50:100)


@variable(model, 0 <= x[1:n], Int)
@variable(model, z[1:n], Bin)
@variable(model, y[i = 1:n, j = 1:m], Bin)

@constraint(model, c0, 2000sum(z) + 12000sum(x) + sum(ϕ .* x) <= 100000)
@constraint(model, c1,  x .<= ν .* z)
@constraint(model, c3, K .* x .<=  λ )
c4 = @constraint(model, [i = 1:n], α[i, :]' * y[i, :] <= K .* x[i])

@objective(model, Max, α[:]' * y[:])
#optimize!(model)
if termination_status(model) == OPTIMAL
    println(value.(x))
    println(value.(z))
end