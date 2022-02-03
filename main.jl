using Pkg
Pkg.activate(".")

using JuMP
using Cbc

model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 8) 

n = 100
m = 10000
ϕ = 1000rand(n) .+ 1500
δ = rand(0:1, n, m)
α = rand(n, m) .* δ
λ = 100rand(n)
ν = rand(1:6, n)
K = rand(50:100)


@time @variable(model, 0 <= x[1:n], Int)
@time @variable(model, z[1:n], Bin)
@time @variable(model, y[i = 1:n, j = 1:m], Bin)

@time @constraint(model, c0, 2000sum(z) + 12000sum(x) + sum(ϕ .* x) <= 100000)
@time @constraint(model, c1,  x .<= ν .* z)
@time @constraint(model, c3, K .* x .<=  λ )
 
@time c4 = @constraint(model, [i = 1:n], α[i, :]' * y[i, :] <= K .* x[i])

@time @objective(model, Max, α[:]' * y[:])
optimize!(model)
if termination_status(model) == OPTIMAL
    println(value.(x))
    println(value.(z))
end