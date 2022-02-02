using Pkg
Pkg.activate(".")

using JuMP
using Cbc

model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 8) 

n = 100;
m = 10000;
ϕ = 1000rand(n) .+ 1500;
δ = rand(0:1, n, m);
α = rand(n, m) .* δ ;
λ = 100rand(n);
ν = rand(1:6, n)



@variable(model, 0 <= x[1:n], Int)
@variable(model, z[1:n], Bin)

@constraint(model, c0, 2000sum(z) + 15000sum(x) + sum(ϕ .* x) <= 1000000)
@constraint(model, c1,  x .<= ν .* z)

@objective(model, Max, sum(sum(α, dims=2) .* x))
optimize!(model)
if termination_status(model) == OPTIMAL
    println(value.(x))
    println(value.(z))
end