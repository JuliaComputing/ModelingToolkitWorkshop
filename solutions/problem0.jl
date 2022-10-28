# https://github.com/JuliaComputing/ModelingToolkitWorkshop
using ModelingToolkit, Plots, DifferentialEquations

@variables t x(t) y(t)
@parameters α β δ γ
D = Differential(t)
eqs = [
           D(x) ~ α*x - β*x*y
           D(y) ~ δ*x*y - γ*y
       ];

@named model = ODESystem(eqs, t);

prob = ODEProblem(model, [x => 0.9, y=>1.8], (0, 20.0),
                        [α => 2/3, β => 4/3, γ => 1, δ => 1])
sol = solve(prob)
plot(sol)
