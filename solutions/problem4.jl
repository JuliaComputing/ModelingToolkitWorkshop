using BuildingModelLibrary
using JuliaSimModelOptimizer
using ModelingToolkit
using OrdinaryDiffEq
using Unitful
using DataFrames
using OptimizationOptimJL
using DataInterpolations

model = initialize_model(datapath=joinpath(@__DIR__, "..", "data", "USA_AZ_Phoenix.722780_TMY2.mos"))

day = ustrip(u"s", 24u"hr")
n_days = 5.0
t0 = 2.0 * day
tend = t0 + n_days

pre_trial = Trial(nothing, model; tspan=(t0, tend), alg=QNDF(autodiff=false))

pre_invprob = InverseProblem([pre_trial], model,
    [
        @nonamespace(model.T_fluids_1₊k) => (270.0, 280.0),
        @nonamespace(model.T_fluids_2₊k) => (270.0, 280.0),
        @nonamespace(model.T_fluids_3₊k) => (270.0, 280.0)
    ]
)

τ = ustrip(u"K", 5.0u"°C")
data = DataFrame(solve_trial(pre_trial, fill(τ, 3), pre_invprob))

trial = Trial(data, model; tspan=(t0, tend), alg=QNDF(autodiff=false))

invprob = InverseProblem([trial], model,
    [
        @nonamespace(model.T_fluids_1₊k) => (270.0, 280.0),
        @nonamespace(model.T_fluids_2₊k) => (270.0, 280.0),
        @nonamespace(model.T_fluids_3₊k) => (270.0, 280.0)
    ]
)

alg = SplineCollocate(maxiters=1000, solver=BFGS(), interp=CubicSpline)
result = calibrate(invprob, alg)
uconvert.(u"°C", result.*u"K")
