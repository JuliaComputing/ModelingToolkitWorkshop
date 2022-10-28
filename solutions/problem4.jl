using BuildingModelLibrary
using JuliaSimModelOptimizer
using ModelingToolkit
using OrdinaryDiffEq
using Unitful
using DataFrames, CSV
using OptimizationOptimJL
using DataInterpolations

model = initialize_model(datapath=joinpath(@__DIR__, "..", "data", "USA_AZ_Phoenix.722780_TMY2.mos"))

day = ustrip(u"s", 24u"hr")
n_days = 5.0
t0 = 2.0 * day
tend = t0 + n_days

data = CSV.read(joinpath(@__DIR__, "..", "data", "building_data.csv"), DataFrame)

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
