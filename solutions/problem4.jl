using BuildingModelLibrary
using JuliaSimModelOptimizer
using ModelingToolkit
using OrdinaryDiffEq
using Unitful
using DataFrames, CSV
using OptimizationOptimJL
using DataInterpolations

datapath = joinpath(@__DIR__, "..", "data", "USA_AZ_Phoenix.722780_TMY2.mos")
model = initialize_model(;datapath)

df = BuildingModelLibrary.getdata(datapath)
idx_range = 4560:4728+1
ts = df[idx_range, :t] .- df[idx_range[1], :t]
tspan = (0.0, ts[end-1])

data = CSV.read(joinpath(@__DIR__, "..", "data", "building_data.csv"), DataFrame)

trial = Trial(data, model; tspan, alg=QNDF(autodiff=false))

invprob = InverseProblem([trial], model,
    [
        @nonamespace(model.T_fluids_1₊k) => (270.0, 280.0),
        @nonamespace(model.T_fluids_2₊k) => (270.0, 280.0),
        @nonamespace(model.T_fluids_3₊k) => (270.0, 280.0)
    ]
)

# alg = SplineCollocate(maxiters=1000, solver=BFGS(), interp=CubicSpline)
alg = StochGlobalOpt(maxiters=100)
result = calibrate(invprob, alg)
uconvert.(u"°C", result .* u"K")
