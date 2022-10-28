module RC

using ModelingToolkit

@variables t
@connector function Pin(;name)
    @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, [v, i], []; name=name)
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t; name=name), g)
end

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    @variables v(t)=1.0 i(t)=1.0
    eqs = [v ~ p.v - n.v, 0 ~ p.i + n.i, i ~ p.i]
    compose(ODESystem(eqs, t, [v, i], []; name=name), [p, n])
end

function ConstantVoltageSource(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    @parameters V = V
    eqs = [v ~ V]
    extend(ODESystem(eqs, t; name=name), oneport)
end

function CustomVoltageSource(;name, V)
    @named oneport = OnePort()
    @unpack v = oneport
    eqs = [v ~ V(t)]
    extend(ODESystem(eqs, t; name=name), oneport)
end

function Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @parameters C = C
    D = Differential(t)
    eqs = [D(v) ~ i / C]
    extend(ODESystem(eqs, t; name=name), oneport)
end
function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @parameters R = R
    eqs = [v ~ i * R]
    extend(ODESystem(eqs, t; name=name), oneport)
end



end # module

using ModelingToolkit, Plots, DifferentialEquations, .RC
t = RC.t
@named source = RC.ConstantVoltageSource(V = 1.0)
@named resistor = RC.Resistor()
@named capacitor = RC.Capacitor()
@named ground = RC.Ground()
eqs = [connect(source.p, resistor.p),
       connect(resistor.n, capacitor.p),
       connect(capacitor.n, ground.g),
       connect(capacitor.n, source.n)]
@named model = ODESystem(eqs, t)
model = compose(model, [source, resistor, capacitor, ground])
sys = structural_simplify(model)
prob = ODEProblem(sys, [capacitor.v => 0.0], (0.0, 10.0))
sol = solve(prob)
plot(sol)

source_fun = t -> sum(x->sin(x * t), (1, 2, 10, 100))

@named source = RC.CustomVoltageSource(V = source_fun)
@named resistor = RC.Resistor()
@named capacitor = RC.Capacitor()
@named ground = RC.Ground()
eqs = [connect(source.p, resistor.p),
       connect(resistor.n, capacitor.p),
       connect(capacitor.n, ground.g),
       connect(capacitor.n, source.n)]
@named model = ODESystem(eqs, t)
model = compose(model, [source, resistor, capacitor, ground])
sys = structural_simplify(model)
prob = ODEProblem(sys, [capacitor.v => 0.0], (0.0, 10.0))
sol = solve(prob)
ts = range(0, 10, length=5000)
p1 = plot(sol, lab = "v(t)")
p2 = plot(ts, source_fun.(ts), label="source")
plot(p1, p2, layout=(2,1))
