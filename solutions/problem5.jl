using JuliaSimControl, ModelingToolkit, OrdinaryDiffEq, JuliaSimControl.MPC, LowLevelParticleFilters, LinearAlgebra, Plots
using LinearAlgebra:cross, inv

rcam_constants = Dict{Symbol, Float64}(:m => 120000,
		:c̄ => 6.6,
		:lt => 24.8,
		:S => 260,
		:St => 64,
		:Xcg => 0.23 * 6.6,
		:Ycg => 0.0,
		:Zcg => 0.10 * 6.6,
		:Xac => 0.12 * 6.6,
		:Yac => 0.0,
		:Zac => 0.0,
		:Xapt1 => 0,
		:Yapt1 => -7.94,
		:Zapt1 => -1.9,
		:Xapt2 => 0,
		:Yapt2 => 7.94,
		:Zapt2 => -1.9,
		:g => 9.81,
		:depsda => 0.25,
		:α_L0 => deg2rad(-11.5),
		:n => 5.5,
		:a3 => -768.5,
		:a2 => 609.2,
		:a1 => -155.2,
		:a0 => 15.212,
		:α_switch => deg2rad(14.5),
		:ρ => 1.225,)

		# Inertia matrices
		
		Ib_mat =  [40.07          0.0         -2.0923
               		0.0            64.0        0.0  
              	  -2.0923        0.0         99.92]*rcam_constants[:m]
		
		Ib = Ib_mat
		
		invIb = inv(Ib_mat)

        ps = @parameters(
		ρ             = rcam_constants[:ρ], [description="kg/m3 - air density"],
		m             = rcam_constants[:m], [description="kg - total mass"],
		c̄             = rcam_constants[:c̄], [description="m - mean aerodynamic 																		chord"],
		lt            = rcam_constants[:lt], [description="m - tail AC distance to 																			CG"],
		S             = rcam_constants[:S], [description="m2 - wing area"],
		St            = rcam_constants[:St], [description="m2 - tail area"],
		Xcg           = rcam_constants[:Xcg], [description="m - x pos of CG in Fm"],
		Ycg           = rcam_constants[:Ycg], [description="m - y pos of CG in Fm"],
		Zcg           = rcam_constants[:Zcg], [description="m - z pos of CG in Fm"],
		Xac           = rcam_constants[:Xac], [description="m - x pos of aerodynamic 															center in Fm"],
		Yac           = rcam_constants[:Yac], [description="m - y pos of aerodynamic 															center in Fm"],
		Zac           = rcam_constants[:Zac], [description="m - z pos of aerodynamic 															center in Fm"],
		Xapt1         = rcam_constants[:Xapt1], [description="m - x position of 																engine 1 in Fm"],
		Yapt1         = rcam_constants[:Yapt1], [description="m - y position of 																engine 1 in Fm"],
		Zapt1         = rcam_constants[:Zapt1], [description="m - z position of 																engine 1 in Fm"],
		Xapt2         = rcam_constants[:Xapt2], [description="m - x position of 																engine 2 in Fm"],
		Yapt2         = rcam_constants[:Yapt2], [description="m - y position of 																engine 2 in Fm"],
		Zapt2         = rcam_constants[:Zapt2], [description="m - z position of 																engine 2 in Fm"],
		g             = rcam_constants[:g], [description="m/s2 - gravity"],
		depsda        = rcam_constants[:depsda], [description="rad/rad - change in 																downwash wrt α"],
		α_L0          = rcam_constants[:α_L0], [description="rad - zero lift AOA"],
		n             = rcam_constants[:n], [description="adm - slope of linear 														region of lift slope"],
		a3            = rcam_constants[:a3], [description="adm - coeff of α^3"],
		a2            = rcam_constants[:a2], [description="adm -  - coeff of α^2"],
		a1            = rcam_constants[:a1], [description="adm -  - coeff of α^1"],
		a0            = rcam_constants[:a0], [description="adm -  - coeff of α^0"],
		α_switch      = rcam_constants[:α_switch], [description="rad - kink point of 															lift slope"],
    )

    @parameters t
    Dt =  Differential(t)

    V_b = @variables(
	    u(t), [description="translational velocity along x-axis [m/s]"],
	    v(t), [description="translational velocity along y-axis [m/s]"],
	    w(t), [description="translational velocity along z-axis [m/s]"]
	    )
	
	wbe_b = @variables(
	    p(t), [description="rotational velocity about x-axis [rad/s]"], 
	    q(t), [description="rotational velocity about y-axis [rad/s]"],
	    r(t), [description="rotational velocity about z-axis [rad/s]"]
	    )
	
	rot = @variables(
	        ϕ(t), [description="rotation angle about x-axis/roll or bank angle 																			[rad]"], 
	        θ(t), [description="rotation angle about y-axis/pitch angle [rad]"], 
	        ψ(t), [description="rotation angle about z-axis/yaw angle [rad]"]
	        )
	
	# Controls
	U = @variables(
	    uA(t), [description="aileron [rad]"],
	    uT(t), [description="tail [rad]"],
	    uR(t), [description="rudder [rad]"],
	    uE_1(t), [description="throttle 1 [rad]"],
	    uE_2(t), [description="throttle 2 [rad]"],
	)

    # Auxiliary Variables to define model.
	Auxiliary_vars = @variables Va(t) α(t) β(t) Q(t) CL_wb(t) ϵ(t) α_t(t) CL_t(t) CL(t) CD(t) CY(t) F1(t) F2(t)
	
	@variables FA_s(t)[1:3] C_bs(t)[1:3,1:3] FA_b(t)[1:3] eta(t)[1:3] dCMdx(t)[1:3, 1:3] dCMdu(t)[1:3, 1:3] CMac_b(t)[1:3] MAac_b(t)[1:3] rcg_b(t)[1:3] rac_b(t)[1:3] MAcg_b(t)[1:3] FE1_b(t)[1:3] FE2_b(t)[1:3] FE_b(t)[1:3] mew1(t)[1:3] mew2(t)[1:3] MEcg1_b(t)[1:3] MEcg2_b(t)[1:3] MEcg_b(t)[1:3] g_b(t)[1:3] Fg_b(t)[1:3] F_b(t)[1:3] Mcg_b(t)[1:3] H_phi(t)[1:3,1:3]

	# Scalarizing all the array variables. 
	FA_s, C_bs, FA_b, eta, dCMdx, dCMdu, CMac_b, MAac_b, rcg_b, rac_b, MAcg_b, FE1_b, FE2_b, FE_b, mew1, mew2, MEcg1_b, MEcg2_b, MEcg_b, g_b, Fg_b, F_b, Mcg_b, H_phi = collect(FA_s), collect(C_bs), collect(FA_b), collect(eta), collect(dCMdx), collect(dCMdu), collect(CMac_b), collect(MAac_b), collect(rcg_b), collect(rac_b), collect(MAcg_b), collect(FE1_b), collect(FE2_b), collect(FE_b), collect(mew1), collect(mew2), collect(MEcg1_b), collect(MEcg2_b), collect(MEcg_b), collect(g_b), collect(Fg_b), collect(F_b), collect(Mcg_b), collect(H_phi)

	array_vars = vcat(vec(FA_s), vec(C_bs), vec(FA_b), vec(eta), vec(dCMdx), vec(dCMdu), vec(CMac_b), vec(MAac_b), vec(rcg_b), vec(rac_b), vec(MAcg_b), vec(FE1_b), vec(FE2_b), vec(FE_b), vec(mew1), vec(mew2), vec(MEcg1_b), vec(MEcg2_b), vec(MEcg_b), vec(g_b), vec(Fg_b), vec(F_b), vec(Mcg_b), vec(H_phi))

    eqns =[
    # Step 1. Intermediate variables 
    # Airspeed
    Va ~ sqrt(u^2 + v^2 + w^2)

    # α and β
    α ~ atan(w,u)
    β ~ asin(v/Va)

    # dynamic pressure
    Q ~ 0.5*ρ*Va^2


    # Step 2. Aerodynamic Force Coefficients
    # CL - wing + body
    CL_wb ~  n*(α - α_L0)

    # CL thrust
    ϵ ~ depsda*(α - α_L0)
    α_t ~ α - ϵ + uT + 1.3*q*lt/Va
    CL_t ~ 3.1*(St/S) * α_t

    # Total CL
    CL ~ CL_wb + CL_t

    # Total CD
    CD ~ 0.13 + 0.07 * (n*α + 0.654)^2

    # Total CY
    CY ~ -1.6*β + 0.24*uR


    # Step 3. Dimensional Aerodynamic Forces
    # Forces in F_s
    FA_s .~ [-CD * Q * S
             CY * Q * S
             -CL * Q * S] 


    # rotate forces to body axis (F_b)  
    vec(C_bs .~ [cos(α)      0.0      -sin(α)
                0.0          1.0      0.0
                sin(α)      0.0      cos(α)])


    FA_b .~ C_bs*FA_s 

    # Step 4. Aerodynamic moment coefficients about AC
    # moments in F_b
    eta .~ [ -1.4 * β 
            -0.59 - (3.1 * (St * lt) / (S * c̄)) * (α - ϵ)
            (1 - α * (180 / (15 * π))) * β
    ]


    vec(dCMdx .~ (c̄ / Va)* [-11.0              0.0                           5.0
                            0.0     (-4.03 * (St * lt^2) / (S * c̄^2))        0.0
                            1.7                 0.0                          -11.5])


    vec(dCMdu .~ [-0.6                   0.0                 0.22
                   0.0   (-3.1 * (St * lt) / (S * c̄))         0.0
                   0.0                    0.0                -0.63])

    # CM about AC in Fb
    CMac_b .~ eta + dCMdx*wbe_b + dCMdu*[uA
                                        uT
                                        uR]
                                        
    # Step 5. Aerodynamic moment about AC 
    # normalize to aerodynamic moment
    MAac_b .~ CMac_b * Q * S * c̄

    # Step 6. Aerodynamic moment about CG
    rcg_b .~    [Xcg
                Ycg
                Zcg]

    rac_b .~ [Xac
             Yac
             Zac]

    MAcg_b .~ MAac_b + cross(FA_b, rcg_b - rac_b)

    # Step 7. Engine force and moment
    # thrust
    F1 ~ uE_1 * m * g
    F2 ~ uE_2 * m * g

    # thrust vectors (assuming aligned with x axis)
    FE1_b .~ [F1
              0
              0]

    FE2_b .~ [F2
              0
              0]

    FE_b .~ FE1_b + FE2_b

    # engine moments
    mew1 .~  [Xcg - Xapt1
              Yapt1 - Ycg
              Zcg - Zapt1]

    mew2 .~ [ Xcg - Xapt2
              Yapt2 - Ycg
              Zcg - Zapt2]

    MEcg1_b .~ cross(mew1, FE1_b)
    MEcg2_b .~ cross(mew2, FE2_b)

    MEcg_b .~ MEcg1_b + MEcg2_b

    # Step 8. Gravity effects
    g_b .~ [-g * sin(θ)
             g * cos(θ) * sin(ϕ)
             g * cos(θ) * cos(ϕ)]

    Fg_b .~ m * g_b

    # Step 9: State derivatives

    # form F_b and calculate u, v, w dot
    F_b .~ Fg_b + FE_b + FA_b

    Dt.(V_b) .~ (1 / m)*F_b - cross(wbe_b, V_b)

    # form Mcg_b and calc p, q r dot
    Mcg_b .~ MAcg_b + MEcg_b

    Dt.(wbe_b) .~ invIb*(Mcg_b - cross(wbe_b, Ib*wbe_b))

    # phi, theta, psi dot
    vec(H_phi .~ [1.0         sin(ϕ)*tan(θ)       cos(ϕ)*tan(θ)
                  0.0         cos(ϕ)              -sin(ϕ)
                  0.0         sin(ϕ)/cos(θ)       cos(ϕ)/cos(θ)])

    Dt.(rot) .~  H_phi*wbe_b        
]

all_vars = vcat(V_b,wbe_b,rot,U, Auxiliary_vars, array_vars)

@named rcam_model = ODESystem(eqns, t, all_vars, ps)

inputs = [uA, uT, uR, uE_1, uE_2]
outputs = [u,v, w, p, q, r, ϕ, θ, ψ]
sys0, diff_idxs0, alge_idxs0, input_idxs0 = ModelingToolkit.io_preprocessing(rcam_model, inputs, outputs)
reduced_states = states(sys0)
sys0

x0 = Dict(
    u => 87.0,
    v => 0.0,
    w => 1.2713,
    p => 0.0,
    q => 0.0,
    r => 0.0,
    ϕ => 0.0,
    θ => 0.01495, # approx 5.73 degrees
    ψ => 0.0
    )

u0 = Dict(
    uA => 0.0,
    uT => -0.1,
    uR => 0.0,
    uE_1 => 0.08,
    uE_2 => 0.08
    )

# Vector form to keep ordering
x0_vec = map(elem -> float(x0[elem]), reduced_states)
u0_vec = map(elem -> float(u0[elem]), inputs)

tspan0 = (0.0, 180)
prob0 = ODEProblem(sys0, x0, tspan0, u0, jac = true)
sol0 = solve(prob0, Tsit5())
plot(sol0, idxs = reduced_states, layout = length(reduced_states))

desired_states = Dict(u => 85, v => 0, ϕ => 0, ψ => 0, uE_1 => 0.1, uE_2 => 0.1)

hard_eq_cons =  [
    Va ~ 85
    θ - atan(w,u) ~ 0.0
]

hard_ineq_cons = [
    - uA + deg2rad(-25)    
    uA - deg2rad(25)
    -uT + deg2rad(-25)
    uT - deg2rad(10)
    -uR + deg2rad(-30)
    uR - deg2rad(30)
    -uE_1 + deg2rad(0.5) 
    uE_1 - deg2rad(10)
    -uE_2 + deg2rad(0.5) 
    uE_2 - deg2rad(10)
]

penalty_multipliers = Dict(:desired_states => 10.0, :trim_cons => 10.0, :sys_eqns => 10.0)

sol_trim, trim_states = trim(rcam_model; penalty_multipliers, desired_states, inputs, hard_eq_cons, hard_ineq_cons)

x0_trim = Dict(state => trim_states[state] for state in reduced_states)

u0_trim = Dict(input => trim_states[input] for input in inputs)

x0_trim_vec = map(elem -> x0_trim[elem], reduced_states)
u0_trim_vec = map(elem -> u0_trim[elem], inputs)

(; A, B, C, D), ssys = ModelingToolkit.linearize(rcam_model, inputs, outputs; op = merge(x0_trim, u0_trim))

linsys = named_ss(ss(A, B, C, D), x=Symbol.(states(ssys)), u=Symbol.(ModelingToolkit.inputs(ssys)), y=Symbol.(ModelingToolkit.outputs(ssys)))

# Discretization:
Ts = 0.02
disc(x) = c2d(x, Ts)
G = disc(linsys)

func_sys = JuliaSimControl.build_controlled_dynamics(rcam_model, inputs, outputs; ps)

pms = ModelingToolkit.varmap_to_vars(ModelingToolkit.defaults(rcam_model), func_sys.p)

discrete_dynamics = JuliaSimControl.rk4(func_sys, Ts)

nx = func_sys.nx
ny = func_sys.ny
nu = func_sys.nu

R1 = Matrix(1e-5*I(nx)) # Dynamics covariance
R2 = Matrix(1e-20*I(ny))  # Measurement covariance
Ru = Matrix(1e-5*I(nu))
d0 = MvNormal(float.(x0_vec), Matrix(R1)) # Initial point of the simulation x0.

ukf = JuliaSimControl.UnscentedKalmanFilter(discrete_dynamics, R1, R2, d0);

op_trim_wrapper = OperatingPoint(x0_trim_vec, u0_trim_vec, x0_trim_vec);

observer = JuliaSimControl.OperatingPointWrapper(ukf, op_trim_wrapper);

# 5. Formulate MPC problem:
umin = [
    deg2rad(-25)
    deg2rad(-25) 
    deg2rad(-30)
    deg2rad(0.5)
    deg2rad(0.5)
]

umax = [
    deg2rad(25)
    deg2rad(10)
    deg2rad(30)
    deg2rad(10)
    deg2rad(10)
]
constraints = MPCConstraints(; umin, umax)

MPC_models = LinearMPCModel(G, observer; constraints, x0 = x0_vec, op=op_trim_wrapper)

solver = OSQPSolver(
    eps_rel = 1e-5,
    eps_abs = 1e-5,
    max_iter = 500,
    check_termination = 5,
    sqp_iters = 1,
    dynamics_interval = 1,
    verbose = false,
    polish = false,
)
	
N = 100 # prediction horizon
Q1 = 10I(nx)
Q2 = 1.0*I(nu)
qs = 100
qs2 = 100000

# The reference is the trim state we want to get to:
MPCProblem = LQMPCProblem(MPC_models; Q1, Q2, qs, qs2, N, r=x0_trim_vec, solver, p = pms)

# Case 1. MPC - Without disturbance:
@time hist_lin_no_disturb = MPC.solve(MPCProblem; x0 = x0_vec, T = 300, verbose = false, noise=0, dyn_actual=discrete_dynamics, reset_observer = false)
plot(hist_lin_no_disturb, plot_title="Linear MPC", legend=:bottomright, layout=length(x0_vec) + length(u0_vec), sp=[repeat(1:9, outer=2); 10:14], title="", size=(1920, 1080))

# Case 2. MPC - With disturbance:
function disturbance_1(u, t)
    (nu, ) = size(u)
    d0 = zeros(nu)
    if t > 100 && t < 140
        d0[5] = deg2rad(-40)
    end 
    d0
end

@time hist_lin_disturbance = MPC.solve(MPCProblem; x0 = x0_vec, T = 300, verbose = false, noise=0, dyn_actual=discrete_dynamics, disturbance = disturbance_1, reset_observer = false)
plot(hist_lin_disturbance, plot_title="Linear MPC", legend=:bottomright, layout=length(x0_vec) + length(u0_vec), sp=[repeat(1:9, outer=2); 10:14], title="", size=(1920, 1080))
