type MpcModel
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    c_Vx::Array{JuMP.NonlinearParameter,1}
    c_Vy::Array{JuMP.NonlinearParameter,1}
    c_Psi::Array{JuMP.NonlinearParameter,1}
    coeffTermConst::Array{JuMP.NonlinearParameter,3}
    coeffTermCost::Array{JuMP.NonlinearParameter,2}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    ParInt::JuMP.Variable

    dsdt::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    laneCost::JuMP.NonlinearExpression
    constZTerm::JuMP.NonlinearExpression
    costZTerm::JuMP.NonlinearExpression
    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression

    uCurr::Array{JuMP.NonlinearParameter,1}

    function MpcModel(mpcParams::MpcParams,mpcCoeff::MpcCoeff,modelParams::ModelParams,trackCoeff::TrackCoeff,z_Init::Array{Float64,1})
        m = new()
        dt   = modelParams.dt
        L_a  = modelParams.l_A
        L_b  = modelParams.l_B
        c0   = modelParams.c0
        u_lb = modelParams.u_lb
        u_ub = modelParams.u_ub
        z_lb = modelParams.z_lb
        z_ub = modelParams.z_ub

        N               = mpcParams.N
        Q_term          = mpcParams.Q_term
        R               = mpcParams.R
        order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
        ey_max          = trackCoeff.width/2

        QderivZ         = mpcParams.QderivZ::Array{Float64,1}
        QderivU         = mpcParams.QderivU::Array{Float64,1}
        Q_term_cost     = mpcParams.Q_term_cost::Float64

        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
        
        mdl = Model(solver = IpoptSolver(print_level=3,max_cpu_time=0.5))#,check_derivatives_for_naninf="yes"))#,linear_solver="ma57",print_user_options="yes"))

        @variable( mdl, z_Ol[1:(N+1),1:6], start = 0)
        @variable( mdl, u_Ol[1:N,1:2],     start = 0)
        @variable( mdl, 0 <= ParInt <= 1,  start = 0)

        z_lb_6s = ones(mpcParams.N+1,1)*[0.1 -Inf -Inf -Inf -Inf -Inf]                      # lower bounds on states
        z_ub_6s = ones(mpcParams.N+1,1)*[3.0  Inf  Inf  Inf  Inf  Inf]                      # upper bounds
        u_lb_6s = ones(mpcParams.N,1) * [-1.0  -pi/6]                                       # lower bounds on steering
        u_ub_6s = ones(mpcParams.N,1) * [3.0    pi/6]                                       # upper bounds

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb_6s[j,i])
                setupperbound(u_Ol[j,i], u_ub_6s[j,i])
            end
        end
        for i=1:6
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb_6s[j,i])
                setupperbound(z_Ol[j,i], z_ub_6s[j,i])
            end
        end

        @NLparameter(mdl, z0[i=1:6] == z_Init[i])
        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == 0);
        @NLparameter(mdl, c_Vx[i=1:4]  == 0)
        @NLparameter(mdl, c_Vy[i=1:4]  == 0)
        @NLparameter(mdl, c_Psi[i=1:3] == 0)
        @NLparameter(mdl, coeffTermConst[i=1:order+1,j=1:2,k=1:5] == 0)
        @NLparameter(mdl, coeffTermCost[i=1:order+1,j=1:2] == 0)
        @NLparameter(mdl, uCurr[1:2] == 0)

        # Conditions for first solve:
        setvalue(c_Vx[4],1)
        setvalue(c_Vy[4],1)
        setvalue(c_Psi[3],1)

        @NLconstraint(mdl, [i=1:6], z_Ol[1,i] == z0[i])

        @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,6]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])
        @NLexpression(mdl, dsdt[i = 1:N], (z_Ol[i,1]*cos(z_Ol[i,4]) - z_Ol[i,2]*sin(z_Ol[i,4]))/(1-z_Ol[i,5]*c[i]))
        
        println("Initializing model...")

        # System dynamics
        for i=1:N
            @NLconstraint(mdl, z_Ol[i+1,1]  == z_Ol[i,1] + c_Vx[1]*z_Ol[i,2] + c_Vx[2]*z_Ol[i,3] + c_Vx[3]*z_Ol[i,1]^2 + c_Vx[4]*u_Ol[i,1])                             # xDot
            @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + c_Vy[1]*z_Ol[i,2]/z_Ol[i,1] + c_Vy[2]*z_Ol[i,1]*z_Ol[i,3] + c_Vy[3]*z_Ol[i,3]/z_Ol[i,1] + c_Vy[4]*u_Ol[i,2]) # yDot
            @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + c_Psi[1]*z_Ol[i,3]/z_Ol[i,1] + c_Psi[2]*z_Ol[i,2]/z_Ol[i,1] + c_Psi[3]*u_Ol[i,2])                            # psiDot
            @NLconstraint(mdl, z_Ol[i+1,4]  == z_Ol[i,4] + dt*(z_Ol[i,3]-dsdt[i]*c[i]))                                                                                 # ePsi
            @NLconstraint(mdl, z_Ol[i+1,5]  == z_Ol[i,5] + dt*(z_Ol[i,1]*sin(z_Ol[i,4])+z_Ol[i,2]*cos(z_Ol[i,4])))                                                      # eY
            @NLconstraint(mdl, z_Ol[i+1,6]  == z_Ol[i,6] + dt*dsdt[i]  )
        end

        # Cost functions

        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:6} +
                                          sum{QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum{(u_Ol[i,j]-u_Ol[i+1,j])^2,i=1:N-1}),j=1:2})

        # Lane cost
        # ---------------------------------
        @NLexpression(mdl, laneCost, 10*sum{z_Ol[i,5]^2*((0.5+0.5*tanh(10*(z_Ol[i,5]-ey_max))) + (0.5-0.5*tanh(10*(z_Ol[i,5]+ey_max)))),i=1:N+1})

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*sum{R[j]*sum{(u_Ol[i,j])^2,i=1:N},j=1:2})

        # Terminal constraints (soft), starting from 2nd lap
        # ---------------------------------

        #println("Q_Term = $Q_term")
        #println("coeffTermConst = $(getvalue(coeffTermConst))")
        #println("coeffTermCost = $(getvalue(coeffTermCost))")
        #println("order = $order")
        @NLexpression(mdl, constZTerm, sum{Q_term[j]*(ParInt*sum{coeffTermConst[i,1,j]*z_Ol[N+1,6]^(order+1-i),i=1:order+1}+
                                            (1-ParInt)*sum{coeffTermConst[i,2,j]*z_Ol[N+1,6]^(order+1-i),i=1:order+1}-z_Ol[N+1,j])^2,j=1:5})
        
        # Terminal cost
        # ---------------------------------
        # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
        
        @NLexpression(mdl, costZTerm, Q_term_cost*(ParInt*sum{coeffTermCost[i,1]*z_Ol[N+1,6]^(order+1-i),i=1:order+1}+
                                      (1-ParInt)*sum{coeffTermCost[i,2]*z_Ol[N+1,6]^(order+1-i),i=1:order+1}))
        
        @NLobjective(mdl,Min,sum{(z_Ol[i,1]-0.5)^2,i=1:N+1})

        # solve model once
        println("solving...")
        solve(mdl)
        @NLobjective(mdl, Min, derivCost + controlCost + laneCost + constZTerm + costZTerm)
        solve(mdl)
        println("finished")
        m.mdl = mdl
        m.z0 = z0
        m.coeff = coeff
        m.z_Ol = z_Ol
        m.u_Ol = u_Ol
        m.dsdt = dsdt
        m.c = c
        m.c_Vx = c_Vx
        m.c_Vy = c_Vy
        m.c_Psi = c_Psi
        m.ParInt = ParInt
        m.uCurr = uCurr

        m.coeffTermCost = coeffTermCost
        m.coeffTermConst = coeffTermConst

        m.derivCost = derivCost
        m.controlCost = controlCost
        m.laneCost = laneCost
        m.constZTerm = constZTerm
        m.costZTerm  = costZTerm
        return m
    end
end

type MpcModel_pF
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}

    dsdt::Array{JuMP.NonlinearExpression,1}
    bta::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    derivCost::JuMP.NonlinearExpression
    costZ::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression

    uCurr::Array{JuMP.NonlinearParameter,1}

    function MpcModel_pF(mpcParams::MpcParams,modelParams::ModelParams,trackCoeff::TrackCoeff,z_Init::Array{Float64,1})
        println("Starting creation of pf model")
        m = new()
        dt          = modelParams.dt
        L_a         = modelParams.l_A
        L_b         = modelParams.l_B
        c0          = modelParams.c0
        u_lb        = modelParams.u_lb
        u_ub        = modelParams.u_ub
        z_lb        = modelParams.z_lb
        z_ub        = modelParams.z_ub

        N           = mpcParams.N
        Q           = mpcParams.Q
        R           = mpcParams.R
        QderivZ     = mpcParams.QderivZ::Array{Float64,1}
        QderivU     = mpcParams.QderivU::Array{Float64,1}

        v_ref       = mpcParams.vPathFollowing

        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation

        # Create function-specific parameters
        z_Ref::Array{Float64,2}
        z_Ref           = cat(2,zeros(N+1,3),v_ref*ones(N+1,1))       # Reference trajectory: path following -> stay on line and keep constant velocity
        u_Ref           = zeros(N,2)

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=3,max_cpu_time=5.0))#,linear_solver="ma57",print_user_options="yes"))

        # Create variables (these are going to be optimized)
        @variable( mdl, z_Ol[1:(N+1),1:4])          # z = s, ey, epsi, v
        @variable( mdl, u_Ol[1:N,1:2])

        # Set bounds (hard constraints)
        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], modelParams.u_lb[j,i])
                setupperbound(u_Ol[j,i], modelParams.u_ub[j,i])
            end
        end
        for i=1:4
            for j=1:N+1
                setlowerbound(z_Ol[j,i], modelParams.z_lb[j,i])
                setupperbound(z_Ol[j,i], modelParams.z_ub[j,i])
            end
        end

        @NLparameter(mdl, z0[i=1:4] == z_Init[i])
        @NLparameter(mdl, uCurr[i=1:2] == 0)

        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == 0)

        @NLexpression(mdl, c[i = 1:N],    sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])
        @NLexpression(mdl, bta[i = 1:N],  atan( L_a / (L_a + L_b) * ( u_Ol[i,2] ) ) )
        @NLexpression(mdl, dsdt[i = 1:N], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))

        # System dynamics
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])         # initial condition
        for i=1:N
            @NLconstraint(mdl, z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )                                                # s
            @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                        # ey
            @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_a*sin(bta[i])-dsdt[i]*c[i])  )               # epsi
            @NLconstraint(mdl, z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1] - modelParams.c_f*abs(z_Ol[i,4]) * z_Ol[i,4])) # v
        end

        # Cost definitions
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:4} +
                                            sum{QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum{(u_Ol[i,j]-u_Ol[i+1,j])^2,i=1:N-1}),j=1:2})

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*sum{R[j]*sum{(u_Ol[i,j])^2,i=1:N},j=1:2})

        # State cost
        # ---------------------------------
        @NLexpression(mdl, costZ, 0.5*sum{Q[i]*sum{(z_Ol[j,i]-z_Ref[j,i])^2,j=2:N+1},i=1:4})    # Follow trajectory

        # Objective function
        @NLobjective(mdl, Min, costZ)# + derivCost + controlCost)

        # create first artificial solution (for warm start)
        setvalue(z_Ol[1,:],z_Init')
        for i=2:N+1
            setvalue(z_Ol[i,:],[0 0 0 v_ref])
        end
        # First solve
        #JuMP.build(mdl)
        #println("Initialized path following controller. Starting first solve...")
        #println("z0 = $(getvalue(z0))")
        #println("Initial guess z: $(getvalue(z_Ol))")
        #println(mdl)
        sol_status = solve(mdl)
        #println("Solved with status $sol_status")
        #println("Solution_z = $(getvalue(z_Ol))")
        #println("Solution_u = $(getvalue(u_Ol))")
        m.mdl = mdl
        m.z0 = z0
        m.coeff = coeff
        m.z_Ol = z_Ol
        m.u_Ol = u_Ol
        m.dsdt = dsdt
        m.bta = bta
        m.c = c
        m.uCurr = uCurr
        m.derivCost = derivCost
        m.costZ = costZ
        m.controlCost = controlCost
        return m
    end
end

# why does the solution change when I start the simulation twice, totally independently? Is not everything initialized equally?
# why do I get a failed restoration phase even though there are *almost* no constraints and the solution should be super clear?

# does it actually make sense to use the first (initial) state as a @variable and combine it with a constraint to fix it to z0 ?
# or better use only a parameter for z0? -> Found out that solutions are different if I set constraints for first state or not!?