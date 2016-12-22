# type MpcModel
#     mdl::JuMP.Model

#     z0::Array{JuMP.NonlinearParameter,1}
#     coeff::Array{JuMP.NonlinearParameter,1}
#     c_Vx::Array{JuMP.NonlinearParameter,1}
#     c_Vy::Array{JuMP.NonlinearParameter,1}
#     c_Psi::Array{JuMP.NonlinearParameter,1}
#     coeffTermConst::Array{JuMP.NonlinearParameter,3}
#     coeffTermCost::Array{JuMP.NonlinearParameter,2}

#     z_Ol::Array{JuMP.Variable,2}
#     u_Ol::Array{JuMP.Variable,2}
#     ParInt::JuMP.Variable

#     laneCost::JuMP.NonlinearExpression
#     constZTerm::JuMP.NonlinearExpression
#     costZTerm::JuMP.NonlinearExpression
#     derivCost::JuMP.NonlinearExpression
#     controlCost::JuMP.NonlinearExpression
#     costZ::JuMP.NonlinearExpression

#     uPrev::Array{JuMP.NonlinearParameter,2}

#     function MpcModel(mpcParams::MpcParams,mpcCoeff::MpcCoeff,modelParams::ModelParams,trackCoeff::TrackCoeff)
#         m = new()
#         dt   = modelParams.dt
#         L_a  = modelParams.l_A
#         L_b  = modelParams.l_B
#         c0   = modelParams.c0
#         u_lb = modelParams.u_lb
#         u_ub = modelParams.u_ub
#         z_lb = modelParams.z_lb
#         z_ub = modelParams.z_ub

#         N               = mpcParams.N
#         Q               = mpcParams.Q
#         Q_term          = mpcParams.Q_term
#         R               = mpcParams.R
#         order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
#         ey_max          = trackCoeff.width/2

#         QderivZ         = mpcParams.QderivZ::Array{Float64,1}
#         QderivU         = mpcParams.QderivU::Array{Float64,1}
#         Q_term_cost     = mpcParams.Q_term_cost::Float64
#         delay_df        = mpcParams.delay_df

#         n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
        
#         # Path following mode:
#         # Create function-specific parameters
#         v_ref       = mpcParams.vPathFollowing
#         z_Ref::Array{Float64,2}
#         z_Ref       = cat(2,v_ref*ones(N+1,1),zeros(N+1,5))       # Reference trajectory: path following -> stay on line and keep constant velocity
#         u_Ref       = zeros(N,2)

#         mdl = Model(solver = IpoptSolver(print_level=0,max_cpu_time=1.0))#,check_derivatives_for_naninf="yes"))#,linear_solver="ma57",print_user_options="yes"))

#         @variable( mdl, z_Ol[1:(N+1),1:6])
#         @variable( mdl, u_Ol[1:N,1:2])
#         @variable( mdl, 0 <= ParInt <= 1)
#         @variable( mdl, eps[1:2] >= 0)                   # eps for soft lane constraints

#         z_lb_6s = ones(mpcParams.N+1,1)*[0.1 -1.0 -Inf -Inf -Inf -Inf]                      # lower bounds on states
#         z_ub_6s = ones(mpcParams.N+1,1)*[5.0  2.0  Inf  Inf  Inf  Inf]                      # upper bounds
#         u_lb_6s = ones(mpcParams.N,1) * [-3.0  -0.3]                                         # lower bounds on steering
#         u_ub_6s = ones(mpcParams.N,1) * [ 6.0   0.3]                                         # upper bounds

#         for i=1:2
#             for j=1:N
#                 setlowerbound(u_Ol[j,i], u_lb_6s[j,i])
#                 setupperbound(u_Ol[j,i], u_ub_6s[j,i])
#             end
#         end
#         for i=1:6
#             for j=1:N+1
#                 setlowerbound(z_Ol[j,i], z_lb_6s[j,i])
#                 setupperbound(z_Ol[j,i], z_ub_6s[j,i])
#             end
#         end

#         @NLparameter(mdl, z0[i=1:6] == 0)
#         @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == 0)
#         @NLparameter(mdl, c_Vx[i=1:4]  == 0)
#         @NLparameter(mdl, c_Vy[i=1:4]  == 0)
#         @NLparameter(mdl, c_Psi[i=1:3] == 0)
#         @NLparameter(mdl, coeffTermConst[i=1:order+1,j=1:2,k=1:5] == 0)
#         @NLparameter(mdl, coeffTermCost[i=1:order+1,j=1:2] == 0)
#         @NLparameter(mdl, uPrev[1:10,1:2] == 0)

#         # Conditions for first solve:
#         setvalue(z0[1],1)
#         setvalue(c_Vx[3],1)

#         @NLconstraint(mdl, [i=1:6], z_Ol[1,i] == z0[i])
#         @NLconstraint(mdl, [i=1:N+1], z_Ol[i,5] <=  ey_max + eps[1])
#         @NLconstraint(mdl, [i=1:N+1], z_Ol[i,5] >= -ey_max - eps[2])

#         @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,6]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])
#         @NLexpression(mdl, dsdt[i = 1:N], (z_Ol[i,1]*cos(z_Ol[i,4]) - z_Ol[i,2]*sin(z_Ol[i,4]))/(1-z_Ol[i,5]*c[i]))
        
#         println("Initializing model...")

#         # System dynamics
#         for i=1:N
#             if i<=delay_df
#                 @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + c_Vy[1]*z_Ol[i,2]/z_Ol[i,1] + c_Vy[2]*z_Ol[i,1]*z_Ol[i,3] + c_Vy[3]*z_Ol[i,3]/z_Ol[i,1] + c_Vy[4]*uPrev[delay_df+1-i,2]) # yDot
#                 @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + c_Psi[1]*z_Ol[i,3]/z_Ol[i,1] + c_Psi[2]*z_Ol[i,2]/z_Ol[i,1] + c_Psi[3]*uPrev[delay_df+1-i,2])                            # psiDot
#             else
#                 @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + c_Vy[1]*z_Ol[i,2]/z_Ol[i,1] + c_Vy[2]*z_Ol[i,1]*z_Ol[i,3] + c_Vy[3]*z_Ol[i,3]/z_Ol[i,1] + c_Vy[4]*u_Ol[i-delay_df,2]) # yDot
#                 @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + c_Psi[1]*z_Ol[i,3]/z_Ol[i,1] + c_Psi[2]*z_Ol[i,2]/z_Ol[i,1] + c_Psi[3]*u_Ol[i-delay_df,2])                            # psiDot
#             end
#             @NLconstraint(mdl, z_Ol[i+1,1]  == z_Ol[i,1] + c_Vx[1]*z_Ol[i,2] + c_Vx[2]*z_Ol[i,3] + c_Vx[3]*z_Ol[i,1] + c_Vx[4]*u_Ol[i,1])                                                   # xDot
#             #@NLconstraint(mdl, z_Ol[i+1,1]  == z_Ol[i,1] + dt*(u_Ol[i,1] - 0.5*z_Ol[i,1]))
#             @NLconstraint(mdl, z_Ol[i+1,4]  == z_Ol[i,4] + dt*(z_Ol[i,3]-dsdt[i]*c[i]))                                                                                 # ePsi
#             @NLconstraint(mdl, z_Ol[i+1,5]  == z_Ol[i,5] + dt*(z_Ol[i,1]*sin(z_Ol[i,4])+z_Ol[i,2]*cos(z_Ol[i,4])))                                                      # eY
#             @NLconstraint(mdl, z_Ol[i+1,6]  == z_Ol[i,6] + dt*dsdt[i]  )                                                                                                # s
#         end

#         # Cost functions

#         # Derivative cost
#         # ---------------------------------
#         @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:6} +
#                                           QderivU[1]*((uPrev[1,1]-u_Ol[1,1])^2+sum{(u_Ol[i,1]-u_Ol[i+1,1])^2,i=1:N-1})+
#                                           QderivU[2]*((uPrev[1,2]-u_Ol[1,2])^2+sum{(u_Ol[i,2]-u_Ol[i+1,2])^2,i=1:N-delay_df-1}))

#         # Lane cost
#         # ---------------------------------
#         @NLexpression(mdl, laneCost, sum{100*eps[i] + 100*eps[i]^2,i=1:2})

#         # Control Input cost
#         # ---------------------------------
#         @NLexpression(mdl, controlCost, R[1]*sum{(u_Ol[i,1])^2,i=1:N}+
#                                         R[2]*sum{(u_Ol[i,2])^2,i=1:N-delay_df})

#         # Terminal constraints (soft), starting from 2nd lap
#         # ---------------------------------
#         @NLexpression(mdl, constZTerm, sum{Q_term[j]*(ParInt*(sum{coeffTermConst[i,1,j]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermConst[order+1,1,j])+
#                                             (1-ParInt)*(sum{coeffTermConst[i,2,j]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermConst[order+1,2,j])-z_Ol[N+1,j])^2,j=1:5})
        
#         # Terminal cost
#         # ---------------------------------
#         # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
#         @NLexpression(mdl, costZTerm, Q_term_cost*(ParInt*(sum{coeffTermCost[i,1]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermCost[order+1,1])+
#                                       (1-ParInt)*(sum{coeffTermCost[i,2]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermCost[order+1,2])))
        
#         # State cost (only for path following mode)
#         # ---------------------------------
#         #@NLexpression(mdl, costZ, 0.5*sum{Q[i]*sum{(z_Ol[j,i]-z_Ref[j,i])^2,j=1:N+1},i=1:6})    # Follow trajectory
#         @NLexpression(mdl, costZ, 0.0)
#         #@NLexpression(mdl, costZ, 0.5*sum{(Q[1]*(sqrt(z_Ol[j,1]^2+z_Ol[j,2]^2)-0.8)^2 + Q[4]*z_Ol[j,4]^2 + Q[5]*z_Ol[j,5]^2),j=2:N+1})

#         # Solve model once
#         @NLobjective(mdl, Min, derivCost + constZTerm + costZTerm + laneCost)
#         #@NLobjective(mdl, Min, derivCost + costZ + laneCost)
#         sol_stat=solve(mdl)
#         println("Finished solve 1: $sol_stat")
#         sol_stat=solve(mdl)
#         println("Finished solve 2: $sol_stat")
#         m.mdl = mdl
#         m.z0 = z0
#         m.coeff = coeff
#         m.z_Ol = z_Ol
#         m.u_Ol = u_Ol
#         m.c_Vx = c_Vx
#         m.c_Vy = c_Vy
#         m.c_Psi = c_Psi
#         m.ParInt = ParInt
#         m.uPrev = uPrev

#         m.coeffTermCost = coeffTermCost
#         m.coeffTermConst = coeffTermConst

#         m.derivCost = derivCost
#         m.controlCost = controlCost
#         m.laneCost = laneCost
#         m.constZTerm = constZTerm
#         m.costZTerm  = costZTerm
#         m.costZ  = costZ
#         return m
#     end
# end


type MpcModel
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    eps::Array{JuMP.Variable,1}

    coeffTermConst::Array{JuMP.NonlinearParameter,3}
    coeffTermCost::Array{JuMP.NonlinearParameter,2}

    ParInt::JuMP.Variable

    laneCost::JuMP.NonlinearExpression
    constZTerm::JuMP.NonlinearExpression
    costZTerm::JuMP.NonlinearExpression
    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression
    costZ::JuMP.NonlinearExpression

    uPrev::Array{JuMP.NonlinearParameter,2}

    function MpcModel(mpcParams::MpcParams,mpcCoeff::MpcCoeff,modelParams::ModelParams,trackCoeff::TrackCoeff)
        println("Starting creation of the LMPC model")
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

        Q_term          = mpcParams.Q_term
        order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
        ey_max          = trackCoeff.width/2

        Q_term_cost     = mpcParams.Q_term_cost::Float64


        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=0,max_cpu_time=1.0))#,linear_solver="ma57",print_user_options="yes"))

        # Create variables (these are going to be optimized)
        @variable( mdl, z_Ol[1:(N+1),1:4], start = 0)          # z = s, ey, epsi, v
        @variable( mdl, u_Ol[1:N,1:2], start = 0)
        @variable( mdl, 0 <= ParInt <= 1)
        @variable( mdl, eps[1:N+1] >= 0)                   # eps for soft lane constraints

        # Set bounds
        z_lb_4s = ones(mpcParams.N+1,1)*[-Inf -Inf -Inf -0.5]                  # lower bounds on states
        z_ub_4s = ones(mpcParams.N+1,1)*[ Inf  Inf  Inf  1.5]                  # upper bounds
        u_lb_4s = ones(mpcParams.N,1) * [-0.2  -0.3]                           # lower bounds on steering
        u_ub_4s = ones(mpcParams.N,1) * [2.0   0.3]                            # upper bounds

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb_4s[j,i])
                setupperbound(u_Ol[j,i], u_ub_4s[j,i])
            end
        end
        for i=1:4
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb_4s[j,i])
                setupperbound(z_Ol[j,i], z_ub_4s[j,i])
            end
        end

        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,2] <=  ey_max + eps[i])
        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,2] >= -ey_max - eps[i])

        @NLparameter(mdl, z0[i=1:4] == 0)
        @NLparameter(mdl, uPrev[1:10,1:2] == 0)

        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == 0)
        @NLparameter(mdl, coeffTermConst[i=1:order+1,j=1:2,k=1:5] == 0)
        @NLparameter(mdl, coeffTermCost[i=1:order+1,j=1:2] == 0)

        @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])

        # System dynamics
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])         # initial condition
        for i=1:N
            @NLexpression(mdl, bta[i],  atan( L_a / (L_a + L_b) * tan( u_Ol[i,2] ) ) )
            @NLexpression(mdl, dsdt[i], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))
            @NLconstraint(mdl, z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )                                                # s
            @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                        # ey
            @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_a*sin(bta[i])-dsdt[i]*c[i])  )               # epsi
            @NLconstraint(mdl, z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1] - 0.5*z_Ol[i,4]))                              # v
        end

        # Cost definitions
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:4} +
                                            QderivU[1]*((uPrev[1,1]-u_Ol[1,1])^2+sum{(u_Ol[i,1]-u_Ol[i+1,1])^2,i=1:N-1})+
                                            QderivU[2]*((uPrev[1,2]-u_Ol[1,2])^2+sum{(u_Ol[i,2]-u_Ol[i+1,2])^2,i=1:N-1}))
        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*R[1]*sum{(u_Ol[i,1])^2,i=1:N}+
                                        0.5*R[2]*sum{(u_Ol[i,2])^2,i=1:N})
        
        # State cost
        # ---------------------------------
        @NLexpression(mdl, costZ, 1)    # Follow trajectory

        # Lane cost
        # ---------------------------------
        @NLexpression(mdl, laneCost, sum{10*eps[i] + 100*eps[i]^2,i=2:N+1})

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, R[1]*sum{(u_Ol[i,1])^2,i=1:N}+
                                        R[2]*sum{(u_Ol[i,2])^2,i=1:N})

        # Terminal constraints (soft), starting from 2nd lap
        # ---------------------------------
        @NLexpression(mdl, constZTerm, sum{Q_term[j]*(ParInt*(sum{coeffTermConst[i,1,j]*z_Ol[N+1,1]^(order+1-i),i=1:order}+coeffTermConst[order+1,1,j])+
                                            (1-ParInt)*(sum{coeffTermConst[i,2,j]*z_Ol[N+1,1]^(order+1-i),i=1:order}+coeffTermConst[order+1,2,j])-z_Ol[N+1,j+1])^2,j=1:3})
        
        # Terminal cost
        # ---------------------------------
        # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
        @NLexpression(mdl, costZTerm, Q_term_cost*(ParInt*(sum{coeffTermCost[i,1]*z_Ol[N+1,1]^(order+1-i),i=1:order}+coeffTermCost[order+1,1])+
                                      (1-ParInt)*(sum{coeffTermCost[i,2]*z_Ol[N+1,1]^(order+1-i),i=1:order}+coeffTermCost[order+1,2])))

        # Objective function
        @NLobjective(mdl, Min, costZ + derivCost + controlCost + constZTerm + costZTerm + laneCost)

        m.mdl = mdl
        m.z0 = z0
        m.coeff = coeff
        m.z_Ol = z_Ol
        m.u_Ol = u_Ol
        m.uPrev = uPrev
        m.derivCost = derivCost
        m.costZ = costZ
        m.controlCost = controlCost

        m.ParInt = ParInt

        m.coeffTermCost = coeffTermCost
        m.coeffTermConst = coeffTermConst

        m.laneCost = laneCost
        m.constZTerm = constZTerm
        m.costZTerm  = costZTerm
        m.eps = eps

        return m
    end
end


type MpcModel_pF
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}

    derivCost::JuMP.NonlinearExpression
    costZ::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression

    uPrev::Array{JuMP.NonlinearParameter,2}

    function MpcModel_pF(mpcParams::MpcParams,modelParams::ModelParams,trackCoeff::TrackCoeff)
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
        z_Ref       = cat(2,zeros(N+1,3),v_ref*ones(N+1,1))       # Reference trajectory: path following -> stay on line and keep constant velocity
        u_Ref       = zeros(N,2)

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=0,max_cpu_time=1.0))#,linear_solver="ma57",print_user_options="yes"))

        # Create variables (these are going to be optimized)
        @variable( mdl, z_Ol[1:(N+1),1:4], start = 0)          # z = s, ey, epsi, v
        @variable( mdl, u_Ol[1:N,1:2], start = 0)

        # Set bounds
        z_lb_4s = ones(mpcParams.N+1,1)*[-Inf -Inf -Inf -0.5]                  # lower bounds on states
        z_ub_4s = ones(mpcParams.N+1,1)*[ Inf  Inf  Inf  1.5]                  # upper bounds
        u_lb_4s = ones(mpcParams.N,1) * [-0.2  -0.3]                           # lower bounds on steering
        u_ub_4s = ones(mpcParams.N,1) * [2.0   0.3]                            # upper bounds

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb_4s[j,i])
                setupperbound(u_Ol[j,i], u_ub_4s[j,i])
            end
        end
        for i=1:4
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb_4s[j,i])
                setupperbound(z_Ol[j,i], z_ub_4s[j,i])
            end
        end

        @NLparameter(mdl, z0[i=1:4] == 0)
        @NLparameter(mdl, uPrev[1:10,1:2] == 0)

        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == 0)

        @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])

        # System dynamics
        setvalue(z0[4],v_ref)
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])         # initial condition
        for i=1:N
            @NLexpression(mdl, bta[i],  atan( L_a / (L_a + L_b) * tan( u_Ol[i,2] ) ) )
            @NLexpression(mdl, dsdt[i], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))
            @NLconstraint(mdl, z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )                                                # s
            @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                        # ey
            @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_a*sin(bta[i])-dsdt[i]*c[i])  )               # epsi
            @NLconstraint(mdl, z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1] - 0.5*z_Ol[i,4]))#modelParams.c_f*abs(z_Ol[i,4]) * z_Ol[i,4])) # v
        end

        # Cost definitions
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:4} +
                                            QderivU[1]*((uPrev[1,1]-u_Ol[1,1])^2+sum{(u_Ol[i,1]-u_Ol[i+1,1])^2,i=1:N-1})+
                                            QderivU[2]*((uPrev[1,2]-u_Ol[1,2])^2+sum{(u_Ol[i,2]-u_Ol[i+1,2])^2,i=1:N-1}))
        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*R[1]*sum{(u_Ol[i,1])^2,i=1:N}+
                                        0.5*R[2]*sum{(u_Ol[i,2])^2,i=1:N})
        # State cost
        # ---------------------------------
        @NLexpression(mdl, costZ, 0.5*sum{Q[i]*sum{(z_Ol[j,i]-z_Ref[j,i])^2,j=2:N+1},i=1:4})    # Follow trajectory

        # Objective function
        @NLobjective(mdl, Min, costZ + derivCost + controlCost)

        m.mdl = mdl
        m.z0 = z0
        m.coeff = coeff
        m.z_Ol = z_Ol
        m.u_Ol = u_Ol
        m.uPrev = uPrev
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
