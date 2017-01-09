# Variable definitions
# mdl.z_Ol[i,j] = z_OpenLoop, open loop prediction of the state, i = state, j = step

# States in path following mode:
# i = 1 -> s
# i = 2 -> ey
# i = 3 -> epsi
# i = 4 -> v

# States in LMPC and system ID mode:
# i = 1 -> xDot
# i = 2 -> yDot
# i = 3 -> psiDot
# i = 4 -> ePsi
# i = 5 -> eY
# i = 6 -> s

function solveMpcProblem_LMPC(mdl::MpcModel,mpcSol::MpcSol,mpcCoeff::MpcCoeff,mpcParams::MpcParams,trackCoeff::TrackCoeff,lapStatus::LapStatus,posInfo::PosInfo,modelParams::ModelParams,zCurr::Array{Float64},uPrev::Array{Float64})

    # Load Parameters
    sol_status::Symbol
    sol_u::Array{Float64,2}
    sol_z::Array{Float64,2}

    # Update current initial condition, curvature and System ID coefficients
    setvalue(mdl.z0,zCurr)
    setvalue(mdl.uPrev,uPrev)
    setvalue(mdl.coeff,trackCoeff.coeffCurvature)       # Track curvature
    setvalue(mdl.coeffTermCost,mpcCoeff.coeffCost)      # Terminal cost
    setvalue(mdl.coeffTermConst,mpcCoeff.coeffConst)    # Terminal constraints

    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)
    k = 1
    while sol_status != :Optimal && k<=10
        println("Not solved optimally. Solving again.")
        sol_status  = solve(mdl.mdl)
        k += 1
    end

    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)

    # export data
    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    #mpcSol.cost = zeros(6)
    mpcSol.cost = [getvalue(mdl.costZ),getvalue(mdl.costZTerm),getvalue(mdl.constZTerm),getvalue(mdl.derivCost),getvalue(mdl.controlCost),getvalue(mdl.laneCost)]
    #mpcSol.cost = [getvalue(mdl.costZ),0,0,getvalue(mdl.derivCost),0,0]

    # Print information
    # println("--------------- MPC START -----------------------------------------------")
    # println("z0             = ", zCurr')
    # println("u              = ", sol_u[1,:])
    # println("ParInt         = ", getvalue(mdl.ParInt))
    # println("Solved, status = ", sol_status)
    # println("Predict. to s  = ", sol_z[end,1])
    # println("costZ          = ", mpcSol.cost[1])
    # println("termCost       = ", mpcSol.cost[2])
    # println("termConst      = ", mpcSol.cost[3])
    # println("derivCost      = ", mpcSol.cost[4])
    # println("controlCost    = ", mpcSol.cost[5])
    # println("laneCost       = ", mpcSol.cost[6])
    # println("eps            = ", getvalue(mdl.eps))
    # println("--------------- MPC END ------------------------------------------------")
    nothing
end

function solveMpcProblem_pathFollow(mdl::MpcModel_pF,mpcSol::MpcSol,mpcParams::MpcParams,trackCoeff::TrackCoeff,posInfo::PosInfo,modelParams::ModelParams,zCurr::Array{Float64},uPrev::Array{Float64})

    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}

    sol_status::Symbol
    sol_u::Array{Float64,2}
    sol_z::Array{Float64,2}

    # Update current initial condition, curvature and previous input
    setvalue(mdl.z0,zCurr)
    setvalue(mdl.uPrev,uPrev)
    setvalue(mdl.coeff,coeffCurvature)

    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)
    k = 1
    while sol_status != :Optimal && k<=10
        println("Not solved optimally. Solving again.")
        sol_status  = solve(mdl.mdl)
        k += 1
    end
    sol_u       = getvalue(mdl.u_Ol)
    sol_z       = getvalue(mdl.z_Ol)

    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    #mpcSol.cost = zeros(6)
    mpcSol.cost = [getvalue(mdl.costZ),0,0,getvalue(mdl.derivCost),getvalue(mdl.controlCost),0]

    # Print information
    # println("--------------- MPC PF START -----------------------------------------------")
    # println("z0             = ", zCurr')
    # println("Solved, status = ", sol_status)
    # println("Predict. to s  = ", sol_z[end,1])
    # println("--------------- MPC PF END ------------------------------------------------")
    nothing
end

# Cost structure: costZ, costZTerm, constZTerm, derivCost, controlCost, laneCost
