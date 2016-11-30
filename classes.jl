# VARIOUS TYPES FOR CALCULATIONS

type LapStatus
    currentLap::Int32       # current lap number
    currentIt::Int32        # current iteration in current lap
    switchLap::Bool
    nextLap::Bool
    s_lapTrigger::Float32
end

# Structure of coeffConst:
# 1st dimension is the polynomial coefficient
# 2nd dimension specifies the state (1 = eY, 2 = ePsi, 3 = v or 1 = xDot, 2 = yDot, 3 = psiDot, 4 = ePsi, 5 = eY)
# 3th dimension specifies one of the two lap numbers between which are iterated

type MpcCoeff           # coefficients for trajectory approximation
    coeffCost::Array{Float32}
    coeffConst::Array{Float32}
    order::Int32
    pLength::Int32      # small values here may lead to numerical problems since the functions are only approximated in a short horizon
                        # "small" values are about 2*N, good values about 4*N
                        # numerical problems occur at the edges (s=0, when v is almost 0 and s does not change fast and at s=s_target)
    c_Vx::Array{Float32,1}
    c_Vy::Array{Float32,1}
    c_Psi::Array{Float32,1}
    MpcCoeff(coeffCost=Float32[], coeffConst=Float32[], order=4, pLength=0,c_Vx=Float32[],c_Vy=Float32[],c_Psi=Float32[]) = new(coeffCost, coeffConst, order, pLength, c_Vx, c_Vy, c_Psi)
end

type OldTrajectory      # information about previous trajectories
    oldTraj::Array{Float32}             # contains all states over all laps
    oldInput::Array{Float32}            # contains all inputs
    oldTimes::Array{Float32}            # contains times related to states and inputs
    oldCost::Array{Int32}               # contains costs of laps
    count::Array{Int32}                 # contains the counter for each lap
    prebuf::Int32
    postbuf::Int32
    idx_start::Array{Int32}             # index of the first measurement with s > 0
    idx_end::Array{Int32}               # index of the last measurement with s < s_target
    OldTrajectory(oldTraj=Float32[],oldInput=Float32[],oldTimes=Float32[],oldCost=Float32[],count=Int32[],prebuf=50,postbuf=50,idx_start=Int32[],idx_end=Int32[]) = new(oldTraj,oldInput,oldTimes,oldCost,count,prebuf,postbuf,idx_start,idx_end)
end

type MpcParams          # parameters for MPC solver
    N::Int32
    nz::Int32
    OrderCostCons::Int32
    Q::Array{Float32,1}
    Q_term::Array{Float32,1}
    R::Array{Float32,1}
    vPathFollowing::Float32
    QderivZ::Array{Float32,1}
    QderivU::Array{Float32,1}
    Q_term_cost::Float32
    delay_df::Int32
    delay_a::Int32
    MpcParams(N=0,nz=0,OrderCostCons=0,Q=Float32[],Q_term=Float32[],R=Float32[],vPathFollowing=1.0,QderivZ=Float32[],QderivU=Float32[],Q_term_cost=1.0,delay_df=0,delay_a=0) = new(N,nz,OrderCostCons,Q,Q_term,R,vPathFollowing,QderivZ,QderivU,Q_term_cost,delay_df,delay_a)
end

type PosInfo            # current position information
    s::Float32
    s_target::Float32
    PosInfo(s=0,s_target=0) = new(s,s_target)
end

type MpcSol             # MPC solution output
    a_x::Float32
    d_f::Float32
    solverStatus::Symbol
    u::Array{Float32}
    z::Array{Float32}
    cost::Array{Float32}
    MpcSol(a_x=0.0,d_f=0.0,solverStatus=Symbol(),u=Float32[],z=Float32[],cost=Float32[]) = new(a_x,d_f,solverStatus,u,z,cost)
end

type TrackCoeff         # coefficients of track
    coeffAngle::Array{Float32,1}
    coeffCurvature::Array{Float32,1}
    nPolyCurvature::Int32      # order of the interpolation polynom
    width::Float32               # lane width -> is used in cost function as soft constraints (to stay on track)
    TrackCoeff(coeffAngle=Float32[],coeffCurvature=Float32[],nPolyCurvature=4,width=1.0) = new(coeffAngle,coeffCurvature,nPolyCurvature,width)
end

type ModelParams
    l_A::Float32
    l_B::Float32
    m::Float32
    I_z::Float32
    dt::Float32
    u_lb::Array{Float32}        # lower bounds for u
    u_ub::Array{Float32}        # upper bounds
    z_lb::Array{Float32}
    z_ub::Array{Float32}
    c0::Array{Float32}
    c_f::Float32
    ModelParams(l_A=0.25,l_B=0.25,m=1.98,I_z=0.24,dt=0.1,u_lb=Float32[],u_ub=Float32[],z_lb=Float32[],z_ub=Float32[],c0=Float32[],c_f=0.0) = new(l_A,l_B,m,I_z,dt,u_lb,u_ub,z_lb,z_ub,c0,c_f)
end
