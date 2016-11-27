# This function evaluates and returns the coefficients for constraints and cost which are used later in the MPC formulation
# Inputs are:
# oldTraj   -> contains information about previous trajectories and Inputs
# mpcCoeff  -> contains information about previous coefficient results
# posInfo   -> contains information about the car's current position along the track
# mpcParams -> contains information about the MPC formulation (e.g. Q, R)
# stateIn   -> current state of the car
# inputIn   -> last input to the system

# structure of oldTrajectory: 1st dimension = state number, 2nd dimension = step number (time equiv.), 3rd dimennsion = lap number

# z[1] = xDot
# z[2] = yDot
# z[3] = psiDot
# z[4] = ePsi
# z[5] = eY
# z[6] = s

function coeffConstraintCost(oldTraj::OldTrajectory, mpcCoeff::MpcCoeff, posInfo::PosInfo, mpcParams::MpcParams,lapStatus::LapStatus)
    # this computes the coefficients for the cost and constraints

    # Outputs: 
    # coeffConst
    # coeffCost

    # Read Inputs
    s               = posInfo.s
    s_target        = posInfo.s_target

    # Parameters
    N               = mpcParams.N
    nz              = mpcParams.nz
    R               = mpcParams.R
    Order           = mpcCoeff.order                # interpolation order for cost and constraints
    pLength         = mpcCoeff.pLength              # interpolation length for polynomials

    n_prev          = 20                            # number of points before current s for System ID
    n_ahead         = 60                            # number of points ahead of current s for System ID

    coeffCost       = zeros(Order+1,2)              # polynomial coefficients for cost
    coeffConst      = zeros(Order+1,2,5)            # nz-1 beacuse no coeff for s

    n_laps_sysID    = 2                             # number of previous laps that are used for sysID
    #n_laps_sysID    = lapStatus.currentLap-2
    #selected_laps = [lapStatus.currentLap-1,lapStatus.currentLap-2]
    selected_laps = collect(lapStatus.currentLap-1:-1:lapStatus.currentLap-n_laps_sysID)
    # Select the old data
    oldxDot         = oldTraj.oldTraj[:,1,selected_laps]::Array{Float64,3}
    oldyDot         = oldTraj.oldTraj[:,2,selected_laps]::Array{Float64,3}
    oldpsiDot       = oldTraj.oldTraj[:,3,selected_laps]::Array{Float64,3}
    oldePsi         = oldTraj.oldTraj[:,4,selected_laps]::Array{Float64,3}
    oldeY           = oldTraj.oldTraj[:,5,selected_laps]::Array{Float64,3}
    oldS            = oldTraj.oldTraj[:,6,selected_laps]::Array{Float64,3}
    oldt            = oldTraj.oldTimes[:,selected_laps]::Array{Float64,2}
    olda            = oldTraj.oldInput[:,1,selected_laps]::Array{Float64,3}
    olddF           = oldTraj.oldInput[:,2,selected_laps]::Array{Float64,3}

    N_points        = size(oldTraj.oldTraj,1)     # second dimension = length (=buffersize)

    s_total::Float64        # initialize
    DistS::Array{Float64}   # initialize
    idx_s::Array{Int64}     # initialize
    idx_s_target        = 0
    dist_to_s_target    = 0
    qLength             = 0
    vec_range::Tuple{UnitRange{Int64},UnitRange{Int64}}
    bS_Vector::Array{Float64}
    s_forinterpy::Array{Float64}

    # Compute the total s (current position along track)
    s_total = s % s_target

    # Compute the index
    DistS = ( s_total - oldS ).^2

    idx_s = findmin(DistS,1)[2]              # contains both indices for the closest distances for both oldS !!

    vec_range = (idx_s[1]:idx_s[1]+pLength,idx_s[2]:idx_s[2]+pLength)

    # Create the vectors used for the interpolation
    # bS_vector contains the s-values for later interpolation
    bS_Vector       = zeros(pLength+1,2)
    for i=1:pLength+1
        bS_Vector[i,1] = oldS[vec_range[1][i]]
        bS_Vector[i,2] = oldS[vec_range[2][i]]
    end
    if norm(bS_Vector[1,1]-s_total) > 0.3 || norm(bS_Vector[1,2]-s_total) > 0.3
        warn("Couldn't find a close point to current s.")
    end

    # println("************************************** COEFFICIENTS **************************************")
    # println("idx_s[1]  = $(idx_s[1]), idx_s[2] = $(idx_s[2])")
    # println("s_total   = $s_total")
    # println("bS_Vector[1,:] = $(bS_Vector[1,:]')")
    # These matrices (above) contain two vectors each (for both old trajectories), stored in the 3rd dimension
    
    # The states are parametrized with resprect to the curvilinear abscissa,
    # so we select the point used for the interpolation. Need to subtract an
    # offset to be coherent with the MPC formulation
    s_forinterpy   = bS_Vector
    if s_total < 0
        s_forinterpy += s_target
    end
    # println("s_forinterpy[:,:]' = $(s_forinterpy[:,:]')")
    # Create the Matrices for the interpolation
    MatrixInterp = zeros(pLength+1,Order+1,2)

    for k = 0:Order
        MatrixInterp[:,Order+1-k,:] = s_forinterpy[:,:].^k
    end
    
    # Compute the coefficients
    coeffConst = zeros(Order+1,2,5)
    for i=1:2
        coeffConst[:,i,1]    = MatrixInterp[:,:,i]\oldxDot[vec_range[i]]
        coeffConst[:,i,2]    = MatrixInterp[:,:,i]\oldyDot[vec_range[i]]
        coeffConst[:,i,3]    = MatrixInterp[:,:,i]\oldpsiDot[vec_range[i]]
        coeffConst[:,i,4]    = MatrixInterp[:,:,i]\oldePsi[vec_range[i]]
        coeffConst[:,i,5]    = MatrixInterp[:,:,i]\oldeY[vec_range[i]]
    end

    # xDotInterp = [s_forinterpy[:,1].^5 s_forinterpy[:,1].^4 s_forinterpy[:,1].^3 s_forinterpy[:,1].^2 s_forinterpy[:,1].^1 s_forinterpy[:,1].^0]*coeffConst[:,1,1]
    # plot(s_forinterpy,oldxDot[vec_range[1]],"*",s_forinterpy,xDotInterp)
    # title("xDotInterp")
    # grid("on")
    # readline()

    # Finished with calculating the constraint coefficients
    
    # Now compute the final cost coefficients

    # The Q-function contains for every point in the sampled safe set the minimum cost-to-go-value
    # These values are calculated for both old trajectories
    # The vector bQfunction_Vector contains the cost at each point in the interpolated area to reach the finish line
    # From this vector, polynomial coefficients coeffCost are calculated to approximate this cost
    for i=1:2   
            dist_to_s_target  = oldTraj.oldCost[i] - (idx_s[i]-N_points*(i-1))  # number of iterations from idx_s to s_target
            bQfunction_Vector = collect(linspace(dist_to_s_target,dist_to_s_target-pLength,pLength+1))    # build a vector that starts at the distance and
                                                                                                    # decreases in equal steps
            coeffCost[:,i]    = MatrixInterp[:,:,i]\bQfunction_Vector           # interpolate this vector with the given s
    end
    # println("coeffCost: $coeffCost")
    # println("coeffConst: $coeffConst")
    # println("Finished coefficients, starting system ID")
    # --------------- SYSTEM IDENTIFICATION --------------- #
    # ----------------------------------------------------- #

    vec_range_ID    = ()
    for i=1:n_laps_sysID
        vec_range_ID    = tuple(vec_range_ID...,idx_s[i]-n_prev:idx_s[i]+n_ahead)     # related index range
    end

    cC = oldTraj.count[lapStatus.currentLap]-1      # current count
    cL = lapStatus.currentLap                       # current lap
    vec_range_ID2   = cC-n_prev:cC-1                # index range for current lap
    
    # TODO:*******************
    # CHECK IF DIFFERENCES ARE ALIGNED PROPERLY
    # ADD 2-step-delay to sysID
    # ************************
    sz_last = n_prev+n_ahead+1      # number of sysID points from last laps (each lap)
    sz_curr = n_prev                # number of sysID points from current lap
    n_sysID = n_laps_sysID*sz_last + sz_curr    # total number of sysID points

    y_psi = zeros(n_sysID)
    A_psi = zeros(n_sysID,3)
    y_xDot = zeros(n_sysID)
    A_xDot = zeros(n_sysID,4)
    y_yDot = zeros(n_sysID)
    A_yDot = zeros(n_sysID,4)

    for i=1:n_laps_sysID
        y_psi[(1:sz_last)+(i-1)*sz_last]    = diff(oldpsiDot[idx_s[i]-n_prev:idx_s[i]+n_ahead+1])
        A_psi[(1:sz_last)+(i-1)*sz_last,:]  = [oldpsiDot[vec_range_ID[i]]./oldxDot[vec_range_ID[i]] oldyDot[vec_range_ID[i]]./oldxDot[vec_range_ID[i]] olddF[vec_range_ID[i]]]
        y_xDot[(1:sz_last)+(i-1)*sz_last]   = diff(oldxDot[idx_s[i]-n_prev:idx_s[i]+n_ahead+1])
        A_xDot[(1:sz_last)+(i-1)*sz_last,:] = [oldyDot[vec_range_ID[i]] oldpsiDot[vec_range_ID[i]] oldxDot[vec_range_ID[i]] olda[vec_range_ID[i]]]
        y_yDot[(1:sz_last)+(i-1)*sz_last]   = diff(oldyDot[idx_s[i]-n_prev:idx_s[i]+n_ahead+1])
        A_yDot[(1:sz_last)+(i-1)*sz_last,:] = [oldyDot[vec_range_ID[i]]./oldxDot[vec_range_ID[i]] oldpsiDot[vec_range_ID[i]].*oldxDot[vec_range_ID[i]] oldpsiDot[vec_range_ID[i]]./oldxDot[vec_range_ID[i]] olddF[vec_range_ID[i]]]
    end

    # psiDot
    y_psi[(1:sz_curr)+n_laps_sysID*sz_last]   = diff(oldTraj.oldTraj[cC-n_prev:cC,3,cL])
    A_psi[(1:sz_curr)+n_laps_sysID*sz_last,:] = [oldTraj.oldTraj[vec_range_ID2,3,cL]./oldTraj.oldTraj[vec_range_ID2,1,cL] oldTraj.oldTraj[vec_range_ID2,2,cL]./oldTraj.oldTraj[vec_range_ID2,1,cL] oldTraj.oldInput[vec_range_ID2,2,cL]]
    
    # xDot
    y_xDot[(1:sz_curr)+n_laps_sysID*sz_last] = diff(oldTraj.oldTraj[cC-n_prev:cC,1,cL])
    A_xDot[(1:sz_curr)+n_laps_sysID*sz_last,:] = [oldTraj.oldTraj[vec_range_ID2,2,cL] oldTraj.oldTraj[vec_range_ID2,3,cL] oldTraj.oldTraj[vec_range_ID2,1,cL] oldTraj.oldInput[vec_range_ID2,1,cL]]
    
    # yDot
    y_yDot[(1:sz_curr)+n_laps_sysID*sz_last] = diff(oldTraj.oldTraj[cC-n_prev:cC,2,cL])
    A_yDot[(1:sz_curr)+n_laps_sysID*sz_last,:] = [oldTraj.oldTraj[vec_range_ID2,2,cL]./oldTraj.oldTraj[vec_range_ID2,1,cL] oldTraj.oldTraj[vec_range_ID2,3,cL].*oldTraj.oldTraj[vec_range_ID2,1,cL] oldTraj.oldTraj[vec_range_ID2,3,cL]./oldTraj.oldTraj[vec_range_ID2,1,cL] oldTraj.oldInput[vec_range_ID2,2]]

    if any(isnan,y_yDot)            # check if any value in the y_yDot value is NaN
        println(y_yDot)
        warn("NaN value detected in y_yDot! Press to continue...")
    end
    if any(isnan,coeffCost)
        println(coeffCost)
        warn("NaN value detected in coeffCost! Press to continue...")
    end
    if any(isnan,coeffConst)
        println(coeffCost)
        warn("NaN value detected in coeffConst! Press to continue...")
    end

    mpcCoeff.c_Psi = zeros(3)
    mpcCoeff.c_Vx = zeros(4)
    mpcCoeff.c_Vy = zeros(4)

    if det(A_psi'*A_psi) != 0
        mpcCoeff.c_Psi = (A_psi'*A_psi)\A_psi'*y_psi
    else
        warn("det y_psi = 0")
    end
    if det(A_xDot'*A_xDot) != 0
        mpcCoeff.c_Vx  = (A_xDot'*A_xDot)\A_xDot'*y_xDot         # the identity matrix is used to scale the coefficients
    else
        warn("det vx = 0")
    end
    if det(A_yDot'*A_yDot) != 0
        mpcCoeff.c_Vy  = (A_yDot'*A_yDot)\A_yDot'*y_yDot
    else
        warn("det vy = 0")
    end
    mpcCoeff.coeffCost  = coeffCost
    mpcCoeff.coeffConst = coeffConst
    nothing
end

# Notes about oldTrajectory:
# oldTrajectory[1:prebuf] = states before finish line (s < 0)
# oldTrajectory[prebuf+1:prebuf+cost] = states between start and finish line (0 <= s < s_target)
# oldTrajectory[prebuf+cost+1:prebuf+cost+postbuf] = states after finish line (s_target <= s)
# once one lap is over, the states are saved (interval prebuf:s_target)
# during the next lap, the states behind the finish line are appended (interval s_target:postbuf)
