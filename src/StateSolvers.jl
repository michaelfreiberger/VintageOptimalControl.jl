
"""
    state_PDE_solver(Stat_conc::Array{Float64,3},
                     Stat_dist::Array{Float64,3},
                     Stat_agg::Array{Float64,3},
                     Con_conc::Array{Float64,3},
                     Con_dist::Array{Float64,3},
                     Para::Dict)

state\\_PDE\\_solver solves the dynamic system of state equations with given initial conditions and control variables.
"""
function state_PDE_solver(Stat_conc::Array{Float64,3},
                          Stat_dist::Array{Float64,3},
                          Stat_agg::Array{Float64,3},
                          Con_conc::Array{Float64,3},
                          Con_dist::Array{Float64,3},
                          Para::Dict)
    
    Stat_conc[1,1,:] = Para["InitStat_conc"]
    Stat_dist[:,1,:] = Para["InitStat_dist"]
    Stat_dist[1,:,:] = Para["BoundStat_dist"]

    Dt_conc = zeros(Para["nStat_conc"])
    Dt_dist = zeros(Para["nStat_dist"])

    for ii = 1:Para["nTime"]-1
        # Aggregation for time ii
        Aggregate(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,ii,Para)
        # Initial values for second stage variable
        InitialDistVariables(Stat_conc,Stat_dist,Stat_agg,Con_conc,ii,Para)
        
        # First part of Heun-Method
        f_ODE_conc(Stat_conc,Stat_agg,Con_conc,ii,Dt_conc,Para)
        Stat_conc[1,ii+1,:] = Stat_conc[1,ii,:] + Para["hstep"]*Dt_conc
        for jj = 1:(Para["nVintage"]-1)
            f_PDE_dist(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,ii,jj,Dt_dist,Para)
            Stat_dist[jj+1,ii+1,:] = Stat_dist[jj,ii,:] + Para["hstep"]*Dt_dist
        end

        # Second part of Heun-Method
        Aggregate(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,ii+1,Para)
        InitialDistVariables(Stat_conc,Stat_dist,Stat_agg,Con_conc,ii+1,Para)
        f_ODE_conc(Stat_conc,Stat_agg,Con_conc,ii+1,Dt_conc,Para)
        Stat_conc[1,ii+1,:] = 0.5 * Stat_conc[1,ii+1,:] + 0.5 * (Stat_conc[1,ii,:] + Para["hstep"]*Dt_conc)
        for jj = 1:(Para["nVintage"]-1)
            f_PDE_dist(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,ii+1,jj+1,Dt_dist,Para)
            Stat_dist[jj+1,ii+1,:] = 0.5 * Stat_dist[jj+1,ii+1,:] + 0.5 * (Stat_dist[jj,ii,:] + Para["hstep"]*Dt_dist)
        end 
    end       
    #   Calculate aggregated variable for last step
    Aggregate(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,Para["nTime"],Para)
    #   Calculate initial value of distributed variable for last step
    InitialDistVariables(Stat_conc,Stat_dist,Stat_agg,Con_conc,Para["nTime"],Para)
end

"""
    costate_PDE_solver(Stat_conc::Array{Float64,3},
                       Stat_dist::Array{Float64,3}, 
                       Stat_agg::Array{Float64,3},
                       Con_conc::Array{Float64,3},
                       Con_dist::Array{Float64,3},
                       CoStat_conc::Array{Float64,3},
                       CoStat_dist::Array{Float64,3},
                       CoStat_agg::Array{Float64,3},
                       Para::Dict)

costate\\_PDE\\_solver solves the dynamic system of co-state equations backwards with given endvalues for given controls and state variables.
For the special case of endconstraints for the STATE-variables, we use the FOC to derive the respective value of CO-STATE variable.
"""
function costate_PDE_solver(Stat_conc::Array{Float64,3},
                            Stat_dist::Array{Float64,3}, 
                            Stat_agg::Array{Float64,3},
                            Con_conc::Array{Float64,3},
                            Con_dist::Array{Float64,3},
                            CoStat_conc::Array{Float64,3},
                            CoStat_dist::Array{Float64,3},
                            CoStat_agg::Array{Float64,3},
                            Para::Dict)
    Dt_conc = zeros(Para["nStat_conc"])
    Dt_dist = zeros(Para["nStat_dist"])

    CoStat_conc .= 0
    CoStat_agg .= 0    
    CoStat_dist .= 0
    
    EndConstraintCoStat_conc(Stat_conc,Stat_dist,CoStat_conc,Para)
    EndConstraintCoStat_dist(Stat_conc,Stat_dist,CoStat_dist,Para)
    CoStat_dist[end,:,:] .= 0.0

    for ii = Para["nTime"]:(-1):2
        if Para["nStat_agg"] > 0
            Aggregate_co(Stat_conc, Stat_dist, Stat_agg, Con_conc, Con_dist, CoStat_conc, CoStat_dist, CoStat_agg, ii, Para)
        end
        if Para["nStat_conc"] > 0
            f_ODE_conc_co(Stat_conc, Stat_dist, Stat_agg, Con_conc, Con_dist, CoStat_conc, CoStat_dist, CoStat_agg, ii, Dt_conc, Para)
            CoStat_conc[1,ii-1,:] =  CoStat_conc[1,ii,:] - Para["hstep"]*Dt_conc
        end
        if Para["nStat_dist"] > 0
            for jj = Para["nVintage"]:(-1):2
                f_PDE_dist_co(Stat_conc, Stat_dist, Stat_agg, Con_conc, Con_dist, CoStat_conc, CoStat_dist, CoStat_agg, ii, jj, Dt_dist, Para)
                CoStat_dist[jj-1,ii-1,:] = CoStat_dist[jj,ii,:] - Para["hstep"]*Dt_dist
            end
        end
        if Para["nStat_agg"] > 0
            Aggregate_co(Stat_conc, Stat_dist, Stat_agg, Con_conc, Con_dist, CoStat_conc, CoStat_dist, CoStat_agg, ii-1, Para)
        end
        if Para["nStat_conc"] > 0        
            f_ODE_conc_co(Stat_conc, Stat_dist, Stat_agg, Con_conc, Con_dist, CoStat_conc, CoStat_dist, CoStat_agg, ii-1, Dt_conc, Para)
            CoStat_conc[1,ii-1,:] =  0.5*CoStat_conc[1,ii-1,:] + 0.5*(CoStat_conc[1,ii,:] - Para["hstep"]*Dt_conc)
        end
        
        if Para["nStat_dist"] > 0
            for jj = Para["nVintage"]:(-1):2
                f_PDE_dist_co(Stat_conc, Stat_dist, Stat_agg, Con_conc, Con_dist, CoStat_conc, CoStat_dist, CoStat_agg, ii-1, jj-1, Dt_dist, Para)
                CoStat_dist[jj-1,ii-1,:] = 0.5*CoStat_dist[jj-1,ii-1,:] + 0.5*(CoStat_dist[jj,ii,:] - Para["hstep"]*Dt_dist)
            end
        end
    end
    if Para["nStat_agg"] > 0
        Aggregate_co(Stat_conc, Stat_dist, Stat_agg, Con_conc, Con_dist, CoStat_conc, CoStat_dist, CoStat_agg, 1, Para)
    end
end
