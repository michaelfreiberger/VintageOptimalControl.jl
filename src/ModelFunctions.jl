

"""
ObjectValue(Stat_conc::Array{Float64,3}, 
            Stat_dist::Array{Float64,3}, 
            Stat_agg::Array{Float64,3}, 
            Con_conc::Array{Float64,3},
            Con_dist::Array{Float64,3},
            Para::Dict )

Calculate the aggregated objective value of the optimal control problem for the current solutions of control and state variables
"""
function ObjectValue(Stat_conc::Array{Float64,3}, 
                     Stat_dist::Array{Float64,3}, 
                     Stat_agg::Array{Float64,3}, 
                     Con_conc::Array{Float64,3},
                     Con_dist::Array{Float64,3}, 
                     Para::Dict )
    F = zeros(Para["nTime"])
    Fhelp = zeros(Para["nVintage"])
    for ii = 1:Para["nTime"]
        Fhelp .= 0
        for jj = 1:Para["nVintage"]
            Fhelp[jj] = ObjectiveIntegrand(Stat_conc[1,ii,:] ,Stat_dist[jj,ii,:],Stat_agg[1,ii,:],
                                            Con_conc[1,ii,:] ,Con_dist[jj,ii,:],            
                                            Para["tmesh"][ii],Para["amesh"][jj],Para)
        end
        F[ii] = exp(-Para["TimeDiscountRate"]*Para["tmesh"][ii])*Integral(Para["nVintage"],Para["hstep"],Fhelp)
    end

    FSalv = zeros(Para["nVintage"])
    for jj = 1:Para["nVintage"]
        FSalv[jj] = SalvageFunction(Stat_conc[1,end,:],Stat_dist[jj,end,:],Para["tmesh"][end],Para["amesh"][jj],Para)
    end
    
    VV = Integral(Para["nTime"],Para["hstep"],F) + Integral(Para["nVintage"],Para["hstep"],FSalv)
    return VV
end

"""
Aggregate(Stat_conc::Array{Float64,3},
            Stat_dist::Array{Float64,3},
            Stat_agg::Array{Float64,3},
            Con_conc::Array{Float64,3},
            Con_dist::Array{Float64,3},
            ii::Int64,
            Para::Dict )

Calculate the aggregated variables over all vintages
"""
function Aggregate(Stat_conc::Array{Float64,3},
                   Stat_dist::Array{Float64,3},
                   Stat_agg::Array{Float64,3},
                   Con_conc::Array{Float64,3},
                   Con_dist::Array{Float64,3},
                   ii::Int64,
                   Para::Dict )
    F = zeros(Para["nVintage"],Para["nStat_agg"])
    for jj = 1:Para["nVintage"]
        F[jj,:] .= AggregationIntegrand(Stat_conc[1,ii,:], Stat_dist[jj,ii,:],
                                        Con_conc[1,ii,:], Con_dist[jj,ii,:],
                                        Para["tmesh"][ii],Para["amesh"][jj],Para)
    end
    for kk = 1:Para["nStat_agg"]
        Stat_agg[1,ii,kk] = Integral(Para["nVintage"],Para["hstep"],F[:,kk])
    end
end

"""
InitialDistVariables(Stat_conc::Array{Float64,3},
                     Stat_dist::Array{Float64,3},
                     Stat_agg::Array{Float64,3},
                     Con_conc::Array{Float64,3},
                     ii::Int64,
                     Para::Dict )

Calculate the starting values for all vintages.
"""
function InitialDistVariables(Stat_conc::Array{Float64,3},
                              Stat_dist::Array{Float64,3},
                              Stat_agg::Array{Float64,3},
                              Con_conc::Array{Float64,3},
                              ii::Int64,
                              Para::Dict )
    Stat_dist[1,ii,:] .= InitDist(Stat_conc[1,ii,:],Stat_agg[1,ii,:],Con_conc[1,ii,:],Para["tmesh"][ii],Para);
end

"""
f_ODE_conc( Stat_conc::Array{Float64,3}, 
            Stat_agg::Array{Float64,3}, 
            Con_conc::Array{Float64,3}, 
            tt::Int64, 
            Dt::Array{Float64,1}, 
            Para::Dict )

Right hand side of the ODE for the concentrated state variables at time index tt.
"""
function f_ODE_conc(Stat_conc::Array{Float64,3}, 
                    Stat_agg::Array{Float64,3}, 
                    Con_conc::Array{Float64,3}, 
                    tt::Int64,  
                    Dt::Array{Float64,1}, 
                    Para::Dict )
    Dt .= StateDynamic_conc(Stat_conc[1,tt,:],Stat_agg[1,tt,:],Con_conc[1,tt,:],Para["tmesh"][tt],Para)
end

"""
f_ODE_interstep(Stat_conc, Stat_agg, Con_conc, t::Float64, Para::Dict)

Right hand side of the ODE for the concentrated state variables. Only used for the 4th order RK-method where the control variables have to be interpolated.
"""
function f_ODE_interstep(Stat_conc, 
                         Stat_agg, 
                         Con_conc, 
                         t::Float64, 
                         Para::Dict)
    k = StateDynamic_conc(Stat_conc,Stat_agg,Con_conc,t,Para)
    return k
end

"""
f_PDE_dist(Stat_conc::Array{Float64,3},
            Stat_dist::Array{Float64,3}, 
            Stat_agg::Array{Float64,3},
            Con_conc::Array{Float64,3}, 
            Con_dist::Array{Float64,3}, 
            tt::Int64, ss::Int64, 
            Dt::Array{Float64,1},
            Para::Dict )

Right hand side of the PDE for the state variables in the second stage at time index tt and vintage ss.
"""
function f_PDE_dist(Stat_conc::Array{Float64,3},
                    Stat_dist::Array{Float64,3}, 
                    Stat_agg::Array{Float64,3},
                    Con_conc::Array{Float64,3}, 
                    Con_dist::Array{Float64,3}, 
                    tt::Int64, 
                    ss::Int64, 
                    Dt::Array{Float64,1},
                    Para::Dict  )
    Dt .= StateDynamic_dist(Stat_conc[1,tt,:],Stat_dist[ss,tt,:],Stat_agg[1,tt,:],Con_conc[1,tt,:],Con_dist[ss,tt,:],Para["tmesh"][tt],Para["amesh"][ss],Para)
end


"""
Hamiltonian(Stat_conc, 
            Stat_dist, 
            Stat_agg, 
            Con_conc, 
            Con_dist, 
            CoStat_conc,
            CoStat_dist,
            CoStat_agg, 
            t::Float64,
            s::Float64,
            Para::Dict )

Definition of the Hamiltonian.
"""
function Hamiltonian(Stat_conc, 
                     Stat_dist, 
                     Stat_agg, 
                     Con_conc, 
                     Con_dist, 
                     CoStat_conc,
                     CoStat_dist,
                     CoStat_agg, 
                     t::Float64,
                     s::Float64,
                     Para::Dict )

    if Para["nStat_conc"] > 0
        Stat_conc_term = 1/Para["ω"] * dot(CoStat_conc,StateDynamic_conc(Stat_conc,Stat_agg,Con_conc,t,Para)) 
    else
        Stat_conc_term = 0.0
    end
    if Para["nStat_agg"] > 0
        Stat_agg_term = dot(CoStat_agg,AggregationIntegrand(Stat_conc, Stat_dist,Con_conc, Con_dist,t,s,Para))
    else
        Stat_agg_term = 0.0
    end
    if Para["nStat_dist"] > 0
        Stat_dist_term = dot(CoStat_dist,StateDynamic_dist(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,t,s,Para)) + 
                        + dot(CoStat_dist,InitDist(Stat_conc,Stat_agg,Con_conc,t,Para))
    else
        Stat_dist_term = 0.0
    end
    
    #exp(-Para["TimeDiscountRate"]*t)*

    return ObjectiveIntegrand(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,t,s,Para) + 
            + Stat_conc_term + Stat_dist_term + Stat_agg_term

    # return ObjectiveIntegrand(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,t,s,Para) + 
    #         + 1/Para["ω"] * dot(CoStat_conc,StateDynamic_conc(Stat_conc,Stat_agg,Con_conc,t,Para)) + 
    #         + dot(CoStat_dist,StateDynamic_dist(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,t,s,Para)) +
    #         + dot(CoStat_agg,AggregationIntegrand(Stat_conc, Stat_dist,Con_conc, Con_dist,t,s,Para)) 
    #         + dot(CoStat_dist,InitDist(Stat_conc,Stat_agg,Con_conc,t,Para))
end


"""
f_ODE_co(Stat_conc::Array{Float64,3}, 
         Stat_dist::Array{Float64,3},
         Stat_agg::Array{Float64,3},
         Con_conc::Array{Float64,3},
         Con_dist::Array{Float64,3},
         CoStat_conc::Array{Float64,3},
         CoStat_dist::Array{Float64,3},
         CoStat_agg::Array{Float64,3},
         tt::Int64,
         Dt::Array{Float64,1},
         Para::Dict )

Right hand side of the PDE for the concentrated costate variables at time index tt.
"""
function f_ODE_conc_co(Stat_conc::Array{Float64,3}, 
                       Stat_dist::Array{Float64,3}, 
                       Stat_agg::Array{Float64,3}, 
                       Con_conc::Array{Float64,3}, 
                       Con_dist::Array{Float64,3}, 
                       CoStat_conc::Array{Float64,3}, 
                       CoStat_dist::Array{Float64,3}, 
                       CoStat_agg::Array{Float64,3}, 
                       tt::Int64, 
                       Dt::Array{Float64,1}, 
                       Para::Dict )
    
    F = zeros(Para["nVintage"],Para["nStat_conc"])
    for jj = 1:Para["nVintage"]
        F[ii,:] = ForwardDiff.gradient(X->Hamiltonian(X, Stat_dist[jj,tt,:], Stat_agg[1,tt,:], Con_conc[1,tt,:], Con_dist[jj,tt,:], 
                                                        CoStat_conc[1,tt,:],CoStat_dist[jj,tt,:],CoStat_agg[1,tt,:], Para["tmesh"][tt],Para["amesh"][jj],Para),
                                                     Stat_conc[1,tt,:])
    end   
            
    Dt .= -ForwardDiff.gradient(X -> dot(CoStat_dist[1,tt,:],InitDist(X,Stat_agg[1,tt,:],Con_conc[1,tt,:],Para["tmesh"][tt],Para)),
                                    Stat_conc[1,tt,:])
    for kk = 1:Para["nStat_conc"]    
        Dt[kk] = Dt[kk] - Integral(Para["nVintage"],Para["hstep"],F[:,kk])
    end 
    Dt .= Para["TimeDiscountRate"]*CoStat_conc[1,tt,:] + Dt
end

"""
f_PDE_dist_co( Stat_conc::Array{Float64,3}, 
                Stat_dist::Array{Float64,3},
                Stat_agg::Array{Float64,3},
                Con_conc::Array{Float64,3},
                Con_dist::Array{Float64,3},
                CoStat_conc::Array{Float64,3},
                CoStat_dist::Array{Float64,3},
                CoStat_agg::Array{Float64,3},
                tt::Int64,
                ss::Int64,
                Dt::Array{Float64,1},
                Para::Dict )

Right hand side of the PDE for the distributed costate variables at time index tt and vintage ss.
"""
function f_PDE_dist_co(Stat_conc::Array{Float64,3}, 
                       Stat_dist::Array{Float64,3}, 
                       Stat_agg::Array{Float64,3}, 
                       Con_conc::Array{Float64,3}, 
                       Con_dist::Array{Float64,3}, 
                       CoStat_conc::Array{Float64,3}, 
                       CoStat_dist::Array{Float64,3}, 
                       CoStat_agg::Array{Float64,3}, 
                       tt::Int64, 
                       ss::Int64, 
                       Dt::Array{Float64,1}, 
                       Para::Dict )
    Dt .= -ForwardDiff.gradient(X->Hamiltonian(Stat_conc[1,tt,:], X, Stat_agg[1,tt,:], Con_conc[1,tt,:], Con_dist[ss,tt,:], 
                                CoStat_conc[1,tt,:],CoStat_dist[ss,tt,:],CoStat_agg[1,tt,:], Para["tmesh"][tt],Para["amesh"][ss],Para),Stat_dist[ss,tt,:])

    Dt .= Para["TimeDiscountRate"]*CoStat_dist[ss,tt,:] + Dt
end


"""
Aggregate_co( Stat_conc::Array{Float64,3}, 
                Stat_dist::Array{Float64,3},
                Stat_agg::Array{Float64,3},
                Con_conc::Array{Float64,3},
                Con_dist::Array{Float64,3},
                CoStat_conc::Array{Float64,3},
                CoStat_dist::Array{Float64,3},
                CoStat_agg::Array{Float64,3},
                tt::Int64,
                Para::Dict )

Calculation of the aggregated costate variables at time index tt.
"""
function Aggregate_co(Stat_conc::Array{Float64,3}, 
                      Stat_dist::Array{Float64,3}, 
                      Stat_agg::Array{Float64,3}, 
                      Con_conc::Array{Float64,3}, 
                      Con_dist::Array{Float64,3}, 
                      CoStat_conc::Array{Float64,3}, 
                      CoStat_dist::Array{Float64,3},
                      CoStat_agg::Array{Float64,3}, 
                      tt::Int64, 
                      Para::Dict)
    
    CoStat_agg[1,tt,:] .= ForwardDiff.gradient(X -> dot(CoStat_dist[1,tt,:],InitDist(Stat_conc[1,tt,:],X,Con_conc[1,tt,:],Para["tmesh"][tt],Para)),
                                                        Stat_agg[1,tt,:])
    
    F = zeros(Para["nVintage"],Para["nStat_agg"])
    for jj = 1:Para["nVintage"]
        F[jj,:] .= ForwardDiff.gradient(X->Hamiltonian(Stat_conc[1,tt,:], Stat_dist[jj,tt,:], X, Con_conc[1,tt,:], Con_dist[jj,tt,:], 
                                                        CoStat_conc[1,tt,:],CoStat_dist[jj,tt,:],CoStat_agg[1,tt,:], Para["tmesh"][tt],Para["amesh"][jj],Para),
                                                      Stat_agg[1,tt,:])
    end   
            
    for kk = 1:Para["nStat_agg"]    
        CoStat_agg[1,tt,kk] = CoStat_agg[1,tt,kk] + Integral(Para["nVintage"],Para["hstep"],F[:,kk])
    end 
end


"""
GradHamiltonian(Stat_conc::Array{Float64,3},
                Stat_dist::Array{Float64,3},
                Stat_agg::Array{Float64,3},
                Con_conc::Array{Float64,3},
                Con_dist::Array{Float64,3},
                CoStat_conc::Array{Float64,3},
                CoStat_dist::Array{Float64,3},
                CoStat_agg::Array{Float64,3},
                dHam_conc::Array{Float64,3},
                dHam_dist::Array{Float64,3},
                Para::Dict)

Calculation of the gradient of the Hamiltonian.
"""
function GradHamiltonian(Stat_conc::Array{Float64,3},
                         Stat_dist::Array{Float64,3},
                         Stat_agg::Array{Float64,3},
                         Con_conc::Array{Float64,3},
                         Con_dist::Array{Float64,3},
                         CoStat_conc::Array{Float64,3},
                         CoStat_dist::Array{Float64,3},
                         CoStat_agg::Array{Float64,3},
                         dHam_conc::Array{Float64,3},
                         dHam_dist::Array{Float64,3},
                         Para::Dict)
    
    dHam_conc .= 0
    dHam_dist .= 0
    
    # if Para["ParallelComputing"] == true
    #     Threads.@threads for ii = 1:Para["nTime"]    
    #         for jj = 1:Para["nVintage"]
    #             Fhelp[jj,ii,:] = ForwardDiff.gradient(X->Hamiltonian(Stat_conc[jj,ii,:], Stat_dist[jj,ii,:], Stat_agg[1,ii,:], X, Con_dist[jj,ii,:], 
    #                                                         CoStat_conc[1,ii,:],CoStat_dist[jj,ii,:],CoStat_agg[1,ii,:], Para["tmesh"][ii],Para["amesh"][jj],Para),Con_conc[1,ii,:])
    #             dHam_dist[jj,ii,:] = ForwardDiff.gradient(X->Hamiltonian(Stat_conc[jj,ii,:], Stat_dist[jj,ii,:], Stat_agg[1,ii,:], Con_conc[1,ii,:], X, 
    #                                                         CoStat_conc[1,ii,:],CoStat_dist[jj,ii,:],CoStat_agg[1,ii,:], Para["tmesh"][ii],Para["amesh"][jj],Para),Con_dist[jj,ii,:])
    #         end
    #         dHam_conc[1,ii,:] = ForwardDiff.gradient(X -> dot(CoStat_dist[1,ii,:],InitDist(Stat_conc[1,ii,:],Stat_agg[1,ii,:],X,Para["tmesh"][ii],Para)),Con_conc[1,ii,:])
    #         for kk = 1:Para["nCon_conc"]
    #             dHam_conc[1,ii,kk] = dHam_conc[1,ii,kk] + Integral(Para["nVintage"],Para["hstep"],Fhelp[:,ii,kk]) 
    #         end
    #     end
    # else
        Fhelp = zeros(Para["nVintage"],Para["nTime"],Para["nCon_conc"])
        if Para["nCon_conc"] > 0
            for ii = 1:Para["nTime"]    
                for jj = 1:Para["nVintage"]
                    Fhelp[jj,ii,:] = ForwardDiff.gradient(X->Hamiltonian(Stat_conc[1,ii,:], Stat_dist[jj,ii,:], Stat_agg[1,ii,:], X, Con_dist[jj,ii,:], 
                                                        CoStat_conc[1,ii,:],CoStat_dist[jj,ii,:],CoStat_agg[1,ii,:], Para["tmesh"][ii],Para["amesh"][jj],Para),Con_conc[1,ii,:])
                end
                dHam_conc[1,ii,:] = ForwardDiff.gradient(X -> dot(CoStat_dist[1,ii,:],InitDist(Stat_conc[1,ii,:],Stat_agg[1,ii,:],X,Para["tmesh"][ii],Para)),Con_conc[1,ii,:])
                for kk = 1:Para["nCon_conc"]
                    dHam_conc[1,ii,kk] = dHam_conc[1,ii,kk] + Integral(Para["nVintage"],Para["hstep"],Fhelp[:,ii,kk]) 
                end
            end
        end

        if Para["nCon_dist"] > 0
            for ii = 1:Para["nTime"]    
                for jj = 1:Para["nVintage"]    
                    dHam_dist[jj,ii,:] = ForwardDiff.gradient(X->Hamiltonian(Stat_conc[1,ii,:], Stat_dist[jj,ii,:], Stat_agg[1,ii,:], Con_conc[1,ii,:], X, 
                                            CoStat_conc[1,ii,:],CoStat_dist[jj,ii,:],CoStat_agg[1,ii,:], Para["tmesh"][ii],Para["amesh"][jj],Para),Con_dist[jj,ii,:])
                end
            end
        end 
    # end
end

"""
NewDirection(Stat_conc::Array{Float64,3},
             Stat_dist::Array{Float64,3},
             Stat_agg::Array{Float64,3},
             Con_conc::Array{Float64,3},
             Con_dist::Array{Float64,3},
             CoStat_conc::Array{Float64,3},
             CoStat_dist::Array{Float64,3},
             CoStat_agg::Array{Float64,3},
             dHam_conc::Array{Float64,3},
             dHam_dist::Array{Float64,3},
             Para::Dict)

Adjustment of the gradient of the Hamiltonian depending on the optimisation type chosen. Available types are:
- "Gradient" -> Uses the gradient without adjustments
- "Newton-Raphson" -> Adjusts the gradient as described in the Newton-Raphson method using the Hessian of the Hamiltonian
"""
function NewDirection(Stat_conc::Array{Float64,3},
                      Stat_dist::Array{Float64,3},
                      Stat_agg::Array{Float64,3},
                      Con_conc::Array{Float64,3},
                      Con_dist::Array{Float64,3},
                      CoStat_conc::Array{Float64,3},
                      CoStat_dist::Array{Float64,3},
                      CoStat_agg::Array{Float64,3},
                      dHam_conc::Array{Float64,3},
                      dHam_dist::Array{Float64,3},
                      Para::Dict)

    if Para["OptiType"] == "Gradient"
        dHam_conc .= dHam_conc
        dHam_dist .= dHam_dist
    elseif Para["OptiType"] == "Newton-Raphson"
        if Para["nCon_dist"] > 0
            HessianDist(ii,jj) = ForwardDiff.hessian(X->Hamiltonian(Stat_conc[1,ii,:], Stat_dist[jj,ii,:], Stat_agg[1,ii,:], Con_conc[1,ii,:], X, 
                                                        CoStat_conc[1,ii,:],CoStat_dist[jj,ii,:],CoStat_agg[1,ii,:], Para["tmesh"][ii],Para["amesh"][jj],Para),Con_dist[jj,ii,:])
            for jj = 1:Para["nVintage"]
                for ii = 1:Para["nTime"]
                    if abs(det(HessianDist(ii,jj))) > 1e-9
                        dHam_dist[jj,ii,:] .= - inv(HessianDist(ii,jj))*dHam_dist[jj,ii,:] 
                    end                    
                end
            end
        end
        if Para["nCon_conc"] > 0
            HessianConc = function(ii)
                Hess = zeros(Para["nVintage"],Para["nCon_conc"],Para["nCon_conc"])
                for jj = 1:Para["nVintage"]
                    Hess[jj,:,:] = ForwardDiff.hessian(X->Hamiltonian(Stat_conc[1,ii,:], Stat_dist[jj,ii,:], Stat_agg[1,ii,:], X, Con_dist[jj,ii,:], 
                                                CoStat_conc[1,ii,:],CoStat_dist[jj,ii,:],CoStat_agg[1,ii,:], Para["tmesh"][ii],Para["amesh"][jj],Para),Con_conc[1,ii,:])
                end

                Hess2 = ForwardDiff.hessian(X -> dot(CoStat_dist[1,ii,:],InitDist(Stat_conc[1,ii,:],Stat_agg[1,ii,:],X,Para["tmesh"][ii],Para)),Con_conc[1,ii,:])                
                for kk1 = 1:Para["nCon_conc"]
                    for kk2 = 1:Para["nCon_conc"]
                        Hess2[kk1,kk2] = Hess2[kk1,kk2] + Integral(Para["nVintage"],Para["hstep"],Hess[:,kk1,kk2])
                    end
                end
                return Hess2
            end
            
            for ii = 1:Para["nTime"]
                Hess = HessianConc(ii)
                if abs(det(Hess)) > 1e-9
                    dHam_conc[1,ii,:] .= -inv(Hess)*dHam_conc[1,ii,:]
                end
            end
        end
    else
        Para["OptiType"] = "Gradient"
    end

    #=----------------------------
        Smoothing
    ----------------------------=#
    if Para["GradSmooth"] == true
        GradSmooth(dHam_dist,dHam_conc,Para)
    end

    #=-------------------------------------------------
        Reducing the model
    -------------------------------------------------=#
    for kk in Para["ReduceConc"]
        dHam_conc[:,:,kk] .= 0
    end
    for kk in Para["ReduceDist"]
        dHam_dist[:,:,kk] .= 0
    end
end

"""
ConMapping_conc(Con_conc::Array{Float64,3},Para)

Map concentrated control variables into the feasible region described by the lower and upper bounds of the controls.
"""
function ConMapping_conc(Con_conc::Array{Float64,3},Para)
    for kk = 1:Para["nCon_conc"]
        if !(kk in Para["ReduceConc"])
            Con_conc[1,:,kk] .= max.(Con_conc[1,:,kk],Para["Con_concMin"][kk])
            Con_conc[1,:,kk] .= min.(Con_conc[1,:,kk],Para["Con_concMax"][kk])
        end
    end
end

"""
ConMapping_dist(Con_dist::Array{Float64,3},Para)

Map distributed control variables into the feasible region described by the lower and upper bounds of the controls.
"""
function ConMapping_dist(Con_dist::Array{Float64,3},Para)
    for kk = 1:Para["nCon_dist"]
        if !(kk in Para["ReduceDist"])
            Con_dist[:,:,kk] = max.(Con_dist[:,:,kk],Para["Con_distMin"][kk])
            Con_dist[:,:,kk] = min.(Con_dist[:,:,kk],Para["Con_distMax"][kk])
        end
    end
end

"""
EndConstraintCoStat_conc(Stat_conc::Array{Float64,3},
                         Stat_dist::Array{Float64,3},
                         CoStat_conc::Array{Float64,3},
                         Para::Dict)

Define the endconstraint for the concentrated costate variables using the derivative of the salvage-function.
"""
function EndConstraintCoStat_conc(Stat_conc::Array{Float64,3},
                                  Stat_dist::Array{Float64,3},
                                  CoStat_conc::Array{Float64,3},
                                  Para::Dict)
    F = zeros(Para["nVintage"],Para["nStat_conc"])
    for jj = 1:Para["nVintage"]
        F[jj,:] = ForwardDiff.gradient(X->SalvageFunction(X,Stat_dist[jj,end,:],Para["tmesh"][end],Para["amesh"][jj],Para),Stat_conc[1,end,:])
    end
    for kk = 1:Para["nStat_conc"]
        CoStat_conc[1,Para["nTime"],kk] = Integral(Para["nVintage"],Para["hstep"],F[:,kk])
    end
end

"""
EndConstraintCoStat_dist(Stat_conc::Array{Float64,3},
                         Stat_dist::Array{Float64,3},
                         CoStat_dist::Array{Float64,3},
                         Para::Dict)

Define the endconstraint for the distributed costate variables using the derivative of the salvage-function.
"""
function EndConstraintCoStat_dist(Stat_conc::Array{Float64,3},
                                  Stat_dist::Array{Float64,3},
                                  CoStat_dist::Array{Float64,3},
                                  Para::Dict)
    CoStat_dist[end,:,:] .= 0.0
    for jj = 1:Para["nVintage"]
        CoStat_dist[jj,Para["nTime"],:] = ForwardDiff.gradient(X->SalvageFunction(Stat_conc[1,end,:],X,Para["tmesh"][end],Para["amesh"][jj],Para),Stat_dist[jj,end,:])
    end
end


