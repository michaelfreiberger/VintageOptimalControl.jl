
"""
LineSearch(Stat_conc::Array{Float64,3}, 
            Stat_dist::Array{Float64,3}, 
            Stat_agg::Array{Float64,3}, 
            Con_conc::Array{Float64,3}, 
            Con_dist::Array{Float64,3},  
            Dir_u_conc::Array{Float64,3}, 
            Dir_u_dist::Array{Float64,3}, 
            ObjValue::Float64, 
            Step, 
            Para::Dict,
            Stat_conc_new::Array{Float64,3}, 
            Stat_dist_new::Array{Float64,3}, 
            Stat_agg_new::Array{Float64,3},
            Con_conc_new::Array{Float64,3}, 
            Con_dist_new::Array{Float64,3}, 
            Stat_conc_best::Array{Float64,3}, 
            Stat_dist_best::Array{Float64,3}, 
            Stat_agg_best::Array{Float64,3},
            Con_conc_best::Array{Float64,3},
            Con_dist_best::Array{Float64,3} 
            )

Searches for an improvement along the gradient using a quadratic approximation for the objective function along the gradient.
"""
function LineSearch(Stat_conc::Array{Float64,3}, 
                    Stat_dist::Array{Float64,3}, 
                    Stat_agg::Array{Float64,3}, 
                    Con_conc::Array{Float64,3}, 
                    Con_dist::Array{Float64,3},  
                    Dir_u_conc::Array{Float64,3}, 
                    Dir_u_dist::Array{Float64,3}, 
                    ObjValue::Float64, 
                    Step, 
                    Para::Dict,
                    Stat_conc_new::Array{Float64,3}, 
                    Stat_dist_new::Array{Float64,3}, 
                    Stat_agg_new::Array{Float64,3},
                    Con_conc_new::Array{Float64,3}, 
                    Con_dist_new::Array{Float64,3}, 
                    Stat_conc_best::Array{Float64,3}, 
                    Stat_dist_best::Array{Float64,3}, 
                    Stat_agg_best::Array{Float64,3},
                    Con_conc_best::Array{Float64,3},
                    Con_dist_best::Array{Float64,3})
            # Perform line search for new optimal control

    impr = 0
    step_increase = Para["LineSearchStepIncrease"]

    Para["LineIter"] = 0
    Para["UpperLineStepTemporary"] = Step*Para["UpperLineStepTemporaryFactor"]
    
    ObjValue_new = 0
    ObjValue_best = 0

    while (impr==0) && (Para["LineIter"]<=Para["MaxLineIter"]) && (Step > Para["stepLowBound"])

        Step = min(Step,Para["UpperLineStep"])

        ObjValue_new, Step = Update(Con_conc, Con_dist, Dir_u_conc, Dir_u_dist, Step, Para,
                             Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new)

        y0 = ObjValue
        y1 = ObjValue_new
        x1 = Step
        if y1 > y0
            impr = 1

            # Save currently best values
            ObjValue_best, step_best = AssignBest(ObjValue_new, Stat_conc_new, Stat_dist_new, Stat_agg_new,
                                                                Con_conc_new, Con_dist_new, Step, 
                                                                Stat_conc_best, Stat_dist_best, Stat_agg_best, 
                                                                Con_conc_best, Con_dist_best)

            Step = (1+step_increase)*Step
            x2 = Step
            ObjValue_new, Step = Update(Con_conc, Con_dist, Dir_u_conc, Dir_u_dist, Step, Para,
                                        Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new)
            y2 = ObjValue_new

            z = x1*(y2-y0) - x2*(y1-y0)
            
            if y2 > y1
                # Save currently best values
                ObjValue_best, step_best = AssignBest(ObjValue_new, Stat_conc_new, Stat_dist_new, Stat_agg_new,
                                                        Con_conc_new, Con_dist_new, Step, 
                                                        Stat_conc_best, Stat_dist_best, Stat_agg_best, 
                                                        Con_conc_best, Con_dist_best)

                if z >= 0       # Check if curve is convex
                    # Curve is convex --> second solution is the best
                else
                    # Curve is concave --> stepsize for maximum of quadratic approximation
                    st = ( (y2-y0)*x1^2 - (y1-y0)*x2^2 ) / ( 2*z )
                    
                    ObjValue_new, st = Update(Con_conc, Con_dist, Dir_u_conc, Dir_u_dist, st, Para,
                                            Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new)
                    
                    if ObjValue_new > y2
                        # Save currently best values
                        ObjValue_best, step_best = AssignBest(ObjValue_new, Stat_conc_new, Stat_dist_new, Stat_agg_new,
                                                                Con_conc_new, Con_dist_new, st, 
                                                                Stat_conc_best, Stat_dist_best, Stat_agg_best, 
                                                                Con_conc_best, Con_dist_best)
                    end
                end
            else    #(y2 <= y1)
                # Curve is concave --> stepsize for maximum of quadratic approximation
                st = ( (y2-y0)*x1^2 - (y1-y0)*x2^2 ) / ( 2*z )
                
                ObjValue_new, st = Update(Con_conc, Con_dist, Dir_u_conc, Dir_u_dist, st, Para,
                                        Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new)
                                
                if ObjValue_new > y1
                    # Save currently best values
                    ObjValue_best, step_best = AssignBest(ObjValue_new, Stat_conc_new, Stat_dist_new, Stat_agg_new,
                                                            Con_conc_new, Con_dist_new, st, 
                                                            Stat_conc_best, Stat_dist_best, Stat_agg_best, 
                                                            Con_conc_best, Con_dist_best)
                end
            end
        else     # (y1 <= y0)
            Step = 0.5*Step
            x2 = Step

            ObjValue_new, Step = Update(Con_conc, Con_dist, Dir_u_conc, Dir_u_dist, Step, Para,
                                        Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new)
            y2 = ObjValue_new

            if y2 > y0
                impr = 1
                # Save currently best values
                ObjValue_best, step_best = AssignBest(ObjValue_new, Stat_conc_new, Stat_dist_new, Stat_agg_new,
                                                        Con_conc_new, Con_dist_new, Step, 
                                                        Stat_conc_best, Stat_dist_best, Stat_agg_best, 
                                                        Con_conc_best, Con_dist_best)
                
                # Curve is concave --> stepsize for maximum of quadratic approximation
                z = x2*(y1-y0) - x1*(y2-y0)
                st = ( x2^2*(y1-y0) -x1^2*(y2-y0) ) / ( 2*z )
                
                # New controls + new states + new objective function()
                ObjValue_new, st = Update(Con_conc, Con_dist, Dir_u_conc, Dir_u_dist, st, Para,
                                            Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new)

                if ObjValue_new > y2
                    # Save currently best values
                    ObjValue_best, step_best = AssignBest(ObjValue_new, Stat_conc_new, Stat_dist_new, Stat_agg_new,
                                                            Con_conc_new, Con_dist_new, st, 
                                                            Stat_conc_best, Stat_dist_best, Stat_agg_best, 
                                                            Con_conc_best, Con_dist_best)
                end
            elseif y2 > 0.5*(y0+y1)     # && y2 < (3*y0+y1)/4
                z = x2*(y1-y0) - x1*(y2-y0)
                st = ( x2^2*(y1-y0) -x1^2*(y2-y0) ) / ( 2*z )
                if st > 0
                    ObjValue_new, st = Update(Con_conc, Con_dist, Dir_u_conc, Dir_u_dist, st, Para,
                                            Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new)

                    if ObjValue_new >= y0
                        impr = 1
                        # Save currently best values
                        ObjValue_best, step_best = AssignBest(ObjValue_new, Stat_conc_new, Stat_dist_new, Stat_agg_new,
                                                                Con_conc_new, Con_dist_new, st, 
                                                                Stat_conc_best, Stat_dist_best, Stat_agg_best, 
                                                                Con_conc_best, Con_dist_best)
                    end
                end
            end
        end
        if Para["SearchDISP"] == true
            if impr==1
                printSuccess(Para, ObjValue_best, step_best)
            else
                printNoSuccess(Para, ObjValue, Step)
            end
        end
    end

    if impr == 0
        # If there is no improvement, return original values
        step_best = Step
        # only step-size gets adjusted --> new try from the beginning with the new step-size
    else
        # Assign the best values to the basic variables without _best or _new extension
        ObjValue, Step = AssignBest(ObjValue_best, Stat_conc_best, Stat_dist_best, Stat_agg_best,
                                    Con_conc_best, Con_dist_best, step_best, 
                                    Stat_conc, Stat_dist, Stat_agg, 
                                    Con_conc, Con_dist)
    end

    return ObjValue,step_best,impr
end


"""
Update(Con_conc::Array{Float64,3},
        Con_dist::Array{Float64,3},
        Dir_u_conc::Array{Float64,3},
        Dir_u_dist::Array{Float64,3},
        Step,
        Para::Dict,
        Stat_conc_new::Array{Float64,3},
        Stat_dist_new::Array{Float64,3}, 
        Stat_agg_new::Array{Float64,3},
        Con_conc_new::Array{Float64,3},
        Con_dist_new::Array{Float64,3})

Update calculates the new values for Controls `Con_new`, States `Stat_new` and Objectiv value given start values for the Controls `Con`, the search direction `Dir_u` and the step-size in the direction `step`
"""
function Update(Con_conc::Array{Float64,3},
                Con_dist::Array{Float64,3},
                Dir_u_conc::Array{Float64,3},
                Dir_u_dist::Array{Float64,3},
                Step,
                Para::Dict,
                Stat_conc_new::Array{Float64,3},
                Stat_dist_new::Array{Float64,3}, 
                Stat_agg_new::Array{Float64,3},
                Con_conc_new::Array{Float64,3},
                Con_dist_new::Array{Float64,3})

    Step = min(Step,Para["UpperLineStep"],Para["UpperLineStepTemporary"])

    Con_conc_new .= Con_conc + Step*Dir_u_conc
    Con_dist_new .= Con_dist + Step*Dir_u_dist
    ConMapping_conc(Con_conc_new,Para)
    ConMapping_dist(Con_dist_new,Para)
    
    state_PDE_solver(Stat_conc_new,Stat_dist_new,Stat_agg_new,Con_conc_new,Con_dist_new,Para)
    ObjValue_new = ObjectValue( Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new,Con_dist_new, Para)
    Para["LineIter"] = Para["LineIter"] + 1

    return ObjValue_new, Step
end

"""
AssignBest(ObjValue::Float64,
            Stat_conc::Array{Float64,3},
            Stat_dist::Array{Float64,3}, 
            Stat_agg::Array{Float64,3}, 
            Con_conc::Array{Float64,3},
            Con_dist::Array{Float64,3},
            step, 
            Stat_conc_best::Array{Float64,3}, 
            Stat_dist_best::Array{Float64,3},  
            Stat_agg_best::Array{Float64,3},
            Con_conc_best::Array{Float64,3}, 
            Con_dist_best::Array{Float64,3})

Assign the current variables to the "best" variables
"""
function AssignBest(ObjValue::Float64,
                    Stat_conc::Array{Float64,3},
                    Stat_dist::Array{Float64,3},
                    Stat_agg::Array{Float64,3}, 
                    Con_conc::Array{Float64,3},
                    Con_dist::Array{Float64,3},
                    step, 
                    Stat_conc_best::Array{Float64,3},
                    Stat_dist_best::Array{Float64,3},
                    Stat_agg_best::Array{Float64,3},
                    Con_conc_best::Array{Float64,3},
                    Con_dist_best::Array{Float64,3})

    Stat_conc_best .= Stat_conc
    Stat_dist_best .= Stat_dist
    Stat_agg_best .= Stat_agg
    Con_conc_best .= Con_conc
    Con_dist_best .= Con_dist
    
    return ObjValue, step
end


"""
maxabsGradient(Con_conc::Array{Float64,3},
                Con_dist::Array{Float64,3},
                dHam_conc::Array{Float64,3},
                dHam_dist::Array{Float64,3},
                Para::Dict)

Calculate the maximum of the absolute value of the gradient of the Hamiltonian accounting for controls which are on the boundary of their feasible regions.
"""
function maxabsGradient(Con_conc::Array{Float64,3},
                        Con_dist::Array{Float64,3},
                        dHam_conc::Array{Float64,3},
                        dHam_dist::Array{Float64,3},
                        Para::Dict)
    dHam_conc_abs = zeros(1,Para["nTime"],Para["nCon_conc"])
    dHam_dist_abs = zeros(Para["nVintage"],Para["nTime"],Para["nCon_dist"])
    
    for ii = 1:Para["nTime"]
        for kk = 1:Para["nCon_conc"]
            if Con_conc[1,ii,kk] >= Para["Con_concMin"][kk]
                dHam_conc_abs[1,ii,kk] = abs(max(dHam_conc[1,ii,kk],0))
            elseif Con_conc[1,ii,kk] <= Para["Con_concMax"][kk]
                dHam_conc_abs[1,ii,kk] = abs(min(dHam_conc[1,ii,kk],0))
            else
                dHam_conc_abs[1,ii,kk] = abs(dHam[1,ii,kk])
            end
        end
    end

    for ii = 1:Para["nTime"]
        for jj = 1:Para["nVintage"]
            for kk = 1:Para["nCon_dist"]
                if Con_dist[jj,ii,kk] <= Para["Con_distMin"][kk]
                    dHam_dist_abs[jj,ii,kk] = abs(max(dHam_dist[jj,ii,kk],0))
                elseif Con_dist[jj,ii,kk] >= Para["Con_distMax"][kk]
                    dHam_dist_abs[jj,ii,kk] = abs(min(dHam_dist[jj,ii,kk],0))
                else
                    dHam_dist_abs[jj,ii,kk] = abs(dHam_dist[jj,ii,kk])
                end
            end
        end
    end 

    if Para["nCon_conc"] > 0
        dHam_conc_absValue = maximum(dHam_conc_abs)
    else
        dHam_conc_absValue = 0.0
    end
    if Para["nCon_dist"] > 0
        dHam_dist_absValue = maximum(dHam_dist_abs)
    else
        dHam_dist_absValue = 0.0
    end

    return max(dHam_conc_absValue,dHam_dist_absValue)
end