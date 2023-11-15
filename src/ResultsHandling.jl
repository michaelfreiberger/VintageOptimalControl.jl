
"""
    AssignResults(Results::Dict, 
                  Stat_conc::Array{Float64,3}, 
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

Assign the control and state variables to the Results dictionary
"""
function AssignResults(Results::Dict,
                       Stat_conc::Array{Float64,3},
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

    Results["Stat_conc"] = Stat_conc
    Results["Stat_dist"] = Stat_dist
    Results["Stat_agg"] = Stat_agg
    Results["Con_conc"] = Con_conc
    Results["Con_dist"] = Con_dist
    Results["CoStat_conc"] = CoStat_conc
    Results["CoStat_dist"] = CoStat_dist
    Results["CoStat_agg"] = CoStat_agg
    Results["dHam_conc"] = dHam_conc
    Results["dHam_dist"] = dHam_dist
    Results["Para"] = Para
    Results["ObjectiveValue"] = ObjectValue( Stat_conc,Stat_dist,Stat_agg,Con_conc, Con_dist,Para)
    return
end

"""
    PlotResults(Results2::Dict;
                SavePlot=false,
                Display=true,
                sizeX=600,sizeY=400)

Basic plots of all controls + gradient and state + costate variables
- SavePlot -> Indicator whether plots should be saved to location specifice in Para["SavePlotPath"]
- Display  -> Should the plots be displayed or supressed.
- sizeX, sizeY -> dimensions of the plots. Default = 600x400
"""
function PlotResults(Results2::Dict;SavePlot=false,Display=true,sizeX=600,sizeY=400)

    Results = deepcopy(Results2)
    Stat_conc = Results["Stat_conc"]
    Stat_dist = Results["Stat_dist"]
    Stat_agg = Results["Stat_agg"]
    Con_conc = Results["Con_conc"]
    Con_dist = Results["Con_dist"]
    CoStat_conc = Results["CoStat_conc"]
    CoStat_dist = Results["CoStat_dist"]
    CoStat_agg = Results["CoStat_agg"]
    dHam_conc = Results["dHam_conc"]
    dHam_dist = Results["dHam_dist"]
    Para = Results["Para"]

    for ii = 1:Para["nTime"]
        for kk = 1:Para["nCon_conc"]
            if Con_conc[1,ii,kk] <= Para["Con_concMin"][kk]
                dHam_conc[1,ii,kk] = max(dHam_conc[1,ii,kk],0)
            elseif Con_conc[1,ii,kk] >= Para["Con_concMax"][kk]
                dHam_conc[1,ii,kk] = min(dHam_conc[1,ii,kk],0)
            end
        end
    end

    for ii = 1:Para["nTime"]
        for jj = 1:Para["nVintage"]
            for kk = 1:Para["nCon_dist"]
                if Con_dist[jj,ii,kk] <= Para["Con_distMin"][kk]
                    dHam_dist[jj,ii,kk] = max(dHam_dist[jj,ii,kk],0)
                elseif Con_dist[jj,ii,kk] >= Para["Con_distMax"][kk]
                    dHam_dist[jj,ii,kk] = min(dHam_dist[jj,ii,kk],0)
                end
            end
        end
    end 

    
    PLOTS = Dict();

    nPlotHelp = min(Para["nVintage"],Para["nVintagePlot"])
    nPlotHelp2 = floor(Int,nPlotHelp * Para["T"]/Para["ω"])
    index = unique(round.(Int,range(1,Para["nVintage"],length=nPlotHelp)))
    index2 = unique(round.(Int,range(-Para["nVintage"],Para["nTime"],length=nPlotHelp2)))
    
    nPlotHelp = length(index)
    nPlotHelp2 = length(index2)
    
    ylab = string("")

    if SavePlot == true
        pyplot();
        StateConcLabels = Para["StateConcLabelsLatex"]
        StateDistLabels = Para["StateDistLabelsLatex"]
        StateAggLabels = Para["StateAggLabelsLatex"]
        CoStateConcLabels = Para["CoStateConcLabelsLatex"]
        CoStateDistLabels = Para["CoStateDistLabelsLatex"]
        CoStateAggLabels = Para["CoStateAggLabelsLatex"]
        ControlConcLabels = Para["ControlConcLabelsLatex"]
        ControlDistLabels = Para["ControlDistLabelsLatex"]
    else
        gr();
        StateConcLabels = Para["StateConcLabelsSimple"]
        StateDistLabels = Para["StateDistLabelsSimple"]
        StateAggLabels = Para["StateAggLabelsSimple"]
        CoStateConcLabels = Para["CoStateConcLabelsSimple"]
        CoStateDistLabels = Para["CoStateDistLabelsSimple"]
        CoStateAggLabels = Para["CoStateAggLabelsSimple"]
        ControlConcLabels = Para["ControlConcLabelsSimple"]
        ControlDistLabels = Para["ControlDistLabelsSimple"]
    end

	colors = get(colorschemes[:rainbow_bgyr_35_85_c72_n256],range(0,stop=1,length=nPlotHelp))
    colors2 = get(colorschemes[:rainbow_bgyr_35_85_c72_n256],range(0,stop=1,length=nPlotHelp2))

    #------------------------------------------------------------------------------------
    #   Plot of non-distributed Controls
    PLOTS["ControlsConc"] = plot(Para["tmesh"],Con_conc[1,:,:],lw=2,title = "Concentrated Controls",legend=:best,label=reshape(ControlConcLabels,1,:),size = (sizeX, sizeY))

    #   Plot of Gradient w.r.t. non-distributed Controls
    PLOTS["GradientsConc"] = plot(Para["tmesh"],dHam_conc[1,:,:],lw=2,title = "Gradients of concentrated controls",legend=:best,label=reshape(ControlConcLabels,1,:),size = (sizeX, sizeY))

    #   Plot of non-distributed States
    PLOTS["StatesConc"] = plot(Para["tmesh"],Stat_conc[1,:,:],lw=2,title = "Concentrated States",legend=:best,label=reshape(StateConcLabels,1,:),size = (sizeX, sizeY))

    #   Plot of non-distributed Co-states
    PLOTS["CoStatesConc"] = plot(Para["tmesh"],CoStat_conc[1,:,:],lw=2,title = "Concentrated CoStates",legend = :best,label=reshape(CoStateConcLabels,1,:),size = (sizeX, sizeY))

    #   Plot of non-distributed States
    PLOTS["StatesAgg"] = plot(Para["tmesh"],Stat_agg[1,:,:],lw=2,title = "Aggregated States",legend=:best,label=reshape(StateAggLabels,1,:),size = (sizeX, sizeY))

    #   Plot of non-distributed Co-states
    PLOTS["CoStatesAgg"] = plot(Para["tmesh"],CoStat_agg[1,:,:],lw=2,title = "Aggregated CoStates",legend = :best,label=reshape(CoStateAggLabels,1,:),size = (sizeX, sizeY))


    #------------------------------------------------------------------------------------
    #   Plot of distributed controls
    for k = 1:Para["nCon_dist"]
        PLOTS[string("ControlsDist",k)] = plot(title = string("Distributed Control ",ControlDistLabels[k], "| Cross-sectional"),size = (sizeX, sizeY),ylabel = ylab)
        for ii=1:nPlotHelp
            PLOTS[string("ControlsDist",k)] = plot!(Para["tmesh"],Con_dist[index[ii],:,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("ControlsDist",k,"_V2")] = plot(title = string("Distributed Control ",ControlDistLabels[k], "| Longitudinal"),size = (sizeX, sizeY),ylabel = ylab)
        for ii=1:nPlotHelp2
            if index2[ii] < 0
                PLOTS[string("ControlsDist",k,"_V2")] = plot!(Para["tmesh"][1:Para["nVintage"]+1+index2[ii]],diag(Con_dist[-index2[ii]:end,1:Para["nVintage"]+1+index2[ii],k]),lw = 2,color=colors2[ii],label=false)
            elseif index2[ii] > Para["nTime"] - Para["nVintage"]
                PLOTS[string("ControlsDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]:end],diag(Con_dist[1:Para["nTime"]+1-index2[ii],index2[ii]:end,k]),lw = 2,color=colors2[ii],label=false)
            else
                PLOTS[string("ControlsDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]+1:index2[ii]+Para["nVintage"]],diag(Con_dist[:,index2[ii]+1:index2[ii]+Para["nVintage"],k]),lw = 2,color=colors2[ii],label=false)
            end
        end
        PLOTS[string("ControlsDist",k,"_V3")] = plot(title = string("Distributed Control ",ControlDistLabels[k], "| Heatmap"),size = (sizeX, sizeY),ylabel = ylab)
        PLOTS[string("ControlsDist",k,"_V3")] = heatmap!(Para["tmesh"],Para["amesh"],Con_dist[:,:,k])
    end

    #------------------------------------------------------------------------------------
    #   Plot of distributed gradient
    for k = 1:Para["nCon_dist"]
        PLOTS[string("GradDist",k)] = plot(title = string("Gradient of distributed Control ",ControlDistLabels[k], "| Cross-sectional"),size = (sizeX, sizeY))
        for ii=1:nPlotHelp
            PLOTS[string("GradDist",k)] = plot!(Para["tmesh"],dHam_dist[index[ii],:,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("GradDist",k,"_V2")] = plot(title = string("Gradient of distributed Control ",ControlDistLabels[k], "| Longitudinal"),size = (sizeX, sizeY),ylabel = ylab)
        for ii=1:nPlotHelp2
            if index2[ii] < 0
                PLOTS[string("GradDist",k,"_V2")] = plot!(Para["tmesh"][1:Para["nVintage"]+1+index2[ii]],diag(dHam_dist[-index2[ii]:end,1:Para["nVintage"]+1+index2[ii],k]),lw = 2,color=colors2[ii],label=false)
            elseif index2[ii] > Para["nTime"] - Para["nVintage"]
                PLOTS[string("GradDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]:end],diag(dHam_dist[1:Para["nTime"]+1-index2[ii],index2[ii]:end,k]),lw = 2,color=colors2[ii],label=false)
            else
                PLOTS[string("GradDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]+1:index2[ii]+Para["nVintage"]],diag(dHam_dist[:,index2[ii]+1:index2[ii]+Para["nVintage"],k]),lw = 2,color=colors2[ii],label=false)
            end
        end
        PLOTS[string("GradDist",k,"_V3")] = plot(title = string("Gradient of distributed Control ",ControlDistLabels[k], "| Heatmap"),size = (sizeX, sizeY),ylabel = ylab)
        PLOTS[string("GradDist",k,"_V3")] = heatmap!(Para["tmesh"],Para["amesh"],dHam_dist[:,:,k])
    end

    #------------------------------------------------------------------------------------
    #   Plot of distributed states
    for k = 1:Para["nStat_dist"]
        PLOTS[string("StatDist",k)] = plot(title = string("Distributed State ",StateDistLabels[k], "| Cross-sectional"),size = (sizeX, sizeY))
        for ii=1:nPlotHelp
            PLOTS[string("StatDist",k)] = plot!(Para["tmesh"],Stat_dist[index[ii],:,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("StatDist",k,"_V2")] = plot(title = string("Distributed State ",StateDistLabels[k], "| Longitudinal"),size = (sizeX, sizeY),ylabel = ylab)
        for ii=1:nPlotHelp2
            if index2[ii] < 0
                PLOTS[string("StatDist",k,"_V2")] = plot!(Para["tmesh"][1:Para["nVintage"]+1+index2[ii]],diag(Stat_dist[-index2[ii]:end,1:Para["nVintage"]+1+index2[ii],k]),lw = 2,color=colors2[ii],label=false)
            elseif index2[ii] > Para["nTime"] - Para["nVintage"]
                PLOTS[string("StatDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]:end],diag(Stat_dist[1:Para["nTime"]+1-index2[ii],index2[ii]:end,k]),lw = 2,color=colors2[ii],label=false)
            else
                PLOTS[string("StatDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]+1:index2[ii]+Para["nVintage"]],diag(Stat_dist[:,index2[ii]+1:index2[ii]+Para["nVintage"],k]),lw = 2,color=colors2[ii],label=false)
            end
        end
        PLOTS[string("StatDist",k,"_V3")] = plot(title = string("Distributed State ",StateDistLabels[k], "| Heatmap"),size = (sizeX, sizeY),ylabel = ylab)
        PLOTS[string("StatDist",k,"_V3")] = heatmap!(Para["tmesh"],Para["amesh"],Stat_dist[:,:,k])
    end

    #------------------------------------------------------------------------------------
    #   Plot of distributed costates
    for k = 1:Para["nStat_dist"]
        PLOTS[string("CoStatDist",k)] = plot(title = string("Distributed CoState ",CoStateDistLabels[k], "| Cross-sectional"),size = (sizeX, sizeY))
        for ii=1:nPlotHelp
            PLOTS[string("CoStatDist",k)] = plot!(Para["tmesh"],CoStat_dist[index[ii],:,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("CoStatDist",k,"_V2")] = plot(title = string("Distributed CoState ",CoStateDistLabels[k], "| Longitudinal"),size = (sizeX, sizeY),ylabel = ylab)
        for ii=1:nPlotHelp2
            if index2[ii] < 0
                PLOTS[string("CoStatDist",k,"_V2")] = plot!(Para["tmesh"][1:Para["nVintage"]+1+index2[ii]],diag(CoStat_dist[-index2[ii]:end,1:Para["nVintage"]+1+index2[ii],k]),lw = 2,color=colors2[ii],label=false)
            elseif index2[ii] > Para["nTime"] - Para["nVintage"]
                PLOTS[string("CoStatDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]:end],diag(CoStat_dist[1:Para["nTime"]+1-index2[ii],index2[ii]:end,k]),lw = 2,color=colors2[ii],label=false)
            else
                PLOTS[string("CoStatDist",k,"_V2")] = plot!(Para["tmesh"][index2[ii]+1:index2[ii]+Para["nVintage"]],diag(CoStat_dist[:,index2[ii]+1:index2[ii]+Para["nVintage"],k]),lw = 2,color=colors2[ii],label=false)
            end
        end
        PLOTS[string("CoStatDist",k,"_V3")] = plot(title = string("Distributed CoState ",CoStateDistLabels[k], "| Heatmap"),size = (sizeX, sizeY),ylabel = ylab)
        PLOTS[string("CoStatDist",k,"_V3")] = heatmap!(Para["tmesh"],Para["amesh"],CoStat_dist[:,:,k])
    end

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if SavePlot == true
        mkpath(string(Para["SavePlotPath"],"/BasePlots"))
        for str in sort(collect(keys(PLOTS)))
            savefig(PLOTS[str],string(Para["SavePlotPath"],"/BasePlots/",str,".pdf"))
        end
    end

    if Display == true
        if Para["nCon_conc"] > 0
            display(plot(PLOTS["ControlsConc"],PLOTS["GradientsConc"],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
        end
        if Para["nStat_conc"] > 0
            display(plot(PLOTS["StatesConc"],PLOTS["CoStatesConc"],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
        end
        if Para["nStat_agg"] > 0
            display(plot(PLOTS["StatesAgg"],PLOTS["CoStatesAgg"],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
        end

        for k=1:Para["nCon_dist"]
            display(plot(PLOTS[string("ControlsDist",k)],PLOTS[string("GradDist",k)],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
            display(plot(PLOTS[string("ControlsDist",k,"_V2")],PLOTS[string("GradDist",k,"_V2")],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
            display(plot(PLOTS[string("ControlsDist",k,"_V3")],PLOTS[string("GradDist",k,"_V3")],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
        end
        for k=1:Para["nStat_dist"]
           display(plot(PLOTS[string("StatDist",k)],PLOTS[string("CoStatDist",k)],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
           display(plot(PLOTS[string("StatDist",k,"_V2")],PLOTS[string("CoStatDist",k,"_V2")],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
           display(plot(PLOTS[string("StatDist",k,"_V3")],PLOTS[string("CoStatDist",k,"_V3")],layout = grid(2, 1, heights=[0.5 ,0.5]),size=[sizeX,sizeY*2]))
        end
    end
    return PLOTS
end


"""
    SaveResults(Results::Dict,filepath)

Save the Results Dictionary to filepath
"""
function SaveResults(Results::Dict,filepath)
    Results["Para"]["SavePlotPath"] = filepath
    save(string(filepath,".jld2"),Results)
end

"""
    LoadResults(filepath)

Load the Results Dictionary from the filepath
"""
function LoadResults(filepath)
    Results = load(string(filepath,".jld2"))
    return Results
end