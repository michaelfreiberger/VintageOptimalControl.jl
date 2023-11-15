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