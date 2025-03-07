
begin
    amat = reshape(log.(a),N,N)
    Lmat = reshape(log.(L),N,N)
    wmat = reshape(log.(w),N,N)
    rmat = reshape(log.(r),N,N)
    Pmat = reshape(log.(P),N,N)
    heatmap(amat,xlabel = "Lon",ylabel = "Lat",title = "Log Productivtiy", fontsize = 12, c = :viridis)
    savefig("Figures/a.pdf")
    p = plot(layout = (2,2))
    heatmap!(p[1],Lmat,xlabel = "Lon",ylabel = "Lat",title = "Panel A: Log Population", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p[2],wmat,xlabel = "Lon",ylabel = "Lat",title = "Panel B: Log Wages",titlefontsize=10,  fontsize = 8, c = :viridis)
    heatmap!(p[3],rmat,xlabel = "Lon",ylabel = "Lat",title = "Panel C: Log Land Prices", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p[4],Pmat,xlabel = "Lon",ylabel = "Lat",title = "Panel D: Log Price Index",titlefontsize=10,  fontsize = 8, c = :viridis)
    savefig("Figures/initial.pdf")
end

begin
    dL = cL./L 
    dw = cw./w
    dr = cr./r 
    dP = cP./P
    dLmat = reshape(log.(dL),N,N)
    dwmat = reshape(log.(dw),N,N)
    drmat = reshape(log.(dr),N,N)
    dPmat = reshape(log.(dP),N,N)
    p2 = plot(layout = (2,2))
    heatmap!(p2[1],dLmat,xlabel = "Lon",ylabel = "Lat",title = "Panel A: Log Relative Population", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p2[2],dwmat,xlabel = "Lon",ylabel = "Lat",title = "Panel B: Log Relative Wages", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p2[3],drmat, xlabel = "Lon",ylabel = "Lat",title = "Panel C: Log Relative Land Prices",titlefontsize=10,  fontsize = 8, c = :viridis)
    heatmap!(p2[4],dPmat,xlabel = "Lon",ylabel = "Lat",title = "Panel D: Log Relative Price Index", titlefontsize=10, fontsize = 8, c = :viridis)
    savefig("Figures/c.pdf")
    # multi-panel figures levels 
    p3 = plot(layout = (2,2))
    heatmap!(p3[1],Lmat + dLmat,xlabel = "Lon",ylabel = "Lat",title = "Panel A: Log Population", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p3[2],wmat + dwmat,xlabel = "Lon",ylabel = "Lat",title = "Panel B: Log Wages", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p3[3],rmat + drmat, xlabel = "Lon",ylabel = "Lat",title = "Panel C: Log Land Prices", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p3[4],Pmat + dPmat,xlabel = "Lon",ylabel = "Lat",title = "Panel D: Log Price Index", titlefontsize=10,  fontsize = 8, c = :viridis)
    savefig("Figures/c_lev.pdf")
end

begin
    ddL = ccL./L 
    ddw = ccw./w
    ddr = ccr./r 
    ddP = ccP./P
    ddLmat = reshape(log.(ddL),N,N)
    ddwmat = reshape(log.(ddw),N,N)
    ddrmat = reshape(log.(ddr),N,N)
    ddPmat = reshape(log.(ddP),N,N)
    p4 = plot(layout = (2,2))
    heatmap!(p4[1], ddLmat, xlabel = "Lon", ylabel = "Lat",title = "Panel A: Log Relative Population", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p4[2],ddwmat, clim = (-0.05,0.05), xlabel = "Lon",ylabel = "Lat",title = "Panel B: Log Relative Wages (Truncated)", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p4[3],ddrmat, xlabel = "Lon", ylabel = "Lat",title = "Panel C: Log Relative Land Prices", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p4[4],ddPmat, clim = (-1.35,-1.15), xlabel = "Lon",ylabel = "Lat",title = "Panel D: Log Relative Price Index (Truncated)", titlefontsize=10, fontsize = 8, c = :viridis)
    savefig("Figures/cc.pdf")
    # multi-panel figures levels 
    p5 = plot(layout = (2,2))
    heatmap!(p5[1],Lmat + ddLmat, xlabel = "Lon",ylabel = "Lat",title = "Panel A: Log Population", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p5[2],wmat + ddwmat, xlabel = "Lon",ylabel = "Lat",title = "Panel B: Log Wages", titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p5[3],rmat + ddrmat, xlabel = "Lon",ylabel = "Lat",title = "Panel C: Log Land Prices",titlefontsize=10, fontsize = 8, c = :viridis)
    heatmap!(p5[4],Pmat + ddPmat, xlabel = "Lon",ylabel = "Lat",title = "Panel D: Log Price Index",titlefontsize=10, fontsize = 8, c = :viridis)
    savefig("Figures/cc_lev.pdf")
end