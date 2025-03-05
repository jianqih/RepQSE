# Author: Jianqi Huang
using Distributions
using LinearAlgebra
using Statistics
using StatsBase
using CSV,DataFrames,DataFramesMeta
using Parameters
using Plots
include("solver.jl")
include("prices_and_welf_calculation.jl")
# fixed parameters
N = 30;
NN = N*N;
ltd = range(0,4,N)';
lgd = range(0,4,N);

# transport costs for each point on grid. 
tt = 1;
tau = zeros(N);
tau .= tt;

# compute weighted distance using transport costs;
dist = Matrix{Float64}(undef, NN, NN)

# #  using ImageMorphology
# function my_dist_transform(X)
#     return distance_transform(feature_transform(X))
# end

## Same done manually
function my_dist_transform(X)
    N = size(X,1)
    dist = Matrix{Float64}(undef, N, N)
    x, y = Tuple(findall(X .== 1)[1])
    for i in 1:N, j in 1:N
        dist[i,j] = sqrt((i-x)^2 + (j-y)^2)
    end
    return dist
end

# Distance weighting does not work yet
for z in 1:NN
    seed = falses(N,N)
    seed[z] = true
    # w = vectotuple(vec(τ))
    temp = my_dist_transform(seed) # graydist(τ,seed,'quasi-euclidean')
    dist[z,:] = temp[:]
end

dist[diagind(dist)] .= 1 
dist = dist.^0.33

Iwest=zeros(N,N);
Ieast=zeros(N,N);
Iwest[:,1:Int(N/2)] .=1;
Ieast[:,Int(N/2)+1:N] .=1;
Iwest=vec(Iwest);
Ieast=vec(Ieast);

# Broader
bord = fill(2.0, (NN,NN))  
bord[diagind(bord)] .= 1

# Broader frictions between two countries 
bordcty = ones(NN,NN)
bordcty[Iwest.==1, Ieast.==1] .= 2
bordcty[Ieast.==1, Iwest.==1] .= 2

# counterfactual broader
cbord = ones(NN,NN)
cbordcty = ones(NN,NN)

@with_kw struct model
    alpha = 0.75
    sigma = 5
    LL = 153889 # US civilian labor force 2010 (Statistical Abstract, millions)
    LLwest = (sum(Iwest)/(sum(Iwest)+sum(Ieast)))*LL 
    LLeast = (sum(Ieast)/(sum(Iwest)+sum(Ieast)))*LL
end
param = model()

a = CSV.read("Code/a.csv", DataFrame; header = false) # Read a.csv
a = exp.(Matrix(a))
a[Iwest.==1] .= a[Iwest.==1] ./ geomean(a[Iwest.==1])
a[Ieast.==1] .= a[Ieast.==1] ./ geomean(a[Ieast.==1]) 

# Land area
H = 100 * ones(NN)   

F = 1 # fixed labor cost
# matrix of fundamentals.
fund = zeros(NN,4)
fund[:,1] = a
fund[:,2] = H
fund[:,3] = Iwest
fund[:,4] = Ieast

println("Start Wage and Population Convergence ")
w,L,tradesh,dtradesh = solveHLwCtyOpen_E6(fund,dist,bord,bordcty,NN)
# Price index 
P = Hpindex(fund,L,w,dtradesh)
realwage = Hrealw(fund,L,tradesh)
r = Hlandpirce(fund,L,w)

# Counterfactual eliminating border frictions between countries 
cw,cL,ctradesh,cdtradesh = solveHLwCtyOpen_E6(fund,dist,bord,cbordcty,NN)
# Price index 
cP = Hpindex(fund,cL,cw,cdtradesh)
crealwage = Hrealw(fund,cL,ctradesh)
cr = Hlandpirce(fund,cL,cw)
cwelfgain = Hwelfaregains(ctradesh,tradesh,cL,L)
cwelfgain = round.(cwelfgain, digits = 4)
unique(cwelfgain)

# Counterfactual eliminating border frictions between grid points
ccw,ccL,cctradesh,ccdtradesh = solveHLwCtyOpen_E6(fund,dist,cbord,bordcty,NN)
# Price index 
ccP = Hpindex(fund,ccL,ccw,ccdtradesh)
ccrealwage = Hrealw(fund,ccL,cctradesh)
ccr = Hlandpirce(fund,ccL,ccw)
ccwelfgain = Hwelfaregains(cctradesh,tradesh,ccL,L)
ccwelfgain = round.(ccwelfgain, digits = 4)
unique(ccwelfgain)

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