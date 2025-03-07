# Author: Jianqi Huang
using Distributions
using LinearAlgebra
using Statistics
using StatsBase
using CSV,DataFrames,DataFramesMeta
using Parameters
using Plots
include("solver2.jl")
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
using BenchmarkTools

w,L,tradesh,dtradesh = solveHLwCtyOpen_v2(fund,dist,bord,bordcty,NN)
# Price index 
P = Hpindex(fund,L,w,dtradesh)
realwage = Hrealw(fund,L,tradesh)
r = Hlandpirce(fund,L,w)

# Counterfactual eliminating border frictions between countries 
cw,cL,ctradesh,cdtradesh = solveHLwCtyOpen_v2(fund,dist,bord,cbordcty,NN)
# Price index 
cP = Hpindex(fund,cL,cw,cdtradesh)
crealwage = Hrealw(fund,cL,ctradesh)
cr = Hlandpirce(fund,cL,cw)
cwelfgain = Hwelfaregains(ctradesh,tradesh,cL,L)
cwelfgain = round.(cwelfgain, digits = 4)
unique(cwelfgain)

# Counterfactual eliminating border frictions between grid points
ccw,ccL,cctradesh,ccdtradesh = solveHLwCtyOpen_v2(fund,dist,cbord,bordcty,NN)
# Price index 
ccP = Hpindex(fund,ccL,ccw,ccdtradesh)
ccrealwage = Hrealw(fund,ccL,cctradesh)
ccr = Hlandpirce(fund,ccL,ccw)
ccwelfgain = Hwelfaregains(cctradesh,tradesh,ccL,L)
ccwelfgain = round.(ccwelfgain, digits = 4)
unique(ccwelfgain)

include("plots.jl")
