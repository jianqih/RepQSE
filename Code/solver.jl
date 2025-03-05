function solveHLwCtyOpen_E6(fund, dist, bord, bordc, NN; param = param)
    @unpack LL,alpha,sigma,LLwest, LLeast = param
    
    # extract location characteristics from fundamentals matrix;
    a=fund[:,1]; 
    H=fund[:,2]; 
    Iwest=fund[:,3]; 
    Ieast=fund[:,4];
    
    # Initialization based on a symmetric allocation;
    L_i = ones(NN) .* (LL/NN);
    w_i = ones(NN);
    
    # trade cost (with 1-sigma)
    dd = (dist.*bord.*bordc).^(1-sigma);
    
    wind = findall(Iwest .==1) # west index 
    eind = findall(Ieast .==1) # east index 
    for _ in range(1,200000)
        # trade share
        pwmat = (L_i .* (a .^ (sigma-1)) .* (w_i .^ (1-sigma))) * ones(1,NN)
        nummat = dd .* pwmat # numerator
        denom = sum(nummat, dims = 1) # donominator
        denommat = ones(NN) * denom
        tradesh = nummat ./ denommat # trade share 
        
        # income equals expenditure equilibrium condition:
        income = w_i .* L_i
        expend = tradesh * income 

        # domestic trade share 
        dtradesh = diag(tradesh)

        # pop mobility equilibrium condition:
        num = ((a .^ alpha) .* (H .^ (1-alpha)) .* (dtradesh .^ (-alpha/(sigma-1)))) .^ ((sigma-1)/((sigma*(1-alpha))-1)) #\lambda_n

        L_e = zeros(NN)
        L_e[wind] = num[wind]./sum(num[wind])
        L_e[eind] = num[eind]./sum(num[eind])
        L_e[wind] = L_e[wind] .* LLwest
        L_e[eind] = L_e[eind] .* LLeast

        # convergence criterion;
        income_r = round.(income .* (10^6))
        expend_r = round.(expend .* (10^6))
        L_i_r = round.(L_i .* (10^6))
        L_e_r = round.(L_e .* (10^6))

        
        if (income_r == expend_r) & (L_i_r == L_e_r)
            println("Convergence Achieved")
            return w_i,L_i,tradesh,dtradesh
        else
            w_e = w_i .* (expend ./ income) .^ (1/(sigma-1))
            w_i = (0.25 .* w_e) + (0.75 .* w_i)
            L_i = (0.25 .* L_e) + (0.75 .* L_i)
            # Normalization;
            # Choose geometric mean wage in West as numeraire;
            w_i[wind] = w_i[wind] ./ geomean(w_i[wind])
        end
    end
end