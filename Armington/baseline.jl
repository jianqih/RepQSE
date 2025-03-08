function solve_armington_eq(sigma, tau, A, L, a ,tol = 1e-5, damp = .1)
    w, lambda = ones(size(A)),zeros(size(a)) #initial guess
    wage_error = 1e5
    while wage_error > tol
        denominator2, wL = zeros(size(A)), zeros(size(A))
        for k in eachindex(A)
            denominator2 .+= a[k,:] .* (tau[k,:]*w[k] / A[k]).^(1 .- sigma)
        end
        for i in eachindex(A)
            wL .+= a[:,i] .* (tau[:,i] .* w ./ A).^(1 .- sigma) .* w[i] .* L[i] ./ denominator2[i]
        end
        wnew = wL ./ L;
        wage_error = maximum(abs.(wnew - w)./ w)
        w = damp .* wnew .+ (1-damp) .* w

    end
    for o in eachindex(A), d in eachindex(A)
        lambda[o,d] = a[o,d] * (tau[o,d]*w[o]/A[o])^(1-sigma)/ sum(a[:,d] .* (tau[:,d] .* w[:] ./ A[:]).^(1 - sigma))
    end
    return w,lambda
end

# Productivity shock to region 1
tau = [1. 5.;5. 1.] # symmetric matrix
a = [1. 1.;1. 1.]
A = [1.0,1.0]
L = [1.0,1.0]
sigma = 2;
w,lambda = solve_armington_eq(sigma,tau,A,L,a)



w, lambda = solve_armington_eq(sigma, tau, A, L, a);
display(w)
display(lambda)

# autarky 
tau_aut = [1. 1e9; 1e9 1.]; 
a = [1. 1.; 1. 1.];
A = [1., 1.];
L = [1., 1.];
sig = 2.;
w, lambda_aut = solve_armington_eq(sig, tau_aut, A, L, a);
display(w)
display(lambda_aut) # never trade


# New trading partner
tau_ntp = [1. 2. 3.; 3. 1. 2.; 1e9 5. 1.]; 
a_ntp = [1. 1. 1.; 1. 1. 1.; 1. 1. 1.];
A_ntp = [1., 1., 1.];
L_ntp = [1., 1., 1.];
w, lambda_ntp = solve_armington_eq(sig,tau_ntp, A_ntp, L_ntp, a_ntp);
display(w)
display(lambda_ntp)


