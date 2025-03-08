w = [6.0, 3.0];
L = [0.3, 0.7];
lambda = [0.8 0.2; 0.2 0.8];
X = (w .* L)' .* lambda;

Ahat = [0.5, 1.];
sigma = 2

function solve_armington_mig_exact_hat(X,lambda,w,L,Ahat,sigma,tol = 1e-5,damp=0.1)
    what, Lhat, labor_error,wage_error = ones(size(Ahat)), ones(size(Ahat)), 1e5, 1e5;
    while max(labor_error,wage_error) > tol 
        denominator, what_new, Lhat_new = zeros(size(Ahat)), zeros(size(Ahat)), zeros(size(Ahat));
        for k in eachindex(Ahat)
            denominator .+= lambda[k,:] .* (what[k]/Ahat[k])^(1-sigma);
        end
        for i in eachindex(Ahat)
            what_new .+= (X[:,i] ./ (w .* L)) .* (what ./ Ahat).^(1-sigma) .* what[i] .* Lhat[i] ./ denominator[i] ./ Lhat;
        end
        Phat = vec((sum(lambda .* (what ./ Ahat).^(1 .- sigma), dims = 1)).^(1 ./(1 .- sigma)));
        Lhat_new = what ./ Phat ./ sum(L .* what ./ Phat);
        labor_error = maximum(abs.(Lhat_new - Lhat) ./ Lhat);
        wage_error = maximum(abs.(what_new - what) ./ what);
        what = damp .* what_new .+ (1-damp) .* what;
        Lhat = damp .* Lhat_new .+ (1-damp) .* Lhat;
    end
    lambda_hat = (what ./ Ahat).^(1 .- sigma) ./ sum(lambda .* (what ./ Ahat).^(1 .- sigma), dims = 1);
    Phat = vec((sum(lambda .* (what ./ Ahat).^(1 .- sigma), dims = 1)).^(1 ./(1 .- sigma)));
    return what, lambda_hat, Phat, Lhat
end

Ahat = [0.5, 1.];
what, lambda_hat, Phat,Lhat = solve_armington_mig_exact_hat(X, lambda, w, L, Ahat, sigma);
what ./ Phat
Lhat