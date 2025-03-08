w = [1.0,1.0]
L = [1.0,1.0]
lambda = [0.8 0.2; 0.2 0.8]
X = (w .* L)' .* lambda
Ahat = [10.,1.]
sigma = 2.0;


function solve_armington_exact_hat(X, lambda, w, L, Ahat, sigma, tol = 1e-5, damp = .1)
    what = ones(size(Ahat))
    wage_error = 1e5
    while wage_error>tol
        denominator2, what_new = zeros(size(Ahat)), zeros(size(Ahat)) # start from zero
        for k in eachindex(Ahat)
            denominator2 .+= lambda[k,:] .* (what[k]/Ahat[k]).^(1-sigma)
        end # Loop over k
        for i in eachindex(Ahat)
            what_new .+= (X[:,i] ./ (w .* L)) .* (what ./ Ahat).^(1-sigma) .* what[i] ./denominator2[i]
        end # Loop over i
        wage_error = maximum(abs.(what_new - what)./what)
        what = damp .*  what_new .+ (1-damp) .* what
    end
    lambda_hat = (what ./ Ahat).^(1 .- sigma) ./ sum(lambda .* (what ./ Ahat).^(1 .- sigma), dims = 1)
    Phat = vec((sum(lambda .* (what ./ Ahat).^(1 .- sigma), dims = 1)).^(1 ./(1 .- sigma)))
    return what,lambda_hat,Phat
end

Ahat = [10., 1.];
what, lambda_hat, Phat = solve_armington_exact_hat(X, lambda, w, L, Ahat, sigma);
what ./ Phat
lambda_hat