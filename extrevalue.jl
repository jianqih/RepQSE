using Distributions,Plots
mu1 = 0.0
beta1 = 1.0
mu2 = 1.0
beta2 = 1.0
mu3 =  0.0
beta3 = 2.0

gumbel_dist1 = Gumbel(mu1, beta1)
gumbel_dist2 = Gumbel(mu2, beta2)
gumbel_dist3 = Gumbel(mu3, beta3)

x = range(-2, 10, length=100)

pdf_values1 = pdf.(gumbel_dist1, x)
pdf_values2 = pdf.(gumbel_dist2, x)
pdf_values3 = pdf.(gumbel_dist3, x)


# Plot the distribution
plot(x, pdf_values1, 
    label="Gumbel PDF (μ=$mu1, β=$beta1)",
    xlabel="x",
    ylabel="Density",
    title="Gumbel Distribution",
    linewidth=2,
    color=:blue,
    legend=:topright)
plot!(x, pdf_values2, 
    label="Gumbel PDF (μ=$mu2, β=$beta2)",
    xlabel="x",
    ylabel="Density",
    title="Gumbel Distribution",
    linewidth=2,
    color=:red,
    legend=:topright)
plot!(x, pdf_values3, 
    label="Gumbel PDF (μ=$mu3, β=$beta3)",
    xlabel="x",
    ylabel="Density",
    title="Gumbel Distribution",
    linewidth=2,
    color=:green,
    legend=:topright)
