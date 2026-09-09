using Plots, Revise, BenchmarkTools

includet("../src/spinglass.jl")
includet("../src/methods.jl")

N = [10, 50, 100, 200, 1000]  # 1000
params = [[1.2, 0.0], [1.2, -0.1]]

sols = Array{NamedTuple}(undef, length(N), length(params))

for (i, n) in enumerate(N)
  sg = SpinGlass(size=n)

  for (j, p) in enumerate(params)
    time, magnetization = MCMC_test(sg, p[1], p[2])
    sols[i, j] = (N=n, t=time, y=magnetization)
  end
end

paramset=2
titley::String = "Spontaneous magnetization at ß and h at $(params[paramset])"
plt = plot(
  title=titley,
  legend=:topright,
  xlabel="time/N",
  ylabel="magnetization",
)

# plot!(plt, sols[5, 1].y)
for i in 1:5
  plot!(plt, sols[i, paramset].t ./ sols[i, paramset].N, sols[i, paramset].y, label=sols[i, paramset].N)
end

savefig(plt, "plots/test.png")
