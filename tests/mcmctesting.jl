# using Base: julia_exename

# A test script to solve EPFL Ch1 1.1.4 Monte-Carlo-Markov-Chain via the Metropolis-Hasting Algorithm
function MCMC_test(sg::SpinGlass, beta::Float64, h::Float64)
  sg.spins[:] = ones(Int64, sg.size)

  t = 1:(sg.size*100)
  y::Array{Float64} = []
  for _ in t
    MCMC(sg, beta, h)
    m = sum(sg.spins)/sg.size
    append!(y, m)
  end
  return t, y
end

N = [10, 50, 100, 200, 100]
params = [[1.2, 0.0], [1.2, -0.1]]

solutions = Array{NamedTuple}(undef, length(N), length(params))

for (i, n) in enumerate(N)
  sg = SpinGlass(size=n)

  for (j, p) in enumerate(params)
    time, magnetization = MCMC_test(sg, p[1], p[2])
    solutions[i, j] = (N=n, t=time, y=magnetization)
  end
end

plt = plot(
  x->-0.1,
  label="h",
  title="Spontaneous Magnetization",
  legend=:topright,
  xlabel="time",
  ylabel="magnetization", xlims=(0, 1000),
)

for i in 2:5
  plot!(plt, sols[i, 2].y[1:5000], label=sols[i, 2].N)
end
