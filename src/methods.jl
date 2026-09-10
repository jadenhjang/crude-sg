include("spinglass.jl")

# Randomize spins
function reinitialize_spins(sg::SpinGlass)
  sg.spins[:] = floor.(Int64, rand(Float64, sg.size)*2)*2 .- 1
  return sg
end

# Hamiltonian w/ external field h
function hamiltonian(sg::SpinGlass)
  x::Float64 = 0

  for i in 1:sg.size  #  O(n^2)
    for j in 1:sg.size
      x -= sg.spins[i] * sg.spins[j] * sg.J[i, j] * 0.5 # note all self-couplings are zero
    end
    x -= sg.h[i]*sg.spins[i]
  end

  return x
end


# Update Rule
function update(sg_original::SpinGlass)
  sg = deepcopy(sg_original)
  t = floor.(Int64, rand(1000)*sg.size) .+ 1
  # println("First 5 indices: ", t[1:5]) # df.head()

  for k in 1:t.size[1]
    aux::Float64 = 0
    for i in 1:sg.size
      aux += sg.spins[t[k]]*sg.spins[i]*sg.J[t[k], i]
    end
    #sg.spins[k]+= floor(Int64,aux)
    if aux < 0
      sg.spins[t[k]] *= -1   # flip spin
    end
  end

  return sg
end

function MCMC(sg::SpinGlass, beta::Float64, h::Float64)
  #=
  Monte-Carlo-Markov-Chain, the Metropolis_Hastings Algorithm

  Statistical Physics Methods in Optimization and Learning
  by F. Krzakala and L. Zdeborova

  Exercise 1.1.4
  =#

  # E_now = hamiltonian(sg, h)

  k = floor(Int64, rand()*sg.size)+1 # random spin selector
  # sg.spins[k] *= -1
  # E_flip = hamiltonian(sg, h)

  f::Float64 = 0

  for j in 1:sg.size
    f += sg.J[k, j]*sg.spins[j]
  end

  ∆E = 2*sg.spins[k]*(f+h)
  r = rand()  # supposed to be [0.1] but is [0,1)

  if r < exp(-beta*(∆E))
    sg.spins[k] *= -1  # flip
  end

end

function meanfield_rfim!(sg::SpinGlass)
  # one "sweep" with random updates
  #TODO: Are random sweeps and ordered sweeps the same?
  for _ in 1:sg.size
    index::Int64 = rand(1:sg.size)
    if sg.spins[index] != sign(sg.h[index]+sum(sg.spins)/sg.size)
      sg.spins[index] *= -1
    end
  end
  return sg
end
"""
Tester Functions
"""
# Bias Checking
function bias_check(sg::SpinGlass)
  return sum(sg.spins) / sqrt(sg.size)  # random walk normalization
end

# Hamiltonian Minimization
function check2(sg::SpinGlass)
  sg = reinitialize(sg)
  println("hamiltonian_i = ", hamiltonian(sg))
  sg2 = update(sg)
  println("hamiltonian_f = ", hamiltonian(sg2))
  println("change in hamiltonian = $(hamiltonian(sg2)-hamiltonian(sg))")
end

# MCMC Metropolis Tester
function MCMC_test(sg::SpinGlass, beta::Float64, h::Float64)
  #Q1
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
