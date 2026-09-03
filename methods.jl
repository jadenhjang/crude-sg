# Randomize spins
function reinitialize(sg::SpinGlass)
  sg.spins[:] = floor.(Int64, rand(Float64, sg.size)*2)*2 .- 1
  return sg
end

# Hamiltonian w/ external field h
function hamiltonian(sg::SpinGlass, h::Any)
  x::Float64 = 0

  for i in 1:sg.size
    for j in 1:sg.size
      x -= sg.spins[i]*sg.spins[j]*sg.config[i, j]/2 # note all self-couplings are zero
    end
    x -= h*sg.spins[i]
  end

  return x
end


# Update Rule
function update(sg::SpinGlass)
  t = floor.(Int64, rand(1000)*sg.size) .+ 1
  # println("First 5 indices: ", t[1:5]) # df.head()

  for k in 1:t.size[1]
    aux::Float64 = 0
    for i in 1:sg.size
      aux += sg.spins[t[k]]*sg.spins[i]*sg.config[t[k], i]
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

  E_now = hamiltonian(sg, h)

  k = floor(Int64, rand()*sg.size)+1 # random spin selector
  sg.spins[k] *= -1
  E_flip = hamiltonian(sg, h)

  r = rand()  # supposed to be [0.1] but is [0,1)
  if r > exp(beta*(E_now-E_flip))
    sg.spins[k] *= -1  # restore spin, flipped inequality
  end

end

"""
Testing Functions
"""
# Bias Checking
function bias_check(sg::SpinGlass)
  return sum(sg.spins) / sqrt(sg.size)  # random walk normalization
end

function check2(sg::SpinGlass)
  sg = reinitialize(sg)
  println("hamiltonian_i = ", hamiltonian(sg, 0))
  sg2 = update(sg)
  println("hamiltonian_f = ", hamiltonian(sg2, 0))
  println("change in hamiltonian = ", (hamiltonian(sg2, 0)-hamiltonian(sg, 0)))
end

function MCMC_test(sg::SpinGlass, beta::Float64, h::Float64)
  #Q1
  sg.spins[:] = ones(Int64, sg.size)
  # beta::Float64 = 1.2
  # h::Float64 = 0

  t = 1:(sg.size*100)
  y::Array{Float64} = []
  for _ in t
    MCMC(sg, beta, h)
    m = sum(sg.spins)/sg.size
    append!(y, m)
  end
  return t, y
end
