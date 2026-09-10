using Plots, Revise, BenchmarkTools

includet("../src/spinglass.jl")
includet("../src/methods.jl")
includet("../src/CrudeSpinGlass.jl")

# Parameters
size = 1000
h = 1.0

# Constructor
# sg = SpinGlass(size=size, J=curie_weiss(size), h=random_field(size, h))

# Test
function test_rfim()
  sweeps::Int = 0
  sg = SpinGlass(size=size, J=curie_weiss(size), h=random_field(size, h))

  e_old = hamiltonian(sg)
  println(e_old)

  while true
    # make use of deterministic greedy algorithm
    meanfield_rfim!(sg)
    sweeps += 1
    e_new = hamiltonian(sg)

    if abs(e_old - e_new) < 1/sg.size || sweeps >= 100
      break
    end
    e_old = e_new
  end
  println(e_old, " with $sweeps sweeps")
  return e_old
end

v = zeros(15)

for i in 1:15
  v[i] = test_rfim()
end

println("")
println("Mean: ", sum(v)/15)
