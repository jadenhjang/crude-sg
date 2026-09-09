@kwdef struct SpinGlass
  size::Int64
  spins::Array{Int64} = floor.(Int8, rand(Float64, size)*2)*2 .- 1
  config::Matrix{Float64} = curie_weiss(size)
end

# helper functions: generate coupling matrix J_{i,j}
function sherrington_kirkpatrick(n::Int64)
  A::Matrix{Float64} = randn(n, n) ./ sqrt(n)  # normalization factor
  J = (A + A') / 2  # symmetric (hermitian)

  for i in 1:n
    J[i, i] = 0  # zero the diagonal
  end

  return J
end

function curie_weiss(n::Int64)
  J = fill(1/n, n, n)
  for i in 1:n
    J[i, i] = 0.0
  end
  return J
end
