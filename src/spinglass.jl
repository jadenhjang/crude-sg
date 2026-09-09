@kwdef struct SpinGlass
  size::Int64
  spins::Array{Int64} = floor.(Int8, rand(Float64, size)*2)*2 .- 1
  config::Matrix{Float64} = curie_weiss(size)
end

# helper functions: generate coupling matrix J_{i,j}
function sherrington_kirkpatrick(n::Int64)
  J::Matrix{Float64} = randn(n, n) ./ sqrt(n)  # normalization factor

  for i in 1:n
    J[i, i] = 0  # zero the diagonal
    for j in 1:(i-1)
      J[i, j] = J[j, i]  # symmetric, may be dropped
    end
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
