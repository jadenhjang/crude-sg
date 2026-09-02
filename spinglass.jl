@kwdef struct SpinGlass
  size::Int64
  spins::Array{Int64} = floor.(Int8, rand(Float64, size)*2)*2 .- 1
  config::Matrix{Float64} = generate_coupling(size)
end

# helper: generate coupling matrix J_{i,j}
function generate_coupling(n::Int64)
  J::Matrix{Float64} = randn(n, n) ./ sqrt(n)  # normalization factor

  for i in 1:n
    J[i, i] = 0  # zero the diagonal
    for j in 1:(i-1)
      J[i, j] = J[j, i]  # symmetric, may be dropped
    end
  end
  return J
end
