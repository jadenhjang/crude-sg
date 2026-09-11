@kwdef struct SpinGlass
  size::Int64
  spins::Array{Int64} = floor.(Int8, rand(Float64, size)*2)*2 .- 1
  h::Array{Float64} = random_field(size)
  J::Matrix{Float64} = curie_weiss(size)
end

# helper functions: generate coupling matrix J_{i,j}
function sherrington_kirkpatrick(n::Int64)
  A::Matrix{Float64} = randn(n, n) ./ sqrt(2*n)  # normalization factor
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

# spare graph function
function render_sparse(J::Matrix{Float64}, p::Float64)
  println("We will kill all connections with probability $p")
  # exception handling
  if p < 0 || p >= 1.0
    throw(DomainError(p, "probabilities must be [0,1)"))
  end

  # kill with probability p all J_ij
  # J .*= (rand(size(J)) .> p)
  sz = size(J, 1)
  for i in 1:sz
    for j in i:sz
      J[i, j] *= floor(rand()+(1-p))
      J[j, i] = J[i, j]
    end
  end
  return J
end

# helper functions: generate random field h_i
function random_field(n::Int64, h0::Float64=1.0)
  return randn(n) .* h0
end

function const_field(n::Int64, h0::Float64=1.0)
  return ones(n) .* h0
end
