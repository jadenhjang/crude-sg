module CrudeSpinGlass

export SpinGlass
export curie_weiss, sherrington_kirkpatrick, random_field, const_field
export reinitialize, update, mcmc, mcmc_test, hamiltonian, meanfield_rfim!

include("spinglass.jl")
include("methods.jl")
include("statphys.jl")

end
