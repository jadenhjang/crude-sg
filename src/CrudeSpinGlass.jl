module CrudeSpinGlass

export SpinGlass
export curie_weiss, sherrington_kirkpatrick
export reinitialize, update, mcmc, mcmc_test

include("spinglass.jl")
include("methods.jl")
include("statphys.jl")

end
