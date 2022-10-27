using ITensors
using Printf
using Random
using Plots

Random.seed!(1234)

let
  N = 1000

  # Create N spin-one degrees of freedom
  sites = siteinds("Fermion", N)
  # Alternatively can make spin-half sites instead
  #sites = siteinds("S=1/2",N)

  # Input operator terms which define a Hamiltonian
  os = OpSum()
  for j in 1:(N - 1)
    os += "Cdag",j,"C",j+1
    os += "C",j,"Cdag",j+1
  end
  # Convert these terms to an MPO tensor network
  H = MPO(os, sites)

  # Create an initial random matrix product state
  psi0 = randomMPS(sites, 10)

  # Plan to do 5 DMRG sweeps:
  nsweeps = 5
  # Set maximum MPS bond dimensions for each sweep
  maxdim = [10, 20, 100, 100, 200]
  # Set maximum truncation error allowed when adapting bond dimensions
  cutoff = [1E-11]

  # Run the DMRG algorithm, returning energy and optimized MPS
  energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
  corr = correlation_matrix(psi,"Cdag","C")
  one_corr = [abs(corr[1,x]) for x in 1:N]
#   print(corr[1,3])
#   print(one_corr)
#   print("\n")
  @printf("Final energy = %.12f\n", energy)
  x_axis=1:N
  scatter(x_axis,one_corr,xaxis=:log,yaxis=:log)
#   xaxis!("Site", :log10)
#   yaxis!("Corr", :log10)

end
