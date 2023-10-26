using ITensors
using Printf
using Random
using Plots

Random.seed!(1234)

let
  # Create 100 spin-one indices
  N = 100
  sites = siteinds("S=1",N)

  # Input operator terms which define 
  # a Hamiltonian matrix, and convert
  # these terms to an MPO tensor network
  # (here we make the 1D Heisenberg model)
  ampo = AutoMPO()
  for j=1:N-1
    ampo +=     ("Sz",j,"Sz",j+1)
    ampo += (0.5,"S+",j,"S-",j+1)
    ampo += (0.5,"S-",j,"S+",j+1)
  end
  H = MPO(ampo,sites)

  # Create an initial random matrix product state
  psi0 = randomMPS(sites)

  # Plan to do 5 passes or 'sweeps' of DMRG,
  # setting maximum MPS internal dimensions 
  # for each sweep and maximum truncation cutoff
  # used when adapting internal dimensions:
  sweeps = Sweeps(5)
  maxdim!(sweeps, 10,20,100,100,200)
  cutoff!(sweeps, 1E-10)
  @show sweeps

  # Run the DMRG algorithm, returning energy 
  # (dominant eigenvalue) and optimized MPS
  energy, psi = dmrg(H,psi0, sweeps)
  println("Final energy = $energy")

  corr = correlation_matrix(psi,"Sz","Sz")
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

using Plots
x=1:10;y = rand(10);
plot(x,y)