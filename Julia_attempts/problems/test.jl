using ITensors
using PyCall
using Plots

N = 10 # number of sites
d = 3  # dimension of each site
chi = 16 # bond dimension of the MPS
s = siteinds("S=1/2",N)
psi = randomMPS(s;linkdims=chi)

ampo = OpSum()
# Then put operators into ampo...
for j in 1:N
    ampo+="Sz",j
end
sites = siteinds(psi) # Get site indices from your MPS
H = MPO(ampo,sites)

# Compute <psi|H|psi>
energy_psi = inner(psi',H,psi)