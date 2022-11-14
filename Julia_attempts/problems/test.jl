using ITensors
using PyCall
using Plots

N = 3 # number of sites
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
H = MPO(ampo,s)
sweeps = Sweeps(5) # number of sweeps is 5


maxdim!(sweeps,10,20,100,100,200) # gradually increase states kept
cutoff!(sweeps,1E-10) # desired truncation error
psi0 = randomMPS(s,5)
energy,psi = dmrg(H,psi0,sweeps)

magz = expect(psi,"Sz")




# Compute <psi|H|psi>
energy_psi = inner(psi',H,psi)

a = contract(H)
Evals,Evecs = eigen(reshape(a.tensor,(2^N,2^N)))















J=2
B=-.01
N=10
sites = siteinds("S=1/2",N)

hterms = OpSum()
for j=1:(N-1) #add iteraction terms
    hterms -= J,"Sz",j,"Sz",j+1
end
hterms -= J,"Sz",N,"Sz",1  # term 'wrapping' around the ring
for j=1:(N) #add magnetic terms
    hterms -= B,"Sz",j
end

H = MPO(hterms,sites)

sweeps = Sweeps(5) # number of sweeps is 5


maxdim!(sweeps,10,20,100,100,200) # gradually increase states kept
cutoff!(sweeps,1E-10) # desired truncation error

psi0 = randomMPS(ComplexF64, sites)
magz = expect(psi0,"Sz")
print(mean(magz))
energy,psi = dmrg(H,psi0,sweeps)
magz2 = expect(psi,"Sz")
print(mean(magz2))


B=1
# Make H for measuring the energy
for (i,B) in enumerate(-5:5)
terms = OpSum()
for j in 1:(N)
terms += B,"Sz", j
terms += 0.1B, "Sx",j
end
H = MPO(terms, s)
Szterms = OpSum()
for j in 1:N
    Szterms += "Sz",j
end
SzH = MPO(Szterms,s)


psi0 = randomMPS(s)
sweeps = Sweeps(5) # number of sweeps is 5

energy,psi = dmrg(H,psi0,sweeps)
en_vec[i]=(mean(expect(psi,"Sz")))
end















function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ)
    # print("B: ",B)
    h = 1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2) +
        B*op("Sz", s1) * op("Sz", s2)
    return exp(τ * h)
  end
  
N=10
cutoff=1E-8
δτ=0.1 
beta_max=1.0
# Make an array of 'site' indices
s = siteinds("S=1/2", N)
B = 0.1
######################################
# Ground state
######################################
terms = OpSum()
for j in 1:(N-1)
terms += "Sz", j,"Sz",j+1
end
for j in 1:N
terms += B,"Sz", j
end
H = MPO(terms, s)

Szterms = OpSum()
for j in 1:N
Szterms += "Sz",j
end
SzH = MPO(Szterms,s)


psi0 = randomMPS(s)
sweeps = Sweeps(5) # number of sweeps is 5

energy,psi = dmrg(H,psi0,sweeps)

######################################
# thermalization
######################################

# Make gates (1,2),(2,3),(3,4),...
gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,)) for n in 1:(N - 1)], s)
# Include gates in reverse order to complete Trotter formula
append!(gates, reverse(gates))

#make rho from groundstate
rho = outer(psi',psi)
magz = inner(rho,SzH)/N

# en_vec=zeros(Int(beta_max/δτ)+1)
for (i,β) in enumerate(0:δτ:beta_max)
energy = inner(rho, H)
magz = -inner(rho,SzH)/N
@printf("β = %.2f energy = %.8f, magz = %.5f\n", β, energy,magz)
rho = apply(gates, rho; cutoff)
rho = rho / tr(rho)
end
