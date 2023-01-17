#Ha = Himp + Hhyb + Hbath

#Himp = ϵ (n_d↑ + n_d↓) + U n_d↑ n_d↓

#n_dσ = d†_σ d_σ

#H_hyb = -V∑(d†_σ c_2σ + H.c)
#H_bath = -t ∑ (c†_rσ c_(r+1)σ + H.c)
# from r= 2,σ to N-1

using ITensors
using Plots

N=10

sites = siteinds("Fermion", N)
psi0 = randomMPS(sites)
for i in 1:N-1
    l = getfirst(x->hastags(x,"Link,l=$(i)"),inds(psi0[i]))
    print(Int(dim(l)),'\n')
end
os = OpSum()
#impurity terms
for j in 1:(1)
  os .+= 0.5, "N", j
end
H = MPO(os, sites)

# Plan to do 5 DMRG sweeps:
nsweeps = 5
# Set maximum MPS bond dimensions for each sweep
maxdim = [10, 20, 100, 100, 200]
# Set maximum truncation error allowed when adapting bond dimensions
cutoff = [1E-11]

# Run the DMRG algorithm, returning energy and optimized MPS
energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
for i in 1:N-1
    l = getfirst(x->hastags(x,"Link,l=$(i)"),inds(psi[i]))
    print(Int(dim(l)),'\n')
end

#%%%%%%%%%%%%%%%%%%%%%%%
#this way is bad and will use huge bond dims
N=100

sites = siteinds("Electron", N)
psi0 = randomMPS(sites)
for i in 1:N-1
    l = getfirst(x->hastags(x,"Link,l=$(i)"),inds(psi0[i]))
    print(Int(dim(l)),'\n')
end

ϵ = 1
U = 1
V = 1
t=1

os = OpSum()
#impurity terms
for j in 1:(1)
  os .+= ϵ, "n↑",j
  os .+= ϵ, "n↓", j
  os .+= U,"Ntot", j 
end
#hyb terms
for j in 1:1
    os .+= V, "Cdagup", j, "Cup", j + 1
    os .+= V, "Cdagdn", j, "Cdn", j + 1
    os .+= V, "Cdagup", j + 1, "Cup", j
    os .+= V, "Cdagdn", j + 1, "Cdn", j
end
#bath terms
for j in 2:(N - 1)
    os .+= t, "Cdagup", j + 1, "Cup", j
    os .+= t, "Cdagdn", j + 1, "Cdn", j
end
H = MPO(os, sites)



# Plan to do 5 DMRG sweeps:
nsweeps = 5
# Set maximum MPS bond dimensions for each sweep
maxdim = [10, 20, 100, 100, 200]
# Set maximum truncation error allowed when adapting bond dimensions
cutoff = [1E-11]

# Run the DMRG algorithm, returning energy and optimized MPS
energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

bond_dims=[]
for i in 1:N-1
    l = getfirst(x->hastags(x,"Link,l=$(i)"),inds(psi[i]))
    push!(bond_dims,dim(l))
end

plot(bond_dims)

# split site representation

N=100*2

sites = siteinds("Fermion", N)
psi0 = randomMPS(sites)
for i in 1:N-1
    l = getfirst(x->hastags(x,"Link,l=$(i)"),inds(psi0[i]))
    print(Int(dim(l),'\n'))
end

ϵ = 1
U = 1
V = 1
t=1

os = OpSum()
#impurity terms
for j in N/2-1:N/2
  os .+= ϵ, "n↑", Int(j)
  os .+= ϵ, "n↓", Int(j)
  os .+= U, "N", Int(j)
end
#hyb terms
for j in N/2-2:2:N/2
    os .+= V, "Cdag", Int(j), "C", Int(j) + 1
    os .+= V, "Cdag", Int(j), "C", Int(j) + 1
    os .+= V, "Cdag", Int(j) + 1, "C", Int(j)
    os .+= V, "Cdag", Int(j) + 1, "C", Int(j)
end
#bath terms
for j in 1:N/2-1
    os .+= t, "Cdag", Int(j) + 1, "C", Int(j)
    os .+= t, "Cdag", Int(j) + 1, "C", Int(j)
end
for j in N/2+1:N/2
    os .+= t, "Cdag", Int(j) + 1, "C", Int(j)
    os .+= t, "Cdag", Int(j) + 1, "C", Int(j)
end
H = MPO(os, sites)



# Plan to do 5 DMRG sweeps:
nsweeps = 5
# Set maximum MPS bond dimensions for each sweep
maxdim = [10, 20, 100, 100, 200]
# Set maximum truncation error allowed when adapting bond dimensions
cutoff = [1E-11]

# Run the DMRG algorithm, returning energy and optimized MPS
energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

bond_dims=[]
for i in 1:N-1
    l = getfirst(x->hastags(x,"Link,l=$(i)"),inds(psi[i]))
    push!(bond_dims,dim(l))
end

plot(bond_dims)