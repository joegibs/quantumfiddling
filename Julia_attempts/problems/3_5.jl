using ITensors
"""
not done
"""
#differnt BC from handwaving
X = [0 1 ; 1 0]
I = [1 0; 0 1]
A = reshape([[1 0 ; 0 1]  [0 1 ; 1 0]],(2,2,2))

function e_mat(n)
    return ones((n for i in 1:n)...)-LinearAlgebra.I(n)
end

#initialize matricies
N = 5
q = 3
a = [i+1==j ? 1 : 0 for i = 1:q, j=1:q]
a[q,1] = 1
A = reshape(reduce(vcat,transpose([a^i for i in 0:q-1])),(q for i in 1:q)...)
I = LinearAlgebra.I(q)

#initialize sites and MPS
sites = siteinds(q,N)
psi = MPS(sites,linkdims=q)
ampo = OpSum()
# Then put operators into ampo...
ampo += ("I",1)
sites = siteinds(psi) # Get site indices from your MPS
H = MPO(ampo,sites)

#assign data in each MPS site
for i in 1:N
    if i==1
        psi[i] = ITensor(I,inds(psi[i]))
    elseif i==N
        psi[i] = ITensor(A[:,:,q],inds(psi[i]))
    else
        psi[i] = ITensor(A,inds(psi[i]))
    end
end

#Evaluate a state for each site and contract with the mps
el = [2,1,2,2,1]
V = ITensor(1.)
for j=1:N
  global V *= (psi[j]*state(sites[j],el[j]))
end
v = scalar(V)