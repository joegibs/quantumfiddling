#make mps
using ITensors
using Random

N = 4
# sites = siteinds(2,N)
s = siteinds("Qubit", N)

psi = MPS(sites,linkdims=2)

#make a row
function make_row(N,eoo,pc)
    if eoo
        lst =[[i,i+1] for i in 1:2:N-1]
    else
        lst = [[i,i+1] for i in 2:2:N-1]
    end
    if pc
        if !eoo && !Bool(N%2)
            append!(lst,[[N,1]])
        end
    end
return lst
end

#make gates
# RXX
# sites = siteinds("S=1/2",N)
# Sz1 = op("Sz",sites[1])
# # Sp3 = op("S+",sites[3])
# s0 = siteinds("Qubit", N)
# os = OpSum()
# os+=1,"Rxy",s0
# H = MPO(os,s0)
function make_rand_opertor()
    T = randn(2,2,2,2)

    k = Index(2,"index_k")
    l = Index(2,"index_l")
    m = Index(2,"index_m")
    n = Index(2,"index_n")
    G = ITensor(T,k,l,m,n)
    return G
end
#apply gates

function apply_gate(psi,G,sites)
    a,b=sites[1],sites[2]
    orthogonalize!(psi,a)
    wf = (psi[a] * psi[b]) * G
    noprime!(wf)
    inds3 = uniqueinds(psi[a],psi[b])
    U,S,V = svd(wf,inds3,cutoff=1E-8)
    psi[a] = U
    psi[b] = S*V
    return psi
end

#compute 

#Measure
function rec_ent(psi,b)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
    end
    return SvN
end
#record

#plot

N = 5
# sites = siteinds(2,N)
s = siteinds("Qubit", N)
psi = randomMPS(s)

row = make_row(N,true,true)
# sites = [1,2]
# G = make_rand_opertor()


# orthogonalize!(psi,3)

# wf = (psi[3] * psi[4]) * G
# noprime!(wf)
# inds3 = uniqueinds(psi[3],psi[4])
# U,S,V = svd(wf,inds3,cutoff=1E-8)
# psi[3] = U
# psi[4] = S*V
D = 40
rec_vN = []
for d in 1:D
    row = make_row(N,Bool(d%2),true)
    for i in row
        gate = make_rand_opertor()
        psi = apply_gate(psi, gate, i)
    end
    append!(rec_vN,rec_ent(psi,Int(round(N/2))))
end
#%%

using ITensors

#meas funct
ITensors.op(::OpName"Pup",::SiteType"Qubit") =
 [1 0
  0 0]
  ITensors.op(::OpName"Pdn",::SiteType"Qubit") =
 [0 0
  0 1]
#apply <Pup> to one site of the mps to get the prob
function samp_mps(psi,s,site)
    cutoff = 1E-8

    magz = expect(psi,"Pup")
    proj="Pup"
    if magz[site]<rand()
        proj ="Pdn"
    end
    G = ITensor[op(proj,s[site])]
    psi = apply(G, psi; cutoff)
    normalize!(psi)
    return psi
end


N = 10
cutoff = 1E-8
tau = 0.1
ttotal = 5.0

# Make an array of 'site' indices
s = siteinds("Qubit", N)
function gen_step(psi,s,step_num,meas_p)
    #apply gates
    row = make_row(N,Bool(step_num%2),true)
    gates = ITensor[]
    for j in row
        s1 = s[j[1]]
        s2 = s[j[2]]
        hj = op("Rxx",[s1,s2],ϕ=0.5)
        Gj=hj
        push!(gates, Gj)
    end
    psi = apply(gates, psi; cutoff)
    #calculate obs
    mesured_vals = [1]
    #sample as needed
    samp_row=gen_samp_row(N,meas_p)
    for (i,x) in enumerate(samp_row)
        if Bool(x)
            psi = samp_mps(psi,s,i)
        end
    end
    return psi,measured_vals
end
# Make gates (1,2),(2,3),(3,4),...
function 
row = make_row(N,true,true)
gates = ITensor[]
for j in row
    s1 = s[j[1]]
    s2 = s[j[2]]
    hj = op("Rxx",[s1,s2],ϕ=0.5)
    Gj=hj
    push!(gates, Gj)
    
end
# Include gates in reverse order too
# (N,N-1),(N-1,N-2),...
#append!(gates, reverse(gates))

# Initialize psi to be a product state (alternating up and down)
psi = productMPS(s, "Up" )
gates = ITensor[]
for j in 1:N
    s1 = s[j[1]]
    hj = op("H",s1)
    Gj=hj
    push!(gates, Gj)
end
psi = apply(gates, psi; cutoff)

psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

c = div(N, 2) # center site

# Compute and print <Sz> at each time step
# then apply the gates to go to the next time
# for t in 0.0:tau:ttotal
Sz = expect(psi, "Sz"; sites=c)
println("$Sz")

psi = apply(gates, psi; cutoff)
normalize!(psi)
# end

hj =op("Sz", s[1]) * op("Sz", s[2])
hj = op("Rxx",[s[1],s[2]],ϕ=0.1)
function gen_samp_row(N,meas_p)
    return [rand()<meas_p ? 1 : 0 for i in 1:N]
end
function preform_step(psi,s,step_num)
    
end


end
#Given an MPS called "psi",
#and assuming j > i
N = 20
s = siteinds("Qubit", N)
# auto state = InitState(sites,"Up");
psi = randomMPS(s)

i = 4;
j = 8;

#'gauge' the MPS to site i
#any 'position' between i and j, inclusive, would work here
orthogonalize!(psi,i)

psidag = dag(psi);
prime(psidag,"Link");

#//index linking i to i-1:
li_1 = linkind(psi,i);

rho = prime(psi[i],li_1)*prime(psidag[i],"Site");

for k in i:j#(int k = i+1; k < j; ++k)
    print(k)
    lk = linkind(psi,k)
    rho *= prime(psi[k],lk);
    rho *= prime(psidag[k],"Site");
end
#index linking j to j+1:
lj = linkind(psi,j+1);

rho *= prime(psi[j],lj);
rho *= prime(psidag[j],"Site");
D,U = eigen(rho)
SvN = 0.0
for n in 1:dim(D, 1)
    print(n)
    p = D[n,n]^2
    if p != 0
    SvN -= p * log(p)
    end
end