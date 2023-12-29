using ITensors
using Random
using Plots
using Statistics
using LaTeXStrings
using LinearAlgebra
using PyCall


function RandomUnitaryMatrix(dim)
    random_matrix=randn(ComplexF64,(dim,dim))
    Q, _ = NDTensors.qr_positive(random_matrix)
    return Q
  end
function IXgate(N::Int)
    theta = rand()*2*pi
    eps=0.5#rand()
    u =exp(1im*theta).*(sqrt((1-eps)).*Matrix(I,4,4) + 1im.*sqrt(eps).*[[0,0,0,1] [0,0,1,0] [0,1,0,0] [1,0,0,0]])
    return u
end

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


function rec_ent(psi,b,s)  
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
    SvN = 0.0
    for n in 1:dim(S, 1)
      p = S[n,n]^2
      if p != 0
        SvN -= p * log2(p)
      end
    end
    return SvN
end
function rec_ent_rho(rho::MPO,b,s)
    #=
    Bipartite entropy across b for itensor mpo
    =#
    n = length(rho)
    # orthogonalize!(rho,b)
    rho_temp = deepcopy(rho)
    # s = siteinds("Qubit",n) 
  
    #contract half   x x x x | | | |
    L = ITensor(1.0)
      for i = 1:b
        L *= tr(rho_temp[i])
      end
      # absorb
      rho_temp[b+1] *= L
      # no longer a proper mpo
      M =MPO(n-b)
      for i in 1:(n-b)
          M[i]=rho_temp[b+i]
      end
      M=M/tr(M)
      #turn this mpo into a single tensor
      T = prod(M)
   
      # @show T
      _,S,_ = svd(T,s)#[inds(T)[i] for i = 1:2:length(inds(T))])
      SvN = 0.0
      for n in 1:dim(S, 1)
        p = S[n,n]
        if p != 0
          SvN -= p * log2(p)
        end
    end
    return real(SvN)
end
function split_ren(psi,b)
    rho=outer(psi',psi)
    n = length(rho)
    rho_temp = deepcopy(rho)
    s = siteinds("Qubit",n) 
  
    #contract half   x x x x | | | |
    L = ITensor(1.0)
    for i = 1:b
      L *= tr(rho_temp[i])
    end
    # absorb
    rho_temp[b+1] *= L
    # no longer a proper mpo
    M =MPO(n-b)
    for i in 1:(n-b)
        M[i]=rho_temp[b+i]
    end
    ren = -log2(tr(apply(M,M)))
    return ren
  end

function entropy_subsys(psi,inds)
    i = inds[1];
    j = inds[end];
    #'gauge' the MPS to site i
    #any 'position' between i and j, inclusive, would work here
    orthogonalize!(psi,i)

    psidag = dag(psi);
    prime(psidag,"Link");

    #//index linking i to i-1:
    li_1 = linkind(psi,i);

    rho = prime(psi[i],li_1)*prime(psidag[i],"Site");

    for k in inds[2:end-1]#(int k = i+1; k < j; ++k)
        # print(k)
        lk = linkind(psi,k)
        rho *= prime(psi[k],lk);
        rho *= prime(psidag[k],"Site");
    end
    for k in setdiff([i:j...],[inds])
        rho *= psi[k];
        rho *= psidag[k];
    end
    #index linking j to j+1:
    lj = linkind(psi,inds[end]+1);
    rho *= prime(psi[j],lj);
    rho *= prime(psidag[j],"Site");

    D,U = eigen(rho)
    SvN = 0.0
    for n in 1:dim(D, 1)
        # print(n,'\n')
        p = D[n,n]^2
        if p != 0
        SvN -= p * log(p)
        end
    end
    return SvN
end

function two_point_entropy(psi,i,j)
    #'gauge' the MPS to site i
    #any 'position' between i and j, inclusive, would work here
    orthogonalize!(psi,i)

    psidag = dag(psi);
    prime(psidag,"Link");

    #//index linking i to i-1:
    li_1 = linkind(psi,i);

    rho = prime(psi[i],li_1)*prime(psidag[i],"Site");

    for k in i+1:j-1
        rho *= psi[k];
        rho *= psidag[k];
    end
    #index linking j to j+1:
    lj = linkind(psi,j+1);
    rho *= prime(psi[j],lj);
    rho *= prime(psidag[j],"Site");
    # print(rho)
    D,U = eigen(rho)
    SvN = 0.0
    for n in 1:dim(D, 1)
        # print(n,'\n')
        p = D[n,n]*conj(D[n,n])
        if p != 0
        SvN -= p * log(p)
        end
    end
    return real(SvN)
end

function tri_part_MI(psi_abc,sysa,sysb,sysc)
    hab = entropy_subsys(psi_abc, vcat(sysa,sysb))
    hac = entropy_subsys(psi_abc, vcat(sysa,sysc))
    hbc = entropy_subsys(psi_abc, vcat(sysc,sysb))
    
    habc = entropy_subsys(psi_abc, vcat(sysa,sysb,sysc))

    ha = entropy_subsys(psi_abc, sysa)
    hb = entropy_subsys(psi_abc, sysb)
    hc = entropy_subsys(psi_abc, sysc)
    return hb + ha + hc - hab - hac - hbc + habc
end

function two_point_MI(psi, index_a, index_b)
    ha = rec_ent(psi, index_a)
    hb = rec_ent(psi, index_b)
    hab = two_point_entropy(psi, index_a,index_b)
    print("entropies" , ha, " " , hb , " ",hab,"\n")
    return ha + hb - hab
end

ITensors.op(::OpName"Pup",::SiteType"Qubit") =
 [1 0
  0 0]
ITensors.op(::OpName"Pdn",::SiteType"Qubit") =
 [0 0
  0 1]
ITensors.op(::OpName"Rand",::SiteType"Qubit") = 
    RandomUnitaryMatrix(4)
ITensors.op(::OpName"IX",::SiteType"Qubit") = 
    IXgate(4)
ITensors.op(::OpName"Ihaar",::SiteType"Qubit") = 
    Ihaargate(4)
ITensors.op(::OpName"Iden",::SiteType"Qubit") = 
[1 0
0 1]
function gen_samp_row(N,meas_p)
    return [rand()<meas_p ? 1 : 0 for i in 1:N]
end
function samp_mps(psi::MPS,s,site)
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

function gen_step(N,psi,s,step_num,meas_p)
    #apply gates
    row = make_row(N,Bool(step_num%2),false)
    gates = ITensor[]
    for j in row
        s1 = s[j[1]]
        s2 = s[j[2]]
        hj = op("Rand",[s1,s2])
        Gj=hj
        push!(gates, Gj)
    end
    cutoff = 1E-15

    psi = apply(gates, psi; cutoff)
    #calculate obs
    rho = outer(psi',psi)
    measured_vals = (rec_ent(psi,Int(round(N/2)),s),rec_ent_rho(rho,Int(round(N/2)),s))
    
    #metts shenanagins
    # psi = metts(psi,s)
    #sample as needed
    samp_row=gen_samp_row(N,meas_p)
    for (i,x) in enumerate(samp_row)
        if Bool(x)
            psi = samp_mps(psi,s,i)
        end
    end
    return psi,measured_vals
end
function rho_to_dense(rho,s)
    Hitensor = ITensor(1.)
    N = length(s)
    for i = 1:N
        Hitensor *= rho[i]
    end
  
    A=Array(Hitensor,prime(s),s)
    return reshape(A,2^N,2^N)
  end

N=2
step_num=1
meas_p=0.5
s = siteinds("Qubit", N)
psi = productMPS(s, "Up" )
rho = outer(psi',psi)


row = make_row(N,Bool(step_num%2),false)
gates = ITensor[]
for j in row
    s1 = s[j[1]]
    s2 = s[j[2]]
    hj = op("Rand",[s1,s2])
    Gj=hj
    push!(gates, Gj)
end
cutoff = 1E-15

psi = apply(gates, psi; cutoff)
rho = apply(gates, rho; apply_dag=true, cutoff)

rhoc=outer(psi',psi)
isapprox(rho_to_dense(rho,s),rho_to_dense(rhoc,s))

#sample as needed
samp_row=gen_samp_row(N,meas_p)

psim=deepcopy(psi)
for (i,x) in enumerate(samp_row)
    if Bool(x)
        psim = samp_mps(psim,s,i)
    end
end

rhom=samp_mps(rho,s,samp_row)

rhocm=outer(psim',psim)
rho_to_dense(rhom,s)
rho_to_dense(rhocm,s)

isapprox(rho_to_dense(rhom,s),rho_to_dense(rhocm,s))

function samp_mps(rho::MPO,s,samp_row)
    #=
    sample an itensor mpo and return the next resulting mps
    =#
    N = length(rho)
    samp =deepcopy(rho)
    samples= sample(samp)
    magz = [x == 1 ? "Pup" : "Pdn" for x in samples]
  
    gates = ITensor[]
  
    for i in 1:N
      if Bool(samp_row[i])
          hj = op(magz[i],s[i])
          push!(gates, hj)
      end
    end
    rho = apply(gates, rho;apply_dag=true)
   
    rho=rho/tr(rho)
    return rho
  end

p = plot([0.0:0.1:0.8...],decays,title=string("Bip_ent Gat: 2haar, varying meas_p"), label=string.(transpose([6:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
p = plot([0.0:0.1:0.8...],real(growths),title=string("Bip_ent Gat: 2haar, varying meas_p"), label=string.(transpose([6:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn