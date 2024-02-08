using ITensors
using Random
using LinearAlgebra
using Statistics
using NDTensors
using Plots
#define swapping function

function generate_indices(N)
    start_ind=Int(N/2)
    end_ind=2
    list = []
    while start_ind < N
        lst =[[i,i+1] for i in start_ind:-1:end_ind]
        list = vcat(list,lst)
        start_ind +=1
        end_ind +=2
    end
    return list
end
function alternating_swap(psi)
    row = generate_indices(length(psi))
    gates = ITensor[]
    for j in row
        s1 = s[j[1]]
        s2 = s[j[2]]
        hj = op("SWAP",[s1,s2])
        Gj=hj
        push!(gates, Gj)
    end
    psi = apply(gates,psi)
    return psi
end

####entanglement measures
function entanglement(psi, b)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)) )
    SVN = 0.0
    for n=1: dim(S, 1)
      p = S[n, n]^2
      if p > 0.0
        SVN = (SVN - p * log(p) )
      end
    end
    SVN= SVN/log(2)
    return SVN
  end

  function rec_ent(rho::MPO,s)
    #=
    Bipartite entropy across b for itensor mpo
        only even length mpo
    =#
    n = length(rho)
    # orthogonalize!(rho,b)
    rho_temp = deepcopy(rho)
    # s = siteinds("Qubit",n) 
  
    #contract half   x x x x | | | |
    
      for i = 2:2:N
        L = ITensor(1.0)
        L *= tr(rho_temp[i])
        rho_temp[i-1] *= L
      end
      # absorb
      
      # no longer a proper mpo
      M =MPO(Int(N/2))
      for i in 1:Int(N/2)
          M[i]=rho_temp[2*i-1]
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
    return SvN
end
############ generate random thingd
function RandomUnitaryMatrix(dim)
    random_matrix=randn(ComplexF64,(dim,dim))
    Q, _ = NDTensors.qr_positive(random_matrix)
    return Q
end
ITensors.op(::OpName"Rand",::SiteType"Qubit") = 
    RandomUnitaryMatrix(4)
ITensors.op(::OpName"Rand1",::SiteType"Qubit") = 
    RandomUnitaryMatrix(2)

function make_row(N,eoo,pc)
    #=
    N: number of sites
    eoo: even or odd step
    pc: periodic
    =#
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
function gen_step(N,psi,s,step_num)
    #=
    perform one step of brickwork
    =#
    #apply gates
  row = make_row(N,Bool(step_num%2),false)
  gates = ITensor[]
  measured_vals=([0],[0])
  for j in row
      s1 = s[j[1]]
      s2 = s[j[2]]
      hj = op("Rand",[s1,s2])
      Gj=hj
      push!(gates, Gj)
  end
  cutoff = 1E-8

  psi = apply(gates,psi)

  #calculate obs
  normalize!(psi)

  return psi
end

function do_exp(N,steps,psi,s)
    for i in 1:steps
        psi= gen_step(N,psi,s,i)
    end
    return psi
end


#### get some random haar state

ents = []
ents_swapped =[]
for i in 1:2
    N=4
    steps = 10
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    psi = do_exp(N,steps,psi,s)
    append!(ents,entanglement(psi,Int(N/2)))
    append!(ents_swapped,entanglement(alternating_swap(psi),Int(N/2)))
end
histogram(ents-ents_swapped)

#need to get alternating boundary
entanglement(alternating_swap(psi),Int(N/2))
rec_ent(outer(psi',psi),s)
