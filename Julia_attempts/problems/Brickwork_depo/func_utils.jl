using ITensors
using Random
using LinearAlgebra
using NDTensors

function RandomUnitaryMatrix(dim)
  random_matrix=randn(ComplexF64,(dim,dim))
  Q, _ = NDTensors.qr_positive(random_matrix)
  return Q
end

function rho_to_dense(rho,s)
    #=
    MPO to Dense Matrix
    rho: itensor mpo
    s: itensor sites
    =#
    Hitensor = ITensor(1.)
    N = length(s)
    for i = 1:N
        Hitensor *= rho[i]
    end
  
    A=Array(Hitensor,prime(s),s)
    return reshape(A,2^N,2^N)
end

function partial_transpose(rho::MPO, sites)
    #=
    wikipedia.org/wiki/Peres–Horodecki_criterion
    =#
    rho1 = copy(rho)
    for n in sites
      rho1[n] = swapinds(rho1[n], siteinds(rho1, n)...)
    end
    return rho1
end

function trace_norm_dense(A)
    #=
    A: dense matrix
    =#
    e,_=eigen(A)
    # _,S,_=svd(A)
    return (sum(abs.(e))-sum(e))/2
end

function purity(rho::MPO)
  pur = tr(apply(rho,rho))
  return pur
end

function negativity(rho::MPO, b, s)
    #=
    negativity of an itensor mpo
    =#
    n = length(rho)
    #orthogonalize!(rho,b)
    rho_temp = deepcopy(rho)
    # normalize!(rho_temp)
    M=partial_transpose(rho_temp,[b:n...])

    #turn this mpo into a single tensor
    T = rho_to_dense(M,s)
    return  (trace_norm_dense(T))
end
function log_negativity(A::MPO, b, s)
  neg = negativity(A, b, s)
  return log2(2*neg+1)
end

# function rec_ent(rho::MPO,b,s)
#   #=
#   Bipartite entropy across b for itensor mpo
#   =#
#   n = length(rho)
#   # orthogonalize!(rho,b)
#   rho_temp = deepcopy(rho)
#   # s = siteinds("Qubit",n) 

#   #contract half   x x x x | | | |
#   L = ITensor(1.0)
#     for i = 1:b
#       L *= tr(rho_temp[i])
#     end
#     # absorb
#     rho_temp[b+1] *= L
#     # no longer a proper mpo
#     M =MPO(n-b)
#     for i in 1:(n-b)
#         M[i]=rho_temp[b+i]
#     end
#     M=M/tr(M)
#     #turn this mpo into a single tensor
#     SvN = -log2(tr(apply(M,M)))
#   #   # @show T
#   #   _,S,_ = svd(T,s)#[inds(T)[i] for i = 1:2:length(inds(T))])
#   #   SvN = 0.0
#   #   for n in 1:dim(S, 1)
#   #     p = S[n,n]
#   #     if p != 0
#   #       SvN -= p * log2(p)
#   #     end
#   # end
#   return SvN
# end

function rec_ent(rho::MPO,b,s)
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
    return SvN
end

#get bond dims of a mpo
function bond_dim_array(rho)
  arr=[]
  for i in rho
      for j in inds(i)
          if occursin("Link",string(tags(j)))
              push!(arr,dim(j))
          end
      end
  end
  return Int.(arr)
end

function max_bond_dim(rho)
  max_dim = maximum(bond_dim_array(rho))
  return max_dim
end

#### one pass alg for mean and var
function welford_mean_var(existing_agg, new_value)
  (cnt, men, M2)  = existing_agg
  cnt= cnt+1
  delta = new_value - men
  men = men + delta/cnt
  delta2 = new_value - men
  M2 = M2 + delta .* delta2
  return (cnt,men,M2)
end

function welford_extract(existing_agg)
  (cnt, men, M2)  = existing_agg
  (men, variance) = (men, M2 ./ cnt)
  return (men, variance)
end



#file shenanagins
function open_csv(flname)
  m = []#Vector{Tuple{Matrix{Int}, Vector{Float64},Int, Float64, Vector{Any}, Vector{Any}}}()
  open(flname, "r") do io
      while !eof(io)
          push!(m, eval(Meta.parse(readline(io))))
      end
  end
  sits = [i for i in m[1]]
  interval = [i for i in m[2]]
  num_samp = m[3]
  noise_val = m[4]
  Int_Svns_Mean = hcat([i for i in m[5]]...)
  Int_Negs_Mean = hcat([i for i in m[6]]...)
  Int_Svns_Var = hcat([i for i in m[7]]...)
  Int_Negs_Var = hcat([i for i in m[8]]...)
  return sits,interval,num_samp,noise,Int_Svns_Mean,Int_Negs_Mean,Int_Svns_Var,Int_Negs_Var
end

function open_csv_purity(flname)
  m = []#Vector{Tuple{Matrix{Int}, Vector{Float64},Int, Float64, Vector{Any}, Vector{Any}}}()
  open(flname, "r") do io
      while !eof(io)
          push!(m, eval(Meta.parse(readline(io))))
      end
  end
  N = m[1]
  depth = m[2]
  num_samp = m[3]
  noise_int = [i for i in m[4]]
  meas_int = [i for i in m[5]]
  purity_arr = [i for i in m[6]]
  
  return N,depth,num_samp,noise_int,meas_int,purity_arr
end
