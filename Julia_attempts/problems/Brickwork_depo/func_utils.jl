using ITensors
using Random
using LinearAlgebra


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
    rho = copy(rho)
    for n in sites
      rho[n] = swapinds(rho[n], siteinds(rho, n)...)
    end
    return rho
end

function trace_norm_dense(A)
    #=
    A: dense matrix
    =#
    e,_=eigen(A)
    return sum(abs.(e))
end

function negativity(rho::MPO, b, s)
    #=
    log negativity of an itensor mpo
    =#
    n = length(rho)
    orthogonalize!(rho,b)
    rho_temp = deepcopy(rho)
    # normalize!(rho_temp)
    M=partial_transpose(rho_temp,[b:n-b...])

    #turn this mpo into a single tensor
    T = rho_to_dense(M,s)
    return  (trace_norm_dense(T)-1)/2
end

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
