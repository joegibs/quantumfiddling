using ITensors
using Random
using LinearAlgebra
using Statistics
using NDTensors

#standard functions

ITensors.op(::OpName"Rand",::SiteType"Qubit") = 
    RandomUnitaryMatrix(4)
ITensors.op(::OpName"Rand1",::SiteType"Qubit") = 
    RandomUnitaryMatrix(2)
ITensors.op(::OpName"PPgate",::SiteType"Qubit") = 
    PPgate()

function PPgate()
    u1=RandomUnitaryMatrix(2)
    u2=RandomUnitaryMatrix(2)
    u1=dephase(u1)
    u2=dephase(u2)
    u = [u1[1,1] 0 0 u1[1,2];
         0 u2[1,1] u2[1,2] 0;
         0 u2[2,1] u2[2,2] 0;
         u1[2,1] 0 0 u1[2,2]]
    return u
end

function RandomUnitaryMatrix(dim)
    random_matrix=randn(ComplexF64,(dim,dim))
    Q, _ = NDTensors.qr_positive(random_matrix)
    return Q
end

function dephase(unitary)
    glob = det(unitary)
    theta = atan(imag(glob)/real(glob))/2
    unitary=unitary*exp(-im*theta)
    if real(det(unitary))<0
        unitary=unitary .* im
    end
    return unitary
end
function calculate_r(psi)
    b = 6
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    schmidt = zeros(1,dim(S,1))
    for n = 1:dim(S, 1)
        schmidt[n] = S[n,n]
    end

    schmidt= schmidt.^2

    d_k = [schmidt[i] - schmidt[i+1] for i in 1:(length(schmidt)-1)]

    r = [min(d_k[i],d_k[i+1])/max(d_k[i],d_k[i+1]) for i in 1:(length(d_k)-1)]
    return mean(r)
end


function sing_vals(psi,b)
    sing = []
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
    for n in 1:dim(S, 1)
        append!(sing, S[n,n]^2)
    end
    return sing
end
function calculate_trpk(psi)
    
    n = length(psi)
    #calculate rhoA
    rhoA = sing_vals(psi,Int(n/2))
    dep = deepcopy(rhoA)
    #exp it and recursively add to array
    tprk= [sum(rhoA.^0)]
    for _ in 1:20
        append!(tprk,sum(rhoA))
        rhoA = rhoA.*dep
    end
    return tprk
end

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
function gen_step(N,psi,s,step_num,gate)
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
      hj = op(gate,[s1,s2])
      Gj=hj
      push!(gates, Gj)
  end
  cutoff = 1E-8

  psi = apply(gates,psi)

  #calculate obs
  normalize!(psi)

  return psi
end

function do_exp(N,steps,psi,s,gate)
    for i in 1:steps
        psi= gen_step(N,psi,s,i,gate)
    end
    return psi
end

function rec_ent(psi,b,s) 
    #s is sites
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

function rec_ren(psi,b,s)  
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
    SvN = 0.0
    for n in 1:dim(S, 1)
      p = S[n,n]^2
      if p != 0
        SvN += p^2
      end
      SvN = -log2(SvN)
    end
    return SvN
end
function rec_ren_vec(psi,b,s)  
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
    SvN = [0.0 for i in 3:7]
    for n in 1:dim(S, 1)
        #singular values
      p = S[n,n]^2
      if p != 0
        SvN += [p^(i) for i in 3:7]
      end
    end

    SvN = log2.(SvN)
    SvN = [1/(1-i) for i in 3:7] .* SvN
    return SvN
end

#random initial state ppgate case
arr_r = []
arr_ent = []
for i in 1:100
    N=12
    steps = 144
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    psi = do_exp(N,steps,psi,s,"PPgate")
    append!(arr_r,calculate_r(psi))
    append!(arr_ent,[rec_ren_vec(psi,Int(N/2),s)])
end
mean(arr_r)
mean(arr_ent,dims=1)

#random initial state haar case
arr_r = []
arr_ent = []
for i in 1:100
    N=12
    steps = 144
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    psi = do_exp(N,steps,psi,s,"Rand")
    append!(arr_r,calculate_r(psi))
    append!(arr_ent,[rec_ren_vec(psi,Int(N/2),s)])
end
mean(arr_r)
mean(arr_ent,dims=1)

#cnot every other thing with ppgate
#random initial state  cnots PPgate case
arr_r = []
arr_ent = []
for i in 1:100
    N=12
    steps = 144
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    gates = ITensor[]
    for i in 1:2:11
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s,"PPgate")
    append!(arr_r,calculate_r(psi))
    append!(arr_ent,[rec_ren_vec(psi,Int(N/2),s)])
end
mean(arr_r)
mean(arr_ent,dims=1)

#random initial state  cnots with conjugation PPgate case
arr_r = []
arr_ent = []
for i in 1:100
    N=12
    steps = 144
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    gates = ITensor[]
    for i in 1:2:11
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s,"PPgate")

    gates=ITensor[]
    for i in 1:2:11
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    append!(arr_r,calculate_r(psi))
    append!(arr_ent,[rec_ren_vec(psi,Int(N/2),s)])
end
mean(arr_r)
mean(arr_ent,dims=1)

#random initial state  cnots with conjugation PPgate case
arr_r = []
arr_ent = []
for i in 1:100
    N=12
    steps = 144
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    gates = ITensor[]
    for i in 1:11
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s,"PPgate")

    gates=ITensor[]
    for i in 1:11
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    append!(arr_r,calculate_r(psi))
    append!(arr_ent,[rec_ren_vec(psi,Int(N/2),s)])
end
mean(arr_r)
mean(arr_ent,dims=1)
#random product state state PPgate case
arr_r = []
arr_ent = []
for i in 1:100
    N=12
    steps = 144
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:2:11
        s1 = s[i]
        s2 = s[i+1]
        hj = op("Rand",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)


    psi = do_exp(N,steps,psi,s,"PPgate")

    append!(arr_r,calculate_r(psi))
    append!(arr_ent,[rec_ren_vec(psi,Int(N/2),s)])
end
mean(arr_r)
mean(arr_ent,dims=1)

#random product state state PPgate case
arr_r = []
arr_ent = []
for i in 1:100
    N=12
    steps = 144
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:11
        s1 = s[i]
        s2 = s[i+1]
        hj = op("Rand",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)


    psi = do_exp(N,steps,psi,s,"PPgate")

    append!(arr_r,calculate_r(psi))
    append!(arr_ent,rec_ent(psi,Int(N/2),s))
end
mean(arr_r)
mean(arr_ent)