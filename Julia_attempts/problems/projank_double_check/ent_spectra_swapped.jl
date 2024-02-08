using ITensors
using Random
using LinearAlgebra
using Statistics
using NDTensors


ITensors.op(::OpName"Rand",::SiteType"Qubit") = 
    RandomUnitaryMatrix(4)
ITensors.op(::OpName"Rand1",::SiteType"Qubit") = 
    RandomUnitaryMatrix(2)
ITensors.op(::OpName"PPgate",::SiteType"Qubit") = 
    PPgate()
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

N=6
s = siteinds("Qubit", N) #+1 for ancilla
psi = productMPS(s, "Up" )
# rho = outer(psi,psi)
tp = calculate_trpk(psi)

arr_r = []
for i in 1:10
    N=12
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = randomMPS(s,20)


    append!(arr_r,calculate_r(psi))
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
      hj = op("PPgate",[s1,s2])
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


#random ppgate case
arr_r = []
arr_swapped = []
for i in 1:20
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

    psi = do_exp(N,steps,psi,s)
    
    append!(arr_r,calculate_r(psi))
    append!(arr_swapped,calculate_r(alternating_swap(psi)))
end
mean(arr_r)
mean(arr_swapped)
#T1 T1 = (CNOT_{12}CNOT_{23}CNOT_{34}...H_1 H_2 H_3 ...)
arr_r = []
for i in 1:10
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
    #set initial gate
    gates = ITensor[]
    for i in 1:N-1
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s1,s2])
        push!(gates, hj)
    end
    for i in 1:N
        s1 = s[i]
        hj = op("H",s1)
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s)
    append!(arr_r,calculate_r(psi))
end
mean(arr_r)

#T2 = (CNOT_{12}CNOT_{34}CNOT_{56}...CNOT_{32}CNOT{54}...)
arr_r = []
for i in 1:10
    N=12
    steps = 100
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    #set initial gate
    gates = ITensor[]
    for i in 1:2:N-1
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s1,s2])
        push!(gates, hj)
    end
    for i in 2:2:N-1
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s2,s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s)
    append!(arr_r,calculate_r(psi))
end
mean(arr_r)
#T3 = H_1 H_2 CNOT_{12}
arr_r = []
for i in 1:100
    N=12
    steps = 100
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    #set initial gate
    gates = ITensor[]
    for i in 1:2
        s1 = s[i]
        hj = op("H",s1)
        push!(gates, hj)
    end
    for i in 1:1
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s1,s2])
        push!(gates, hj)
    end

    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s)
    append!(arr_r,calculate_r(psi))
end
mean(arr_r)
#T4 = H_1 H_2 H_3 CNOT_{12} CNOT_{23}
arr_r = []
for i in 1:100
    N=12
    steps = 100
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)
    #set initial gate
    gates = ITensor[]
    for i in 1:3
        s1 = s[i]
        hj = op("H",s1)
        push!(gates, hj)
    end
    for i in 1:2
        s1 = s[i]
        s2 = s[i+1]
        hj = op("CNOT",[s1,s2])
        push!(gates, hj)
    end

    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s)
    append!(arr_r,calculate_r(psi))
end
mean(arr_r)

# let
#random ppgate case with trpk
arr_r = zeros(21)
for i in 1:10
    N=12
    steps = 60
    s = siteinds("Qubit", N)
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:N
        s1 = s[i]
        hj = op("Rand1",[s1])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s)
    arr_r += calculate_trpk(psi)
end
show(arr_r./10)
#random ppgate case with trpk
arr_r2 = zeros(21)
for i in 1:10
    N=12
    steps = 60
    s = siteinds("Qubit", N)
    psi = productMPS(s, "Up" )
    gates=ITensor[]
    for i in 1:2:N-1
        s1 = s[i]
        s2 = s[i+1]
        hj = op("Rand",[s1,s2])
        push!(gates, hj)
    end
    psi = apply(gates,psi)

    psi = do_exp(N,steps,psi,s)
    arr_r2 += calculate_trpk(psi)
end
show(arr_r2./10)
show(arr_r-arr_r2)

# end