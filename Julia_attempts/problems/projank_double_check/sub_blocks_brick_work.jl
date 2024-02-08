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
function calculate_r(psi,b)
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

function make_row(N,eoo,pc,start,length)
    #=
    N: number of sites
    eoo: even or odd step
    pc: periodic
    =#
    if start+length>N-1
        if eoo
            lst =[[i,i+1] for i in start:2:N-1]
        else
            lst = [[i,i+1] for i in start+1:2:N-1]
        end
        if pc
            if !eoo && !Bool(N%2)
                append!(lst,[[N,1]])
            end
        end
    else
        if eoo
            lst =[[i,i+1] for i in start:2:start+length-2]
        else
            lst = [[i,i+1] for i in start+1:2:start+length-2]
        end
        if pc
            if !eoo && !Bool(N%2)
                append!(lst,[[N,1]])
            end
        end
    end
    return lst
end
function gen_step(N,psi,s,step_num,block_size)
    #=
    perform one step of brickwork
    =#
    #apply gates
    for start in 1:block_size:N
        row = make_row(N,Bool(step_num%2),false,start,block_size)
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
    end

  return psi
end

function do_exp(N,steps,psi,s,block_size)
    for i in 1:steps
        psi= gen_step(N,psi,s,i,block_size)
    end
    return psi
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


N=10
block_size=10
s = siteinds("Qubit", N) #+1 for ancilla
psi = productMPS(s, "Up" )
psi = do_exp(N,block_size*2,psi,s,block_size)
calculate_r(psi,5)

arr_r=[]
arr_swapped = []
for i in  1:100
    N=10
    block_size=10
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    psi = do_exp(N,block_size*2,psi,s,block_size)

    append!(arr_r,calculate_r(psi,5))
    append!(arr_swapped,calculate_r(alternating_swap(psi),5))
end
mean(arr_r)
mean(arr_swapped)

arr_r=[]
arr_swapped = []
for i in  1:100
    N=10
    block_size=2
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    psi = do_exp(N,block_size*2,psi,s,block_size)

    append!(arr_r,calculate_r(psi,5))
    append!(arr_swapped,calculate_r(alternating_swap(psi),5))
end
mean(arr_r)
mean(arr_swapped)


arr_r=[]
arr_swapped = []
for i in  1:100
    N=20
    block_size=5
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    psi = do_exp(N,block_size*2,psi,s,block_size)

    append!(arr_r,calculate_r(psi,8))
    append!(arr_swapped,calculate_r(alternating_swap(psi),5))
end
mean(arr_r)
mean(arr_swapped)
arr_r = []
#T4 = H_1 H_2 H_3 CNOT_{12} CNOT_{23}
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