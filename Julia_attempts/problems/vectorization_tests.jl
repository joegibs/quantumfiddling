using ITensors
using Plots
using Statistics



function rho_to_dense(rho,s)
    Hitensor = ITensor(1.)
    N = length(s)
    for i = 1:N
        Hitensor *= rho[i]
    end
  
    A=Array(Hitensor,prime(s),s)
    return reshape(A,2^N,2^N)
  end

function rec_ent(rho,b,s)
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
    #turn this mpo into a single tensor
    T = prod(M)
 
    _,S,_ = svd(T,s)
    SvN = 0.0
    for n in 1:dim(S, 1)
      p = S[n,n]
      if p != 0
        SvN -= p * log2(p)
      end
    end
    return SvN
end

function rec_ent_mps(psi,b)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    if length(psi)==2
        if dim(U,1) != dim(U,2)
            for n in diag(reshape(U.tensor,(dim(U,1),dim(U,2))))
                p = n^2
                if p>0.
                    SvN -= p * log2(p)
                end
            end
        else
        for n=1:dim(U, 1)
            p = U[n,n]^2
                if p>0.
                    SvN -= p * log2(p)
                end
            end
        end
    else
        for n=1:dim(S, 1)
            
        p = S[n,n]^2
            if p>0.
                SvN -= p * log2(p)
            end
        end
    end
    return SvN
end
  #
let #up
    #should return 0 for each case
    N=2
    s = siteinds("Qubit", N)
    psi = productMPS(s, "Up" )
    rho = outer(psi',psi)

    @show "UP",rec_ent(rho,1,s)
    @show "UP",rec_ent_mps(psi,2)
    nothing
end

let #bell
    #should return 1 for each case
    N=2
    s = siteinds("Qubit", N)
    psi1 = productMPS(s, "Up" )
    rho = outer(psi1',psi1)
    gates = ITensor[]

    s1 = s[1]
    s2 = s[2]
    hj = op("H",s1)
    push!(gates, hj)
    hj = op("CNOT",[s1,s2])
    push!(gates, hj)

    rho = apply(gates, rho; apply_dag=true, cutoff = 1E-8)
    psi = apply(gates,psi1)

    chk_psi = contract(psi)
    reshape(chk_psi.tensor,(2^N,1))

    @show "bell", rec_ent(rho,1,s)
    @show "bell",rec_ent_mps(psi,2)
    nothing
end


let#rand this is weird, but checked manually and rho is doing right thing
    # does depend on link dim, if dim =1 it works fine
    N=2
    s = siteinds("Qubit", N)
    psi = randomMPS(s,1)
    rho = outer(psi',psi)
    @show reshape(contract(psi).tensor,(4,1))
    @show rho_to_dense(rho,s)

    @show "rand", rec_ent(rho,1,s)
    @show "rand", rec_ent_mps(psi,2)
    nothing
end
let #just another state
    N = 4
    sites = siteinds("S=1/2", N)
    states = [isodd(n) ? "Up" : "Dn" for n=1:N]
    psi = MPS(ComplexF64, sites, states)
    rho = outer(psi',psi)
    
    @show "rand", rec_ent(rho,1,sites)
    @show "rand", rec_ent_mps(psi,2)
    
end
let 
    N = 6
    s = siteinds("S=1/2",N)
    rho = MPO(s, "Id")
    rho=rho/tr(rho)

    return rec_ent(rho,3,s)
end
let
    N = 4
    s = siteinds("S=1/2",N)
    A = ITensor(s...)
    A[1,1] = 1/√2
    A[2,2] = 1/√2
    @show A
    psiA = MPS(A,s)
    @show reshape(contract(psiA).tensor,(2^N,1))
    @show norm(psiA)
    @show rec_ent_mps(psiA,2)
end

let #nothing test
    N = 4
    s = siteinds("S=1/2",N)
    A = ITensor(s...)

    A[1,1,1,1] = 1/√2
    A[1,2,2,1] = 1/√2
    psiA = MPS(A,s)
    @show reshape(contract(psiA).tensor,(2^N,1))
    @show norm(psiA)
    @show rec_ent_mps(psiA,2)
end

let #check entropy functions for bell state
    N=4
    s = siteinds("Qubit", N)
    psi1 = productMPS(s, "Up" )
    rho = outer(psi1',psi1)
    gates = ITensor[]

    s1 = s[2]
    s2 = s[3]
    hj = op("H",s1)
    push!(gates, hj)
    hj = op("CNOT",[s1,s2])
    push!(gates, hj)

    rho = apply(gates, rho; apply_dag=true, cutoff = 1E-8)
    psi = apply(gates,psi1)

    chk_psi = contract(psi)
    reshape(chk_psi.tensor,(2^N,1))

    @show("rho ent: ", rec_ent(rho,2,s))
    @show("psi ent: ",rec_ent_mps(psi,2))
    nothing
end

let # this no work psi isn't nicely normalized
    N=2
    s = siteinds("Qubit", N)
    psi = randomMPS(s)
    rho = outer(psi',psi)
    
    @show("rho ent: ", rec_ent(rho,1,s))
    @show tr(rho)
    @show("psi ent: ",rec_ent_mps(psi,2))
    nothing
end

N=2
s = siteinds("Qubit", N)
psi = randomMPS(s,1)
rho = outer(psi',psi)
a=rho_to_dense(rho,s)
@show reshape(contract(psi).tensor,(4,1))
@show rho_to_dense(rho,s)

@show "rand", rec_ent(rho,1,s)
@show "rand", rec_ent_mps(psi,1)

function ptr(a)

    return[[a[1,1]+a[2,2] a[3,1]+a[4,2]] ;[a[3,1]+a[4,2] a[3,3]+a[4,4]]]
end

let #computational state
    N = 4
    sites = siteinds("S=1/2", N)
    states = [isodd(n) ? "Up" : "Dn" for n=1:N]
    psi = MPS(ComplexF64, sites, states)
    rho = outer(psi',psi)
    @show reshape(contract(psi).tensor,(16,1))

    @show "rand", rec_ent(rho,1,sites)
    @show "rand", rec_ent_mps(psi,2)
    nothing
end
#vectorization

#start with bell state
N=2
s = siteinds("Qubit", N)
psi1 = productMPS(s, "Up" )
rho = outer(psi1',psi1)
gates = ITensor[]

s1 = s[1]
s2 = s[2]
hj = op("H",s1)
push!(gates, hj)
hj = op("CNOT",[s1,s2])
push!(gates, hj)

rho = apply(gates, rho; apply_dag=true, cutoff = 1E-8)
psi = apply(gates,psi1)

#get inds and combine them
# test_tens = rho[1]
# (a,b,c)=inds(test_tens)
# C=combiner(a,b)

# test_tensa= C*test_tens

sites = siteinds("Qudit", N; dim=4)
#get list of contracted rho[i]s
arr=[]
for i in 1:length(rho)
    test_tens = rho[i]
    site_inds=inds(test_tens)
    (a,b)=site_inds[[occursin("Site,n=",string(tags(i))) for i in site_inds]]
    @show (a,b)
    @show tags(sites[i])
    C=combiner(a,b; tags=tags(sites[i]))

    append!(arr,[ITensor(C*test_tens)])
end

N = 2
# sites = siteinds("Qudit", N; dim=4)
fin =MPS(sites)
for i in 1:length(fin)
    fin[i]=arr[i]
end

function vectorize(rho,s)
    arr=[]
    for i in 1:length(rho)
        test_tens = rho[i]
        site_inds=inds(test_tens)
        (a,b)=site_inds[[occursin("Site,n=",string(tags(i))) for i in site_inds]]
        C=combiner(a,b; tags=tags(s[i]))
    
        append!(arr,[ITensor(C*test_tens)])
    end
    
    fin =MPS(s)
    for i in 1:length(fin)
        fin[i]=arr[i]
        ind=inds(fin[i])[[occursin("Site,n=",string(tags(i))) for i in inds(fin[i])]]
        swapinds!(fin[i],ind[1],s[i])
    end
    return fin
end
vectorize(rho,sites)



arr=[]
for i in 1:length(rho)
    test_tens = rho[i]
    site_inds=inds(test_tens)
    (a,b)=site_inds[[occursin("Site,n=",string(tags(i))) for i in site_inds]]
    C=combiner(a,b; tags=tags(s[i]))

    append!(arr,[ITensor(C*test_tens)])
end

# sites = new_sites#siteinds("Qudit", N; dim=4)
fin =MPS(sites)
for i in 1:length(fin)
    #find size of link needed
    linksize=dim([arr[i]...]/4)
    #get link index

    fin[i]=ITensor([arr[i]...],inds(fin[i]))
end
return fin
show(reshape(prod(fin).tensor,(16,1)))
show(reshape(prod(vectorize(rho)).tensor,(16,1)))
## SIIIIIICK ok lets get log negativity


function partial_transpose(A::MPO, sites)
    A = copy(A)
    for n in sites
      A[n] = swapinds(A[n], siteinds(A, n)...)
    end
    return A
  end
function trace_norm_dense(rho)
    e,_=eigen(rho)
    return sum(abs.(e))
end
function negativity(A::MPO, b, s)
    n = length(A)
    orthogonalize!(A,b)
    rho_temp = deepcopy(A)
    # s = siteinds("Qubit",n) 

    # #contract half   x x x x | | | |
    # L = ITensor(1.0)
    # for i = 1:b
    #   L *= tr(rho_temp[i])
    # end
    # # absorb
    # rho_temp[b+1] *= L
    # # no longer a proper mpo
    # M =MPO(n-b)
    # for i in 1:(n-b)
    #     M[i]=rho_temp[b+i]
    # end
    
    M=partial_transpose(rho_temp,[b:n-b...])

    #turn this mpo into a single tensor

    T = rho_to_dense(M,s)
    
    
    return  (trace_norm_dense(T)-1)/2
end

function log_negativity(A::MPO, b, s)
    neg = negativity(A, b, s)
    return log2(2*neg+1)
end



N=2
s = siteinds("Qubit", N)
psi1 = productMPS(s, "Up" )
rho = outer(psi1',psi1)
gates = ITensor[]
negativity(rho,1,s)
log_negativity(rho,1,s)

s1 = s[1]
s2 = s[2]
hj = op("H",s1)
push!(gates, hj)
hj = op("CNOT",[s1,s2])
push!(gates, hj)

rho = apply(gates, rho; apply_dag=true, cutoff = 1E-8)
psi = apply(gates,psi1)

negativity(rho,1,s)
log_negativity(rho,1,s)

# ok that kinda works

N=8
s = siteinds("Qubit", N)
psi1 = productMPS(s, "Up" )
rho = outer(psi1',psi1)
gates = ITensor[]

negativity(rho,4,s)
log_negativity(rho,4,s)

h = op("H",[s[1]])
h1 = op("H",[s[3]])
h2 = op("H",[s[5]])
h3 = op("H",[s[7]])
push!(gates, h)
push!(gates, h1)
push!(gates, h2)
push!(gates, h3)
h = op("CNOT",[s[1],s[2]])
h1 = op("CNOT",[s[3],s[4]])
h2 = op("CNOT",[s[5],s[6]])
h3 = op("CNOT",[s[7],s[8]])
push!(gates, h)
push!(gates, h1)
push!(gates, h2)
push!(gates, h3)

cutoff = 1E-8

rho = apply(gates, rho; apply_dag=true)
psi = apply(gates,psi1)

negativity(rho,4,s)
log_negativity(rho,4,s)




###################################
# now do superoperator

let #base plot 
    N = 6
    cutoff = 1E-8
    tau = 0.5
    ttotal = 25.0
  
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)
  
    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        op("Sx", s1)*op("I",s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  
    c = div(N, 2) # center site
  
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    ent_z=[]
    for t in 0.0:tau:ttotal
      append!(ent_z,rec_ent_mps(psi,3))
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
    return plot([0.0:tau:ttotal...],ent_z, linewidth=3)
end

builder = qtn.SpinHam1D(S=3/2)
#hmmmmmm
j=1
builder += 1/4, kron(pauli('X'),np.identity(2)),kron(pauli('X'),np.identity(2))
builder += -1/4, kron(np.identity(2),pauli('X').T),kron(np.identity(2),pauli('x').T)

builder += 1/2, kron(pauli('z'),np.identity(2))
builder += -1/2, kron(np.identity(2),pauli('z').T)

let #vectorized plot 
    N = 6
    cutoff = 1E-8
    tau = 0.5
    ttotal = 25.0
  
    # Make an array of 'site' indices
    s = siteinds("Qudit", N; dim=4)
    sm = siteinds("S=1/2", N)

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = sm[j]
      s2 = sm[j + 1]
      hj =
      1op(kron(op("Sz",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("Sz",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sz",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,transpose(op("Sz",s2).tensor)),s[j+1]) +
      1*op(kron(op("Sx",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sx",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1])
    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")
    rho = outer(psi1',psi1)
    psi = vectorize(rho,s)
  
    c = div(N, 2) # center site
  
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    ent_z=[]
    for t in 0.0:tau:ttotal
      append!(ent_z,rec_ent_mps(psi,3))
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
    return plot([0.0:tau:ttotal...],ent_z, linewidth=4)
end

#ok thats running need to pull out vectorized entropy in hopefully a smart way.

#create mpo
asd=MPO(sm)
#need to split indicies
#get data
#create mpo
#fill things

function vect_to_rho(vect,sm,s)
    arr=[]
    rho = MPO(sm)
    for i in 1:length(rho)
        test_tens = rho[i]
        site_inds=inds(test_tens)
        (a,b)=site_inds[[occursin("Site,n=",string(tags(i))) for i in site_inds]]
        C=combiner(a,b; tags=tags(s[i]))
        swapinds!(C,inds(C)[1],s[i])
        
        append!(arr,[ITensor(dag(C)*vect[i])])
    end
        # sites = new_sites#siteinds("Qudit", N; dim=4)
    for i in 1:length(rho)
        rho[i]=arr[i]
    end
    return rho
end

let #up state check vect un vect check
    N=2
    s = siteinds("Qubit", N)
    sl = siteinds("Qudit", N,dim=4)

    psi1 = productMPS(s, "Up" )
    rho_t = outer(psi1',psi1)

    hyup=vectorize(rho_t,sl)
    hydn=vect_to_rho(hyup,s,sl)
    @show rho_to_dense(rho_t,s)
    @show rho_to_dense(rho_t,s)
    nothing
end
let #random state check vect un vect check
    N=2
    s = siteinds("Qubit", N)
    sl = siteinds("Qudit", N,dim=4)

    psi1 = randomMPS(s )
    rho_t = outer(psi1',psi1)

    hyup=vectorize(rho_t,sl)
    hydn=vect_to_rho(hyup,s,sl)
    @show rho_to_dense(rho_t,s)
    @show rho_to_dense(rho_t,s)
    nothing
end
let #bell state check vect un vect check
    N=2
    s = siteinds("Qubit", N)
    psi1 = productMPS(s, "Up" )
    rho_t = outer(psi1',psi1)
    gates = ITensor[]

    s1 = s[1]
    s2 = s[2]
    hj = op("H",s1)
    push!(gates, hj)
    hj = op("CNOT",[s1,s2])
    push!(gates, hj)

    rho_t = apply(gates, rho_t; apply_dag=true, cutoff = 1E-8)
    sl = siteinds("Qudit", N,dim=4)

    hyup=vectorize(rho_t,sl)
    hydn=vect_to_rho(hyup,s,sl)
    @show rho_to_dense(rho_t,s)
    @show rho_to_dense(rho_t,s)
    nothing
end

let #test vectorize then entropy
    N=2
    s = siteinds("Qubit", N)
    psi1 = productMPS(s, "Up" )
    rho_t = outer(psi1',psi1)
    gates = ITensor[]

    s1 = s[1]
    s2 = s[2]
    hj = op("H",s1)
    push!(gates, hj)
    hj = op("CNOT",[s1,s2])
    push!(gates, hj)

    rho_t = apply(gates, rho_t; apply_dag=true, cutoff = 1E-8)
    psi = apply(gates,psi1)

    sl = siteinds("Qudit", N,dim=4)

    hyup=vectorize(rho_t,sl)
    hydn=vect_to_rho(hyup,s,sl)
    # @show rho_to_dense(rho_t,s)
    # @show rho_to_dense(hydn,s)

    @show("psi ent: ",rec_ent_mps(psi,2))
    @show("rho ent: ", rec_ent(rho_t,1,s))
    @show("rho ent: ", rec_ent(hydn,1,s))

    
    nothing
end

let #check ent values after some evolution
    N = 6
    cutoff = 1E-8
    tau = 1.
    ttotal = 10.0
  
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)
  
    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        op("Sx", s1)*op("I",s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  
    c = div(N, 2) # center site
  
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    ent_base=[]
    for t in 0.0:tau:ttotal
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
    psi_base = psi
    @show psi_base
    @show rec_ent_mps(psi_base,Int(N/2))
    # @show display(round.(rho_to_dense(outer(psi_base',psi_base),s),digits=3))
    @show rec_ent(outer(psi_base',psi_base),Int(N/2),s)

    # Make an array of 'site' indices
    sl = siteinds("Qudit", N; dim=4)
    # sm = siteinds("S=1/2", N)

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
      1op(kron(op("Sz",s1).tensor,op("I",s1).tensor),sl[j])* op(kron(op("Sz",s2).tensor,op("I",s2).tensor),sl[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sz",s1).tensor)),sl[j]) * op(kron(op("I",s2).tensor,transpose(op("Sz",s2).tensor)),sl[j+1]) +
      1*op(kron(op("Sx",s1).tensor,op("I",s1).tensor),sl[j])* op(kron(op("I",s2).tensor,op("I",s2).tensor),sl[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sx",s1).tensor)),sl[j]) * op(kron(op("I",s2).tensor,op("I",s2).tensor),sl[j+1])
    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi1 = MPS(s, n -> isodd(n) ? "Up" : "Dn")
    rho = outer(psi1',psi1)
    psi = vectorize(rho,sl)
    
    c = div(N, 2) # center site
  
    for t in 0.0:tau:ttotal
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end

    rho_end_vec = rho_to_dense(vect_to_rho(psi,s,sl),s)
    # @show display(round.(rho_end_vec,digits=3))
    @show rec_ent_mps(psi,Int(N/2))/2
    @show rec_ent(vect_to_rho(psi,s,sl),Int(N/2)-1,s)
    @show tr(rho_end_vec)
    @show "purity: ", tr(rho_end_vec^2)

    
    return nothing
end

let #comparison plot 
    N = 6
    cutoff = 1E-8
    tau = 0.5
    ttotal = 25.0
  
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)
  
    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        op("Sx", s1)*op("I",s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  
    c = div(N, 2) # center site
  
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    ent_base=[]
    for t in 0.0:tau:ttotal
      append!(ent_base,rec_ent_mps(psi,3))
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
  
    # Make an array of 'site' indices
    s = siteinds("Qudit", N; dim=4)
    sm = siteinds("S=1/2", N)

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = sm[j]
      s2 = sm[j + 1]
      hj =
      1op(kron(op("Sz",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("Sz",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sz",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,transpose(op("Sz",s2).tensor)),s[j+1]) +
      1*op(kron(op("Sx",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sx",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1])
    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")
    rho = outer(psi1',psi1)
    psi = vectorize(rho,s)
    @show psi
    c = div(N, 2) # center site
  
    ent_z=[]
    ent_vec=[]
    for t in 0.0:tau:ttotal
      append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),3,sm))
      append!(ent_vec,rec_ent_mps(psi,3))
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
    return plot([0.0:tau:ttotal...],[ent_z,ent_vec,ent_base], linewidth=3)
end

#############
#great thats working now just need to get dephase going.

let#plot of different parts of the vect form of rho
    N=2
    cutoff = 1E-8
    tau = 0.1
    ttotal = 10

    s = siteinds("Qubit", N)
    psi1 = productMPS(s, "Up" )
    rho_t = outer(psi1',psi1)
    gates = ITensor[]

    s1 = s[1]
    s2 = s[2]
    hj = op("H",s1)
    push!(gates, hj)
    hj = op("CNOT",[s1,s2])
    push!(gates, hj)

    rho_t = apply(gates, rho_t; apply_dag=true, cutoff = 1E-8)
    sl = siteinds("Qudit", N,dim=4)
    hyup=vectorize(rho_t,sl)

    gates = ITensor[]
    for j in 1:(N)
    gamma = √(1/2)*im
    s1 = s[j]
    pz = op("Sz",s1).tensor
    pI=op("I",s1).tensor
    hj =gamma*op(kron(pz,pz),sl[j])+
    -1/2*gamma*op(kron(pz*pz,pI),sl[j])+
    -1/2*gamma*op(kron(pI,pz*pz),sl[j])

    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
    end
    append!(gates, reverse(gates))


    @show rho_to_dense(vect_to_rho(hyup,s,sl),s)

    recc =[]
    recl =[]
    for t in 0.0:tau:ttotal
    #   println("$t $Sz")
        append!(recc,rho_to_dense(vect_to_rho(hyup,s,sl),s)[1,4])
        append!(recl,rho_to_dense(vect_to_rho(hyup,s,sl),s)[1,1])

        t≈ttotal && break

        hyup = apply(gates, hyup; cutoff)
        normalize!(hyup)
    end
    @show rho_to_dense(vect_to_rho(hyup,s,sl),s)

    return plot(0.0:tau:ttotal,[recc,recl])

    # tr(rho_to_dense(vect_to_rho(hyup,s,sl),s))

end

let #comparison plot with dephase
    N = 6
    cutoff = 1E-8
    tau = 0.5
    ttotal = 50.0
  
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)
  
    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        op("Sx", s1)*op("I",s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  
    c = div(N, 2) # center site
  
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    ent_base=[]
    for t in 0.0:tau:ttotal
      append!(ent_base,rec_ent_mps(psi,3))
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
  
    # Make an array of 'site' indices
    s = siteinds("Qudit", N; dim=4)
    sm = siteinds("S=1/2", N)

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = sm[j]
      s2 = sm[j + 1]
      hj =
      1op(kron(op("Sz",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("Sz",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sz",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,transpose(op("Sz",s2).tensor)),s[j+1]) +
      1*op(kron(op("Sx",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sx",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1])
    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
    end
    for j in 1:(N)
        gamma = √(.001/2)*im
        s1 = sm[j]
        pz = op("Sz",s1).tensor
        pI=op("I",s1).tensor
        hj =gamma*op(kron(pz,pz),s[j])+
        -1/2*gamma*op(kron(pz*pz,pI),s[j])+
        -1/2*gamma*op(kron(pI,pz*pz),s[j])
    
        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))



    # Initialize psi to be a product state (alternating up and down)
    psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")
    rho = outer(psi1',psi1)
    psi = vectorize(rho,s)
  
    ent_z=[]
    ent_vec=[]
    for t in 0.0:tau:ttotal
      append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),3,sm))
      append!(ent_vec,rec_ent_mps(psi,3))
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
    #   normalize!(psi)
    end
    return plot([0.0:tau:ttotal...],[ent_z,ent_vec,ent_base], linewidth=3)
end

let #comparison plot with logneg
    N = 8
    cutoff = 1E-8
    tau = 0.1
    ttotal = 50.0
  
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)
  
    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        op("Sx", s1)*op("I",s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  
    c = div(N, 2) # center site
  
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    ent_base=[]
    neg_base=[]
    for t in 0.0:tau:ttotal
      append!(ent_base,rec_ent_mps(psi,Int(N/2)))
      append!(neg_base,log_negativity(outer(psi',psi),Int(N/2),s))

    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
    @show("base is done")
    # Make an array of 'site' indices
    s = siteinds("Qudit", N; dim=4)
    sm = siteinds("S=1/2", N)

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = sm[j]
      s2 = sm[j + 1]
      hj =
      1op(kron(op("Sz",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("Sz",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sz",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,transpose(op("Sz",s2).tensor)),s[j+1]) +
      1*op(kron(op("Sx",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1]) +
      -1*op(kron(op("I",s1).tensor,transpose(op("Sx",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1])
    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
    end
    for j in 1:(N)
        gamma = √(.00001/2)*im
        s1 = sm[j]
        pz = op("Sz",s1).tensor
        pI=op("I",s1).tensor
        hj =gamma*op(kron(pz,pz),s[j])+
        -1/2*gamma*op(kron(pz*pz,pI),s[j])+
        -1/2*gamma*op(kron(pI,pz*pz),s[j])
    
        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))



    # Initialize psi to be a product state (alternating up and down)
    psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")
    rho = outer(psi1',psi1)
    psi = vectorize(rho,s)
  
    ent_z=[]
    neg_z=[]
    for t in 0.0:tau:ttotal
      append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
      append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
    #   println("$t $Sz")
  
      t≈ttotal && break
  
      psi = apply(gates, psi; cutoff)
    #   normalize!(psi)
    end
    return plot([0.0:tau:ttotal...],[ent_z,neg_z,ent_base,neg_base], linewidth=3)
end

let #randomintial states
    arr_base_ent =[]
    arr_base_neg=[]
    arr_deph_ent=[]
    arr_deph_neg=[]
    N = 6
    cutoff = 1E-8
    tau = 0.2
    ttotal = 40.0
    for _ in 0:60
    
        # Make an array of 'site' indices
        s = siteinds("S=1/2", N)
        
        # Make gates (1,2),(2,3),(3,4),...
        gates = ITensor[]
        for j in 1:(N - 1)
        s1 = s[j]
        s2 = s[j + 1]
        hj =
            op("Sz", s1) * op("Sz", s2) +
            op("Sx", s1)*op("I",s2)
        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
        end
        # Include gates in reverse order too
        # (N,N-1),(N-1,N-2),...
        append!(gates, reverse(gates))
    
        # Initialize psi to be a product state (alternating up and down)
        psi = randomMPS(s)#MPS(s, n -> isodd(rand(Int)) ? "Up" : "Dn")
    
    
        # Compute and print <Sz> at each time step
        # then apply the gates to go to the next time
        ent_base=[]
        neg_base=[]
        for t in 0.0:tau:ttotal
        append!(ent_base,rec_ent_mps(psi,Int(N/2)))
        append!(neg_base,log_negativity(outer(psi',psi),Int(N/2),s))

        #   println("$t $Sz")
    
        t≈ttotal && break
    
        psi = apply(gates, psi; cutoff)
        normalize!(psi)
        end
        @show("base is done")
        # Make an array of 'site' indices
        append!(arr_base_ent,[ent_base])
        append!(arr_base_neg,[neg_base])

    #     s = siteinds("Qudit", N; dim=4)
    #     sm = siteinds("S=1/2", N)

    # # Make gates (1,2),(2,3),(3,4),...

    #     gates = ITensor[]
    #     for j in 1:(N - 1)
    #     s1 = sm[j]
    #     s2 = sm[j + 1]
    #     hj =
    #     1op(kron(op("Sz",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("Sz",s2).tensor,op("I",s2).tensor),s[j+1]) +
    #     -1*op(kron(op("I",s1).tensor,transpose(op("Sz",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,transpose(op("Sz",s2).tensor)),s[j+1]) +
    #     1*op(kron(op("Sx",s1).tensor,op("I",s1).tensor),s[j])* op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1]) +
    #     -1*op(kron(op("I",s1).tensor,transpose(op("Sx",s1).tensor)),s[j]) * op(kron(op("I",s2).tensor,op("I",s2).tensor),s[j+1])
    #     Gj = exp(-im * tau / 2 * hj)
    #     push!(gates, Gj)
    #     end
    #     for j in 1:(N)
    #         gamma = √(.01/2)*im
    #         s1 = sm[j]
    #         pz = op("Sz",s1).tensor
    #         pI=op("I",s1).tensor
    #         hj =gamma*op(kron(pz,pz),s[j])+
    #         -1/2*gamma*op(kron(pz*pz,pI),s[j])+
    #         -1/2*gamma*op(kron(pI,pz*pz),s[j])
        
    #         Gj = exp(-im * tau / 2 * hj)
    #         push!(gates, Gj)
    #     end
    #     # Include gates in reverse order too
    #     # (N,N-1),(N-1,N-2),...
    #     append!(gates, reverse(gates))



    #     # Initialize psi to be a product state (alternating up and down)
    #     psi1 = MPS(sm, n -> isodd(rand(Int)) ? "Up" : "Dn")
    #     rho = outer(psi1',psi1)
    #     psi = vectorize(rho,s)
    
    #     ent_z=[]
    #     neg_z=[]
    #     for t in 0.0:tau:ttotal
    #     append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
    #     append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
    #     #   println("$t $Sz")
    
    #     t≈ttotal && break
    
    #     psi = apply(gates, psi; cutoff)
    #     #   normalize!(psi)
    #     end
    #     append!(arr_deph_ent,[ent_z])
    #     append!(arr_deph_neg,[neg_z])
    end
    @show( arr_base_ent)
    return plot([0.0:tau:ttotal...],[mean(arr_base_ent,dims=1)], linewidth=3)
end


let #random state gen
    sm = siteinds("S=1/2", 8)
    psi=MPS(sm, n -> isodd(rand(Int)) ? "Up" : "Dn")
    magz = expect(psi,"Sz")
    arr =[]
    for (j,mz) in enumerate(magz)
        append!(arr,mz)
    end
    @show arr
    nothing
end

let #randomintial states
    arr_base_ent =[]
    arr_base_neg=[]
    arr_deph_ent=[]
    arr_deph_neg=[]
    N = 4
    cutoff = 1E-8
    tau = 0.1
    ttotal = 10.0
    for _ in 0:1
    
        # # Make an array of 'site' indices
        # s = siteinds("S=1/2", N)
        
        # # Make gates (1,2),(2,3),(3,4),...
        # gates = ITensor[]
        # for j in 1:(N - 1)
        # s1 = s[j]
        # s2 = s[j + 1]
        # hj =
        #     op("Sz", s1) * op("Sz", s2) +
        #     op("Sx", s1)*op("I",s2)
        # Gj = exp(-im * tau / 2 * hj)
        # push!(gates, Gj)
        # end
        # # Include gates in reverse order too
        # # (N,N-1),(N-1,N-2),...
        # append!(gates, reverse(gates))
    
        # # Initialize psi to be a product state (alternating up and down)
        # psi = randomMPS(s)#MPS(s, n -> isodd(rand(Int)) ? "Up" : "Dn")
    
    
        # # Compute and print <Sz> at each time step
        # # then apply the gates to go to the next time
        # ent_base=[]
        # neg_base=[]
        # for t in 0.0:tau:ttotal
        # append!(ent_base,rec_ent_mps(psi,Int(N/2)))
        # append!(neg_base,log_negativity(outer(psi',psi),Int(N/2),s))

        # #   println("$t $Sz")
    
        # t≈ttotal && break
    
        # psi = apply(gates, psi; cutoff)
        # normalize!(psi)
        # end
        # @show("base is done")
        # # Make an array of 'site' indices
        # append!(arr_base_ent,[ent_base])
        # append!(arr_base_neg,[neg_base])

        s = siteinds("Qudit", N; dim=4)
        sm = siteinds("S=1/2", N)

    # Make gates (1,2),(2,3),(3,4),...

        gates = ITensor[]
        for j in 1:(N - 1)
        s1 = sm[j]
        s2 = sm[j + 1]
        Sz1=op("Sz",s1).tensor
        Sz2 = op("Sz",s2).tensor
        Sx1=op("Sx",s1).tensor
        Sx2 = op("Sx",s2).tensor
        SI1 = op("I",s1).tensor
        SI2 = op("I",s2).tensor
        gamma = √(0.0/2)*im

        hj =
        1*op(kron(Sz1,SI1),s[j])* op(kron(Sz2,SI2),s[j+1]) +
        -1*op(kron(SI1,transpose(Sz1)),s[j]) * op(kron(SI2,transpose(Sz2)),s[j+1]) +
        1*op(kron(Sx1,SI1),s[j])* op(kron(SI2,SI2),s[j+1]) +
        -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,SI2),s[j+1])+
        gamma*op(kron(Sz1,Sz1),s[j]) * op(kron(SI2,SI2),s[j+1])+
            -1/2*gamma*op(kron(Sz1*Sz1,SI1),s[j])* op(kron(SI2,SI2),s[j+1])+
            -1/2*gamma*op(kron(SI1,Sz1*Sz1),s[j])* op(kron(SI2,SI2),s[j+1])

        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
        end
        # for j in 1:(N)
        #     gamma = √(2/2)*im
        #     s1 = sm[j]
        #     pz = op("Sz",s1).tensor
        #     pI=op("I",s1).tensor
        #     hj =gamma*op(kron(pz,pz),s[j])+
        #     -1/2*gamma*op(kron(pz*pz,pI),s[j])+
        #     -1/2*gamma*op(kron(pI,pz*pz),s[j])
        
        #     Gj = exp(-im * tau / 2 * hj)
        #     push!(gates, Gj)
        # end
        # Include gates in reverse order too
        # (N,N-1),(N-1,N-2),...
        append!(gates, reverse(gates))



        # Initialize psi to be a product state (alternating up and down)
        psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")#rand(Int)) ? "Up" : "Dn")
        rho = outer(psi1',psi1)
        psi = vectorize(rho,s)
    
        ent_z=[]
        neg_z=[]
        for t in 0.0:tau:ttotal
        append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
        append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
        #   println("$t $Sz")
    
        t≈ttotal && break
    
        psi = apply(gates, psi; cutoff)
        #   normalize!(psi)
        end
        append!(arr_deph_ent,[ent_z])
        append!(arr_deph_neg,[neg_z])
    end
    @show( arr_base_ent)
    return plot([0.0:tau:ttotal...],[mean(arr_deph_ent,dims=1),mean(arr_deph_neg,dims=1)], linewidth=3)
end
let #randomintial states
    arr_base_ent =[]
    arr_base_neg=[]
    arr_deph_ent=[]
    arr_deph_neg=[]
    N = 8
    cutoff = 1E-8
    tau = 0.1
    ttotal = 20.0
    for _ in 1:1

        # Make an array of 'site' indices
        s = siteinds("S=1/2", N)
        
        # Make gates (1,2),(2,3),(3,4),...
        gates = ITensor[]
        for j in 1:(N - 1)
            s1 = s[j]
            s2 = s[j + 1]
            hj =
                op("Sz", s1) * op("Sz", s2) +
                op("Sx", s1)*op("I",s2)
            Gj = exp(-im * tau / 2 * hj)
            push!(gates, Gj)
        end
        # Include gates in reverse order too
        # (N,N-1),(N-1,N-2),...
        append!(gates, reverse(gates))

        # Initialize psi to be a product state (alternating up and down)
        psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")


        # Compute and print <Sz> at each time step
        # then apply the gates to go to the next time
        ent_base=[]
        neg_base=[]
        for t in 0.0:tau:ttotal
            append!(ent_base,rec_ent_mps(psi,Int(N/2)))
            append!(neg_base,log_negativity(outer(psi',psi),Int(N/2),s))

            #   println("$t $Sz")

            t≈ttotal && break

            psi = apply(gates, psi; cutoff)
            normalize!(psi)
        end
        @show("base is done")
        # Make an array of 'site' indices
        append!(arr_base_ent,[ent_base])
        append!(arr_base_neg,[neg_base])

        s = siteinds("Qudit", N; dim=4)
        sm = siteinds("S=1/2", N)

        # Make gates (1,2),(2,3),(3,4),...

        gates = ITensor[]
        for j in 1:(N - 1)
            s1 = sm[j]
            s2 = sm[j + 1]
            Sz1=op("Sz",s1).tensor
            Sz2 = op("Sz",s2).tensor
            Sx1=op("Sx",s1).tensor
            Sx2 = op("Sx",s2).tensor
            SI1 = op("I",s1).tensor
            SI2 = op("I",s2).tensor
            gamma = √(0.1/2)*im

            hj =
            1*op(kron(Sz1,SI1),s[j])* op(kron(Sz2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sz1)),s[j]) * op(kron(SI2,transpose(Sz2)),s[j+1]) +
            1*op(kron(Sx1,SI1),s[j])* op(kron(SI2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,SI2),s[j+1])+

            gamma*op(kron(Sz1,Sz1),s[j]) * op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(Sz1*Sz1,SI1),s[j])* op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(SI1,Sz1*Sz1),s[j])* op(kron(SI2,SI2),s[j+1])

            Gj = exp(-im * tau / 2 * hj)
            push!(gates, Gj)
        end
        # for j in 1:(N)
        #     gamma = √(2/2)*im
        #     s1 = sm[j]
        #     pz = op("Sz",s1).tensor
        #     pI=op("I",s1).tensor
        #     hj =gamma*op(kron(pz,pz),s[j])+
        #     -1/2*gamma*op(kron(pz*pz,pI),s[j])+
        #     -1/2*gamma*op(kron(pI,pz*pz),s[j])
        
        #     Gj = exp(-im * tau / 2 * hj)
        #     push!(gates, Gj)
        # end
        # Include gates in reverse order too
        # (N,N-1),(N-1,N-2),...
        append!(gates, reverse(gates))



        # Initialize psi to be a product state (alternating up and down)
        psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")#rand(Int)) ? "Up" : "Dn")
        rho = outer(psi1',psi1)
        psi = vectorize(rho,s)

        ent_z=[]
        neg_z=[]
        for t in 0.0:tau:ttotal
            append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
            append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
            #   println("$t $Sz")

            t≈ttotal && break

            psi = apply(gates, psi; cutoff)
            #   normalize!(psi)
        end
        append!(arr_base_ent,[ent_base])
        append!(arr_base_neg,[neg_base])
        append!(arr_deph_ent,[ent_z])
        append!(arr_deph_neg,[neg_z])
    end
    # @show( arr_base_ent)
    labels= ["base ent" "base neg" "dephase ent" "dephase neg"]
    return plot([0.0:tau:ttotal...],[arr_base_ent[1],arr_base_neg[1],arr_deph_ent[1],arr_deph_neg[1]],
     linewidth=3,
     legend=true,
     label=labels,
     xlabel="Time",
     ylabel="Entropy / Negativity")

    # return plot([0.0:tau:ttotal...],[mean(arr_deph_ent,dims=1),mean(arr_deph_neg,dims=1),mean(arr_base_ent,dims=1),mean(arr_base_neg,dims=1)], linewidth=3)
end