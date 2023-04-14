#This is bad need to rework with mpos




using ITensors
using Random
using Plots
using Statistics
using LaTeXStrings
using LinearAlgebra

function RandomUnitaryMatrix(N::Int)
    x = (rand(N,N) + rand(N,N)*im) / sqrt(2)
    f = qr(x)
    diagR = sign.(real(diag(f.R)))
    diagR[diagR.==0] .= 1
    diagRm = diagm(diagR)
    u = f.Q * diagRm
    
    return u
end
function IXgate(N::Int)
    theta = rand()*2*pi
    eps=0.25#rand()
    u =exp(1im*theta).*(sqrt((1-eps)).*Matrix(I,4,4) + 1im.*sqrt(eps).*[[0,0,0,1] [0,0,1,0] [0,1,0,0] [1,0,0,0]])
    return u
end

function make_row(N,eoo,pc)
    N=N #ignore ancilla from gen step
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

function trace_out_anc(psi)
    d=length(psi)
    s = siteinds("Qubit",d)

    # for i = length(psi)-1:length(psi)  #k is the site I stopped. or i = m:length(rho)
    psi[d] = psi[d]*delta(s[d],s[d]')
    # end
    m =MPS(length(psi)-1)
    for i=1:length(psi)-1
        m[i] = psi[i]
    end
    return m
end

function rec_ent(psi,b)
    psi_temp = deepcopy(psi)#trace_out_anc(psi)
    psi_temp = trace_out_anc(psi_temp)
    s = siteinds(psi)
    orthogonalize!(psi_temp, b)
    _,S = svd(psi_temp[b], (linkind(psi_temp, b-1), s[b]))
    SvN = 0.0
    for n in 1:dim(S, 1)
      p = S[n,n]^2
      if p != 0
        SvN -= p * log(p)
      end
    end
    return SvN
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

ITensors.op(::OpName"Iden",::SiteType"Qubit") = 
[1 0
0 1]

function ITensors.op(::OpName"expSWAP", ::SiteType"Qubit")

    h = 
    [1 0 0 0 
     0 0 1 0
     0 1 0 0
     0 0 0 1]
     #op("S+", s1) * op("S-", s2) +
    #     op("S-", s1) * op("S+", s2)
        # op("Iden", s1) * op("Iden", s2)
        # B*(op("Sz", s1) * op("I", s2)+op("I", s1) * op("Sz", s2))
    return exp(0.001* h)
end

function gen_samp_row(N,meas_p)
    return [rand()<meas_p ? 1 : 0 for i in 1:N]
end
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


# function ancilla_evolv(psi,N,sites)
#     # sites = (N-1,N)
#     s = sites
#     # Make gates (1,2),(2,3),(3,4),...
#     gates = ops([("expτSS", (n, n + 1), (τ=1,)) for n in [N]], s)
#     # append!(gates, reverse(gates))

#     cutoff = 1E-8
#     psi = apply(gates, psi; cutoff)
#     normalize!(psi)
#     return psi
# end

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
    cutoff = 1E-8

    psi = apply(gates, psi; cutoff)
    #calculate obs
    measured_vals = rec_ent(psi,Int(round(N/2)))
    
    #metts shenanagins
    gates = ITensor[]
    hj = op("SWAP", [s[N],s[N+1]])
    Gj=hj
    push!(gates,Gj)
    # psi = apply(gates,psi; cutoff)
    #psi = ancilla_evolv(psi,N,s)
    #sample as needed
    samp_row=gen_samp_row(N,meas_p)
    for (i,x) in enumerate(samp_row)
        if Bool(x)
            psi = samp_mps(psi,s,i)
        end
    end
    return psi,measured_vals
end

function do_exp(N,steps,meas_p)
    s = siteinds("Qubit", N+1) #+1 for ancilla
    psi = productMPS(s, "Up" )

    svn =[]
    for i in 1:steps
        psi ,meas_svn= gen_step(N,psi,s,i,meas_p)
        append!(svn,meas_svn)
    end
    #tri_mut = tri_part_MI(psi,[1,2],[3,4],[5,6])
    tri_mut = []
    # for i in 3:length(psi)-1
    #     arr = two_point_MI(psi,2,i)
    #     append!(tri_mut,arr)
    # end
    return svn,tri_mut
end

function do_trials(N,steps,meas_p,trials)
    svn_trials,tri_trials = do_exp(N,steps,meas_p)
    for i in 2:trials
        print(i)
        nSvn,ntm = do_exp(N,steps,meas_p)
        svn_trials = 2*mean([(i-1)/i*svn_trials,1/i*nSvn])
        tri_trials = 2*mean([(i-1)/i*tri_trials,1/i*ntm])

        
    end
    return svn_trials,tri_trials
end

decays=[]
# for n in [8,10]
N = 8
# cutoff = 1E-8
steps = 2*N
meas_p=0.
svns=[]
mut = []
for i in [0 0.1 0.3 0.8 1]#0.0:0.4:1
    print("\n meas_p $i \n")
    svn,tri_mut =do_trials(N,steps,i,50)
    avgsvn = [(svn[x]+svn[x+1])/2 for x in 1:2:(size(svn)[1]-1)]
    append!(svns,[avgsvn])
    append!(mut,tri_mut)
end
# decay = [svns[i][end] for i in 1:size(svns)[1]]
# append!(decays,[decay])
# end

# N=10
p = plot(svns,title=string("Gate Rand", ", ", N, " qubit sites, varying meas_p"), label=string.([0 0.1 0.3 0.8 1]), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
# p = plot([0.1:0.2:1...],decays,title=string("Bip_ent Gat: IX0.9", ", ", N, " qubit sites, varying meas_p"), label=string.(transpose([6:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
# m = plot(real(mut))
display(p)