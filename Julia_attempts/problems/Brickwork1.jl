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
function Ihaargate(N::Int)
    theta = rand()*2*pi
    eps=0.9#rand()
    u =(sqrt((1-eps)).*Matrix(I,4,4) + 1im.*sqrt(eps).*RandomUnitaryMatrix(4))
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

function rec_ent(psi,b)
    s = siteinds(psi)  
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
    SvN = 0.0
    for n in 1:dim(S, 1)
      p = S[n,n]^2
      if p != 0
        SvN -= p * log(p)
      end
    end
    return SvN
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
function ITensors.op(::OpName"expτSS", ::SiteType"Qubit", s1::Index, s2::Index; τ,B)
    h = op("Iden",s1)*op("Iden",s2)
        # 1 / 2 * op("S+", s1) * op("S-", s2) +
        # 1 / 2 * op("S-", s1) * op("S+", s2) +
        # op("Sz", s1) * op("Sz", s2)+
        # B*(op("Sz", s1) * op("I", s2)+op("I", s1) * op("Sz", s2))
        return exp(τ * h)
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


function metts(psi,sites)
    δτ=1.
    beta=2.
    NMETTS=1
    Nwarm=0
    b=0
    s = sites
    # Make gates (1,2),(2,3),(3,4),...
    gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,B=b)) for n in 1:(N - 1)], s)
    # Include gates in reverse order to complete Trotter formula
    append!(gates, reverse(gates))

    # Make y-rotation gates to use in METTS collapses
    Ry_gates = ops([("Ry", n, (θ=π / 2,)) for n in 1:N], s)


    # Make τ_range and check δτ is commensurate
    τ_range = δτ:δτ:(beta / 2)
    if norm(length(τ_range) * δτ - beta / 2) > 1E-10
        error("Time step δτ=$δτ not commensurate with beta/2=$(beta/2)")
    end

    for step in 1:(Nwarm + NMETTS)
        # if step <= Nwarm
        #   println("Making warmup METTS number $step")
        # else
        #   println("Making METTS number $(step-Nwarm)")
        # end

        # Do the time evolution by applying the gates
        for τ in τ_range
        cutoff = 1E-8
        psi = apply(gates, psi; cutoff)
        normalize!(psi)
        end

        # Measure properties after >= Nwarm 
        # METTS have been made
        # if step > Nwarm
        #   energy = inner(psi', H, psi)
        #   sz = inner(psi',Sz,psi)
        #   ent = entrp(psi,Int(N/2))
        #   push!(energies, energy)
        #   push!(magz,sz)
        #   push!(ents,ent)
        #   @printf("  Energy of METTS %d = %.4f\n", step - Nwarm, energy)
        # end

        # Measure in X or Z basis on alternating steps
        # if step % 2 == 1
        # psi = apply(Ry_gates, psi)
        # samp = sample!(psi)
        # new_state = [samp[j] == 1 ? "X+" : "X-" for j in 1:N]
        # else
        # samp = sample!(psi)
        # new_state = [samp[j] == 1 ? "Z+" : "Z-" for j in 1:N]
        # end
        # psi = productMPS(s, new_state)
    end

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
    cutoff = 1E-8

    psi = apply(gates, psi; cutoff)
    #calculate obs
    measured_vals = rec_ent(psi,Int(round(N/2)))
    
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

function do_exp(N,steps,meas_p)
    s = siteinds("Qubit", N)
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
    for i in [0,0.1,0.4]#0.0:0.4:1
        print("\n meas_p $i \n")
        svn,tri_mut =do_trials(N,steps,i,10)
        avgsvn = [(svn[x]+svn[x+1])/2 for x in 1:2:(size(svn)[1]-1)]
        append!(svns,[avgsvn])
        append!(mut,tri_mut)
    end
decay = [svns[i][end] for i in 1:size(svns)[1]]
append!(decays,[decay])
# end

N=10
p = plot(svns,title=string("Gate Rand", ", ", N, " qubit sites, varying meas_p"), label=string.(transpose([0.00:0.3:1...])), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
# p = plot([0.1:0.2:1...],decays,title=string("Bip_ent Gat: IX0.9", ", ", N, " qubit sites, varying meas_p"), label=string.(transpose([6:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
# m = plot(real(mut))
display(p)