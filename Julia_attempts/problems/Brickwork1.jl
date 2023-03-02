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

function tri_part_MI(psi_abc,sysa,sysb,sysc)
    hab = entropy_subsys(psi_abc, vcat(sysa,sysb))
    hac = entropy_subsys(psi_abc, vcat(sysa,sysc))
    hbc = entropy_subsys(psi_abc, vcat(sysc,sysb))
    
    habc = entropy_subsys(psi_abc, vcat(sysa,sysb,sysc))

    ha = entropy_subsys(psi_abc, sysa)
    hb = entropy_subsys(psi_abc, sysb)
    hc = entropy_subsys(psi_abc, sysb)
    return hb + ha - hab - hac - hbc + habc
end
ITensors.op(::OpName"Pup",::SiteType"Qubit") =
 [1 0
  0 0]
ITensors.op(::OpName"Pdn",::SiteType"Qubit") =
 [0 0
  0 1]
ITensors.op(::OpName"Rand",::SiteType"Qubit") = 
    RandomUnitaryMatrix(4)

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

function gen_step(psi,s,step_num,meas_p)
    #apply gates
    row = make_row(N,Bool(step_num%2),true)
    gates = ITensor[]
    for j in row
        s1 = s[j[1]]
        s2 = s[j[2]]
        hj = op("Rand",[s1,s2])
        Gj=hj
        push!(gates, Gj)
    end
    psi = apply(gates, psi; cutoff)
    #calculate obs
    measured_vals = rec_ent(psi,Int(round(N/2)))
    
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
        psi ,meas_svn= gen_step(psi,s,i,meas_p)
        append!(svn,meas_svn)
    end
    tri_mut = tri_part_MI(psi,[1,2],[3,4],[5,6])
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

N = 8
cutoff = 1E-8
steps = 20
meas_p=0.3
svns=[]
mut = []
for i in 0:1:0.3
    print("\n meas_p $i \n")
    svn,tri_mut =do_trials(N,steps,i,5)
    append!(svns,[svn])
    append!(mut,tri_mut)
end
p = plot(svns,title=string("Gate Rand", ", ", N, " qubit sites, varying meas_p"), label=string.(transpose([0.:0.1:0.3...])), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
display(p)