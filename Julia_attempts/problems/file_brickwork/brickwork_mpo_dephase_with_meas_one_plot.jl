
using ITensors
using Random
using Plots
using Statistics
using LaTeXStrings
using LinearAlgebra
using PyCall

function rho_to_dense(rho,s)
    Hitensor = ITensor(1.)
    N = length(s)
    for i = 1:N
        Hitensor *= rho[i]
    end
  
    A=Array(Hitensor,prime(s),s)
    return reshape(A,2^N,2^N)
  end
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
    # normalize!(rho_temp)
    M=partial_transpose(rho_temp,[b:n-b...])

    #turn this mpo into a single tensor

    T = rho_to_dense(M,s)
    
    
    return  (trace_norm_dense(T)-1)/2
end

function log_negativity(A::MPO, b, s)
    neg = negativity(A, b, s)
    return log2(2*neg+1)
end
function rho_to_dense(rho,s)
  Hitensor = ITensor(1.)
  N = length(s)
  for i = 1:N
      Hitensor *= rho[i]
  end

  A=Array(Hitensor,prime(s),s)
  return reshape(A,2^N,2^N)
end

function kraus_dephase(rho,s,p)
    #define the two operators
    #(1-p)ρ + pZρZ
    N=length(rho)
    gates = ITensor[]
    for i in 1:N
      hj = op("Z", s[i])
      push!(gates, hj)
    end
    #apply the operators
    rho = (1-p)*rho + p*apply(gates,rho;apply_dag=true)
    return rho
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


function ren(rho)
  return -log(tr(apply(rho,rho)))
end
function split_ren(rho,b)
  n = length(rho)
  rho_temp = deepcopy(rho)
  s = siteinds("Qubit",n) 

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
  ren = -log2(tr(apply(M,M)))
  return ren
end

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
ITensors.op(::OpName"Iden2",::SiteType"Qubit") = 
  [1 0 0 0
  0 1 0 0
  0 0 1 0
  0 0 0 1]



function gen_samp_row(N,meas_p)
  return [rand()<meas_p ? 1 : 0 for i in 1:N]
end
function samp_mps(rho,s,samp_row)
    cutoff = 1E-8
    N = length(rho)
    samp =deepcopy(rho)
    samp = samp/tr(samp)
    samples= sample(samp)
    magz = [x == 1 ? "Pup" : "Pdn" for x in samples]
  
    gates = ITensor[]
  
    for i in 1:N
      if Bool(samp_row[i])
          hj = op(magz[i],s[i])
          push!(gates, hj)
      end
    end
    rho = apply(gates, rho;apply_dag=true)
   
    rho=rho/tr(rho)
    return rho
  end


function gen_step(N,rho,s,step_num,meas_p,noise)
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
  measured_vals = (rec_ent(rho,Int(round(N/2)),s),log_negativity(rho,Int(round(N/2)),s))
  rho = apply(gates, rho;apply_dag=true,cutoff=1E-8)

  #calculate obs
  rho=rho/tr(rho)
  #kraus_dephase
  rho = kraus_dephase(rho,s,noise)
  rho=rho/tr(rho)
  
  #sample as needed
    

  samp_row=gen_samp_row(N,meas_p)
  rho = samp_mps(rho,s,samp_row)
  rho=rho/tr(rho)
  
  
  return rho,measured_vals
end

function do_exp(N,steps,meas_p,noise)
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    rho=outer(psi',psi)
  
    svn =[]
    neg= []
    for i in 1:steps
        rho ,(meas_svn,meas_neg)= gen_step(N,rho,s,i,meas_p,noise)
        append!(svn,meas_svn)
        append!(neg,meas_neg)
      #   @show(tr(rho))
    end
    #tri_mut = tri_part_MI(psi,[1,2],[3,4],[5,6])
    
    # for i in 3:length(psi)-1
    #     arr = two_point_MI(psi,2,i)
    #     append!(tri_mut,arr)
    # end
    return svn,neg
  end

function do_trials(N,steps,meas_p,trials,noise)
    svn_trials,tri_trials = do_exp(N,steps,meas_p,noise)
    for i in 2:trials
        print(i)
        nSvn,ntm = do_exp(N,steps,meas_p,noise)
        svn_trials = 2*mean([(i-1)/i*svn_trials,1/i*nSvn])
        tri_trials = 2*mean([(i-1)/i*tri_trials,1/i*ntm])
  
        
    end
    return svn_trials,tri_trials
end

function main(meas_ps=[0.0:0.2:1...],trials=10,noise=0.0)
  decays=[]
  svns=[]
  negs=[]
  for n in [4]
  # N = 6
  # cutoff = 1E-8
  steps = 8*n
  
  mut = []
  for i in meas_ps
      print("\n meas_p $i \n")
      svn,neg =do_trials(n,steps,i,trials,noise)
      avgsvn = [(svn[x]+svn[x+1])/2 for x in 1:2:(size(svn)[1]-1)]
      avgneg = [(neg[x]+neg[x+1])/2 for x in 1:2:(size(neg)[1]-1)]

      append!(svns,[avgsvn])
      append!(negs,[avgneg])
    #   append!(mut,tri_mut)
  end
  decay = [svns[i][end] for i in 1:size(svns)[1]]
  append!(decays,[decay])
  end
  p = plot(svns,title=string("MPO Gate Rand qubit sites, varying meas_p"), label=string.(transpose(meas_ps)), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
  p2= plot(negs,title=string("MPO Gate Rand qubit sites, varying meas_p"), label=string.(transpose(meas_ps)), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{N_{vn}}(L/2)$")
  # p = plot(meas_ps,decays[end-2:end],title=string("Bip_ent Gat: 2Haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
  # m = plot(real(mut))
  display(p)
  display(p2)
end

function main_noise(meas_ps=0,trials=20,noise=[0.0:0.1:0.5...])
    decays=[]
    svns=[]
    negs=[]
    for n in [4]
    # N = 6
    # cutoff = 1E-8
    steps = 8*n
    
    mut = []
    for i in noise
        print("\n noise_p $i \n")
        svn,neg =do_trials(n,steps,0.0,trials,i)
        avgsvn = [(svn[x]+svn[x+1])/2 for x in 1:2:(size(svn)[1]-1)]
        avgneg = [(neg[x]+neg[x+1])/2 for x in 1:2:(size(neg)[1]-1)]
  
        append!(svns,[avgsvn])
        append!(negs,[avgneg])
      #   append!(mut,tri_mut)
    end
    decay = [svns[i][end] for i in 1:size(svns)[1]]
    append!(decays,[decay])
    end
    p = plot(svns,title=string("MPO Gate Rand qubit sites, varying noise"), legend=:outertopright,label=string.(transpose(noise)), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
    p2= plot(negs,title=string("MPO Gate Rand qubit sites, varying noise"), legend=:outertopright,label=string.(transpose(noise)), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{N_{vn}}(L/2)$")
    # p = plot(meas_ps,decays[end-2:end],title=string("Bip_ent Gat: 2Haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
    # m = plot(real(mut))
    display(p)
    display(p2)
  end
function main2(meas_ps=[0.1:0.1:1...],trials=5,noise=[0.01:0.05:0.5...])

    decays=[]
    growths=[]
    svns=[]
  muts = []
  ratio=[]
#   labels=[]
#   labels2=[]
  n=4
  steps=4*n
  for meas_p in meas_ps
  # N = 6
  # cutoff = 1E-8
  
    dec=[]
    gro=[]
    for i in noise
        print("\n meas_p $meas_p noise $i \n")
        svn,neg =do_trials(n,steps,meas_p,trials,i)
        #   @show(svn)
        #   @show(neg)
        avgsvn = [(svn[x]+svn[x+1])/2 for x in 1:2:(size(svn)[1]-1)]
        avgneg = [(neg[x]+neg[x+1])/2 for x in 1:2:(size(neg)[1]-1)]
        @show(meas_p/i)
        append!(ratio, i/meas_p)
        append!(svns,[avgsvn])
        append!(muts,[avgneg])
        #   push!(labels,"base ent $i")
        #   push!(labels,"base neg $i")
        append!(decays,avgsvn[end])
        append!(growths,avgneg[end])
    end
  end
#   @show(ratio)
#   @show(growths)
#   @show(decays)
  labels=["neg" "ent"]
    p2 = plot(ratio,[growths decays],
        title=string("MPO Gate Rand qubit sites, varying meas_p"),
        labels=labels,
        linewidth=3,
        seriestype=:scatter,
        xlabel = "ratio γ/p",
        xaxis=:log,

        ylabel = L"$\textbf{S_{vn}}(L/2)$",
        )

  display(p2)
end

function main_proj(meas_ps=[0:0.2:1...],trials=10,noise=[0.0:0.05:0.5...])

    decays=[]
    growths=[]
    svns=[]
  muts = []
  labels=[]
  labels2=[]
  labels3=[]
  n=6
  for meas_p in meas_ps
    steps=60
  # N = 6
  # cutoff = 1E-8
  
    dec=[]
    gro=[]
  for i in noise
      print("\n meas_p $i \n")
      svn,neg =do_trials(n,steps,meas_p,trials,i)
    #   @show(svn)
    #   @show(neg)
      avgsvn = [(svn[x]+svn[x+1])/2 for x in 1:2:(size(svn)[1]-1)]
      avgneg = [(neg[x]+neg[x+1])/2 for x in 1:2:(size(neg)[1]-1)]

      append!(svns,[avgsvn])
      append!(muts,[avgneg])
      push!(labels,"base ent $i")
      push!(labels,"base neg $i")
      append!(dec,avgsvn[end])
      append!(gro,avgneg[end])
  end

  append!(growths,[gro])
  append!(decays,[dec])
  push!(labels3, "$meas_p")
  push!(labels2,"base ent $meas_p")
  push!(labels2,"base neg $meas_p")
#   append!(decays,[decay])
#   append!(growths,[growth])
  #ok i want last value of ent/neg for different noise values
  end

#   svns=[svns[end]]
#   muts=[muts[end]]
#   @show(svns)
  
#   plot_data = collect(Iterators.flatten(zip(svns,muts)))
#   plot_data2 = collect(Iterators.flatten(zip(decays,growths)))
# #   p = plot(1:16,[plot_data[end-1:end]...],
# #             title=string("MPO Gate Rand qubit sites, varying meas_p"),
# #             legend=:outertopright,
# #             label=permutedims(labels),
# #             linewidth=3,
# #             xlabel = "Steps",
# #             ylabel = L"$\textbf{S_{vn}}(L/2)$",
# #             )
#     @show svns
#     p2 = plot(noise,decays,
#         title=string("MPO Gate Rand qubit sites, varying meas_p"),
#         legend=:outertopright,
#         label=permutedims(labels3),
#         linewidth=3,
#         xaxis=:log,
#         xlabel = "γ",
#         ylabel = L"$\textbf{S_{vn}}(L/2)$",
#         palette=:batlow10
#         )
#     p3 = plot(noise,growths,
#         title=string("MPO Gate Rand qubit sites, varying meas_p"),
#         legend=:outertopright,
#         label=permutedims(labels3),
#         linewidth=3,
#         xaxis=:log,
#         xlabel = "γ",
#         ylabel = L"$\textbf{N_{vn}}(L/2)$",
#         palette=:batlow10
#         )
# #   p2 = plot(muts,title=string("MPO Gate Rand qubit sites, varying meas_p"), label=string.(transpose(noise)), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{N_{vn}}(L/2)$")
#   # p = plot(meas_ps,decays[end-2:end],title=string("Bip_ent Gat: 2Haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
#   # m = plot(real(mut))
  
# #   display(p)
#   display(p2)
#   display(p3)
  return decays, growths, noise
end

