using ITensors
using Random
using LinearAlgebra
using Statistics

include("func_utils.jl")

#=
Some ITensor ops
=#
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
ITensors.op(::OpName"K0",::SiteType"Qubit"; p::Number=0) =
  [1 0
   0 √(1-p)]

function kraus_dephase(rho,s,p)
    #define the two operators
    #(1-p)ρ + pZρZ
    N=length(rho)
    
    for i in 1:N
        gates = ITensor[]
        hj = op("Z", s[i])
        push!(gates, hj)
        rho = (1-p)*rho + p*apply(gates,rho;apply_dag=true)
        rho = rho/tr(rho)
    end
    truncate!(rho)
    return rho
end

function kraus_amp_damp(rho,s,p)
    N=length(rho)


    for i in 1:N
      gates1 = ITensor[]
      gates2 = ITensor[]
      hj = op("S+", s[i])
      push!(gates1, hj)

      hj = op("K0", s[i];p)
      push!(gates2, hj)
      
      #apply the operators
      rho = apply(gates2,rho;apply_dag=true) +p*apply(gates1,rho;apply_dag=true)
    end
    return rho
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


function gen_samp_row(N,meas_p,eoo)
    #=
    generate array of bools to sample or not
    =#
    lst= [rand()<meas_p ? 1 : 0 for i in 1:N]
    if !eoo
      #ensure only one chance of measuring after unitary
      lst[1]=0
      lst[end]=0
    end
    return lst
end

function samp_mps(rho,s,samp_row)
  #=
  sample an itensor mpo and return the next resulting mps
  =#
  N = length(rho)
  samp =deepcopy(rho)
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

function gen_step(N,rho,s,step_num,meas_p,noise,noise_type,meas_during)
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
      hj = op("Rand",[s1,s2])
      Gj=hj
      push!(gates, Gj)
  end
  cutoff = 1E-8
  rho = apply(gates, rho;apply_dag=true,cutoff=1E-15)
  #calculate obs
  if meas_during
    measured_vals = (rec_ent(rho,Int(round(N/2)),s),log_negativity(rho,Int(round(N/2)),s))
  end

  #noise channel
  if noise_type=="damp"
    rho = kraus_amp_damp(rho,s,noise)
  elseif noise_type=="deph"
    rho = kraus_dephase(rho,s,noise)
  end
  #sample as needed
    

  samp_row=gen_samp_row(N,meas_p,Bool(step_num%2))
  rho = samp_mps(rho,s,samp_row)
  rho=rho/tr(rho)

  return rho,measured_vals
end

function do_exp(N,steps,meas_p,noise,noise_type)
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    rho=outer(psi',psi)
  
    svn =[]
    neg= []
    for i in 1:steps
        rho ,(meas_svn,meas_neg)= gen_step(N,rho,s,i,meas_p,noise,noise_type)
        append!(svn,meas_svn)
        append!(neg,meas_neg)
    end
    return svn,neg
  end
function do_exp(N,steps,meas_p,noise,noise_type,meas_during)
    s = siteinds("Qubit", N) #+1 for ancilla
    psi = productMPS(s, "Up" )
    rho=outer(psi',psi)
  
    svn =[]
    neg= []
    for i in 1:steps
        rho ,(meas_svn,meas_neg)= gen_step(N,rho,s,i,meas_p,noise,noise_type,meas_during)
        if meas_during
            append!(svn,meas_svn)
            append!(neg,meas_neg)
        end
      #   @show(tr(rho))
    end
    if !meas_during
        (meas_svn,meas_neg)=(rec_ent(rho,Int(round(N/2)),s),log_negativity(rho,Int(round(N/2)),s))
        append!(svn,meas_svn)
        append!(neg,meas_neg)
    end
    return svn,neg
  end

function do_trials(N,steps,meas_p,trials,noise,noise_type)
    svn_trials,tri_trials = do_exp(N,steps,meas_p,noise,noise_type,true)
    for i in 2:trials
        print(i)
        nSvn,ntm = do_exp(N,steps,meas_p,noise,noise_type,true)
        svn_trials = 2*mean([(i-1)/i*svn_trials,1/i*nSvn])
        tri_trials = 2*mean([(i-1)/i*tri_trials,1/i*ntm])
    end
    return svn_trials,tri_trials
end

function do_trials(N,steps,meas_p,trials,noise,noise_type,meas_during)
    svn_trials,neg_trials = do_exp(N,steps,meas_p,noise,noise_type,meas_during)
    svn_egg = (0,zeros(length(svn_trials)),zeros(length(svn_trials)))
    neg_egg = (0,zeros(length(neg_trials)),zeros(length(neg_trials)))
    for i in 2:trials
        print(i)
        nSvn,ntm = do_exp(N,steps,meas_p,noise,noise_type,meas_during)
        svn_egg = welford_mean_var(svn_egg,nSvn)
        neg_egg = welford_mean_var(neg_egg,ntm)
    end
    svn_mean, svn_var = welford_extract(svn_egg)
    neg_mean, neg_var = welford_extract(svn_egg)
    return (svn_mean, svn_var, neg_mean, neg_var) 
end