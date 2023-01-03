using ITensors
using Printf
using Plots
using Statistics
using FiniteDifferences
#=

This example code implements the minimally entangled typical thermal state (METTS).
For more information on METTS, see the following references:
- "Minimally entangled typical quantum states at finite temperature", Steven R. White,
  Phys. Rev. Lett. 102, 190601 (2009)
  and arxiv:0902.4475 (https://arxiv.org/abs/0902.4475)
- "Minimally entangled typical thermal state algorithms", E M Stoudenmire and Steven R White,
  New Journal of Physics 12, 055026 (2010) https://doi.org/10.1088/1367-2630/12/5/055026

=#
function finite_diff(interval,data)
    N= length(data)
    data_diff = Float64[]
    intervals = Float64[]
    for i in 1:(N-1)
        push!(data_diff,(data[i]-data[i+1])/(interval[i]-interval[i+1]))
        push!(intervals,(interval[i]+interval[i+1])/2)
    end
    return intervals,data_diff
end
function central_finite_diff(interval,data)
    N= length(data)
    data_diff = Float64[]
    intervals = Float64[]
    for i in 2:(N-1)
        push!(data_diff,(data[i-1]-data[i+1])/(2*(interval[i-1]-interval[i+1])))
        push!(intervals,(interval[i]))
    end
    return intervals,data_diff
end
function entrp(psi,b)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
    end
    return SvN
end
function ITensors.op(::OpName"expτSS", ::SiteType"S=1", s1::Index, s2::Index; τ,B)
  h =
    op("Sx", s1) * op("Sx", s2) +
    op("Sy", s1) * op("Sy", s2) +
    op("Sz", s1) * op("Sz", s2)+
    B*(op("Sz", s1) * op("I", s2)+op("I", s1) * op("Sz", s2))
    return exp(τ * h)
end

function main(; N=36, cutoff=1E-8, δτ=0.1, beta=2.0, NMETTS=3000, Nwarm=10,b=1)
  # Make an array of 'site' indices
  s = siteinds("S=1", N)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,B=b)) for n in 1:(N - 1)], s)
  # Include gates in reverse order to complete Trotter formula
  append!(gates, reverse(gates))

  # Make y-rotation gates to use in METTS collapses
#   Ry_gates = ops([("Sy",n) for n in 1:N], s)

  # Arbitrary initial state
  psi = randomMPS(s)

  # Make H for measuring the energy
  terms = OpSum()
  for j in 1:(N - 1)
    terms += "Sx", j, "Sx", j + 1
    terms += "Sy", j, "Sy", j + 1
    terms += "Sz", j,"Sz", j+1
  end
  for j in 1:(N)
    terms += b,"Sz", j
  end

  H = MPO(terms, s)
  terms2 = OpSum()
  for j in 1:(N)
    terms2 += "Sz", j
  end
  Sz = MPO(terms2, s)


  # Make τ_range and check δτ is commensurate
  τ_range = δτ:δτ:(beta / 2)
  if norm(length(τ_range) * δτ - beta / 2) > 1E-10
    error("Time step δτ=$δτ not commensurate with beta/2=$(beta/2)")
  end

  energies = Float64[]
  magz = Float64[]
  ents = Float64[]
  for step in 1:(Nwarm + NMETTS)
    if step <= Nwarm
      println("Making warmup METTS number $step")
    else
      println("Making METTS number $(step-Nwarm)")
    end

    # Do the time evolution by applying the gates
    for τ in τ_range
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end

    # Measure properties after >= Nwarm 
    # METTS have been made
    if step > Nwarm
      energy = inner(psi', H, psi)
      sz = inner(psi',Sz,psi)
      ent = entrp(psi,Int(N/2))
      push!(energies, energy)
      push!(magz,sz)
      push!(ents,ent)
      @printf("  Energy of METTS %d = %.4f\n", step - Nwarm, energy)
    end

    # Measure in X or Z basis on alternating steps
    if step % 2 == 1
    #   psi = apply(Ry_gates, psi)
      samp = sample!(psi)
      new_state = []
      for i in samp
        if i == 1
            push!(new_state,"X+")
        elseif i == 2
            push!(new_state,"X0")
        else push!(new_state,"X-")
        end
      end
    else
      samp = sample!(psi)
      new_state = []
      for i in samp
        if i == 1
            push!(new_state,"Z+")
        elseif i == 2
            push!(new_state,"Z0")
        else
            push!(new_state,"Z-")
        end
      end
    end
    psi = productMPS(s, new_state)
  end
  return energies,magz,ents
end

data1 = []
data_ent1=[]
interval = -2:0.5:2
N=18
for i in interval
  eng,mag,ent = main(N=N,b=i,NMETTS=100,δτ=.2, beta=2.0,)
  push!(data1,mean(mag)/N)
  push!(data_ent1,mean(ent))
  p = histogram(eng)
  # display(p)
end
plot(interval,[data1,data_ent1])

data_eng = []
data_var = []
interval =.1:0.1:3.1#10 .^ range(0, stop=1.5, length=40)
N=100
for i in interval
  beta = i
  δτ=i/10
  eng,mag,ent = main(N=N,b=0,NMETTS=1000,δτ=δτ, beta=beta)
#   push!(data,mean(mag)/N)
  push!(data_eng,mean(eng))
  push!(data_var,var(eng))
#   p = histogram(eng)
#   display(p)
end
plot(interval,data_eng)
plot(interval,data_var.*interval.*interval./N)
ints,datas=finite_diff(interval,data_eng)
plot(ints,datas)
erv= [interval...][1:7]
plot(erv,data_var.*erv.*erv./N)
#need to do better point diferrentiation method.....
#structure factor

#1d analouge of cuperates, or fustration, scalin of specific heat
#https://www.nature.com/articles/s41467-018-06800-2#MOESM1
#mag field for thermo quantities

#todo clean up code, check magnetic things , ham from the paper
#talk to kyle 0 temp xxz dft to static structure factor\
#schollwock book


# bisected system entropy scaling
#site-site scaling