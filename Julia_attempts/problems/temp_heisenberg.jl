using ITensors
using Printf
using Plots
using Statistics
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
function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ,B)
  h =
    1 / 2 * op("S+", s1) * op("S-", s2) +
    1 / 2 * op("S-", s1) * op("S+", s2) +
    op("Sz", s1) * op("Sz", s2)+
    B*(op("Sz", s1) * op("I", s2)+op("I", s1) * op("Sz", s2))
    return exp(τ * h)
end

function main(; N=36, cutoff=1E-8, δτ=0.1, beta=2.0, NMETTS=3000, Nwarm=10,b=1)
  # Make an array of 'site' indices
  s = siteinds("S=1/2", N)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,B=b)) for n in 1:(N - 1)], s)
  # Include gates in reverse order to complete Trotter formula
  append!(gates, reverse(gates))

  # Make y-rotation gates to use in METTS collapses
  Ry_gates = ops([("Ry", n, (θ=π / 2,)) for n in 1:N], s)

  # Arbitrary initial state
  psi = randomMPS(s)

  # Make H for measuring the energy
  terms = OpSum()
  for j in 1:(N - 1)
    terms += 1 / 2, "S+", j, "S-", j + 1
    terms += 1 / 2, "S-", j, "S+", j + 1
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
      psi = apply(Ry_gates, psi)
      samp = sample!(psi)
      new_state = [samp[j] == 1 ? "X+" : "X-" for j in 1:N]
    else
      samp = sample!(psi)
      new_state = [samp[j] == 1 ? "Z+" : "Z-" for j in 1:N]
    end
    psi = productMPS(s, new_state)
  end
  return energies,magz,ents
end

data = []
data_ent=[]
interval = -2:0.2:2
N=18
for i in interval
  eng,mag,ent = main(N=N,b=i,NMETTS=100,δτ=.2, beta=2.0,)
  push!(data,mean(mag)/N)
  push!(data_ent,mean(ent))
  p = histogram(eng)
  # display(p)
end
plot(interval,[data,data_ent])
data1=data
plot(interval,[data,data1])

# Need Cv = dE/dT loop over beta, get energy thn do some discrete drivatives... fun
# data = []
data_eng=[]
interval =1:1:31#10 .^ range(0, stop=1.5, length=40)
N=36
for i in interval
    beta = i
    δτ=beta/10
  eng,mag,ent = main(N=N,b=0,NMETTS=50,δτ=δτ, beta=beta)
#   push!(data,mean(mag)/N)
  push!(data_eng,mean(eng))
#   p = histogram(eng)
#   display(p)
end
plot(interval,data_eng)
ints,datas=finite_diff(interval,data_eng)
plot(ints,datas.*ints.*ints)
