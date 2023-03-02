using ITensors
using Printf
using Plots


function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ)
    h =
      1 / 2 * op("S+", s1) * op("S-", s2) +
      1 / 2 * op("S-", s1) * op("S+", s2) +
      op("Sz", s1) * op("Sz", s2)
    return exp(τ * h)
  end
  
  """
  Given a Vector of numbers, returns
  the average and the standard error
  (= the width of distribution of the numbers)
  """
  function avg_err(v::Vector)
    N = length(v)
    avg = v[1] / N
    avg2 = v[1]^2 / N
    for j in 2:N
      avg += v[j] / N
      avg2 += v[j]^2 / N
    end
    return avg, √((avg2 - avg^2) / N)
  end
  
#   function main(; N=10, cutoff=1E-8, δτ=0.1, beta=2.0, NMETTS=3000, Nwarm=10)

function entangle_cluster_sites(psi,sites,s1::Int,s2::Int)
"""
bleh this is not right

psi: state
s1: start of cluster
s2: end of cluster
"""
    H_gates = ops([isodd(n) ? ("H",n) : ("I",n) for n in s1:s2],sites)
    CNOT_gates = ops([("CNOT",(n,n+1)) for n in s1:2:s2],sites)
    psi = apply(H_gates,psi)
    psi = apply(CNOT_gates,psi)
    return psi
end

function imag_time_gates(sites,s1::Int,s2::Int,N::Int)
    s = [1:s1-1;s1:2:s2-1;s2+1:N-1]
    gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,)) for n in s], sites)
    return gates
end 

N=50
s1=2
s2=48
cutoff=1E-8
δτ=0.1
beta=2.0
NMETTS=300
Nwarm=10
esites = siteinds("S=1/2",2)
# Make an array of 'site' indices
s = siteinds("S=1/2", N)

# Make gates (1,2),(2,3),(3,4),...
gates = imag_time_gates(s,s1,s2,N)

# Include gates in reverse order to complete Trotter formula
append!(gates, reverse(gates))

# Arbitrary initial state
psi = randomMPS(s)
psi = entangle_cluster_sites(psi,s,s1,s2)
# Make H for measuring the energy
terms = OpSum()
for j in 1:(N-1)#[1:s1-1;s1:2:s2-1;s2+1:N-1]
  terms += 1 / 2, "S+", j, "S-", j + 1
  terms += 1 / 2, "S-", j, "S+", j + 1
  terms += "Sz", j, "Sz", j + 1
end
H = MPO(terms, s)

# Make τ_range and check δτ is commensurate
τ_range = δτ:δτ:(beta / 2)
if norm(length(τ_range) * δτ - beta / 2) > 1E-10
  error("Time step δτ=$δτ not commensurate with beta/2=$(beta/2)")
end

energies = Float64[]
eng = Float64[]

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
    push!(energies, energy)
    @printf("  Energy of METTS %d = %.4f\n", step - Nwarm, energy)
    a_E, err_E = avg_err(energies)
    push!(eng,a_E)
    @printf(
      "  Estimated Energy = %.4f +- %.4f  [%.4f,%.4f]\n",
      a_E,
      err_E,
      a_E - err_E,
      a_E + err_E
    )
  end

  # Measure in X or Z basis on alternating steps
  if step % 2 == 1
    samp = sample!(psi)
    new_state = [samp[j] == 1 ? "X+" : "X-" for j in 1:N]
  else
    samp = sample!(psi)
    new_state = [samp[j] == 1 ? "Z+" : "Z-" for j in 1:N]
  end
  psi = productMPS(s, new_state)
  psi = entangle_cluster_sites(psi,s,s1,s2)

end

# return nothing
#   end
