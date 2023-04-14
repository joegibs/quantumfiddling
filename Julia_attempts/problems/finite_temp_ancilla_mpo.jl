using ITensors
using Printf

#=

This example code implements the purification or "ancilla" method for 
finite temperature quantum systems.

For more information see the following references:
- "Finite-temperature density matrix renormalization using an enlarged Hilbert space",
  Adrian E. Feiguin and Steven R. White, Phys. Rev. B 72, 220401(R)
  and arxiv:cond-mat/0510124 (https://arxiv.org/abs/cond-mat/0510124)

=#

function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ)
  h =
    1 / 2 * op("S+", s1) * op("S-", s2) +
    1 / 2 * op("S-", s1) * op("S+", s2) +
    op("Sz", s1) * op("Sz", s2)
  return exp(τ * h)
end

function main(; N=64, cutoff=1E-8, δτ=0.1, beta_max=10.0)
#N=50, cutoff=1E-8, δτ=0.1, beta_max=2.0
  # Make an array of 'site' indices
  s = siteinds("S=1/2", N; conserve_qns=true)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,)) for n in 1:(N - 1)], s)
  # Include gates in reverse order to complete Trotter formula
  append!(gates, reverse(gates))

  # Initial state is infinite-temperature mixed state
  rho = MPO(s, "Id") ./ √2

  # Make H for measuring the energy
  terms = OpSum()
  for j in 1:(N - 1)
    terms += 1 / 2, "S+", j, "S-", j + 1
    terms += 1 / 2, "S-", j, "S+", j + 1
    terms += "Sz", j, "Sz", j + 1
  end
  H = MPO(terms, s)

  # Do the time evolution by applying the gates
  # for Nsteps steps
  eng = []
  for β in 0.5:δτ:beta_max
    energy = inner(rho, H)
    @printf("β = %.2f energy = %.8f\n", β, energy)
    push!(eng,energy)
    rho = apply(gates, rho; cutoff)
    rho = rho / tr(rho)
  end

  return eng
end

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

interval = 0.5:0.1:10
testarr=main()
ints,datas = finite_diff(interval,testarr)

plot(ints.^-1,-datas/64)