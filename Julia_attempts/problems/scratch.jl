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

function rec_ent(rho,b)
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
  #turn this mpo into a single tensor
  T = prod(M)

  _,S,_ = svd(T,s)
  SvN = 0.0
  for n in 1:dim(S, 1)
    p = S[n,n]^2
    if p != 0
      SvN -= p * log(p)
    end
  end
  return SvN
end
function ren(rho)
  return log(tr(apply(rho,rho)))
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

  ren = log(tr(apply(M,M)))
  return ren
end
function main(; N=10, cutoff=1E-8, δτ=0.1, beta_max=2.0)

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
  for j in 2:(N - 2)
    terms += 1 / 2, "S+", j, "S-", j + 1
    terms += 1 / 2, "S-", j, "S+", j + 1
    terms += "Sz", j, "Sz", j + 1
  end
  H = MPO(terms, s)

  # Do the time evolution by applying the gates
  # for Nsteps steps
  for β in 0:δτ:beta_max
    energy = inner(rho, H)
    entr = rec_ent(rho,5)#3,8,5)
    @printf("β = %.2f energy = %.8f ent = %.4f\n", β, energy,entr)
    rho = apply(gates, rho; cutoff)
    rho = rho / tr(rho)
  end

  return nothing
end



N=10
s = siteinds("S=1/2", N)

# # Make gates (1,2),(2,3),(3,4),...
# gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,)) for n in 1:(N - 1)], s)
# # Include gates in reverse order to complete Trotter formula
# append!(gates, reverse(gates))

# Initial state is infinite-temperature mixed state
rho = MPO(s, "Id") ./ √2
#contract half   x x x x | | | |
L = ITensor(1.0)
for i = 1:5
  L *= tr(rho[i])
  @show inds(L)
end
# absorb
rho[5+1] *= L
# no longer a proper mpo
M =MPO(10-5)
for i in 1:(10-5)
    M[i]=rho[5+i]
end

re = log(tr(apply(M,M)))
n = length(rho)

rho=rho/tr(rho)
rho_temp = deepcopy(rho)
s = siteinds("Qubit",n) 

#contract half   x x x x | | | |
L = ITensor(1.0)
for i = 1:5
  L *= tr(rho_temp[i])
end
# absorb
rho_temp[5+1] *= L
# no longer a proper mpo
M =MPO(n-5)
for i in 1:(n-5)
    M[i]=rho_temp[5+i]
end
#turn this mpo into a single tensor
T = prod(M)

_,S,_ = svd(T,s)
SvN = 0.0
for n in 1:dim(S, 1)
  p = S[n,n]^2
  if p != 0
    SvN -= p * log(p)
  end
end
SvN