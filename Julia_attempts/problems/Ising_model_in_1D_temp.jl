using ITensors
using Printf

#
# Set up
#
function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ)
  h = J * op("Sz", s1) * op("Sz", s2) +
      B*op("Sz", s1) * op("I", s2)
  return exp(τ * h)
end

N=10
cutoff=1E-8
δτ=0.1 
beta_max=1.0
# Make an array of 'site' indices
s = siteinds("S=1/2", N)
B = 2.5
J=0
n=100
iter = -B:2*B/n:B
mg_vec3 = zeros(n+1)
for (i,B) in enumerate(iter)
  ######################################
  # Ground state
  ######################################
  terms = OpSum()
  for j in 1:(N-1)
    terms += "Sz", j,"Sz",j+1
  end
  for j in 1:N
    terms += B,"Sz", j
  end
  H = MPO(terms, s)
  Szterms = OpSum()
  for j in 1:N
      Szterms += "Sz",j
  end
  SzH = MPO(Szterms,s)


  psi0 = randomMPS(s)
  sweeps = Sweeps(5) # number of sweeps is 5

  energy,psi = dmrg(H,psi0,sweeps)

  ######################################
  # thermalization
  ######################################

  # Make gates (1,2),(2,3),(3,4),...
  gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,)) for n in 1:(N - 1)], s)
  # Include gates in reverse order to complete Trotter formula
  append!(gates, reverse(gates))

  #make rho from groundstate
  rho = outer(psi',psi)
  magz = inner(rho,SzH)/N

  # en_vec=zeros(Int(beta_max/δτ)+1)
  for (i,β) in enumerate(0:δτ:beta_max)
  energy = inner(rho, H)
  magz = -inner(rho,SzH)/N
  @printf("β = %.2f energy = %.8f, magz = %.5f\n", β, energy,magz)
  rho = apply(gates, rho; cutoff)
  rho = rho / tr(rho)
  end
  mg_vec3[i] = magz

#   return nothing
end
plot(iter,[mg_vec mg_vec2 mg_vec3])
