using ITensors

a = zeros(10)
Threads.@threads for i = 1:10
    a[i] = Threads.threadid()
end
print(a)

a = zeros(10)
Threads.@threads for i = 1:10
    a[i] = Threads.threadid()+1
end
print(a)

a = [[] for i in 1:1000]
@time begin
Threads.@threads for i = 1:1000
    a[i] = [rand((1,2,3,4,5)) for i in 1:10000]
end
end




function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ,B)
    h =
      1 / 2 * op("S+", s1) * op("S-", s2) +
      1 / 2 * op("S-", s1) * op("S+", s2) +
      op("Sz", s1) * op("Sz", s2)+
      B*(op("Sz", s1) * op("I", s2)+op("I", s1) * op("Sz", s2))
      return exp(τ * h)
  end

steps = 20
a = [[] for i in 1:steps]
@time begin
Threads.@threads for i = 1:steps
    N = 10
    s = siteinds("S=1/2", N)
    gates = ops([("expτSS", (n, n + 1), (τ=-0.1 / 2,B=1)) for n in 1:(N - 1)], s)
    # Include gates in reverse order to complete Trotter formula
    append!(gates, reverse(gates))
  
    a[i] = [rand((1,2,3,4,5)) for i in 1:100]
end
end

steps = 20
a = [[] for i in 1:steps]
@time begin
Threads.@threads for i = 1:steps
    N = 10
    s = siteinds("S=1/2", N)
    gates = ops([("expτSS", (n, n + 1), (τ=-0.1 / 2,B=1)) for n in 1:(N - 1)], s)
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
    terms += 1,"Sz", j
  end

  H = MPO(terms, s)
  terms2 = OpSum()
  for j in 1:(N)
    terms2 += "Sz", j
  end
  Sz = MPO(terms2, s)

    a[i] = [rand((1,2,3,4,5)) for i in 1:100]
end
end

steps = 20
a = [[] for i in 1:steps]
@time begin
Threads.@threads for i = 1:steps
    N = 10
    s = siteinds("S=1/2", N)
    gates = ops([("expτSS", (n, n + 1), (τ=-0.1 / 2,B=1)) for n in 1:(N - 1)], s)
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
    terms += 1,"Sz", j
  end

  H = MPO(terms, s)
  terms2 = OpSum()
  for j in 1:(N)
    terms2 += "Sz", j
  end
  Sz = MPO(terms2, s)
    # Make τ_range and check δτ is commensurate
    τ_range = 0.1:0.1:(0.6 / 2)
    if norm(length(τ_range) * 0.1 - 0.6 / 2) > 1E-10
      error("Time step δτ=0.1 not commensurate with 2=(0.5/2)")
    end
  
    energies = Float64[]
    magz = Float64[]
    ents = Float64[]

    a[i] = [rand((1,2,3,4,5)) for i in 1:100]
end
end

###########################################
using ITensors
using BenchmarkTools

function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ,B)
    h =
      1 / 2 * op("S+", s1) * op("S-", s2) +
      1 / 2 * op("S-", s1) * op("S+", s2) +
      op("Sz", s1) * op("Sz", s2)+
      B*(op("Sz", s1) * op("I", s2)+op("I", s1) * op("Sz", s2))
      return exp(τ * h)
  end
steps = 5
a = [[] for i in 1:steps]
@benchmark begin
Threads.@threads  for i = 1:steps
    N = 10
    s = siteinds("S=1/2", N)
    gates = ops([("expτSS", (n, n + 1), (τ=-0.1 / 2,B=1)) for n in 1:(N - 1)], s)
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
    terms += 1,"Sz", j
  end

  H = MPO(terms, s)
  terms2 = OpSum()
  for j in 1:(N)
    terms2 += "Sz", j
  end
  Sz = MPO(terms2, s)
    # Make τ_range and check δτ is commensurate
    τ_range = 0.1:0.1:(0.6 / 2)
    if norm(length(τ_range) * 0.1 - 0.6 / 2) > 1E-10
      error("Time step δτ=0.1 not commensurate with 2=(0.5/2)")
    end
  
    energies = Float64[]
    magz = Float64[]
    ents = Float64[]
    for step in 1:(10 + 100)
        # if step <= Nwarm
        #   println("Making warmup METTS number $step")
        # else
        #   println("Making METTS number $(step-Nwarm)")
        # end
    
        # Do the time evolution by applying the gates
        for τ in τ_range
            a[i] = [rand((1,2,3,4,5)) for i in 1:100]
        end
        # if step > 10
        #     energy = inner(psi', H, psi)
        #     sz = inner(psi',Sz,psi)
        #     push!(energies, energy)
        #     push!(magz,sz)

        #     print("  Energy of METTS ",Threads.threadid(), "\n")
        #   end
    end

end
end

###########################################
using ITensors
using LinearAlgebra
using BenchmarkTools
BLAS.set_num_threads(1)

function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ,B)
    h =
      1 / 2 * op("S+", s1) * op("S-", s2) +
      1 / 2 * op("S-", s1) * op("S+", s2) +
      op("Sz", s1) * op("Sz", s2)+
      B*(op("Sz", s1) * op("I", s2)+op("I", s1) * op("Sz", s2))
      return exp(τ * h)
  end
steps = 3
a = [[] for i in 1:steps]
@benchmark begin
Threads.@threads  for i = 1:steps
    N = 10
    s = siteinds("S=1/2", N)
    gates = ops([("expτSS", (n, n + 1), (τ=-0.1 / 2,B=1)) for n in 1:(N - 1)], s)
    # Include gates in reverse order to complete Trotter formula
    append!(gates, reverse(gates))
    # Make y-rotation gates to use in METTS collapses
  Ry_gates = ops([("Ry", n, (θ=π / 2,)) for n in 1:N], s)

  # Arbitrary initial state
  psi = randomMPS(s)
   test =1

  # Make H for measuring the energy
  terms = OpSum()
  for j in 1:(N - 1)
    terms += 1 / 2, "S+", j, "S-", j + 1
    terms += 1 / 2, "S-", j, "S+", j + 1
    terms += "Sz", j,"Sz", j+1
  end
  for j in 1:(N)
    terms += 1,"Sz", j
  end

  H = MPO(terms, s)
  terms2 = OpSum()
  for j in 1:(N)
    terms2 += "Sz", j
  end
  Sz = MPO(terms2, s)
    # Make τ_range and check δτ is commensurate
    τ_range = 0.1:0.1:0.3
  
    energies = Float64[]
    magz = Float64[]
    ents = Float64[]
    for step in 1:(10 + 100)
        # if step <= Nwarm
        #   println("Making warmup METTS number $step")
        # else
        #   println("Making METTS number $(step-Nwarm)")
        # end
    
        # Do the time evolution by applying the gates
        for τ in τ_range
            test=1+Threads.thre/adid()
            # psi = apply(gates, psi,cutoff=1e-8)
            # normalize!(psi)
        end
        if step > 10
            energy = inner(psi', H, psi)
        #     sz = inner(psi',Sz,psi)
        #     push!(energies, energy)
        #     push!(magz,sz)

        #     print("  Energy of METTS ",Threads.threadid(), "\n")
          end
    end

end
end