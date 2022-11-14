using ITensors
"""
def broke
"""
N = 100
cutoff = 1E-8
tau = 0.1
ttotal = .5
J=2
B=1

s = siteinds("S=1/2", N; conserve_qns=true)

gates = ITensor[]
for j in 1:(N - 1)
    s1 = s[j]
    s2 = s[j + 1]
    hj = J*op("Sz", s1) * op("Sz", s2)
    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
end
append!(gates, reverse(gates))

psi = productMPS(s, n -> n!=1||n!=50 ? "Up" : "Dn")

c = div(N, 2) # center site
Sz = zeros((length(0.0:tau:ttotal),N))
# Compute and print <Sz> at each time step
# then apply the gates to go to the next time
for (ii,t) in enumerate(0.0:tau:ttotal)
    global Sz[1+100*(ii-1):100*ii] = expect(psi, "Sz")
    #println("$t $Sz")

    t≈ttotal && break

    global psi = apply(gates, psi; cutoff)
    normalize!(psi)
end

heatmap(Sz)