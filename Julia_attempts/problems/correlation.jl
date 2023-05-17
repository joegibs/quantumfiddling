using ITensors
using Printf
N=4
s = siteinds("Qubit", N)
psi1 = productMPS(s, "Up" )
rho = outer(psi1',psi1)
gates = ITensor[]

s1 = s[1]
s2 = s[2]
s3 = s[3]
s4 = s[4]
hj = op("H",s1)
push!(gates, hj)
hj = op("CNOT",[s1,s2])
push!(gates, hj)
hj = op("H",s3)
push!(gates, hj)
hj = op("CNOT",[s3,s4])
push!(gates, hj)
# rho = apply(gates, rho; apply_dag=true, cutoff = 1E-8)
# gates = ITensor[]
hj = op("H",s2)
push!(gates, hj)
hj = op("CNOT",[s2,s3])
push!(gates, hj)
# hj = op("H",s3)
# push!(gates, hj)
# hj = op("CNOT",[s3,s4])
# push!(gates, hj)


hj = op("Rand",[s1,s2])
push!(gates,hj)
hj = op("Rand",[s3,s4])
push!(gates,hj)
hj = op("Rand",[s2,s3])
push!(gates,hj)
hj = op("Rand",[s1,s2])
push!(gates,hj)
hj = op("Rand",[s3,s4])
push!(gates,hj)
hj = op("Rand",[s2,s3])
push!(gates,hj)


cutoff = 1E-8

rho = apply(gates, rho; apply_dag=true, cutoff = 1E-8)
psi = apply(gates,psi1)
Czz=correlation_matrix(psi,"Sz","Sz")