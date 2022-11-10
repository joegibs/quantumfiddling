J=0.1
B=0
N=10
sites = siteinds("S=1/2",N)

hterms = OpSum()
for j=1:(N-1) #add iteraction terms
    global hterms -= J,"Sz",j,"Sz",j+1
end
hterms -= J,"Sz",N,"Sz",1  # term 'wrapping' around the ring
for j=1:(N) #add magnetic terms
    global hterms -= B,"Sz",j
end

H = MPO(hterms,sites)

sweeps = Sweeps(5) # number of sweeps is 5


maxdim!(sweeps,10,20,100,100,200) # gradually increase states kept
cutoff!(sweeps,1E-10) # desired truncation error

psi0 = randomMPS(ComplexF64, sites)
magz = expect(psi0,"Sz")
print(mean(magz))
energy,psi = dmrg(H,psi0,sweeps)
magz2 = expect(psi,"Sz")
print(mean(magz2),'\n')