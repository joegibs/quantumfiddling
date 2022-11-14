"""
maybe working not sure about temperature 
"""

using ITensors
using PyCall
using Plots
using Statistics

N= 100
n=10
B = .01
J=0

inter = -B:2*B/n:B
en_vec= zeros(length(inter))
for (i,B) in enumerate(inter)
    sites = siteinds("S=1/2",N)

    hterms = OpSum()
    for j=1:(N-1) #add iteraction terms
        hterms -= J,"Sz",j,"Sz",j+1
    end
    hterms -= J,"Sz",N,"Sz",1  # term 'wrapping' around the ring
    for j=1:(N) #add magnetic terms
        hterms -= B,"Sz",j
        hterms -= 0.1B,"Sx",j
    end
    
    H = MPO(hterms,sites)

    sweeps = Sweeps(3) # number of sweeps is 5


    maxdim!(sweeps,10,20,100,100,200) # gradually increase states kept
    cutoff!(sweeps,1E-10) # desired truncation error
    val = zeros(3)
    for i in 1:3
        # psi0 = productMPS(sites, n -> isodd(n) ? "Dn" : "Dn")
        psi0 = randomMPS(sites, 3)
        energy,psi = dmrg(H,psi0,sweeps)
        magz = expect(psi,"Sz")
        global val[i] = 2*mean(magz)
    end
    global en_vec[i] = mean(val)
end

plot(inter,en_vec)