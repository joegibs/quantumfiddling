using ITensors

#this kinda works 
#need to tweek operators a bit to get the exact result with bond dim4

J=1
h=0.5
N=5

sites = siteinds("S=1",N)

hterms = OpSum()
for j=2:(N-1) #add iteraction terms
    hterms -= J,"Sz",j-1,"Sx",j,"Sz",j+1
end
for j=1:(N) #add magnetic terms
    hterms -= h,"Sx",j
end

H = MPO(hterms,sites)

Hitensor = ITensor(1.)
for i = 1:N
    Hitensor *= H[i]
end
