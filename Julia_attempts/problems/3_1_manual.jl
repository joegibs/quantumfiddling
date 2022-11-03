using ITensors

X = [0 1 ; 1 0]
I = [1 0; 0 1]
A = reshape([[1 0 ; 0 1]  [0 1 ; 0 0]],(2,2,2))

#initialize sites and index arrays
N = 10
s_vec="00000010000"
p_ind =Dict('0'=>[1 0], '1' => [0 1])
s_string=[p_ind[x] for x in s_vec]
inds_i = [Index(2, "i$i") for i in 1:N+1]
inds_p = [Index(2, "p$i") for i in 1:N]

#Combine into tensors 
tens = [ITensor(A, inds_i[i],inds_i[i+1],inds_p[i]) for i in 1:N]
push!(tens,ITensor(X,inds_i[N+1],inds_i[1]))

#include physical tensors
sites = [ITensor(s_string[i], inds_p[i]) for i in 1:N]

#contract
cont = accumulate(*,vcat(tens,sites))[end]
print(cont.tensor)