using ITensors

X = [0 1 ; 1 0]
I = [1 0; 0 1]
A = reshape([[1 0 ; 0 1]  [0 1 ; 0 0]],(2,2,2))

N = 10
s_vec=zeros(10)
s_vec[5]=1
# s_vec[6] = 1
p_ind =(a=[1 0], b = [0 1])
s_string=[p_ind[Symbol(Char(x+97))] for x in s_vec]
inds_i = [Index(2, "i$i") for i in 1:N+1]
inds_p = [Index(2, "p$i") for i in 1:N]
tens = [ITensor(A, inds_i[i],inds_i[i+1],inds_p[i]) for i in 1:N]
push!(tens,ITensor(X,inds_i[N+1],inds_i[1]))
# push!(tens,ITensor(X,inds_i[N+1],inds_i[1]))
sites = [ITensor(s_string[i], inds_p[i]) for i in 1:N]

cont = accumulate(*,vcat(tens,sites))[end]
# fin = cont * delta((inds_p[i] for i in 1:N)...)
print(cont.tensor)