using ITensors

#init indicies
α = Index(3,"α")
β = Index(3,"β")
γ = Index(3,"γ")
δ = Index(3,"δ")
ϵ = Index(3,"ϵ")
#init tensors
A=randomITensor(α,β)
B=randomITensor(α,γ)
C=randomITensor(γ,δ,ϵ)
D=randomITensor(β,δ)
E=randomITensor(ϵ)

#step by step contraction and show indicies at each step
print([tags(x) for x in inds(A)],"\n")
first_con = A*B;
print([tags(x) for x in inds(first_con)],"\n")
second_con = first_con*C;
print([tags(x) for x in inds(A)],"\n")
third_con = second_con*D;
print([tags(x) for x in inds(third_con)],"\n")
fourth_con = third_con*E;
print([tags(x) for x in inds(fourth_con)],"\n")

