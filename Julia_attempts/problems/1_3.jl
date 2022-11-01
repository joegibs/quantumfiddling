using ITensors

function e_mat(n , bl)
    mat = zeros(n,n,n)
    for i in 1:n
        mat[i,i,i] = 1
    end
    if bl
        return mat
    else
        return ones(n,n,n)-mat
    end
end

dims = (3 for i in 1:14)
tgs = tuple([string(i) for i in 1:length(dims)]...)
l = [Index([dims...][i],j) for (i,j) in enumerate(tgs)]

Ae = ITensor(LinearAlgebra.I(3),  l[1],l[2])
An1 = ITensor(ones(3,3) -LinearAlgebra.I(3),l[1],l[3])
An2 = ITensor(ones(3,3) -LinearAlgebra.I(3),l[2],l[4])

Be = ITensor(LinearAlgebra.I(3),  l[3],l[5])
Bn = ITensor(ones(3,3) -LinearAlgebra.I(3),l[5],l[7])

Ce = ITensor(LinearAlgebra.I(3),  l[4],l[6])
Cn = ITensor(ones(3,3) -LinearAlgebra.I(3),l[6],l[8])

De = ITensor(e_mat(3,true), l[7],l[10],l[11])
Dn1 = ITensor(ones(3,3) -LinearAlgebra.I(3),l[10],l[9])
Dn2 = ITensor(ones(3,3) -LinearAlgebra.I(3),l[11],l[13])

Ee = ITensor(e_mat(3,true), l[8],l[9],l[12])
En1 = ITensor(ones(3,3) -LinearAlgebra.I(3),l[12],l[14])

Fe = ITensor(LinearAlgebra.I(3),  l[13],l[14])

cont = *(Ae,An1,An2,Be,Bn,Ce,Cn,De,Dn1,Dn2,Ee,En1,Fe)
print(cont.tensor)