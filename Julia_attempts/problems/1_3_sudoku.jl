using ITensors

function e_mat(n , dim, bl)
    mat = zeros((dim for i in 1:n)...)
    for i in 1:dim
        mat[[i for j in 1:n]...] = 1
    end
    if bl
        return mat
    else
        return ones((dim for i in 1:n)...)-mat
    end
end

dim = 4
dims = (dim for i in 1:20)
tgs = tuple([string(i) for i in 1:length(dims)]...)
l = [Index([dims...][i],j) for (i,j) in enumerate(tgs)]

A = ITensor(e_mat(3,dim,true),  l[1],l[2],l[3])
B = ITensor(e_mat(3,dim,true),  l[9],l[4],l[5])
C = ITensor(e_mat(3,dim,true),  l[8],l[11],l[6])
D = ITensor(e_mat(3,dim,true),  l[7],l[10],l[12])
An1 = ITensor(e_mat(2,dim,false),  l[1],l[7])
An2 = ITensor(e_mat(2,dim,false),  l[2],l[8])
An3 = ITensor(e_mat(2,dim,false),  l[3],l[9])
Bn1 = ITensor(e_mat(2,dim,false),  l[4],l[10])
Bn2 = ITensor(e_mat(2,dim,false),  l[5],l[11])
Cn1 = ITensor(e_mat(2,dim,false),  l[6],l[12])

c = A*An1*An2*An3*B*Bn1*Bn2*C*Cn1*D
print(c[], '\n')