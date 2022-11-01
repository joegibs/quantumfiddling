using ITensors

function e_mat(n , bl)
    mat = zeros((3 for i in 1:n)...)
    for i in 1:3
        mat[[i for j in 1:n]...] = 1
    end
    if bl
        return mat
    else
        return ones((3 for i in 1:n)...)-mat
    end
end

dims = (3 for i in 1:20)
tgs = tuple([string(i) for i in 1:length(dims)]...)
l = [Index([dims...][i],j) for (i,j) in enumerate(tgs)]

A = ITensor(e_mat(4,true),  l[1],l[2],l[3],l[4])
B = ITensor(e_mat(4,true),  l[1],l[5],l[9],l[8])
C = ITensor(e_mat(3,true),  l[5],l[6],l[7])
D = ITensor(e_mat(4,true),  l[4],l[7],l[12],l[13])
E = ITensor(e_mat(3,true),  l[2],l[10],l[14])
F = ITensor(e_mat(3,true),  l[3],l[11],l[17])
G = ITensor(e_mat(5,true),  l[6],l[10],l[11],l[15],l[16])
H = ITensor(e_mat(3,true),  l[9],l[15],l[19])
I = ITensor(e_mat(3,true),  l[12],l[16],l[18])
J = ITensor(e_mat(4,true),  l[8],l[14],l[18],l[20])
K = ITensor(e_mat(4,true),  l[13],l[17],l[19],l[20])

cont = *(A,B,C,D,E,F,G,H,I,J,K)
print(cont.tensor)