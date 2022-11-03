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
dims = (dim for i in 1:40)
tgs = tuple([string(i) for i in 1:length(dims)]...)
l = [Index([dims...][i],j) for (i,j) in enumerate(tgs)]

A = ITensor(e_mat(4,dim,true),  l[1],l[2],l[3],l[4])
An1 = ITensor(e_mat(2,dim,false),l[1],l[21])
An2 = ITensor(e_mat(2,dim,false),l[2],l[22])
An3 = ITensor(e_mat(2,dim,false),l[3],l[23])
An4 = ITensor(e_mat(2,dim,false),l[4],l[24])
B = ITensor(e_mat(4,dim,true),  l[21],l[5],l[9],l[8])
Bn1 = ITensor(e_mat(2,dim,false),l[5],l[25])
Bn2 = ITensor(e_mat(2,dim,false),l[9],l[29])
Bn3 = ITensor(e_mat(2,dim,false),l[8],l[28])
C = ITensor(e_mat(3,dim,true),  l[25],l[6],l[7])
Cn1 = ITensor(e_mat(2,dim,false),l[6],l[26])
Cn2 = ITensor(e_mat(2,dim,false),l[7],l[27])
D = ITensor(e_mat(4,dim,true),  l[24],l[27],l[12],l[13])
Dn1 = ITensor(e_mat(2,dim,false),l[12],l[32])
Dn2 = ITensor(e_mat(2,dim,false),l[13],l[33])
E = ITensor(e_mat(3,dim,true),  l[22],l[10],l[14])
En1 = ITensor(e_mat(2,dim,false),l[10],l[30])
En2 = ITensor(e_mat(2,dim,false),l[14],l[34])
F = ITensor(e_mat(3,dim,true),  l[23],l[11],l[17])
Fn1 = ITensor(e_mat(2,dim,false),l[11],l[31])
Fn2 = ITensor(e_mat(2,dim,false),l[17],l[37])
G = ITensor(e_mat(5,dim,true),  l[26],l[30],l[31],l[15],l[16])
Gn1 = ITensor(e_mat(2,dim,false),l[15],l[35])
Gn2 = ITensor(e_mat(2,dim,false),l[16],l[36])
H = ITensor(e_mat(3,dim,true),  l[29],l[35],l[19])
Hn1 = ITensor(e_mat(2,dim,false),l[19],l[39])
I = ITensor(e_mat(3,dim,true),  l[32],l[36],l[18])
In1 = ITensor(e_mat(2,dim,false),l[18],l[38])
J = ITensor(e_mat(4,dim,true),  l[28],l[34],l[38],l[20])
Jn1 = ITensor(e_mat(2,dim,false),l[20],l[40])
K = ITensor(e_mat(4,dim,true),  l[33],l[37],l[39],l[40])

cont = *(A,An1,An2,An3,An4,B,Bn1,Bn2,Bn3,C,Cn1,Cn2,D,Dn1,Dn2,E,En1,En2,F,Fn1,Fn2,G,Gn1,Gn2,H,I,J,K,Hn1,In1,Jn1)
print(cont[])