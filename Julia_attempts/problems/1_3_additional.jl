using ITensors
import Combinatorics: permutations

function Levi_cevita_tensor(L)
    arr = zeros((L for i in 1:L)...)
    for x in permutations(1:L)
        mat = zeros((L,L))
        for (i,j) in zip(1:L, x)
            mat[i,j] = 1
        end 
        arr[x...]=LinearAlgebra.det(mat)
    end
    return arr
end

dat2=Levi_cevita_tensor(2)
dat3 = Levi_cevita_tensor(4)
N=22

i=Index(4,"i")
is = IndexSet(n -> settags(i, "i$n"), N)
A=ITensor(dat3,is[1],is[2],is[3],is[4])
B=ITensor(dat3,is[4],is[5],is[6],is[7])
C=ITensor(dat3,is[3],is[7],is[8],is[9])
D=ITensor(dat3,is[2],is[6],is[9],is[10])
E=ITensor(dat3,is[1],is[5],is[8],is[10])

chk = A*B*C*D*E

dat3 = Levi_cevita_tensor(3)
N=22

i=Index(3,"i")
is = IndexSet(n -> settags(i, "i$n"), N)
A=ITensor(dat3,is[1],is[2],is[3])
B=ITensor(dat3,is[3],is[4],is[5])
C=ITensor(dat3,is[2],is[5],is[6])
D=ITensor(dat3,is[1],is[4],is[6])

chk = A*B*C*D