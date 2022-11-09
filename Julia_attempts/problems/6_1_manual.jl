using ITensors, ITensorGPU

Nx = 2
Ny = 2
data =[1 0  0 0; 0 0 0 1]
data = reshape(data,(2,2,2))

a1 = Index(2,"i1")
a2 = Index(2,"i2")
a3 = Index(2,"i3")
a4 = Index(2,"i4")
p1 = Index(2,"p1")
p2 = Index(2,"p2")
p3 = Index(2,"p3")
p4 = Index(2,"p4")

A = ITensor(data,a1,a2,p1)
B = ITensor(data,a2,a3,p2)
C = ITensor(data,a3,a4,p3)
D = ITensor(data,a4,a1,p4)

chk = A*B*C*D
pdata= [1 0]
pdata2=[0 1]
pA = ITensor(pdata2,p1)
pB = ITensor(pdata,p2)
pC = ITensor(pdata,p3)
pD = ITensor(pdata,p4)

chk2 = chk*pA*pB*pC*pD
print(chk2.tensor)