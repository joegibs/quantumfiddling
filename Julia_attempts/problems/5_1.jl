using ITensors
using PyCall
n=5
en_vec= zeros(n)
for N in 50:100:50+100*(n-1)
    sites = siteinds("S=1/2",N)

    ampo = OpSum()
    for j=1:N-1
        ampo -= "Sx",j,"Sx",j+1
    end
    for j=1:N
        ampo -= "Sz",j
    end
    H = MPO(ampo,sites)

    sweeps = Sweeps(5) # number of sweeps is 5
    maxdim!(sweeps,10,20,100,100,200) # gradually increase states kept
    cutoff!(sweeps,1E-10) # desired truncation error

    psi0 = randomMPS(sites,2)
    energy,psi = dmrg(H,psi0,sweeps)
    global en_vec[Int((N-50)/100)+1] = energy
end


np = pyimport("numpy")
scipy  = pyimport("scipy")

py"""
def func(x, a, b):
    return 1- 1/np.sin(np.pi/(a* x+b))
"""  
               
xdata = np.linspace(50,450,n)
ydata = en_vec

popt, pcov = scipy.optimize.curve_fit(py"func", xdata, ydata)
print("α : ",popt[1], "\nβ : ",popt[2])
