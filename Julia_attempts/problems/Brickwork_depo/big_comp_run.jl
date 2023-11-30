include("pt_mpo.jl")
include("func_utils.jl")
ITensors.set_warn_order(20)

using DelimitedFiles 
using Dates
using Plots
using BenchmarkTools
using LaTeXStrings


cd(@__DIR__)
pwd()

#create data
let
    sits = [4 6];
    interval = 0.0:0.1:1;
    num_samp=10;
    noise_val=0.0;
    fld = Dates.format(now(), "YYYY_mm_dd")
    mkdir(fld)
    tme = Dates.format(now(), "HH_MM")
    flname="pt_mpo_$time.csv"
    touch(fld*"/"*flname)

    @show num_samp
    data= @time main(sits,interval,num_samp,noise_val)
    sits,interval,num_samp,noise_val,decays,growths = data
    #Writing contents to the file 
    
    open(fld*"/"*flname, "w") do io 
        writedlm(io,data,',')
    end; 

    #readlm() method to read the DelimitedFiles 
    # tst = readdlm(flname,',')

    # plot(tst)
    
end 



# plotting
# let 
#     fln = "/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/2023_11_10_11_32/pt_mpo.csv"
#     sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fln)
#     grw2 = real(growths2)
#     p = plot(interval2,grw2,title=string("Bip_ent Gat: 2haar, varying meas_p"), label=string.(transpose(sits2)), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
#     display(p)
# end

let 
    fldr = "/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/2023_11_20_40_damp/"
    growths=[]
    decays =[]
    sits2=0;interval2=0;num_samp2=0;noise_val2=0;
    for fle in readdir(fldr)
        i=1
        if i==1
            sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fldr*'/'*fle)
            growths = growths2
            decays = decays2
        else
            sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fldr*'/'*fle)
            append!(growths,2*mean([(i-1)/i*growths,1/i*growths2]))
            append!(decays, 2*mean([(i-1)/i*decays,1/i*decays2]))
        end
        i=i+1
    end
    # popfirst!(interval2)
    popfirst!(sits2)
    popfirst!(decays)
    popfirst!(growths)
    # fln = "/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/2023_11_10_11_32/pt_mpo.csv"
    # sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fln)
    # @show real(growths)
    # @show size(growths)
    
    # grw2 = real(growths2)
    p = plot(interval2,real(decays),title=string("Bip_ent Gate: 2haar, varying meas_p"), label=string.(transpose(sits2)), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
    display(p)
    pg = plot(interval2,real(growths),title=string("negativity Gate: 2haar, varying meas_p"), label=string.(transpose(sits2)), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
    display(pg)

    py"""
    import numpy as np
    import scipy
    L=$sits2
    interval = $interval2#[x/10 for x in range(9)]
    tot_vonq = np.real($decays)
    def xfunc(p,l,pc,v):
        return (p-pc)*l**(1/v)

    def Spc(pc,l):
        spot, = np.where(np.array(L)==l)
        return np.interp(pc,interval,tot_vonq[spot[0]])

    def yfunc(p,l,pc):
        spot, = np.where(np.array(L)==l)

        a=np.interp(p,interval,tot_vonq[spot[0]])
        b = Spc(pc,l)
        return  a-b 

    def mean_yfunc(p,pc):
        return np.mean([yfunc(p,l,pc) for l in L])
        from scipy.optimize import minimize

    def R(params):
        pc,v = params
        #sum over all the square differences
        x_vals = [[xfunc(p,l,pc,v) for p in interval] for l in L]
        y_vals = [[yfunc(p,l,pc) for p in interval] for l in L]

        min_x = np.max([x[0] for x in x_vals]) #max for smallest value st all overlap
        max_x = np.min([x[-1] for x in x_vals]) # min again to take overlap
        xi = np.linspace(min_x,max_x)
        mean_x_vals = np.mean(x_vals,axis=0)
        mean_y_vals = [mean_yfunc(p,pc) for p in interval]
        
        def mean_y(x):
            return np.interp(x,mean_x_vals,mean_y_vals)
        
        return np.sum([[(np.interp(x,x_vals[i],y_vals[i]) - mean_y(x))**2 for x in xi] for i in range(len(L))]) 
    initial_guess = [0.2,3]
    res = scipy.optimize.minimize(R, initial_guess)
        
    """


    ppc,vv=py"res.x"

    @show(ppc)
    @show(vv)
    py"""
    ppc,vv=res.x
    x_vals = [[xfunc(p,l,ppc,vv) for p in interval] for l in L]
    y_vals = [[yfunc(p,l,ppc) for p in interval] for l in L]
    # mean_y_vals = [mean_yfunc(p,0.26) for p in interval]
    """
    p2 = plot(transpose(py"x_vals"),transpose(py"y_vals"),linewidth=3)
    display(p2)
end