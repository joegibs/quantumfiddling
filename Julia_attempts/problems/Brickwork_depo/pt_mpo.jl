include("func_utils.jl")
include("brickwork_utils.jl")
ITensors.set_warn_order(20)

using Plots
using PyCall
using LaTeXStrings


function main(sits = [6],interval = 0.0:0.1:0.6,num_samp=1,noise=0.00)
    decays=[]
    growths = []
    for n in sits#[6,8,10]
        @show("number of sites: ",sits)
        N = n
        # cutoff = 1E-8
        steps = 8*N
        svns=[]
        negs = []
            for i in interval
                print("\n meas_p $i \n")
                svn,neg =do_trials(N,steps,i,num_samp,noise,"deph",false)

                append!(svns,svn)
                append!(negs,neg)
            end
        append!(decays,[svns])
        append!(growths,[negs])

        # append!(decays,[decay])
    end
    return [sits,interval,num_samp,noise,decays,growths]
end
function main_noise(sits = [6],interval = 0.0,num_samp=1,noise=0.0:0.01:0.1)
    decays=[]
    growths = []
    for n in sits#[6,8,10]
        @show("number of sites: ",sits)
        N = n
        # cutoff = 1E-8
        steps = 8*N
        svns=[]
        negs = []
            for i in noise
                print("\n meas_p $i \n")
                svn,neg =do_trials(N,steps,interval,num_samp,i,"deph",false)

                append!(svns,svn)
                append!(negs,neg)
            end
        append!(decays,[svns])
        append!(growths,[negs])

        # append!(decays,[decay])
    end
    return [sits,interval,num_samp,noise,decays,growths]
end

new_dat = main()
sits,interval,num_samp,noise_val,decays,growths = new_dat
p = plot([0.0:0.1:0.6...],real(growths),title=string("Bip_Neg Gat: 2haar, varying meas_p"), label=string.(transpose([6:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{N_{vn}}(L/2)$")

# p = plot(real(svns),title=string("Gate Rand", ", ", N, " qubit sites, varying meas_p"), label=string.(transpose([interval...])), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
# p = plot([0.0:0.05:1...],decays,title=string("Bip_ent Gat: 2haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
# p = plot([0.0:0.05:1...],real(growths),title=string("Bip_Neg Gat: 2haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{N_{vn}}(L/2)$")

# sits=[4,6,8]
# py"""
# import numpy as np
# import scipy
# L=$sits
# interval = $interval#[x/10 for x in range(9)]
# tot_vonq = $growths
# def xfunc(p,l,pc,v):
#     return (p-pc)*l**(1/v)

# def Spc(pc,l):
#     spot, = np.where(np.array(L)==l)
#     return np.interp(pc,interval,tot_vonq[spot[0]])

# def yfunc(p,l,pc):
#     spot, = np.where(np.array(L)==l)

#     a=np.interp(p,interval,tot_vonq[spot[0]])
#     b = Spc(pc,l)
#     return  a-b 

# def mean_yfunc(p,pc):
#     return np.mean([yfunc(p,l,pc) for l in L])
#     from scipy.optimize import minimize

# def R(params):
#     pc,v = params
#     #sum over all the square differences
#     x_vals = [[xfunc(p,l,pc,v) for p in interval] for l in L]
#     y_vals = [[yfunc(p,l,pc) for p in interval] for l in L]

#     min_x = np.max([x[0] for x in x_vals]) #max for smallest value st all overlap
#     max_x = np.min([x[-1] for x in x_vals]) # min again to take overlap
#     xi = np.linspace(min_x,max_x)
#     mean_x_vals = np.mean(x_vals,axis=0)
#     mean_y_vals = [mean_yfunc(p,pc) for p in interval]
    
#     def mean_y(x):
#         return np.interp(x,mean_x_vals,mean_y_vals)
    
#     return np.sum([[(np.interp(x,x_vals[i],y_vals[i]) - mean_y(x))**2 for x in xi] for i in range(len(L))]) 
# initial_guess = [0.0,0.1]
# res = scipy.optimize.minimize(R, initial_guess)
    
# """


# ppc,vv=py"res.x"

# py"""
# # L=[4,6]
# ppc,vv=res.x
# x_vals = [[xfunc(p,l,ppc,vv) for p in interval] for l in L]
# y_vals = [[yfunc(p,l,ppc) for p in interval] for l in L]
# # mean_y_vals = [mean_yfunc(p,0.26) for p in interval]
# """
# plot(transpose(py"x_vals"),transpose(py"y_vals"),linewidth=3)
