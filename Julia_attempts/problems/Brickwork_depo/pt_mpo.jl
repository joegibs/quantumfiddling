include("func_utils.jl")
include("brickwork_utils.jl")
ITensors.set_warn_order(20)

using Plots
using PyCall
using LaTeXStrings


function main_pt(sits = [4,6,8],interval = 0.0:0.2:0.8,num_samp=10,noise=0.0,noise_type="deph")
    Int_Svns_Mean=[]
    Int_Negs_Mean = []
    Int_Svns_Var=[]
    Int_Negs_Var = []
    for n in sits#[6,8,10]
        @show("number of sites: ",sits)
        N = n
        # cutoff = 1E-8
        steps = 8*N
        meas_p=0.
        svns_mean=[]
        svns_var=[]
        negs_mean = []
        negs_var=[]
            for i in interval
                print("\n meas_p $i \n")
                svn_mean,svn_var,neg_mean,neg_var =do_trials(N,steps,i,num_samp,noise,noise_type,false)

                append!(svns_mean,svn_mean)
                append!(svns_var,svn_var)
                append!(negs_mean,neg_mean)
                append!(negs_var,neg_var)
            end
        append!(Int_Svns_Mean,[svns_mean])
        append!(Int_Negs_Mean,[negs_mean])
        append!(Int_Svns_Var,[svns_var])
        append!(Int_Negs_Var,[negs_var])

        # append!(decays,[decay])
    end
    return [sits,interval,num_samp,noise,Int_Svns_Mean,Int_Negs_Mean,Int_Svns_Var,Int_Negs_Var]
end

# sits,interval,num_samp,noise,Int_Svns_Mean,Int_Negs_Mean,Int_Svns_Var,Int_Negs_Var = main_pt()
# svn_errs = sqrt.(hcat(Int_Svns_Var...) ./ num_samp)
# p = plot÷(real(svns),title=string("Gate Rand", ", ", N, " qubit sites, varying meas_p"), label=string.(transpose([interval...])), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
# p = plot([0.0:0.2:0.8...],real(Int_Svns_Mean),yerr = svn_errs,title=string("Bip_ent Gat: 2haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$",markerstrokecolor = :auto)
# p = plot([0.0:0.1:0.8...],real(growths),title=string("Bip_Neg Gat: 2haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{N_{vn}}(L/2)$")

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
