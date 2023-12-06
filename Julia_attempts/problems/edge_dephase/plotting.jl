

function plot_fldr(fldr)
    growths=[]
    decays =[]
    sits2=0;interval2=0;num_samp2=0;noise_val2=0;
    i=1
    for fle in readdir(fldr)
        # @show(fldr*'/'*fle)
        if i==1
            sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fldr*'/'*fle)
            growths = growths2
            decays = decays2
        else
            sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fldr*'/'*fle)
            # show(length(growths))
            growths = 2*mean([(i-1)/i*growths,1/i*growths2])
            decays = 2*mean([(i-1)/i*decays,1/i*decays2])
        end
        @show(noise_val2)
        i=i+1
        
    end
    # popfirst!(interval2)
    popfirst!(sits2)
    popfirst!(decays)
    popfirst!(growths)
    # popfirst!(sits2)
    # popfirst!(decays)
    # popfirst!(growths)
    # fln = "/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/2023_11_10_11_32/pt_mpo.csv"
    # sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fln)
    # @show real(growths)
    # @show size(growths)
    
    # grw2 = real(growths2)
    @show(string.(transpose(sits2)))
    for i in real(growths)
        @show(maximum(i))
    end
    hyup = basename(fldr)
    p = plot(interval2,real(decays),title=string("Bip_ent Gate: 2haar, varying meas_p $hyup"), label=string.(transpose(sits2)), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{N_{vn}}(L/2)$")
    display(p)
    pg = plot(interval2,real(growths),title=string("negativity $hyup"), label=string.(transpose(sits2)), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
    display(pg)

    # py"""
    # import numpy as np
    # import scipy
    # L=$sits2
    # interval = $interval2#[x/10 for x in range(9)]
    # tot_vonq = np.real($decays)
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
    # initial_guess = [0.2,3]
    # res = scipy.optimize.minimize(R, initial_guess)
        
    # """


    # ppc,vv=py"res.x"

    # @show(ppc)
    # @show(vv)
    # py"""
    # ppc,vv=res.x
    # x_vals = [[xfunc(p,l,ppc,vv) for p in interval] for l in L]
    # y_vals = [[yfunc(p,l,ppc) for p in interval] for l in L]
    # # mean_y_vals = [mean_yfunc(p,0.26) for p in interval]
    # """
    # p2 = plot(transpose(py"x_vals"),transpose(py"y_vals"),linewidth=3)
    # # display(p2)
end

let 
    # fldr = "/Users/joeg/Documents/GitHub/quan/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/trial_run/deph/tumfiddling/Julia_attempts/problems/Brickwork_depo/trial_run/damp/2023_11_16"
    path ="/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/trial_run/deph"
    fulldirpaths=filter(isdir,readdir(path,join=true))
    dirnames=basename.(fulldirpaths)
    for fldr in fulldirpaths
        # @show(fldr,"")
        show(basename(fldr))
        plot_fldr(fldr)
    end
end



let 
    fldr = "/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/trial_run/deph/2023_11_25_10_deph"
    plot_fldr(fldr)
end

arr = []
for i in [1:3...]
    append!(arr, [i])
end