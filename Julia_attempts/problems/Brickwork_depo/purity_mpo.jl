include("func_utils.jl")
include("brickwork_utils.jl")
ITensors.set_warn_order(20)

using Plots
using PyCall
using LaTeXStrings 


function main(N= 6, depth = 5, num_samp=2, noise_int = 0:0.1:0.5, meas_int = 0:0.1:0.3)
 #sweep noise and meas_p look at purity at end only
    purity_arr = zeros((length(noise_int),length(meas_int)))
    tot_arr=[]
    #loops for noise/meas_p
    for (noise_i,noise) in enumerate(noise_int)
        for (meas_i,meas_p) in enumerate(meas_int)
            @show(noise,meas_p)
            fin_arr = []
            #loop to generate average for samps
            for i in 1:num_samp
                s = siteinds("Qubit", N)
                psi2=randomMPS(s,3)
                rho = outer(psi2',psi2)

                #loop to perform the steps
                for _ in 1:depth
                    samp_row = gen_samp_row(N,meas_p)
                    rho = samp_mps(rho,s,samp_row)
                    rho = kraus_dephase(rho,s,noise)
                end
                #get data at end
                pur_arr = [purity(rho)]

                if length(fin_arr) == 0
                    fin_arr = pur_arr
                    # @show("butts")
                else
                    fin_arr = 2*mean([(i-1)/i*fin_arr,1/i*pur_arr])
                end
            end
        purity_arr[noise_i,meas_i] = fin_arr[1]
        end
    end
    return [N,depth,num_samp,noise_int,meas_int,purity_arr]
end

# dat = main()
# (N,depth,averages,noise_int,meas_int,purity_arr)=dat
# contour(meas_int,noise_int,purity_arr,xlabel = "meas", ylabel = "noise",title=string("Purity sites:$N, Depth:$depth"),fill=true,levels=20, color=:turbo)
