using ITensors
include("func_utils.jl")
include("brickwork_utils.jl")

N=4
s = siteinds("Qubit", N)
psi1 = productMPS(s, "Up" )
psi2=randomMPS(s,2)
rho = outer(psi2',psi2)

rho_deph = kraus_dephase(rho,s,0.5)

purity(rho)
purity(rho_deph)

using Plots
using PyCall
using LaTeXStrings


sits = [4 6]
interval = 0.28
num_samp=5
noise=0.0:0.1:0.3
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
            print("\n noise $i \n")
            svn,neg =do_trials(N,steps,interval,num_samp,i,"deph",true)

            append!(svns,svn)
            append!(negs,neg)
        end
    append!(decays,[svns])
    append!(growths,[negs])

    # append!(decays,[decay])
end
# return [sits,interval,num_samp,noise,decays,growths]

# new_dat = main_noise()
# sits,interval,num_samp,noise_val,decays,growths = new_dat
p = plot([1:128...],real(decays[1]),title=string("Bip_Neg Gat: 2haar, varying meas_p"), label=string.(transpose([4:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = "Purity")

p = plot(noise,real(growths),title=string("Bip_Neg Gat: 2haar, varying meas_p"), label=string.(transpose([6:2:14...])), linewidth=3,xlabel = "Meas_P", ylabel = "Purity")

#not sure if this is working right

let #look at purity for fixed noise xor meas_p and dlook at purity over time
    tot_arr=[]
    for nm in 0:0.1:0.5
        fin_arr = []
        for i in 1:20
            N=6
            meas_p=0.1
            noise = nm
            s = siteinds("Qubit", N)
            psi2=randomMPS(s,2)
            rho = outer(psi2',psi2)

            pur_arr = [purity(rho)]




            for _ in 1:20
                samp_row = gen_samp_row(N,meas_p)
                rho = samp_mps(rho,s,samp_row)
                rho = rho/tr(rho)
                rho = kraus_dephase(rho,s,noise)
                rho = rho/tr(rho)
                append!(pur_arr, purity(rho))

            end
            if length(fin_arr) == 0
                fin_arr = pur_arr
                @show("butts")
            else
                fin_arr = 2*mean([(i-1)/i*fin_arr,1/i*pur_arr])
            end
            
        end
        append!(tot_arr,[fin_arr])
    end
    plot(tot_arr,linewidth=3,xlabel = "time_step", ylabel = "Purity",title=string("Purity"))
end

let #sweep noise and meas_p look at purity at end only
    noise_int = 0:0.01:0.5
    meas_int = 0:0.01:0.3
    purity_arr = zeros((length(noise_int),length(meas_int)))
    tot_arr=[]
    #loops for noise/meas_p
    for (noise_i,noise) in enumerate(noise_int)
        for (meas_i,meas_p) in enumerate(meas_int)
            @show(noise,meas_p)
            fin_arr = []
            #loop to generate average for samps
            for i in 1:20
                N=6
                s = siteinds("Qubit", N)
                psi2=randomMPS(s,3)
                rho = outer(psi2',psi2)

                #loop to perform the steps
                for _ in 1:5
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
    contour(meas_int,noise_int,purity_arr,xlabel = "meas", ylabel = "noise",title=string("Purity"),fill=true,levels=20, color=:turbo)
end

