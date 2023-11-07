include("pt_mpo.jl")
include("func_utils.jl")
ITensors.set_warn_order(20)

using DelimitedFiles 
using Dates
using Plots
using BenchmarkTools

cd(@__DIR__)
pwd()

#create data
let
    sits = [4 6];
    interval = 0.0:0.1:1;
    num_samp=10;
    noise_val=0.0;
    fld = Dates.format(now(), "YYYY_mm_dd_HH_MM")
    mkdir(fld)
    flname="pt_mpo.csv"
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



#plotting
let 
    fln = "/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/2023_11_07_14_00/pt_mpo.csv"
    sits2,interval2,num_samp2,noise_val2,decays2,growths2 = open_csv(fln)

    p = plot(interval2,decays2,title=string("Bip_ent Gat: 2haar, varying meas_p"), label=string.(transpose(sits2)), linewidth=3,xlabel = "Meas_P", ylabel = L"$\textbf{S_{vn}}(L/2)$")
    display(p)
end