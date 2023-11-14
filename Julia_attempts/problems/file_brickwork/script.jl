using Plots
using DelimitedFiles 
using Dates
cd(@__DIR__)
include("brickwork_mpo_dephase_with_meas_amp_damp.jl")
include("brickwork_mpo_dephase_with_meas_one_plot.jl")
# using .brickwork_mpo_dephase_with_meas_amp_damp


# function arr_gen(noise,)
fld = Dates.format(now(), "YYYY_mm_dd_HH_MM")

meas = [0:0.1:0.9...]
noise=[0.0:0.1:1...]
decays =[[] for i=0:2]
growths=[[] for i=0:2]
Threads.@threads for i in 0:2
    meas=i*0.3.+[0:0.1:.3...]
    if i <2
    pop!(meas)
    end
    @show(meas)
    decs, gros, noi = main_proj(meas,50,noise)
    decays[i+1] = decs
    growths[i+1] = gros
end


p3 = plot(noise,[growths...],
    title=string("MPO Gate Rand qubit sites, varying meas_p"),
    # legend=:outertopright,
    # label=permutedims(labels3),
    linewidth=3,
    # xaxis=:log,
    xlabel = "γ",
    ylabel = L"$\textbf{N_{vn}}(L/2)$",
    palette=:batlow10
    )


data= [growths...]
arr = []
for i in data
    append!(arr,i)
end
heatmap(1:size(data,1),
    1:size(data,2), data,
    xlabel="x values", ylabel="y values",
    title="My title")

contourf(noise, meas, reduce(hcat,arr)', levels=20, color=:turbo,xlabel="noise", ylabel="y values",)


mkdir(fld)
flname="test_1.csv"
touch(fld*"/"*flname)

A =[[0:5...] [5:10...]]
open(fld*"/"*flname, "w") do io 
    writedlm( fld*"/"*flname,  A, ',')
end; 
tst = readdlm(fld*"/"*flname,',')

plot(tst)
flname="test_2.csv"
touch(fld*"/"*flname)

B =[[10:15...] [15:20...]]
#Writing contents to the file 
open(fld*"/"*flname, "w") do io 
    writedlm( fld*"/"*flname,  B, ',')
end; 
  
#readlm() method to read the DelimitedFiles 
tst = readdlm(fld*"/"*flname,',')

plot(tst)

files = readdir(fld)

arr1 = []
arr2 = []
for fl in files
    tst = transpose(readdlm(fld*"/"*fl,','))
    append!(arr1,tst[1,:])
    append!(arr2,tst[2,:])
end

