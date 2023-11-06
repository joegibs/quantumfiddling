include("brickwork_utils.jl")
include("func_utils.jl")
using Plots

function main(meas_ps=[0.0:0.2:1...],trials=10,noise=0.0)
    decays=[]
    svns=[]
    negs=[]
    for n in [4]
    steps = 8*n
    
    mut = []
    for i in meas_ps
        print("\n meas_p $i \n")
        svn,neg =do_trials(n,steps,i,trials,noise,"deph")
        avgsvn = [(svn[x]+svn[x+1])/2 for x in 1:2:(size(svn)[1]-1)]
        avgneg = [(neg[x]+neg[x+1])/2 for x in 1:2:(size(neg)[1]-1)]
  
        append!(svns,[avgsvn])
        append!(negs,[avgneg])
    end
    decay = [svns[i][end] for i in 1:size(svns)[1]]
    append!(decays,[decay])
    end
    p = plot(svns,title=string("MPO Gate Rand qubit sites, varying meas_p"), label=string.(transpose(meas_ps)), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{S_{vn}}(L/2)$")
    p2= plot(negs,title=string("MPO Gate Rand qubit sites, varying meas_p"), label=string.(transpose(meas_ps)), linewidth=3,xlabel = "Steps", ylabel = L"$\textbf{N_{vn}}(L/2)$")

    display(p)
    display(p2)
  end