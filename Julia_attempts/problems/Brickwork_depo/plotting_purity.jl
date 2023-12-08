

function plot_fldr(fldr)
    growths=[]
    decays =[]
    i=1
    for fle in readdir(fldr)
        # @show(fldr*'/'*fle)
        if i==1
            N,depth,num_samp,noise_int,meas_int,purity_arr2  = open_csv_purity(fldr*'/'*fle)
            purity_arr = purity_arr2
        else
            N2,depth2,num_samp2,noise_int2,meas_int2,purity_arr2 = open_csv_purity(fldr*'/'*fle)
            # show(length(growths))
            purity_arr = 2*mean([(i-1)/i*purity_arr,1/i*purity_arr2])
        end
        i=i+1
    end
    purity_arr = reshape(purity_arr,(length(noise_int),length(meas_int)))
    hyup = basename(fldr)
    p = contour(meas_int,noise_int,purity_arr,xlabel = "meas", ylabel = "noise",title=string("Purity sites:$N, Depth:$depth"),fill=true,levels=20, color=:turbo)
    display(p)

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
    fldr = "/Users/joeg/Documents/GitHub/quantumfiddling/Julia_attempts/problems/Brickwork_depo/2023_12_07_purity"
    plot_fldr(fldr)
end

arr = []
for i in [1:3...]
    append!(arr, [i])
end