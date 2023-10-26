using ITensors
using Plots
using Statistics



function rho_to_dense(rho,s)
    Hitensor = ITensor(1.)
    N = length(s)
    for i = 1:N
        Hitensor *= rho[i]
    end
  
    A=Array(Hitensor,prime(s),s)
    return reshape(A,2^N,2^N)
  end

function rec_ent(rho,b,s)
    n = length(rho)
    # orthogonalize!(rho,b)
    rho_temp = deepcopy(rho)
    # s = siteinds("Qubit",n) 

    #contract half   x x x x | | | |
    L = ITensor(1.0)
    for i = 1:b
      L *= tr(rho_temp[i])
    end
    # absorb
    rho_temp[b+1] *= L
    # no longer a proper mpo
    M =MPO(n-b)
    for i in 1:(n-b)
        M[i]=rho_temp[b+i]
    end
    #turn this mpo into a single tensor
    T = prod(M)
 
    _,S,_ = svd(T,s)
    SvN = 0.0
    for n in 1:dim(S, 1)
      p = S[n,n]
      if p != 0
        SvN -= p * log2(p)
      end
    end
    return SvN
end

function rec_ent_mps(psi,b)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    if length(psi)==2
        if dim(U,1) != dim(U,2)
            for n in diag(reshape(U.tensor,(dim(U,1),dim(U,2))))
                p = n^2
                if p>0.
                    SvN -= p * log2(p)
                end
            end
        else
        for n=1:dim(U, 1)
            p = U[n,n]^2
                if p>0.
                    SvN -= p * log2(p)
                end
            end
        end
    else
        for n=1:dim(S, 1)
            
        p = S[n,n]^2
            if p>0.
                SvN -= p * log2(p)
            end
        end
    end
    return SvN
end
function vectorize(rho,s)
    arr=[]
    for i in 1:length(rho)
        test_tens = rho[i]
        site_inds=inds(test_tens)
        (a,b)=site_inds[[occursin("Site,n=",string(tags(i))) for i in site_inds]]
        C=combiner(a,b; tags=tags(s[i]))
    
        append!(arr,[ITensor(C*test_tens)])
    end
    
    fin =MPS(s)
    for i in 1:length(fin)
        fin[i]=arr[i]
        ind=inds(fin[i])[[occursin("Site,n=",string(tags(i))) for i in inds(fin[i])]]
        swapinds!(fin[i],ind[1],s[i])
    end
    return fin
end

function partial_transpose(A::MPO, sites)
    A = copy(A)
    for n in sites
      A[n] = swapinds(A[n], siteinds(A, n)...)
    end
    return A
  end
function trace_norm_dense(rho)
    e,_=eigen(rho)
    return sum(abs.(e))
end
function negativity(A::MPO, b, s)
    n = length(A)
    orthogonalize!(A,b)
    rho_temp = deepcopy(A)
    
    M=partial_transpose(rho_temp,[b:n-b...])

    #turn this mpo into a single tensor

    T = rho_to_dense(M,s)
    
    
    return  (trace_norm_dense(T)-1)/2
end

function log_negativity(A::MPO, b, s)
    neg = negativity(A, b, s)
    return log2(2*neg+1)
end
function vect_to_rho(vect,sm,s)
    arr=[]
    rho = MPO(sm)
    for i in 1:length(rho)
        test_tens = rho[i]
        site_inds=inds(test_tens)
        (a,b)=site_inds[[occursin("Site,n=",string(tags(i))) for i in site_inds]]
        C=combiner(a,b; tags=tags(s[i]))
        swapinds!(C,inds(C)[1],s[i])
        
        append!(arr,[ITensor(dag(C)*vect[i])])
    end
        # sites = new_sites#siteinds("Qudit", N; dim=4)
    for i in 1:length(rho)
        rho[i]=arr[i]
    end
    return rho
end

let #entropy vs time
    arr_base_ent =[]
    arr_base_neg=[]
    arr_deph_ent=[]
    arr_deph_neg=[]
    N = 8
    cutoff = 1E-8
    tau = 0.1
    ttotal = 20.0
    for _ in 1:1

        # Make an array of 'site' indices
        s = siteinds("S=1/2", N)
        
        # Make gates (1,2),(2,3),(3,4),...
        gates = ITensor[]
        for j in 1:(N - 1)
            s1 = s[j]
            s2 = s[j + 1]
            hj =
                op("Sz", s1) * op("Sz", s2) +
                op("Sx", s1) * op("Sx", s2) +
                op("Sy", s1) * op("Sy", s2) 
                # op("Sx", s1)*op("I",s2)
            Gj = exp(-im * tau / 2 * hj)
            push!(gates, Gj)
        end
        # Include gates in reverse order too
        # (N,N-1),(N-1,N-2),...
        append!(gates, reverse(gates))

        # Initialize psi to be a product state (alternating up and down)
        psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")


        # Compute and print <Sz> at each time step
        # then apply the gates to go to the next time
        ent_base=[]
        neg_base=[]
        for t in 0.0:tau:ttotal
            append!(ent_base,rec_ent_mps(psi,Int(N/2)))
            append!(neg_base,log_negativity(outer(psi',psi),Int(N/2),s))

            #   println("$t $Sz")

            t≈ttotal && break

            psi = apply(gates, psi; cutoff)
            normalize!(psi)
        end
        @show("base is done")
        # Make an array of 'site' indices
        append!(arr_base_ent,[ent_base])
        append!(arr_base_neg,[neg_base])
    # end
        s = siteinds("Qudit", N; dim=4)
        sm = siteinds("S=1/2", N)

        # Make gates (1,2),(2,3),(3,4),...

        gates = ITensor[]
        for j in 1:(N - 1)
            s1 = sm[j]
            s2 = sm[j + 1]
            Sz1=op("Sz",s1).tensor
            Sz2 = op("Sz",s2).tensor
            Sy1=op("Sy",s1).tensor
            Sy2 = op("Sy",s2).tensor
            Sx1=op("Sx",s1).tensor
            Sx2 = op("Sx",s2).tensor
            SI1 = op("I",s1).tensor
            SI2 = op("I",s2).tensor
            gamma = √(.0001/2)*im

            hj =
            1*op(kron(Sz1,SI1),s[j])* op(kron(Sz2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sz1)),s[j]) * op(kron(SI2,transpose(Sz2)),s[j+1]) +
            1*op(kron(Sy1,SI1),s[j])* op(kron(Sy2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sy1)),s[j]) * op(kron(SI2,transpose(Sy2)),s[j+1]) +
            1*op(kron(Sx1,SI1),s[j])* op(kron(Sx2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,transpose(Sx2)),s[j+1])+ 
            # 1*op(kron(Sx1,SI1),s[j])* op(kron(SI2,SI2),s[j+1]) +
            # -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,SI2),s[j+1])+

            gamma*op(kron(Sz1,Sz1),s[j]) * op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(Sz1*Sz1,SI1),s[j])* op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(SI1,Sz1*Sz1),s[j])* op(kron(SI2,SI2),s[j+1])

            Gj = exp(-im * tau / 2 * hj)
            push!(gates, Gj)
        end
        # for j in 1:(N)
        #     gamma = √(2/2)*im
        #     s1 = sm[j]
        #     pz = op("Sz",s1).tensor
        #     pI=op("I",s1).tensor
        #     hj =gamma*op(kron(pz,pz),s[j])+
        #     -1/2*gamma*op(kron(pz*pz,pI),s[j])+
        #     -1/2*gamma*op(kron(pI,pz*pz),s[j])
        
        #     Gj = exp(-im * tau / 2 * hj)
        #     push!(gates, Gj)
        # end
        # Include gates in reverse order too
        # (N,N-1),(N-1,N-2),...
        append!(gates, reverse(gates))



        # Initialize psi to be a product state (alternating up and down)
        psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")#rand(Int)) ? "Up" : "Dn")
        rho = outer(psi1',psi1)
        psi = vectorize(rho,s)

        ent_z=[]
        neg_z=[]
        for t in 0.0:tau:ttotal
            append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
            append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
            #   println("$t $Sz")

            t≈ttotal && break

            psi = apply(gates, psi; cutoff)
            #   normalize!(psi)
        end
        append!(arr_base_ent,[ent_base])
        append!(arr_base_neg,[neg_base])
        append!(arr_deph_ent,[ent_z])
        append!(arr_deph_neg,[neg_z])
    end
    # @show( arr_base_ent)
    labels= ["base ent" "base neg" "dephase ent" "dephase neg"]
    return plot([0.0:tau:ttotal...],[arr_base_ent[1],arr_base_neg[1],arr_deph_ent[1],arr_deph_neg[1]],
     linewidth=3,
     legend=true,
     label=labels,
     xlabel="Time",
     ylabel="Entropy / Negativity")

    # return plot([0.0:tau:ttotal...],[arr_base_ent,arr_base_neg], linewidth=3)
end

let #ent vs gamma
    arr_base_ent =[]
    arr_base_neg=[]
    arr_deph_ent=[]
    arr_deph_neg=[]
    N = 4
    cutoff = 1E-8
    tau = 0.1
    ttotal = 2*N
    for gam in 1 ./(exp10.([0:0.25:4...]))

        # Make an array of 'site' indices
        # s = siteinds("S=1/2", N)
        
        # # Make gates (1,2),(2,3),(3,4),...
        # gates = ITensor[]
        # for j in 1:(N - 1)
        #     s1 = s[j]
        #     s2 = s[j + 1]
        #     hj =
        #         op("Sz", s1) * op("Sz", s2) +
        #         op("Sx", s1) * op("Sx", s2) +
        #         op("Sy", s1) * op("Sy", s2) 
        #         # op("Sx", s1)*op("I",s2)
        #     Gj = exp(-im * tau / 2 * hj)
        #     push!(gates, Gj)
        # end
        # # Include gates in reverse order too
        # # (N,N-1),(N-1,N-2),...
        # append!(gates, reverse(gates))

        # # Initialize psi to be a product state (alternating up and down)
        # psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")


        # # Compute and print <Sz> at each time step
        # # then apply the gates to go to the next time
        # ent_base=[]
        # neg_base=[]
        # for t in 0.0:tau:ttotal
        #     append!(ent_base,rec_ent_mps(psi,Int(N/2)))
        #     append!(neg_base,log_negativity(outer(psi',psi),Int(N/2),s))

        #     #   println("$t $Sz")

        #     t≈ttotal && break

        #     psi = apply(gates, psi; cutoff)
        #     normalize!(psi)
        # end
        # @show("base is done")
        # # Make an array of 'site' indices
        # append!(arr_base_ent,[ent_base])
        # append!(arr_base_neg,[neg_base])
    # end
        s = siteinds("Qudit", N; dim=4)
        sm = siteinds("S=1/2", N)

        # Make gates (1,2),(2,3),(3,4),...

        gates = ITensor[]
        for j in 1:(N - 1)
            s1 = sm[j]
            s2 = sm[j + 1]
            Sz1=op("Sz",s1).tensor
            Sz2 = op("Sz",s2).tensor
            Sy1=op("Sy",s1).tensor
            Sy2 = op("Sy",s2).tensor
            Sx1=op("Sx",s1).tensor
            Sx2 = op("Sx",s2).tensor
            SI1 = op("I",s1).tensor
            SI2 = op("I",s2).tensor
            gamma = √(gam/2)*im

            hj =
            1*op(kron(Sz1,SI1),s[j])* op(kron(Sz2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sz1)),s[j]) * op(kron(SI2,transpose(Sz2)),s[j+1]) +
            1*op(kron(Sy1,SI1),s[j])* op(kron(Sy2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sy1)),s[j]) * op(kron(SI2,transpose(Sy2)),s[j+1]) +
            1*op(kron(Sx1,SI1),s[j])* op(kron(Sx2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,transpose(Sx2)),s[j+1])+ 
            # 1*op(kron(Sx1,SI1),s[j])* op(kron(SI2,SI2),s[j+1]) +
            # -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,SI2),s[j+1])+

            gamma*op(kron(Sz1,Sz1),s[j]) * op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(Sz1*Sz1,SI1),s[j])* op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(SI1,Sz1*Sz1),s[j])* op(kron(SI2,SI2),s[j+1])

            Gj = exp(-im * tau / 2 * hj)
            push!(gates, Gj)
        end
        # for j in 1:(N)
        #     gamma = √(2/2)*im
        #     s1 = sm[j]
        #     pz = op("Sz",s1).tensor
        #     pI=op("I",s1).tensor
        #     hj =gamma*op(kron(pz,pz),s[j])+
        #     -1/2*gamma*op(kron(pz*pz,pI),s[j])+
        #     -1/2*gamma*op(kron(pI,pz*pz),s[j])
        
        #     Gj = exp(-im * tau / 2 * hj)
        #     push!(gates, Gj)
        # end
        # Include gates in reverse order too
        # (N,N-1),(N-1,N-2),...
        append!(gates, reverse(gates))



        # Initialize psi to be a product state (alternating up and down)
        psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")#rand(Int)) ? "Up" : "Dn")
        rho = outer(psi1',psi1)
        psi = vectorize(rho,s)

        ent_z=[]
        neg_z=[]
        for t in 0.0:tau:ttotal
            # append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
            # append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
            # #   println("$t $Sz")

            t≈ttotal && break

            psi = apply(gates, psi; cutoff)
            #   normalize!(psi)
        end

        append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
        append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
            # #   println("$t $Sz")
        # @show(ent_z)
        append!(arr_deph_ent,ent_z[1])
        append!(arr_deph_neg,neg_z[1])
        @show(gam)
    end
    # @show()
    labels= ["base ent" "base neg" "dephase ent" "dephase neg"]
    return plot(1 ./(exp10.([0:0.25:4...])),[arr_deph_ent,arr_deph_neg],
     linewidth=3,
     legend=true,
     label=labels,
     xaxis=:log,
     xlabel="γ",
     ylabel="Entropy / Negativity")

    # return plot([0.0:tau:ttotal...],[arr_base_ent,arr_base_neg], linewidth=3)
end

let #multiplot

    arr_deph_ent=[]
    arr_deph_neg=[]
    N = 8
    cutoff = 1E-8
    tau = 0.1
    labels=[]
    gamint = 0:.5:4
    for ttotal in (1:2:4).*N
        arr_temp_ent =[]
        arr_temp_neg=[]
    for gam in 1 ./(exp10.([gamint...]))

        
        s = siteinds("Qudit", N; dim=4)
        sm = siteinds("S=1/2", N)

        # Make gates (1,2),(2,3),(3,4),...

        gates = ITensor[]
        for j in 1:(N - 1)
            s1 = sm[j]
            s2 = sm[j + 1]
            Sz1=op("Sz",s1).tensor
            Sz2 = op("Sz",s2).tensor
            Sy1=op("Sy",s1).tensor
            Sy2 = op("Sy",s2).tensor
            Sx1=op("Sx",s1).tensor
            Sx2 = op("Sx",s2).tensor
            SI1 = op("I",s1).tensor
            SI2 = op("I",s2).tensor
            gamma = √(gam/2)*im

            hj =
            1*op(kron(Sz1,SI1),s[j])* op(kron(Sz2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sz1)),s[j]) * op(kron(SI2,transpose(Sz2)),s[j+1]) +
            1*op(kron(Sy1,SI1),s[j])* op(kron(Sy2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sy1)),s[j]) * op(kron(SI2,transpose(Sy2)),s[j+1]) +
            1*op(kron(Sx1,SI1),s[j])* op(kron(Sx2,SI2),s[j+1]) +
            -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,transpose(Sx2)),s[j+1])+ 
            # 1*op(kron(Sx1,SI1),s[j])* op(kron(SI2,SI2),s[j+1]) +
            # -1*op(kron(SI1,transpose(Sx1)),s[j]) * op(kron(SI2,SI2),s[j+1])+

            gamma*op(kron(Sz1,Sz1),s[j]) * op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(Sz1*Sz1,SI1),s[j])* op(kron(SI2,SI2),s[j+1])+
                -1/2*gamma*op(kron(SI1,Sz1*Sz1),s[j])* op(kron(SI2,SI2),s[j+1])

            Gj = exp(-im * tau / 2 * hj)
            push!(gates, Gj)
        end
 
        append!(gates, reverse(gates))



        # Initialize psi to be a product state (alternating up and down)
        psi1 = MPS(sm, n -> isodd(n) ? "Up" : "Dn")#rand(Int)) ? "Up" : "Dn")
        rho = outer(psi1',psi1)
        psi = vectorize(rho,s)

        ent_z=[]
        neg_z=[]
        for t in 0.0:tau:ttotal
            # append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
            # append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
            # #   println("$t $Sz")

            t≈ttotal && break

            psi = apply(gates, psi; cutoff)
            #   normalize!(psi)
        end

        append!(ent_z,rec_ent(vect_to_rho(psi,sm,s),Int(N/2),sm))
        append!(neg_z,log_negativity(vect_to_rho(psi,sm,s), Int(N/2), sm))
            # #   println("$t $Sz")
        # @show(ent_z)
        append!(arr_temp_ent,ent_z[1])
        append!(arr_temp_neg,neg_z[1])
        
    end
    append!(arr_deph_ent,[arr_temp_ent])
    append!(arr_deph_neg,[arr_temp_neg])
    push!(labels,"base ent $ttotal")
    push!(labels,"base neg $ttotal")
    @show(ttotal)
    end
    plot_data = collect(Iterators.flatten(zip(arr_deph_ent,arr_deph_neg)))

    return plot(1 ./(exp10.([gamint...])),[plot_data...],
     linewidth=3,
     legend=:outertopright,
     label=permutedims(labels),
     xaxis=:log,
     xlabel="γ",
     ylabel="Entropy / Negativity",
     palette=:batlow25)

    # return plot([0.0:tau:ttotal...],[arr_base_ent,arr_base_neg], linewidth=3)
end