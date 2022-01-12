function Draw_P_vs_samples(samples,P_samples,file_name::String)
    @show n_chains = length(samples)
    @show n_samples = length(samples[1])
    @show n_features = length(samples[1][1])    

    sample_mat = Array{Float64,2}(undef,n_chains*n_samples,n_features)
    P_mat =  Array{Float64,1}(undef,n_chains*n_samples)

    for i = 1:n_chains
        for j = 1:n_samples
            sample_mat[j+(i-1)*n_samples,:] = samples[i][j]
            P_mat[j+(i-1)*n_samples] = P_samples[i][j]
        end
    end

    plot(;layout=n_features,legends=false)
    for i = 1:n_features
        plot!(sample_mat[:,i],P_mat,subplot=i)
    end
    
    savefig("./plots/"*file_name*"_P_vs_Sample.svg") 
end

function Draw_trace_plot(samples,file_name::String)
    @show n_chains = length(samples)
    @show n_samples = length(samples[1])
    @show n_parameters = length(samples[1][1])  

    #plot(layout = n_parameters,legends=false)
    for i = 1:n_parameters
        plot(;legends=false)
        for j = 1:n_chains
            x_val = 1:n_samples 
            y_val = [samples[j][k][i] for k = 1:n_samples]

            plot!(x_val,y_val)
        end
        savefig("./plots/trace_plot/"*file_name*"_$(i).svg")
    end
end

function Draw_histogram(samples,file_name::String)
    @show n_chains = length(samples)
    @show n_samples = length(samples[1])
    @show n_features = length(samples[1][1])    

    sample_mat = Array{Float64,2}(undef,n_chains*n_samples,n_features)

    for i = 1:n_chains
        for j = 1:n_samples
            sample_mat[j+(i-1)*n_samples,:] = samples[i][j]
        end
    end

    histogram(sample_mat,layout=n_features,legend=false)
    savefig("./plots/"*file_name*".svg") 
end

function Find_max_sample(samples,P_samples)
    @show n_chains = length(samples)
    @show n_samples = length(samples[1])
    @show n_features = length(samples[1][1]) 

    max_P_vec = Array{Float64,1}(undef,n_chains)
    max_ind_vec = Array{Int,1}(undef,n_chains)

    for i = 1:n_chains
        P_vex = [P_samples[i][j] for j = 1:n_samples]
        max_P_vec[i] = maximum(P_vex)
        max_ind_vec[i] = argmax(P_vex)
    end

    max_chains = argmax(max_P_vec)
    max_sample = max_ind_vec[max_chains]
    max_features = copy(samples[max_chains][max_sample])

    return max_features
end