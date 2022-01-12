function mcmc_run!(x_current,P_in,p,q_vec,q_mat,u_mat,ind_row::T,ind_chain::T) where {T<:Integer}
    accept_val = 0.0
    P_current = P_in    
    x_new = x_current.+(q_mat[ind_row,:,ind_chain].*q_vec)
    p_new = p(x_new)
    A = min(0,p_new-P_current) #calcualte acceptance ratio
    u = u_mat[ind_row,ind_chain]
        if log(u)<=A
            x_current[1:end] = x_new
            P_current = p_new
            accept_val = 1.0
        end  
    
    return P_current,accept_val
end

function n_mcmc_runs!(x_current,P_in,p,q_vec,q_mat,u_mat,ind_row_start::T,ind_chain::T,n_runs::T) where {T<:Integer}
    n_accepts = 0.0
    P_current = P_in
    for i = 1:n_runs
        ind_row=ind_row_start+i
        P_current,accept_val = mcmc_run!(x_current,P_current,p,q_vec,q_mat,u_mat,ind_row,ind_chain)
        n_accepts += accept_val
    end
    return P_current,n_accepts
end

function metropolis_mcmc(x_current_vec,P_current_vec,p,q_vec::Array{Float64,1},file_name::String
    ;n_chains::Integer=7,n_samples::Integer=100,burn_in::Integer=100,sample_freq::Integer=10)    

    @show n_draws = n_samples*sample_freq+burn_in #number of random number per chain
    n_features = length(q_vec)

    u_mat = rand(n_draws,n_chains)
    q_mat = rand(Normal(),(n_draws,n_features,n_chains))

    @show size(q_mat)

    samples_container = Array{Array{typeof(x_current_vec[1]),1},1}(undef,n_chains)
    P_samples_container = Array{Array{Float64,1},1}(undef,n_chains)

    Threads.@threads for j = 1:n_chains        
        core_id = j 

        x_current,P_current = copy(x_current_vec[j]), copy(P_current_vec[j])

        #Run burn-in
        println("Core id $(core_id): Start Burn-in")
        P_current,na = n_mcmc_runs!(x_current,P_current,p,q_vec,q_mat,u_mat,0,j,burn_in)
        println("Core id $(core_id): Done Burn-in")

        #Create container for samples
        samples = Array{typeof(x_current),1}(undef,n_samples)
        #Propability value of samples
        P_samples = zeros(n_samples)
        accept_rate = 0.0
    
        println("Core id $(core_id): Start Sampling")
        #Run samples
        for i = 1:n_samples
            println("Core id $(core_id): ---Sample: $(i)/$(n_samples)") 
            ind_row_start = burn_in+(i-1)*sample_freq        
            P_current,na = n_mcmc_runs!(x_current,P_current,p,q_vec,q_mat,u_mat,ind_row_start,j,sample_freq)
            accept_rate += na
            #Save sample
            samples[i] = copy(x_current)
            #Save propability values
            P_samples[i] = P_current                
        end
        println("Core id $(core_id): Done Sampling")
        accept_rate /= n_samples*sample_freq 
        println("Core id $(core_id): accept_rate = $(accept_rate)") 
        
        samples_container[j] = samples
        P_samples_container[j] = P_samples
    end
    save("./output/"*file_name*".jld","samples_container",samples_container,"P_samples_container",P_samples_container)
    return nothing
end