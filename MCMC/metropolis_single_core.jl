function proposal(x::Array{Float64,1},q_vec::Array{Float64,1},n_features::Integer)
    y = x.+rand(Normal(),n_features).*q_vec
    return y
end

function mcmc_run!(log_p::Function,x::Array{Float64,1},p_x::Float64,q_vec::Array{Float64,1},n_features::Integer,n_runs::Integer)
    n_accepts = 0
    for i = 1:n_runs
        y = proposal(x,q_vec,n_features)
        p_y = log_p(y)
        alpha = min(one(p_x), exp(p_y - p_x)) # The Metropolis acceptance probability

        if rand() <= alpha
            # This 'accepts' y
            p_x = p_y            
            x[1:end] = y
            n_accepts += 1
        end
    end
    return n_accepts
end

function metropolis_mcmc(log_p::Function,x0::Array{Float64,1},q_vec::Array{Float64,1}
    ;n_samples::Integer=100,burn_in::Integer=100,sample_freq::Integer=5) 

    n_features = length(q_vec)
    x = copy(x0)
    p_x  = log_p(x)   # = log_p(x0); the initial log target

    X = zeros(eltype(x0), length(x0), n_samples) # Allocate output storage 
    p_X = zeros(typeof(p_x), n_samples) # Allocate log target storage

    #Run burn-in
    println("MCMC Sampler: Start Burn-in")
    burn_in_accepts = mcmc_run!(log_p,x,p_x,q_vec,n_features,burn_in)
    println("MCMC Sampler: Done Burn-in, accept rate: $(burn_in_accepts/burn_in)")

    #Run sampling
    println("MCMC Sampler: Start Sampling")
    sampling_accepts = 0
    for i = 1:n_samples
        sampling_accepts += mcmc_run!(log_p,x,p_x,q_vec,n_features,sample_freq)
        #Save sample
        @show X[:,i] = x
        p_X[i] = p_x
    end    
    println("MCMC Sampler: Done Sampling, accept rate: $(sampling_accepts/(n_samples*sample_freq))")
    return X,p_X
end