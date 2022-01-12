#Haario et al. 2001 An  adaptive  Metropolis  algorithm

#proposal function
function mcmc_proposal_func_ad(par,sig)
    return par.+rand(MvNormal(sig))
end

function mcmc_update_Σ_n_μ(x,Σ,μ,ϵ,s_d,T)
    
    μ_new = T*μ/(T+1)+x/(T+1)
    #Σ_new = Hermitian((T-1)/T*Σ+s_d/T*(T*μ*transpose(μ)-(T+1)*μ_new*transpose(μ_new)+x*transpose(x)+UniformScaling(ϵ)))
    Σ_new = (T-1)/T*Σ+s_d/T*(T*(μ*transpose(μ))-(T+1)*(μ_new*transpose(μ_new))+(x*transpose(x))+UniformScaling(ϵ))   

    return Σ_new,μ_new
end

#n runs of the Metropolis algorithm
function mcmc_run_ad!(x_current,P_current,p,sig,n::Integer,n_chains::Integer)
    n_accept = zeros(n_chains)

    Threads.@threads for j = 1:n_chains
    #for j = 1:n_chains
        for i = 1:n        
            x_new = mcmc_proposal_func_ad(x_current[j],sig[j]) #ceate new candidate
            p_new = p(x_new) #Probability of candidate 
            A = min(1,p_new/P_current[j]) #calcualte acceptance ratio
            u = rand() #Probability of accepting candidate 
            if u<=A
                x_current[j],P_current[j]= x_new,p_new
                n_accept[j] += 1
            end            
        end
    end
    return n_accept
end

#n runs of the Metropolis algorithm
function mcmc_run_ad!(x_current,P_current,p,Σ,μ,ϵ,s_d,T,n::Integer,n_chains::Integer)
    n_accept = zeros(n_chains)

    Threads.@threads for j = 1:n_chains
    #for j = 1:n_chains    
        T_local = T
        for i = 1:n        
            x_new = mcmc_proposal_func_ad(x_current[j],Σ[j]) #ceate new candidate
            p_new = p(x_new) #Probability of candidate 
            A = min(1,p_new/P_current[j]) #calcualte acceptance ratio
            u = rand() #Probability of accepting candidate 
            if u<=A
                x_current[j],P_current[j]= x_new,p_new
                n_accept[j] += 1
            end 
            Σ[j],μ[j] = mcmc_update_Σ_n_μ(reshape(x_current[j],(length(x_current[j]),1)),Σ[j],μ[j],ϵ,s_d,T_local)
            T_local += 1           
        end
    end
    T += n 
    return n_accept
end

function mcmc_create_Σ_n_μ!(x_current,P_current,p,sig,n::Integer,n_chains::Integer)
    
    μ = Array{Array{Float64,2},1}(undef,n_chains) #Sample mean
    Σ = Array{Array{Float64,2},1}(undef,n_chains) #Sample Covariance matrix

    Threads.@threads for j = 1:n_chains
    #for j = 1:n_chains    

        samples = Array{Float64,2}(undef,n,length(x_current[j]))

        for i = 1:n        
            x_new = mcmc_proposal_func_ad(x_current[j],sig[j]) #ceate new candidate
            p_new = p(x_new) #Probability of candidate 
            A = min(1,p_new/P_current[j]) #calcualte acceptance ratio
            u = rand() #Probability of accepting candidate 
            if u<=A
                x_current[j],P_current[j]= x_new,p_new                
            end
            samples[i,:] = x_current[j] 
        end
        μ[j] = transpose(mean(samples, dims=1))
        Σ[j] = cov(samples)
    end    
    return Σ,μ
end

function metropolis_mcmc_ad!(x_current,P_current,sig,p,s_d;
    n_samples::Integer=100,burn_in_init::Integer=100,burn_in::Integer=100,
    sample_freq::Integer=5,init_sample::Integer=300,ϵ::Float64=1.0e-9,every_n_sample_accept_check::Integer=10)    
    
    n_chains = length(x_current)    

    #Run initial burn-in
    println("Start: Init Burn-in")
    mcmc_run_ad!(x_current,P_current,p,sig,burn_in,n_chains)
    println("Done: Init Burn-in")

    #Create sample Covariance matrix matrix and mean
    println("Start: Create sample Σ and μ")
    Σ,μ = mcmc_create_Σ_n_μ!(x_current,P_current,p,sig,init_sample,n_chains)
    T = init_sample #Total number of iteration after initial burn-in
    println("Done: Create sample Σ and μ")

    #Run burn-in
    println("Start: Burn-in")
    mcmc_run_ad!(x_current,P_current,p,Σ,μ,ϵ,s_d,T,burn_in,n_chains)    
    println("Done: Burn-in")

    #Create container for samples
    samples = Array{typeof(x_current[1]),2}(undef,n_samples,n_chains)
    accept_rate = zeros(n_chains)
    accept_rate_monitor = zeros(n_chains)
   
    println("Start: Sampling")
    #Run samples
    for i = 1:n_samples
        println("---Sample: ",i,"/",n_samples)
        #Run mcmc
        n_accept = mcmc_run_ad!(x_current,P_current,p,Σ,μ,ϵ,s_d,T,sample_freq,n_chains)         
        accept_rate = accept_rate.+n_accept
        accept_rate_monitor = accept_rate_monitor.+n_accept
        #Save sample
        for j = 1:n_chains
            samples[i,j] = x_current[j]
        end

        if mod(i,every_n_sample_accept_check)==0
            println("Acceptance rate (short term): ",accept_rate_monitor/(every_n_sample_accept_check*sample_freq))    
            println("Acceptance rate (long term): ",accept_rate/(i*sample_freq)) 
            accept_rate_monitor = zeros(n_chains)        
        end
    end
    println("Done: Sampling")
    accept_rate /= n_samples*sample_freq

    return samples,accept_rate
end