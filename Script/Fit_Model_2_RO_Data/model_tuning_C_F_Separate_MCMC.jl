# Sampler in R^d
function RoSampler(log_p::Function, n_steps::Integer, x0::Array{AbstractFloat,1};r_burn_in::Float64=0.3)

    n_burn_in = floor(Integer,n*r_burn_in)
    n_samples = n_steps-n_burn_in

    # Initialise random walk sampler state: r.x current state, r.y proposal
    r = AdaptiveMCMC.RWMState(x0)

    # Initialise Adaptive Metropolis state (with default parameters)
    s = AdaptiveMCMC.AdaptiveMetropolis(x0)
    # Other adaptations are: AdaptiveScalingMetropolis,
    # AdaptiveScalingWithinAdaptiveMetropolis, and RobustAdaptiveMetropolis

    X = zeros(eltype(x0), length(x0), n_samples) # Allocate output storage    
    p_x = log_p(r.x)                     # = log_p(x0); the initial log target
    p_X = zeros(typeof(p_x), n_samples) # Allocate log target storage
    for k = 1:n_steps

        # Draw new proposal r.x -> r.y:
        AdaptiveMCMC.draw!(r, s)

        p_y = log_p(r.y)                      # Calculate log target at proposal
        alpha = min(one(p_x), exp(p_y - p_x)) # The Metropolis acceptance probability

        if rand() <= alpha
            p_x = p_y            
            # This 'accepts', or interchanges r.x <-> r.y:
            # (NB: do not do r.x = r.y; these are (pointers to) vectors!)
            AdaptiveMCMC.accept!(r)
        end

        # Do the adaptation update:
        AdaptiveMCMC.adapt!(s, r, alpha, k)

        if  k>n_burn_in
            #Save Sample and log target  
            X[:,k-n_burn_in] = r.x  # Save the current sample
            p_X[k-n_burn_in] = p_x  #Save the current log target
        end
     end
    return X,p_X
end

function log_post_distri_RO_CCPH(para::Array{Float64,1},parasym::Array{Symbol,1},ranges::Array{Tuple{Float64,Float64},1},
    RO_data::RO_raw_data;stand_type::String="Fertilized",Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit=nothing)
    if any(x->x==false,first.(ranges).<para.<last.(ranges))       
        post = -Inf
        return post
    else        
        try 
            ParaDict = CreateParaDict(parasym,para;ParaDictInit=ParaDictInit)

            model,weatherts,data = Get_Result_RO_CCPH(ParaDict,RO_data;stand_type=stand_type)
                
            logP = Calc_logP(model,data,ParaDict) 
            return logP       
        catch err
            println("Parameters Error: ", err)
            post = -Inf
            return post
        end
    end           
end

function run_sampler(file_name::String,ranges::Array{Tuple{Float64, Float64},1},parasym::Array{Symbol,1},x0::Array{AbstractFloat,1};
    N_samples::Integer=7000,N_BurnIn::Integer=3000,Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,ParaDictInit_C=nothing)    
    file_name_F = file_name*"_F"
    file_name_C = file_name*"_C"

    RO_data = Load_RO_data()  

    X,p_X = RoSampler(x->log_post_distri_RO_CCPH(x,parasym,ranges,RO_data;
    stand_type="Fertilized",Calc_logP::Function=Calc_logP_GPP_Ec,ParaDictInit=ParaDictInit_F),
    n_steps::Integer, x0;r_burn_in::Float64=0.3)
end
