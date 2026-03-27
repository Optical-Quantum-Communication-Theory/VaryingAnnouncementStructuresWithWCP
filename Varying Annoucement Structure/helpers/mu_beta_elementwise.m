function muArray = mu_beta_elementwise(N,F,epsilon, options)
    %N: Total number of rounds
    %F: Table of expectation values
    %epsilon: Epsilon used for parameter estimation
    %t: Acceptance parameter array. Can pass either single value or
    %array
    %
    arguments
        N (1,1) double {mustBeGreaterThanOrEqual(N, 0)}
        F (:, :, :) double {mustBeGreaterThanOrEqual(F, 0)}
        epsilon (1,1) double {mustBeGreaterThan(epsilon, 0)}
        options.t (1,1) double = nan
        options.tau (:, :, :) double = nan
    end

    Fflat = reshape(F,[],1);

    %logic for t/tau vals, ensure defined and satisfy constraints
    if isnan(options.t) == isnan(options.tau)
        %check both defined or neither defined
        throw(MException('mu_beta_func:tVals', ...
        'Must specify exactly one of t or tau.'));
    elseif ~isnan(options.t)
        %if t defined
        try
            mustBeGreaterThanOrEqual(options.t, 0);
            %check >0
        catch error
            rethrow(error);
        end
        tFlat = ones(size(Fflat))*options.t;
    elseif ~isnan(options.tau)
        %if tau defined
        try
            mustBeGreaterThanOrEqual(options.tau, 0);
            mustBeCompatibleSize(options.tau, F);
            %check >0 and compatible size
        catch error
            rethrow(error);
        end
        tFlat = reshape(options.tau, [], 1); 
    end

    mu_init = ones(size(Fflat));
    
    for index=1:numel(Fflat)        
        %Calculate initial mu
        mu_init(index) = mu_beta_expec_floor(N,Fflat(index),tFlat(index),epsilon);
    end

    %Chose smallest mu
    muArray = reshape(mu_init, size(F)); 
end

function mu = mu_beta_expec_floor(N,F,t,epsilon)
    %N: Total number of rounds
    %F: expectation value for frequency
    %t: Acceptance parameter
    %epsilon: Epsilon used for parameter estimation

    %Lower range of the accepted parameters
    l1 = helper_minus(floor(N*(F-t)))-1;

    %Upper range of the accepted parameters
    l2 = floor(N*(F+t));
    
    %Function used for root finding
    fun = @(x) max(betainc(helper_01(1-(F-t-x)), N-l1, l1+1,"upper") , betainc(helper_01(1-(F+t+x)), helper_minus(N-l2), l2+1)) -epsilon;
    
    %Initial guesses
    x0 = 1/sqrt(N); % [0, 1-F-t];%muold(N,epsilon);
    
    %Mu given by root
    mu = fzero(fun,x0);
end

function value = helper_minus(x)
    %Function converts array to an array with non-negative entries
    value = zeros(size(x));
    for index = 1:length(x)
        if x(index) < 0
            value(index) = 0;
        else
            value(index) = x(index);
        end
    end
end

function value = helper_01(x)
    %Function converts array to an array with entries between 0 and 1
    value = zeros(size(x));
    for index = 1:length(x)
        if x(index) < 0
            value(index) = 0;
        elseif x(index) >1
            value(index) = 1;
        else
            value(index) = x(index);
        end
    end
end

function mustBeCompatibleSize(x,y)
    if  size(x) ~= size(y)
        throw(MException('mu_beta_func:FrequenciesTauMustMatch', ...
        'Dimension of the Frequency Array and t Values Array must match.'));
    end
end