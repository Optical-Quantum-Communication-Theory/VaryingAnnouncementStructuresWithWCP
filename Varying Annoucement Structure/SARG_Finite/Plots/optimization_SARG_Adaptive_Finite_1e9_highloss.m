%in case some of the points with high loss fail
lossdB = linspace(20,30,11);
lossEta = 10.^(-lossdB/10);
N = 1e9;
delta = 0.2; %range of the intensity you want to search around last optimal intensity
%initializing

keyrate = [];
mu = [];
loss = [];
full_results = {};
%initial value and initial bounds for lamb
lamb = 0.2;
lamb_up = 0.3;
lamb_low = 0.1;

for i = 1: numel(lossdB)
    eta = lossEta(i);
    qkdInput = SARG04_Finite_Adaptive_WCPPreset();
    qkdInput.addFixedParameter("N", N);
    qkdInput.addFixedParameter("eta",eta);
    qkdInput.addOptimizeParameter("lamb",struct("lowerBound",lamb_low, ...
        "initVal",lamb,"upperBound",lamb_up));
    qkdInput.addFixedParameter("misalignmentAngle",0);
    qkdInput.addFixedParameter("visibilityAngle",0);
    qkdInput.addFixedParameter("pz",0.5);

    results = MainIteration(qkdInput);
    
    %get the intensity and key
    intensity = results.currentParams.lamb;
    key = [results.keyRate];

    %store the results
    keyrate = [keyrate;key];
    mu = [mu,intensity];
    full_results = [full_results,{results}];
    loss = [loss,lossdB(i)];
    %updating the bounds and initial point for next itration
    lamb = intensity;
    lamb_up = intensity*(1+ delta);
    lamb_low = max(0,intensity*(1- delta)); 

    %check consecutive zero keyrate
    if sum(keyrate<0)>2
        break
    end
end

%%
save("Archived_SARG04_Loss_Finite_1e9_highloss20dB-30dB.mat","keyrate","mu","loss","full_results")