lossdB = linspace(0,40,41);
lossEta = 10.^(-lossdB/10);
N = 1e9;
delta = 0.3; %range of the intensity you want to search around last optimal intensity
%initializing

keyrate = [];
mu = [];
loss = [];
full_results = {};
%initial value and initial bounds for lamb
% lamb = 10.^(-11/10);
lamb = 1;
lamb_up = lamb+delta;
lamb_low = 0.7;

for i = 1: numel(lossdB)
    eta = lossEta(i);
    qkdInput = BB84_Finite_WCPPreset();
    qkdInput.addFixedParameter("N", N);
    qkdInput.addFixedParameter("eta",eta);
    qkdInput.addOptimizeParameter("lamb",struct("lowerBound",lamb_low, ...
        "initVal",lamb,"upperBound",lamb_up));
    results = MainIteration(qkdInput);  
    %get the intensity and key
    intensity = results.currentParams.lamb;
    key = [results.keyRate];

    %store the results
    keyrate = [keyrate;key];
    mu = [mu,intensity];
    full_results = [full_results,{results}];
    loss = [loss,lossdB(i)];
    %updating the bounds and initial point for next iteration
    lamb = intensity;
    lamb_up = intensity + delta;
    lamb_low = max(0,intensity - delta);

    %check consecutive zero keyrate
    if sum(keyrate<0)>2
        break
    end
end

save("BB84_LossOnly_Finite_0_40dB_1e9.mat","keyrate","mu","loss","full_results")