clearvars

lossdB = linspace(0,30,31);
lossEta = 10.^(-lossdB/10);
N = 1e6;
delta = 0.3; %range of the intensity you want to search around last optimal intensity
%initializing

keyrate = [];
mu = [];
loss = [];
full_results = {};
%initial value and initial bounds for lamb
% lamb = 10.^(-11/10);
lamb = .9;
lamb_up = lamb+delta;
lamb_low = 0.7;

for i = 1: numel(lossdB)
    eta = lossEta(i);
    qkdInput = NoPABFiniteWCP_Preset();
    qkdInput.addFixedParameter("N", N);
    qkdInput.addFixedParameter("eta",eta);
    qkdInput.addOptimizeParameter("lamb",struct("lowerBound",lamb_low, ...
        "initVal",lamb,"upperBound",lamb_up));
    try
    results = MainIteration(qkdInput);
    catch ME
    end 
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
    lamb_up = intensity + delta;
    lamb_low = max(0,intensity - delta);

    %check consecutive zero keyrate
    if sum(keyrate<0)>2
        break
    end
end

save("ScanLossNoPAB_1e6.mat","keyrate","mu","loss","full_results")