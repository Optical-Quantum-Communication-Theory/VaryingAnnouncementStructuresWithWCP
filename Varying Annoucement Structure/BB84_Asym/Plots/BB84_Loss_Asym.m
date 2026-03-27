lossdB = linspace(0,30,31);
lossEta = 10.^(-lossdB/10);
delta = 0.3; %range of the intensity you want to search around last optimal intensity
%initializing
keyrate = [];
mu = [];
loss = [];
full_results = {};
%initial value and initial bounds for lamb
% lamb = 10.^(-11/10);
lamb = 1.1;
lamb_up = 1.3;
lamb_low = 0.8;

for i = 1: numel(lossdB)
    eta = lossEta(i);
    qkdInput = BB84WCPPreset();
    qkdInput.addFixedParameter("pz",0.1);
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
    lamb_up = intensity*(1+ delta);
    lamb_low = max(0,intensity*(1- delta)); 

    %check consecutive zero keyrate
    if sum(keyrate<0)>2
        break
    end
end

save("BB84_Loss_Asym.mat","keyrate","mu","loss","full_results");

