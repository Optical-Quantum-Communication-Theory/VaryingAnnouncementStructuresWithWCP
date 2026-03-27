function [gains,joints] = Produce_gains_joints_squashing(rhoXY)
% gains: this is the probablity of getting one pair of the annoucements combination
% joints: the probability distributions Prob(RY|AB) of one pair of
% announcements.
gain = zeros(2,4);
joint = cell(2,4);

for beta = 1:2
    for alpha= 1:4
        [x,y] = Prob_RY_given_AB_squashing(alpha,beta,rhoXY);
        gain(beta,alpha) = x;
        joint{beta,alpha} = y;
    end
end

gains = reshape(gain',[1,8])';
joints = reshape(joint',[1,8]);
end