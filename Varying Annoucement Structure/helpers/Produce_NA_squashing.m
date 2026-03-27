function [NA] = Produce_NA_squashing()
%generate the map that map the correct probability of A given X
%A: Alice's announcement
%X: Alice;s measurement

% The matrix has the form NA = sum p(a|x) S(y|y_unsquashed) |a>|x>|y><x|<y_unsquashed|
% where S is the squashing map, mapping the raw click patterns y_unsquashed
% to y.

% For which ever signal Alice sent (e.g H in Z basis), she randomly (p=1/2) pick a state from
% the other basis (X basis) and make an annoucement of these two states.
% There are 4 possible annoucements (HD, HA, VD, VA)
prob_a_given_x = [1/2,0,1/2,0;
                  0,1/2,1/2,0;
                  0,1/2,0,1/2;
                  1/2,0,0,1/2];

NA = 0;
for alpha = 1:4
    for x = 1:4
        NA = NA + kron(zket(4,alpha),zket(4,x)*zket(4,x)') * prob_a_given_x(alpha,x);
    end
end
squashing = Produce_Squashing(0);%squashing map needs to be applied on the statistics first (for Bob's sake)
NA = kron(NA,squashing);
end