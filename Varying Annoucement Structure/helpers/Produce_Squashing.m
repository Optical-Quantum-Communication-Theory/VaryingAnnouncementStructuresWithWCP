function [squashing] = Produce_Squashing(a)
% A function to the squashing map
% Double clicks are randomly assigned to one of the results.
% different a stands for different ordering of clicks
% a =0 : Cols: vac Hsc Vsc HVdc vac Dsc Asc DAdc ; Rows:Hsc Vsc Dsc Asc vac
% a =1 : Cols: Hsc Vsc HVdc vac Dsc Asc DAdc vac ; Rows:Hsc Vsc Dsc Asc vac
if a == 0
    squashing = [[0;0;0;0;1],[1;0;0;0;0],[0;1;0;0;0],[1/2;1/2;0;0;0],[0;0;0;0;1],[0;0;1;0;0],[0;0;0;1;0],[0;0;1/2;1/2;0]];
end

if a == 1
    squashing = [[1;0;0;0;0],[0;1;0;0;0],[1/2;1/2;0;0;0],[0;0;0;0;1],[0;0;1;0;0],[0;0;0;1;0],[0;0;1/2;1/2;0],[0;0;0;0;1]];
end