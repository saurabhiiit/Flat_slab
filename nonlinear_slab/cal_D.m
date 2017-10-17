function [ D ] = cal_D(mod_of_elas )
%CAL_D Summary of this function goes here
%   Detailed explanation goes here
pois_ratio = 0.3;
% mod_of_elas = 2 * 10^11;
a = mod_of_elas * (1-pois_ratio) / ((1- 2 * pois_ratio) * (1 + pois_ratio));
b = mod_of_elas * pois_ratio / ((1- 2 * pois_ratio) * (1 + pois_ratio));
G = mod_of_elas / (2 * (1 + pois_ratio));
D = [ a b b 0 0 0;
      b a b 0 0 0;
      b b a 0 0 0;
      0 0 0 G 0 0;
      0 0 0 0 G 0;
      0 0 0 0 0 G;
    ];


end

