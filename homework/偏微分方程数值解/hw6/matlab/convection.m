function [ret] = convection(h, N, order)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

E = eye(N + 1);

Q{1}(:, :) = (Dzero(E) * h^-1);
Q{2}(:, :) = Q{1}(:, :) - (h^-1 / 6 * Dzero(Dplus(Dminus(E))));
Q{3}(:, :) = Q{2}(:, :) + (h^-1 / 30 .* Dzero(Dplus(Dplus(Dminus(Dminus(E))))));

ret = Q{order};
end