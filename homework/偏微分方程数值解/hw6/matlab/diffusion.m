function [ret] = diffusion(h, N, order)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

E = eye(N + 1);

Q{1}(:, :) = (Dplus(Dminus(E)) * h^-2);
Q{2}(:, :) = Q{1}(:, :) - (h^-2 / 6 * Dplus(Dminus(Dplus(Dminus(E)))));
Q{3}(:, :) = Q{2}(:, :) + (h^-2 / 30 .* Dplus(Dminus(Dplus(Dplus(Dminus(Dminus(E)))))));

ret = Q{order};
end