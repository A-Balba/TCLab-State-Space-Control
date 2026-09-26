function G = gram_matrix(x_lib, T_amb)
%GRAM_MATRIX Build the symmetric Gram matrix for the selected kernel.

N = size(x_lib, 2);
G = zeros(N, N);
for i = 1:N
    G(i, :) = kernel(x_lib(:, i), x_lib, T_amb);
end
end
