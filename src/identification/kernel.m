function K = kernel(x_query, x_lib, T_amb)
%KERNEL Degree-2 polynomial kernel with ambient-temperature centering.
%   K(j) = (1 + phi(x_query)' * phi(x_lib(:,j)))^2
%   where temperature components are shifted by T_amb.

x_query = x_query(:);
x_lib_shifted = x_lib;
x_query_shifted = x_query;

% The first four components are temperature history; the final two are heater inputs.
n_temp = min(4, size(x_query_shifted, 1));
x_query_shifted(1:n_temp) = x_query_shifted(1:n_temp) - T_amb;
x_lib_shifted(1:n_temp, :) = x_lib_shifted(1:n_temp, :) - T_amb;

K = zeros(1, size(x_lib_shifted, 2));
for j = 1:size(x_lib_shifted, 2)
    K(j) = (1 + x_query_shifted' * x_lib_shifted(:, j))^2;
end
end
