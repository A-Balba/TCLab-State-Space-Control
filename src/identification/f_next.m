function y_next = f_next(x_k, x_lib, T_amb, alpha)
%F_NEXT One-step prediction of TCLab outputs using kernel regression.

k = kernel(x_k, x_lib, T_amb);
y_next = (k * alpha).';
end
