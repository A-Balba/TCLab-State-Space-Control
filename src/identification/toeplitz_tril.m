function T_mat = toeplitz_tril(A, B, C, D, L)

    T_mat = kron(eye(L), D);

    for i = 1: (L-1)
        T_mat = T_mat + kron(diag(ones(L-i,1), -i), C*(A^(i-1))*B);
    end
end

