function O = observability_matrix(A,C,L)
    [m,n] = size(C);
    O = zeros( m * L , n );
    for i = 1 : L
        O(m*(i-1)+1:m*i,:) = C * A^(i-1);
    end
end