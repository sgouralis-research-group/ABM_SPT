function x=forwardSubstitution(A,b)
    n = length(A);
    x = NaN(n,1);

    if cond(A) > 1/eps
        error('Matrix is singular')
    end

    for j=1:n
        x(j)=b(j)/A(j,j);
        b(j+1:n)=b(j+1:n)-A(j+1:n,j)*x(j);
    end
end