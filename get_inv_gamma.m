function p = get_inv_gamma(t,A,B)
    t = max(t,realmin);
    p = exp( A*log(B)-gammaln(A)-(A+1)*log(t)-B./t );
end