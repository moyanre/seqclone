function [pdf_binom] = my_binopdf(r,n,p)

ade = gammaln(n+1) - gammaln(n-r+1) - gammaln(r+1) + r*log(p) + (n - r)*log(1-p);
pdf_binom = exp(vpa(ade));

end

