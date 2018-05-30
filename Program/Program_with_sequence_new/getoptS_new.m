function out = getoptS(X,Y,M_E,E,lam1,lam2,lam3,SubSeqSM,SubSeqSMsera)

[n, r] = size(X);
[m, r] = size(Y);
C = X' * ( M_E ) * Y ; C = C(:) ;

LB=diag(sum(SubSeqSM))-SubSeqSM;
LC=diag(sum(SubSeqSMsera))-SubSeqSMsera;

for i = 1:r
        for j = 1:r
                ind = (j-1)*r + i ;
                temp1 = X' * (  (X(:,i) * Y(:,j)').*E ) * Y ;
                temp2 = (lam2)*X' * LB* (X(:,i) * Y(:,j)')  * Y ;
                temp3 = (lam3)*Y' * LC* (Y(:,i) * X(:,j)')  * X ;
                temp=temp1+temp2+temp3;
                A(:,ind) = temp(:) ;
			%	A(ind, ind) = A(ind, ind)+ lam1*m*n ;
        end
end

S = A\C ;
out = reshape(S,r,r) ;
end