function out = F_t(X,Y,S,M_E,E,m0,lam1, lam2,SubSeqSM,lam3,SubSeqSMsera)
[n,r] = size(X) ;
[m,r] = size(Y) ;
% SeqSM is the HA pro sequence similary matrix
% load SeqSM;
s = SubSeqSM;
a = X*S*Y';
[m,n] = size(s);
Sum = 0;
for i = 2:m
     for j = 1:i-1
         Sum = Sum+s(i,j)*norm((a(i,:)-a(j,:)),2).^2;
     end
end

% SeqSMsera is the HA of antisera pro similary matrix 
ss=SubSeqSMsera;
[mm,nn] = size(ss);
Sumsera=0;
for i = 2:mm
     for j = 1:i-1
         Sumsera = Sumsera+ss(i,j)*norm((a(:,i)-a(:,j)),2).^2;
     end
end
out5=lam3*Sumsera;
out4 =lam2*Sum;
out1 = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;
out2 =  lam1*G(Y,m0,r) ;
out3 =  lam1*G(X,m0,r) ;
out = out5+out1+out2+out3+out4 + lam1*m*n*norm(S, 'fro')^2;
end
function out = G(X,m0,r)
z = sum(X.^2,2)/(2*m0*r) ;
y = exp( (z-1).^2 ) - 1 ;
y( find(z < 1) ) = 0 ;
out = sum(y) ;
end