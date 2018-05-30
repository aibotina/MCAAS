function [W, Z] = gradF_t(X,Y,S,M_E,E,m0,lam1,lam2,SubSeqSM,lam3,SubSeqSMsera)
[n, r] = size(X);
[m, r] = size(Y);
XS = X*S ;
YS = Y*S' ;
XSY = XS*Y' ;
% load SeqSM;
B = GetseqM(SubSeqSM);
C = GetseqM(SubSeqSMsera);
G3x = B*XSY*YS;
G3y = YS*X'*B*XS;
Qx = X'* ( (M_E - XSY).*E )*YS /n;
Qy = Y'* ( (M_E - XSY).*E )'*XS /m;
W = ( (XSY - M_E).*E )*YS + X*Qx + lam1*Gp(X,m0,r)+4*lam2*G3x+4*lam3*X*S*Y'*C*Y*S';
Z = ( (XSY - M_E).*E )'*XS + Y*Qy + lam1*Gp(Y,m0,r)+4*lam2*G3y+4*lam3*C*Y*S'*X'*X*S;
end

function out = Gp(X,m0,r)
z = sum(X.^2,2) /(2*m0*r) ;
z = 2*exp( (z-1).^2 ).*(z-1) ;
z( find(z<0) ) = 0;
out = X.*repmat(z,1,r) / (m0*r) ;
end

