function out = getoptT(X,W,Y,Z,S,M_E,E,m0,lam1,lam2,SubSeqSM,lam3,SubSeqSMsera)
norm2WZ = norm(W,'fro')^2 + norm(Z,'fro')^2;
%f(1) = F_t(X, Y,S,M_E,E,m0,lam1,lam2,SubSeqSM) ;
f(1) = F_t(X, Y,S,M_E,E,m0,lam1,lam2,SubSeqSM,lam3,SubSeqSMsera) ;

t = -1e-1 ;
for i = 1:20
     %   f(i+1) = F_t(X+t*W,Y+t*Z,S,M_E,E,m0,lam1,lam2,SubSeqSM) ;
         f(i+1) = F_t(X+t*W,Y+t*Z,S,M_E,E,m0,lam1,lam2,SubSeqSM,lam3,SubSeqSMsera) ;

        if( f(i+1) - f(1) <= .5*(t)*norm2WZ )
            out = t ;
            return;
        end
        t = t/2 ;
end
out = t ;
end