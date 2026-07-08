function R = rate_cal(para,A,W)

Aeq = squeeze(A).';
R = log(det( eye(para.K) + (para.Pt/(para.M_t*para.noise)) * ...
            (Aeq*W) * (Aeq*W)' ));
end
