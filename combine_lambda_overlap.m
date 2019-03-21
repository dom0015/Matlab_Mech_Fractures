function [CROSSED_L,CROSSED_R,changed] = combine_lambda_overlap(CROSSED_L,CROSSED_R,CROSSED_L_old,CROSSED_R_old,lambda)
%COMBINE_LAMBDA_OVERLAP Summary of this function goes here
%   Detailed explanation goes here
n=length(CROSSED_L);
idx=1;
changed=0;
for f=1:n
    crossed_L_old=CROSSED_L_old{f};
    crossed_R_old=CROSSED_R_old{f};
    crossed_L=CROSSED_L{f};
    crossed_R=CROSSED_R{f};
    
    crossed_L(crossed_L_old)=(lambda(idx:idx+sum(crossed_L_old)-1)>0);
    idx=idx+sum(crossed_L_old);
    if sum(crossed_L~=crossed_L_old)>0 
        changed=1;
    end
    
    crossed_R(crossed_R_old)=(lambda(idx:idx+sum(crossed_R_old)-1)>0);
    idx=idx+sum(crossed_R_old);
    if sum(crossed_R~=crossed_R_old)>0 
        changed=1;
    end
    
    CROSSED_L{f}=crossed_L;
    CROSSED_R{f}=crossed_R;
end

