function [CROSSED_L,CROSSED_R,changed] = combine_lambda_overlap(CROSSED_L,CROSSED_R,CROSSED_L_old,CROSSED_R_old,lambda)
%COMBINE_LAMBDA_OVERLAP Summary of this function goes here
%   Detailed explanation goes here

global switching

n=length(CROSSED_L);
changed=0;

% do not switch on new nodes if some lambda is <0
% if sum(lambda<0)>0 % some lambda is <0, switch off that nodes
if switching==1
    idx=1;
    changed=1;
    for f=1:n
        crossed_L_old=CROSSED_L_old{f};
        crossed_R_old=CROSSED_R_old{f};

        crossed_L_old(crossed_L_old)=(lambda(idx:idx+sum(crossed_L_old)-1)>0);
        idx=idx+sum(crossed_L_old);
        
        crossed_R_old(crossed_R_old)=(lambda(idx:idx+sum(crossed_R_old)-1)>0);
        idx=idx+sum(crossed_R_old);

        CROSSED_L{f}=crossed_L_old;
        CROSSED_R{f}=crossed_R_old;
    end
else % no lambda is <0, switch on new contact nodes
    for f=1:n
        crossed_L_old=CROSSED_L_old{f};
        crossed_R_old=CROSSED_R_old{f};
        crossed_L=CROSSED_L{f};
        crossed_R=CROSSED_R{f};

        crossed_L(crossed_L_old)=1;
        if sum(crossed_L~=crossed_L_old)>0 
            changed=1;
        end

        crossed_R(crossed_R_old)=1;
        if sum(crossed_R~=crossed_R_old)>0 
            changed=1;
        end

        CROSSED_L{f}=crossed_L;
        CROSSED_R{f}=crossed_R;
    end
end

