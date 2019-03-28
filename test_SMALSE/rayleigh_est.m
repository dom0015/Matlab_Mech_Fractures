function rayl_c = rayleigh_est(A, it,epsr)
x = ones(size(A,1),1);
x = x / norm(x);
% the backslash operator in Octave solves a linear system
x_old=x;
rayl_c_old=inf;
for i=1:it
    x=A*x_old;
    rayl_c=dot(x_old,x);
    x_old=x/norm(x);
    if(norm(rayl_c_old-rayl_c)<epsr)
        fprintf('Converged: it %d precision %d\n',i,norm(rayl_c_old-rayl_c));
        break;
    end
    if i==it 
        break
    end
    rayl_c_old=rayl_c;
end
if i==it
    fprintf('Maximum iterations reached: it %d precision %d\n',i,norm(rayl_c_old-rayl_c));
end
end