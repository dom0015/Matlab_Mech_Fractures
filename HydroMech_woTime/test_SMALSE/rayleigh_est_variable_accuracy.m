function rayl_c = rayleigh_est_variable_accuracy(A, it,epsr,vypis)
x = ones(size(A,1),1);
x = x / norm(x);
% the backslash operator in Octave solves a linear system
x_old=x;
rayl_c_old=inf;
for i=1:it
    x=A*x_old;
    rayl_c=dot(x_old,x);
    x_old=x/norm(x);
    if(norm(rayl_c_old-rayl_c)/rayl_c<epsr)
        if nargin<4
            fprintf('Converged: it %d precision %d, value = %d\n',i,norm(rayl_c_old-rayl_c),rayl_c);
        end
        break;
    end
    if i==it
        break
    end
    rayl_c_old=rayl_c;
end
if i==it
    if nargin<4
        fprintf('Maximum iterations reached: it %d precision %d, value = %d\n',i,norm(rayl_c_old-rayl_c),rayl_c);
    end
end
end