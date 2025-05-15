% bw_r_oneout
% ==========
% calculates optimal bandwidth for isotropic Gaussian smoothing with
% damping by leading r term

convno1 = 1000; % number of random starting points to try

exitflag0 = 0;
while exitflag0 < convno1
    x_01 = 1 + 10*rand(1);  % starting points
    x_02 = rand(1);
    x_0 = [x_01 x_02]; % transform for positivity
    x_0 = log(x_0);
    [x,fval,exitflag] = fminsearch(@oneout_iG_r,[x_0]);
    if exitflag > 0
        exitflag0 = exitflag0 + 1;
        if exitflag0 == 1
            fval0 = fval + 1;
        end
        if fval < fval0
            fval0 = fval;
            xx = x;
        end
    end
end
x = exp(xx);
lam = x(1)
lsq = lam^2;
alp = x(2)
logL = -fval0

% calculate density
for i = 1:fineness_x+1
    for j = 1:fineness_y+1
        X(i,j) = xcalc(i);
        Y(i,j) = ycalc(j);
        xy0 = [xcalc(i) ycalc(j)];
        temp = 0;
        for k = 1:n
            dij = (xy0 - xy(k,:))';
            dist = sqrt(sum(dij.^2));
            lam_k = lam*(volume(k)^alp);
            lsq_k = lam_k^2;
            temp = temp + (1/(sqrt(2*pi)*lsq_k*lam_k*pi*n))*exp(-0.5*( dist^2 /lsq_k))*dist;
        end
        V(i,j) = temp;
    end
end

% ===========================
function f = oneout_iG_r(x)

global xy volume

n = length(xy);

lam = exp(x(1));
alp = exp(x(2));
%lsq = lam^2;
S = 0;
for i = 1:n
    K = 0;
    for j = 1:n
        if j == i
        else
            Xi = xy(i,:);
            Xj = xy(j,:);
            lam_j = lam*(volume(j)^alp);
            lsq_j = lam_j^2;
            K = K + (1/(sqrt(2*pi)*lsq_j*lam_j*pi*(n-1)))*exp(-0.5*( sum((Xi - Xj).^2)) /lsq_j)*sqrt(sum((Xi - Xj).^2));
        end
    end
    S = S + log(K);
end
f = -S;

end

