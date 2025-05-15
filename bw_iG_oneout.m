% bw_iG_oneout.m
% ==========
% calculates optimal bandwidth for isotropic Gaussian kernel using
% leave-one-out cross validation

convno1 = 1000;  % number of random starting points to try

exitflag0 = 0;
while exitflag0 < convno1
    x_0 = 1 + 10*rand(1);  % starting point
    x_0 = log(x_0);  % transform for positivity
    [x,fval,exitflag] = fminsearch(@oneout_iG,[x_0]);
    if exitflag > 0
        exitflag0 = exitflag0 + 1;
        if exitflag0 == 1
            fval0 = fval + 1;
        end
        if fval < fval0
            fval0 = fval;
            xx = x;   % update best solution
        end
    end
end
lam = exp(xx)  % backtransform to bandwidth
logL = -fval0  % log-Likelihood
lsq = lam^2;  # mel faffing around

% calculate density
for i = 1:fineness_x+1
    for j = 1:fineness_y+1
        X(i,j) = xcalc(i);
        Y(i,j) = ycalc(j);
        xy0 = [xcalc(i) ycalc(j)];
        temp0 = 0;
        for k = 1:n
            dij = (xy0 - xy(k,:))';
            dist = sqrt(sum(dij.^2));
            temp0 = temp0 + (1/(2*lsq*pi*n))*exp(-0.5*( dist^2 /lsq));
        end
        V(i,j) = temp0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = oneout_iG(x)

global xy

n = length(xy);

lam = exp(x);
lsq = lam^2;
S = 0;
for i = 1:n
    K = 0;
    for j = 1:n
        if j == i
        else
            Xi = xy(i,:);
            Xj = xy(j,:);
            K = K + (1/(2*lsq*pi*(n-1)))*exp(-0.5*( sum((Xi - Xj).^2)) /lsq);
        end
    end
    S = S + log(K);
end
f = -S;

end




