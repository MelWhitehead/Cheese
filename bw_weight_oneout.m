% bw_weight_oneout.m
% ==================
% calculates optimal bandwidth for isotropic Gaussian smoothing with
% damping by erupted volume weight

convno1 = 10; % number of random starting points to try

% ================================
function f = oneout_ig_weight(x)

global xy volume

n = length(xy);

lam = exp(x(1));
alp = exp(x(2));
lsq = lam^2;
S = 0;
%diagnost = 1:0;
for i = 1:n
    K = 0;
    normpdf = 0;
    for j = 1:n
        if j == i
        else
            Xi = xy(i,:);
            Xj = xy(j,:);
            K = K + (1/(2*lsq*pi*(n-1)))*exp(-0.5*( sum((Xi - Xj).^2)) /lsq)*(volume(j)^(-alp));
            normpdf = normpdf + volume(j)^(-alp);
        end
    end
    %diagnost = [diagnost; [K normpdf]];
    K = (n-1)*K/normpdf;
    S = S + log(K);
end
f = -S;

end



exitflag0 = 0;
while exitflag0 < convno1
    x_01 = 1 + 10*rand(1); % starting points
    x_02 = rand(1);
    x_0 = [x_01 x_02];
    x_0 = log(x_0);  % transform for positivity
    [x,fval,exitflag] = fminsearch(@oneout_ig_weight,[x_0]);
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
normpdf = sum(volume.^(-alp));
for i = 1:fineness_x+1
    for j = 1:fineness_y+1
        X(i,j) = xcalc(i);
        Y(i,j) = ycalc(j);
        xy0 = [xcalc(i) ycalc(j)];
        temp = 0;
        for k = 1:n
            dij = (xy0 - xy(k,:))';
            dist = sqrt(sum(dij.^2));
            temp = temp + (1/(2*lsq*pi*n))*exp(-0.5*( dist^2 /lsq))*(volume(k)^(-alp));
        end
        V(i,j) = temp;
    end
end
V = n*V/normpdf;


