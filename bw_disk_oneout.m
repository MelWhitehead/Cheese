% bw_disk_oneout.m
% ==========
% calculates optimal bandwidth for isotropic Gaussian smoothing with
% damping by erupted volume (disk) using leave-one-out cross validation

convno1 = 1000; % number of random starting points to try

exitflag0 = 0;
while exitflag0 < convno1
    x_01 = 1 + 10*rand(1); % starting points
    x_02 = rand(1);
    x_03 = -log(rand(1))/10;
    x_0 = [x_01 x_02 x_03];
    x_0 = log(x_0);  % transform for positivity
    [x,fval,exitflag] = fminsearch(@oneout_ig_disk,[x_0]);
    if exitflag > 0
        exitflag0 = exitflag0 + 1;
        if exitflag0 == 1
            fval0 = fval + 1;
        end
        if fval < fval0
            fval0 = fval;
            xx = x;
            %exp(x)
        end
    end
end
x = exp(xx);
lam = x(1)
lsq = lam^2;
alp = x(2)
vol_h = sqrt(x(3))
rad_temp = vol_rad/vol_h;
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
            Kknorm = exp(-0.5*(rad_temp(k)^2)/lsq) + (1 - exp(-0.5*(rad_temp(k)^2)/lsq))*exp(-alp);
            temp = temp + (1/(2*lsq*pi*n))*exp(-0.5*( dist^2 /lsq))*exp(-alp*(dist < rad_temp(k)))/Kknorm;
        end
        V(i,j) = temp;
    end
end

% ===================================================
function f = oneout_ig_disk(x)

global xy distances vol_rad

n = length(xy);

lam = exp(x(1));
alp = exp(x(2));
vol_h = sqrt(exp(x(3)));
rad_temp = vol_rad/vol_h;
lsq = lam^2;
S = 0;
for i = 1:n
    K = 0;
    for j = 1:n
        if j == i
        else
            Xi = xy(i,:);
            Xj = xy(j,:);
            Kjnorm = exp(-0.5*(rad_temp(j)^2)/lsq) + (1 - exp(-0.5*(rad_temp(j)^2)/lsq))*exp(-alp);
            K = K + (1/(2*lsq*pi*(n-1)))*exp(-0.5*( sum((Xi - Xj).^2)) /lsq)*exp(-alp*(distances(i,j) < rad_temp(j)))/Kjnorm;
        end
    end
    S = S + log(K);
end
f = -S;

end
