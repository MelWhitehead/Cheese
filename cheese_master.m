% cheese_master.m
% ================
%
% fits a (an)istropic Gaussian kernel damped by a circular magma lens from other volumes

global xy volume rangitotox rangitotoy coast finenessx finenessy Xplot Yplot distances vol_rad
% xy = volcano locations, ranigtoto = xy reference, distances = matrix of
% interpoint distances, vol_rad = volume->radius relation, rest are plotting
% data

% FLAGS
omitlast = 0; % if true, leave out Rangitoto (to visualize how it might have been forecast)
Ktype = 1; % 1 = standard Gaussian, 2 = disk, 3 = weight, 4 = r, all isotropic

% PLOTTING SETUP
load 'coast.dat';
xlim0 = 0; % map limits
xlim1 = 25;
ylim0 = 0;
ylim1 = 35;
fineness_y = 350; % 100m grid spacing
fineness_x = 250;
temp_x = (0:fineness_x);
xcalc = xlim0 + (xlim1 - xlim0)*temp_x/fineness_x;
temp_y = (0:fineness_y);
ycalc = ylim0 + (ylim1 - ylim0)*temp_y/fineness_y;
X = zeros(fineness_x+1,fineness_y+1);
Y = zeros(fineness_x+1,fineness_y+1);
V = zeros(fineness_x+1,fineness_y+1);

% =========================================================================
% DATA
load 'coords_vol2.txt';  % location data etc.
coords_vol = coords_vol2;
xy_all = coords_vol(:,2:3);
rangitotox = coords_vol(49,2);
rangitotoy = coords_vol(49,3);

if omitlast > 0
    temp = [coords_vol(1:48,:); coords_vol(50:51,:)];
else
    temp = coords_vol;
end

id = temp(:,1);
xcoord = temp(:,2);
ycoord = temp(:,3);
volume1 = temp(:,4); % magmatic, in 10^6 m^3
volume2 = temp(:,5); % tephra, in 10^6 m^3
volume = volume1 + volume2;
volume = volume/1000; % total, in km^3
vol_rad = sqrt(volume/pi); %

n = length(id);
xy = [xcoord ycoord];

% =========================================================================
% intervent distances
distances = zeros(n,n);
for k = 1:n
    distances(k,:) = sqrt((xy(k,1) - xy(:,1)).^2 + (xy(k,2) - xy(:,2)).^2);
end

% =========================================================================
% kernel
if Ktype == 1 % no damping
    bw_iG_oneout; % fit isotropic Gaussian kernel bandwidth by leave-one out cross-validation
end
if Ktype == 2
    bw_disk_oneout;
end
if Ktype == 3
    bw_weight_oneout;
end
if Ktype == 4
    bw_r_oneout;
end

% ====================================
% PLOT RESULTS
% convert to lat long for plotting
% Rangitoto location
rx = xy(49,1);
ry = xy(49,2);
rx2 = 174.86;
ry2 = -36.79;
% centre locations
xk2 = rx2 + (xy(:,1) - rx)/90.87;
yk2 = ry2 + (xy(:,2) - ry)/111.12;
% map limits
xlim02 = 174.7; %rx2 + (xlim0 - rx)/90.87;
ylim02 = -37.051; %ry2 + (ylim0 - ry)/111.12;
xlim12 = 174.94; %rx2 + (xlim1 - rx)/90.87;
ylim12 = -36.75; %ry2 + (ylim1 - ry)/111.12;
% reunit X and Y
X2 = rx2 + (X - rx)/90.87;
Y2 = ry2 + (Y - ry)/111.12;

symbsize = 2+25*sqrt(volume);  % plot larger volumes using larger symbols

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(Ktype)
contour(X2,Y2,V,100);
hold on;
if Ktype == 1
      plot(xk2,yk2,'^k');
else
    for i = 1:length(xk2)
        plot(xk2(i),yk2(i),'^k','MarkerSize',symbsize(i)); hold on;
    end
end
hold on;
plot(coast(:,1),coast(:,2),'-k');
hold off;
xlabel('Longitude');
ylabel('Latitude');
xlim([xlim02 xlim12]);
ylim([ylim02 ylim12]);
