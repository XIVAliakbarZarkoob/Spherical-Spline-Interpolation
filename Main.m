%% Aliakbar Zarkoob, AKA "XIV"
%  Gmail: XIV.Aliakbar.Zarkoob@gmail.com
%  Telegram: @XIVAliakbar

clc, clear, close all, beep off, format long g
set(0,'defaultTextInterpreter','latex')

%#ok<*MINV>
%#ok<*UNRCH>

%% INPUT **********************************************************************

SAVE_OUTPUTS = false;

KERNEL = questdlg('Select The Desired Kernel Function', 'Kernel Funtion', ...
    'Abel-Poisson', 'Singularity', 'Logarithmic', 'Abel-Poisson');
if isempty(KERNEL)
    return
end
METHOD_LIST = {'Cholesky', 'TSVD', 'GCV', 'L-Curve', 'VCE'};
METHOD = listdlg('PromptString','Select Regularization Method', ...
    'ListString',METHOD_LIST,'SelectionMode','single', ...
    'ListSize', [300,100]);
METHOD = string(METHOD_LIST(METHOD));
if isempty(METHOD)
    return
end
NOISE = questdlg('Add White Noise to Input Data?', 'White Noise', ...
    'Yes', 'No', 'No');
if isempty(NOISE)
    NOISE = 'No';
end

% *****************************************************************************

%% Load Data

data.main = readtable('Global_XGM2019_Potential_6-0Deg.dat', 'NumHeaderLines', 32, 'ReadVariableNames', false, 'FileType', 'text');
data.main.Properties.VariableNames = {'idx', 'lon', 'lat', 'h', 'U'};
data.valid = readtable('Global_XGM2019_Potential_3-0Deg.dat', 'NumHeaderLines', 32, 'ReadVariableNames', false, 'FileType', 'text');
data.valid.Properties.VariableNames = {'idx', 'lon', 'lat', 'h', 'U'};

%% Create Data200

% function U = Potential(lat, lon, GM, R, min_l, max_l, C_lm, S_lm)
% 
%     U = zeros(length(lat), 1);
%     m = kron(ones(max_l+1, 1), (0:max_l));
%     l = m';
%     r = ones(length(lat), 1)*6378136.3;
%     C_lm = C_lm(1:max_l+1, 1:max_l+1);
%     S_lm = S_lm(1:max_l+1, 1:max_l+1);
%     for i = 1:length(lat)
%         P_lm = zeros(max_l+1, max_l+1);
%         for k = 0:max_l
%             P_lm(k+1, 1:k+1) = legendre(k, cosd(lat(i)), 'sch');
%         end
%             R_r_ratio = (R/r(i)).^(l + 1);
%             K = R_r_ratio.*P_lm.*(C_lm.*cosd(m*lon(i))+S_lm.*sind(m*lon(i)));
%             U(i,1) = (GM/R)*sum(sum(K(min_l+1:max_l,:)));
%     end
% 
% end
% 
% CS = gfc_read('WHU-SWPU-GOGR2022S.gfc');
% STEP_INPUT = 6; % degrees
% STEP_INT = 3; % degrees
% min_l = 10;
% max_l = 300;
% 
% lat = (1:STEP_INPUT:179)';
% lon = (1:STEP_INPUT:359)'; 
% [lon, lat] = meshgrid(lon, lat);
% lon = lon(:); lat = lat(:);
% data.main = table(lon, lat);
% data.main.U = Potential(lat, lon, CS.GM, CS.R, min_l, max_l, CS.C_lm, CS.S_lm);
% 
% lat = (1:STEP_INT:179)';
% lon = (1:STEP_INT:359)'; 
% [lon, lat] = meshgrid(lon, lat);
% lon = lon(:); lat = lat(:);
% data.valid = table(lon, lat);
% data.valid.U = Potential(lat, lon, CS.GM, CS.R, min_l, max_l, CS.C_lm, CS.S_lm);

%% Add White Noise

if strcmp(NOISE, 'Yes')
    wn = str2double(cell2mat(inputdlg('Enter The Standard Deviation of White Noise', 'White Noise STD', [1,45])));
    % wn = 400;
    data.main.U = data.main.U + randn(height(data.main), 1)*wn;
end

%% Covert Coordinates From Geodetic To Cartesian

% data.main.X = CS.R.*sind(data.main.lat).*cosd(data.main.lon);
% data.main.Y = CS.R.*sind(data.main.lat).*sind(data.main.lon);
% data.main.Z = CS.R.*cosd(data.main.lat);
% 
% data.valid.X = CS.R.*sind(data.valid.lat).*cosd(data.valid.lon);
% data.valid.Y = CS.R.*sind(data.valid.lat).*sind(data.valid.lon);
% data.valid.Z = CS.R.*cosd(data.valid.lat);

function [X, Y, Z] = geo2cart(lat, lon, h)
    wgs84 = wgs84Ellipsoid;
    a = wgs84.SemimajorAxis;
    b = wgs84.SemiminorAxis;
    e = wgs84.Eccentricity;
    N = a./sqrt(1-e^2*sind(lat).^2);
    X = (N+h).*cosd(lat).*cosd(lon);
    Y = (N+h).*cosd(lat).*sind(lon);
    Z = (b^2*N/a^2+h).*sind(lat);
end

[data.main.X, data.main.Y, data.main.Z] = geo2cart(data.main.lat, data.main.lon, data.main.h);
[data.valid.X, data.valid.Y, data.valid.Z] = geo2cart(data.valid.lat, data.valid.lon, data.valid.h);

%% Normalizing

function [Xn, Yn, Zn] = normalized(X, Y, Z)
    NORM = sqrt(X.^2 + Y.^2 + Z.^2);
    Xn = X./NORM; Yn = Y./NORM; Zn = Z./NORM;
end

[data.main.X, data.main.Y, data.main.Z] = normalized(data.main.X, data.main.Y, data.main.Z);
[data.valid.X, data.valid.Y, data.valid.Z] = normalized(data.valid.X, data.valid.Y, data.valid.Z);

%% Kernel Functions & Compution of Design Matrices

function [A, A_int] = create_design(data_main, data_valid, h, kernel)

    n = height(data_main);
    u = height(data_valid);

    XYZ1 = table2array(data_main(:,{'X','Y','Z'}));
    inner_product = sum(kron(XYZ1, ones(n,1)) .* kron(ones(n,1), XYZ1), 2);
    inner_product = reshape(inner_product, n, n);
    
    XYZ2 = table2array(data_valid(:,{'X','Y','Z'}));
    inner_product_int = sum(kron(XYZ1, ones(u,1)) .* kron(ones(n,1), XYZ2), 2);
    inner_product_int = reshape(inner_product_int, u, n);
    
    if strcmp(kernel, 'Abel-Poisson')
        L_h = 1 + h^2 - 2*h*inner_product;
        A = 0.25/pi * (1 - h^2)./(L_h).^(3/2);
        L_h = 1 + h^2 - 2*h*inner_product_int;
        A_int = 0.25/pi * (1 - h^2)./(L_h).^(3/2);
    elseif strcmp(kernel, 'Singularity')
        L_h = 1 + h^2 - 2*h*inner_product;
        A = 0.5/pi./L_h.^0.5;
        L_h = 1 + h^2 - 2*h*inner_product_int;
        A_int = 0.5/pi./L_h.^0.5;
    elseif strcmp(kernel, 'Logarithmic')
        L_h = 1 + h^2 - 2*h*inner_product; 
        A = 0.5/pi/h * log(1 + 2*h./(L_h.^0.5 + 1 - h));
        L_h = 1 + h^2 - 2*h*inner_product_int; 
        A_int = 0.5/pi/h * log(1 + 2*h./(L_h.^0.5 + 1 - h));
    else
        error('The Input Kernel Is Not Valid!')
    end

end

n = height(data.main);
u = height(data.valid);

%% "Signal to Noise Ratio (SNR)" Method to find the best "h"

function [h_values, SNR_values] = SNR(data, method, kernel)

    function ratio = SNR_cal(S, S_hat)
        ratio = 10*log10(sum(S.^2,'all')/sum((S-S_hat).^2,'all'));
    end

    h_values = (0.01:0.01:0.99)';
    n = height(data.main);
    SNR_values = zeros(length(h_values), 1);
    for j = 1:length(h_values)
        [A, A_int] = create_design(data.main, data.valid, h_values(j), kernel);
        if strcmp(method, 'Cholesky')
            coeff_test = lscov(A, data.main.U, eye(n), 'chol');
            SNR_values(j) = SNR_cal(data.main.U, A*coeff_test);
            % SNR_values(j) = SNR_cal(data.valid.U, A_int*coeff_test);
        elseif strcmp(method, 'TSVD')
            k = 5:5:n;
            [U, S, V] = svd(A);
            [coeff_test, ~, ~] = tsvd(U, diag(S), V, data.main.U, k);
            int_test = A_int*coeff_test;
            diff_test = data.valid.U - int_test;
            [~, idx] = min(vecnorm(diff_test));
            SNR_values(j) = SNR_cal(data.main.U, A*coeff_test(:,idx));
            % SNR_values(j) = SNR_cal(data.valid.U, A_int*coeff_test(:,idx));
        elseif strcmp(method, 'VCE')
            [coeff_test, ~, ~, ~, ~, ~] = ...
            vce_iter_opt(A, data.main.U, zeros(n,1), eye(n), eye(n), 1, 5);
            coeff_test = coeff_test(:,end);
            SNR_values(j) = SNR_cal(data.main.U, A*coeff_test);
            % SNR_values(j) = SNR_cal(data.valid.U, A_int*coeff_test);
        elseif strcmp(method, 'GCV')
            [Uq, Sq, ~, ~] = cgsvd(A, eye(n));
            [reg_min, ~, ~] = gcv(Uq(:,1:n), Sq, data.main.U, 'Tikh');
            coeff_test = inv(A'*A + reg_min*eye(n))*A'*data.main.U;
            SNR_values(j) = SNR_cal(data.main.U, A*coeff_test);
            % SNR_values(j) = SNR_cal(data.valid.U, A_int*coeff_test);
        else
            error('The Input Method Is Not Valid!')
        end
    end

end

% [h_values, SNR_values] = SNR(data, METHOD, KERNEL);
% idx = max(SNR_values) == SNR_values;
% h = h_values(idx); h = h(1);
% [A, A_int] = create_design(data.main, data.valid, h, KERNEL);
% 
% figure()
% hold on, grid on, box on
% plot(h_values, SNR_values, '.-k', 'LineWidth', 1.5, 'MarkerSize', 15)
% plot(h, SNR_values(idx), '.r', 'MarkerSize', 35)
% xlabel('$h$', 'FontSize', 15)
% ylabel('$SNR$', 'FontSize', 15)

%% "Cross-Validation" Method to find the best "h"

function cv = CV_h(data, method, kernel, cv_percent, h_decimals)

    for i = 1:h_decimals

        if i == 1
            h_values = (10^-1 : 10^-1 : 9*10^-1)';
        else
            h_values = (best_h - 9*10^-i : 10^-i : best_h + 9*10^-i)';
        end

        n = height(data.main);
        cv_num = ceil(n*cv_percent);
        n = n - cv_num;
        cv_idx = sort(randperm(n, cv_num))';
        data_main = data.main; data_main(cv_idx,:) = [];
        data_cv = data.main(cv_idx, :);
        cv.(sprintf('dec%g', i)) = table(h_values,zeros(length(h_values),1));
        cv.(sprintf('dec%g', i)).Properties.VariableNames = {'h','norm'};

        for j = 1:length(h_values)
            [A, A_int] = create_design(data_main, data_cv, h_values(j), kernel);
            if strcmp(method, 'Cholesky')
                coeff_test = lscov(A, data_main.U, eye(n), 'chol');
                cv.(sprintf('dec%g', i)).norm(j) = norm(data_cv.U - A_int*coeff_test);
            elseif strcmp(method, 'TSVD')
                k = 5:5:n;
                [U, S, V] = svd(A);
                [coeff_test, ~, ~] = tsvd(U, diag(S), V, data_main.U, k);
                int_test = A_int*coeff_test;
                diff_test = data_cv.U - int_test;
                [~, idx] = min(vecnorm(diff_test));
                cv.(sprintf('dec%g', i)).norm(j) = norm(data_cv.U - A_int*coeff_test(:,idx));
            elseif strcmp(method, 'VCE')
                [coeff_test, ~, ~, ~, ~, ~] = ...
                vce_iter_opt(A, data_main.U, zeros(n,1), eye(n), eye(n), 1, 5);
                coeff_test = coeff_test(:,end);
                cv.(sprintf('dec%g', i)).norm(j) = norm(data_cv.U - A_int*coeff_test);
            elseif strcmp(method, 'GCV')
                [Uq, Sq, ~, ~] = cgsvd(A, eye(n));
                [reg_min, ~, ~] = gcv(Uq(:,1:n), Sq, data_main.U, 'Tikh');
                coeff_test = inv(A'*A + reg_min*eye(n))*A'*data_main.U;
                cv.(sprintf('dec%g', i)).norm(j) = norm(data_cv.U - A_int*coeff_test);
            elseif strcmp(method, 'L-Curve')
                [Uq, Sq, ~, ~] = cgsvd(A, eye(n));
                [reg_corner, ~, ~, ~] = l_curve(Uq(:,1:n), Sq, data_main.U, 'Tikh');
                coeff_test = inv(A'*A + reg_corner*eye(n))*A'*data_main.U; 
                cv.(sprintf('dec%g', i)).norm(j) = norm(data_cv.U - A_int*coeff_test);
            end
        end
        best_norm = min(cv.(sprintf('dec%g', i)).norm);
        best_h = best_norm == cv.(sprintf('dec%g', i)).norm;
        best_h = cv.(sprintf('dec%g', i)).h(best_h); best_h = best_h(1);
    end
    cv.best_h = best_h;
    cv.best_norm = best_norm;
end

max_dec = 3;
cv = CV_h(data, METHOD, KERNEL, 0.05, max_dec);
h = cv.best_h;
[A, A_int] = create_design(data.main, data.valid, h, KERNEL);

figure()
for i = 1:max_dec
    subplot(max_dec,1,i)
    hold on, grid on, box on
    plot(cv.(sprintf('dec%g', i)).h, cv.(sprintf('dec%g', i)).norm, '.-k', 'LineWidth', 1.5, 'MarkerSize', 15)
    idx = min(cv.(sprintf('dec%g', i)).norm) == cv.(sprintf('dec%g', i)).norm;
    plot(cv.(sprintf('dec%g', i)).h(idx), cv.(sprintf('dec%g', i)).norm(idx), '.r', 'MarkerSize', 35)
    xlabel('$h$', 'FontSize', 15)
    ylabel('$||S-\hat{S}||$', 'FontSize', 15)
end

%% Chech For Ill-Conditionaliy

figure()
imagesc(A)
title('Design Matrix')
axis equal tight
colormap('turbo')
colorbar

[U, S, V] = svd(A);
figure()
PicardPlot(U, diag(S), data.main.U)

%% Compute Coefficients and Interpolated Values

coeff.A = inv(A)*data.main.U;
int.A = A_int*coeff.A;
diff.A = data.valid.U - int.A;

coeff.pinv = pinv(A)*data.main.U;
int.pinv = A_int*coeff.pinv;
diff.pinv = data.valid.U - int.pinv;

if strcmp(METHOD, 'Cholesky')
    coeff.chol = lscov(A, data.main.U, eye(n), 'chol');
    int.chol = A_int*coeff.chol;
    diff.chol = data.valid.U - int.chol;
elseif strcmp(METHOD, 'TSVD')
    k = 5:5:n;
    [coeff.tsvd, ~, ~] = tsvd(U, diag(S), V, data.main.U, k);
    int.tsvd = A_int*coeff.tsvd;
    diff.tsvd = data.valid.U - int.tsvd;
    [~, idx] = min(vecnorm(diff.tsvd));
    coeff.tsvd = coeff.tsvd(:, idx);
    int.tsvd = int.tsvd(:, idx);
    diff.tsvd = diff.tsvd(:, idx);
    k = k(idx);
    clear idx
elseif strcmp(METHOD, 'VCE')
    [coeff.vce, ~, ~, sigma2_y, sigma2_x, alpha_vce] = ...
        vce_iter_opt(A, data.main.U, zeros(n,1), eye(n), eye(n), 1, 5);
    coeff.vce = coeff.vce(:,end);
    int.vce = A_int*coeff.vce;
    diff.vce = data.valid.U - int.vce;
elseif strcmp(METHOD, 'GCV')
    [Uq, Sq, ~, ~] = cgsvd(A, eye(n));
    figure()
    [reg_min, ~, ~] = gcv(Uq(:,1:n), Sq, data.main.U, 'Tikh');
    coeff.tikh_gcv = inv(A'*A + reg_min*eye(n))*A'*data.main.U; 
    int.tikh_gcv = A_int*coeff.tikh_gcv;
    diff.tikh_gcv = data.valid.U - int.tikh_gcv;
elseif strcmp(METHOD, 'L-Curve')
    [Uq, Sq, ~, ~] = cgsvd(A, eye(n));
    [reg_corner, ~, ~, ~] = l_curve(Uq(:,1:n), Sq, data.main.U, 'Tikh');
    coeff.tikh_lcurve = inv(A'*A + reg_corner*eye(n))*A'*data.main.U; 
    int.tikh_lcurve = A_int*coeff.tikh_lcurve;
    diff.tikh_lcurve = data.valid.U - int.tikh_lcurve;
end


%% Plot Main Data

space = 0;
worldMap = readgeotable('landareas.shp');
figure()
subplot(1, 2, 1)
hold on
scatter(data.main.lon, data.main.lat, 20, data.main.U, 'filled')
geoshow(worldMap, 'FaceColor', 'none');
axis equal 
xlim([min(data.main.lon) max(data.main.lon)]+space*[-1 1])
ylim([min(data.main.lat) max(data.main.lat)]+space*[-1 1])
colormap('turbo')
hhh = colorbar; set(get(hhh,'ylabel'),'String','$\frac{m^2}{s^2}$','FontSize',16,'Interpreter','latex');
xlabel('Longitude')
ylabel('Latitude')
title('Grid 06min')
box on
grid on
set(gca, 'Layer', 'top')
subplot(1, 2, 2)
hold on
scatter(data.valid.lon, data.valid.lat, 20, data.valid.U, 'filled')
geoshow(worldMap, 'FaceColor', 'none');
axis equal 
xlim([min(data.main.lon) max(data.main.lon)]+space*[-1 1])
ylim([min(data.main.lat) max(data.main.lat)]+space*[-1 1])
colormap('turbo')
hhh = colorbar; set(get(hhh,'ylabel'),'String','$\frac{m^2}{s^2}$','FontSize',16,'Interpreter','latex');
xlabel('Longitude')
ylabel('Latitude')
title('Grid 03min')
box on
grid on
set(gca, 'Layer', 'top')
sgtitle('Gravity Potential Data From XGM2019 Model')

clear space hhh

%% Plot Results 

function plot_result(lon, lat, interpolation, difference, marker_size, map_space, method)
    worldMap = readgeotable('landareas.shp');
    figure()
    subplot(1, 2, 1)
    hold on
    scatter(lon, lat, marker_size, interpolation, 'filled')
    geoshow(worldMap, 'FaceColor', 'none');
    axis equal 
    xlim([min(lon) max(lon)]+map_space*[-1 1])
    ylim([min(lat) max(lat)]+map_space*[-1 1])
    colormap('turbo')
    hhh = colorbar; set(get(hhh,'ylabel'),'String','$\frac{m^2}{s^2}$','FontSize',16,'Interpreter','latex');
    xlabel('Longitude')
    ylabel('Latitude')
    title('Interpolated Values')
    box on
    grid on
    set(gca, 'Layer', 'top')
    subplot(1, 2, 2)
    hold on
    scatter(lon, lat, marker_size, difference, 'filled')
    geoshow(worldMap, 'FaceColor', 'none');
    axis equal 
    xlim([min(lon) max(lon)]+map_space*[-1 1])
    ylim([min(lat) max(lat)]+map_space*[-1 1])
    colormap('turbo')
    hhh = colorbar; set(get(hhh,'ylabel'),'String','$\frac{m^2}{s^2}$','FontSize',16,'Interpreter','latex');
    xlabel('Longitude')
    ylabel('Latitude')
    title('Difference From True Values')
    box on
    grid on
    set(gca, 'Layer', 'top')
    sgtitle(method)
    fprintf('\nNorm of Error Vector for %s Method:\n%f\n', method, round(norm(difference),4))
    fprintf('\nMean of Error Vector for %s Method:\n%f\n', method, round(mean(difference),4))
end

space = 0;
% plot_result(data.valid.lon, data.valid.lat, int.A, diff.A, 20, space, '$A^{-1}l$')
% plot_result(data.valid.lon, data.valid.lat, int.pinv, diff.pinv, 20, space, 'PINV')
if strcmp(METHOD, 'Cholesky')
    plot_result(data.valid.lon, data.valid.lat, int.chol, diff.chol, 20, space, 'Cholesky')
elseif strcmp(METHOD, 'TSVD')
        plot_result(data.valid.lon, data.valid.lat, int.tsvd, diff.tsvd, 20, space, 'TSVD')
elseif strcmp(METHOD, 'VCE')
    plot_result(data.valid.lon, data.valid.lat, int.vce, diff.vce, 20, space, 'VCE')
elseif strcmp(METHOD, 'GCV')
    plot_result(data.valid.lon, data.valid.lat, int.tikh_gcv, diff.tikh_gcv, 20, space, 'Tikhonov (GCV)')
elseif strcmp(METHOD, 'L-Curve')
    plot_result(data.valid.lon, data.valid.lat, int.tikh_lcurve, diff.tikh_lcurve, 20, space, 'Tikhonov (L-Curve)')
end

%% Export for Plots in Python

if SAVE_OUTPUTS

    main_data = [data.main.lon data.main.lat data.main.U]; 
    valid_data = [data.valid.lon data.valid.lat data.valid.U];
    
    h_dec1 = table2array(cv.dec1);
    h_dec2 = table2array(cv.dec2);
    h_dec3 = table2array(cv.dec3);

    if strcmp(NOISE, 'Yes')
        noise_str = "wn_";
    else
        noise_str = "";
    end
    
    if strcmp(METHOD, 'Cholesky')
        int_chol = int.chol;
        diff_chol = diff.chol;
        save("Output_"+noise_str+KERNEL+"_"+METHOD+".mat", 'main_data', 'valid_data', 'int_chol', 'diff_chol', ...
            'h_dec1', 'h_dec2', 'h_dec3', 'h')
    elseif strcmp(METHOD, 'TSVD')
        int_tsvd = int.tsvd;
        diff_tsvd = diff.tsvd;
        save("Output_"+noise_str+KERNEL+"_"+METHOD+".mat", 'main_data', 'valid_data', 'int_tsvd', 'diff_tsvd', ...
            'h_dec1', 'h_dec2', 'h_dec3','k', 'h')
    elseif strcmp(METHOD, 'VCE')
        int_vce = int.vce;
        diff_vce = diff.vce;
        lambda = alpha_vce;
        save("Output_"+noise_str+KERNEL+"_"+METHOD+".mat", 'main_data', 'valid_data', 'int_vce', 'diff_vce', ...
            'h_dec1', 'h_dec2', 'h_dec3', 'lambda', 'h')
    elseif strcmp(METHOD, 'GCV')
        int_gcv = int.tikh_gcv;
        diff_gcv = diff.tikh_gcv;
        lambda = reg_min;
        save("Output_"+noise_str+KERNEL+"_"+METHOD+".mat", 'main_data', 'valid_data', 'int_gcv', 'diff_gcv', ...
            'h_dec1', 'h_dec2', 'h_dec3', 'lambda', 'h')
    elseif strcmp(METHOD, 'L-Curve')
        int_lcurve = int.tikh_lcurve;
        diff_lcurve = diff.tikh_lcurve;
        lambda = reg_corner;
        save("Output_"+noise_str+KERNEL+"_"+METHOD+".mat", 'main_data', 'valid_data', 'int_lcurve', 'diff_lcurve', ...
            'h_dec1', 'h_dec2', 'h_dec3', 'lambda')
    end

end

int_A = int.A;
diff_A = diff.A;
main_data = [data.main.lon data.main.lat data.main.U]; 
valid_data = [data.valid.lon data.valid.lat data.valid.U];
save('NoRegResult.mat', 'main_data', 'valid_data', 'int_A', 'diff_A')

%% Finish Sound

% fs = 44100;                    
% t = 0:1/fs:1;                   
% frequency = 440;                 
% signal = sin(2*pi*frequency*t);
% sound(signal, fs);