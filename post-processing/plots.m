%% PLOTS
%
%   Reads the .h5 data file and plots the results
%
% (c) Matteo M. 2022

clear;


%% Settings
plot_profiles = 1;      % Plot profiles along the column
plot_bm_profiles = 1;   % Plot profiles along the biofilm
units = 2;              % 1: ppb, 2: ug/m3
ds = 10;                % Plot every ds samples

%% Open GUI
[FileName,PathName] = uigetfile('../results/*.h5','Select the .h5 data file');
filepath = [PathName FileName];

%% Read data
fid = H5F.open(filepath);
dset_id = H5D.open(fid,'/parameters/use_ts');
use_ts = boolean(H5D.read(dset_id));

dset_id = H5D.open(fid,'/parameters/mixing_model');
mixing_model = boolean(H5D.read(dset_id));

dset_id = H5D.open(fid,'/parameters/mw');
mw = double(H5D.read(dset_id));

dset_id = H5D.open(fid,'/parameters/filename');
yaml = H5D.read(dset_id);

dset_id = H5D.open(fid,'/parameters/export_profiles');
export_profiles = boolean(H5D.read(dset_id));

dset_id = H5D.open(fid,'/parameters/export_bm_profiles');
export_bm_profiles = boolean(H5D.read(dset_id));

% Conversion factor ppb -> ug/m3
if units == 2
    cf = mw/24.45;
else
    cf = 1;
end

if use_ts
    % Inlet
    inlet.data = h5read(filepath,'/data/inlet/c');
    inlet.t = h5read(filepath,'/data/inlet/t');

    % Outlet
    outlet.data = h5read(filepath,'/data/outlet/c');
    outlet.t = h5read(filepath,'/data/outlet/t');
end

% Model
model.c = h5read(filepath,'/results/output/c');
model.t = h5read(filepath,'/results/output/t');

if mixing_model
    model.mix = h5read(filepath,'/results/output/c_mix');
end

%% PLOTS
if use_ts
    % Output concentration
    figure;
    plot(double(inlet.t)/3600,inlet.data*cf,'-r')
    hold all
    plot(double(outlet.t)/3600,outlet.data*cf,'-b')
    plot(double(model.t(1:ds:end))/3600,model.data(1:ds:end)*cf, 'ok')
    xlabel('Time (h)')
    if units == 2
        ylabel('Concentration (\mug/m^3)')
    else
        ylabel('Concentration (ppb)')
    end
    legend('Inlet','Outlet','Model')
end

% Mixing model
if mixing_model
    figure;
    plot(double(model.t)/3600,model.c*cf, 'k')
    hold all
    plot(double(model.t)/3600,model.mix*cf)
    if units == 2
        ylabel('Concentration (\mug/m^3)')
    else
        ylabel('Concentration (ppb)')
    end
    legend('Column outlet','Inlet air - mixed')
end


% Column profiles
if plot_profiles && export_profiles
    profiles.t = h5read(filepath,'/results/profiles/t');
    profiles.z = h5read(filepath,'/results/profiles/z');
    profiles.o = h5read(filepath,'/results/profiles/o');
    profiles.c = h5read(filepath,'/results/profiles/c');
    
    % Oxygen
    figure;
    l = 1;
    for i = 1:length(profiles.t)
        plot(profiles.z,profiles.o(:,i))
        if i == 1
            hold all
        end
        lgnd_o{l} = ['t = ' num2str(profiles.t(i)/3600,'%.2f') ' h'];
        l = l+1;
    end
    xlabel('Distance from inlet (m)')
    ylabel('Oxygen concentration (ppm)')
    legend(lgnd_o)
    
    % VOC
    figure;
    l = 1;
    for i = 1:length(profiles.t)
        plot(profiles.z,profiles.c(:,i)*cf)
        if i == 1
            hold all
        end
        lgnd_i{l} = ['t = ' num2str(profiles.t(i)/3600,'%.2f') ' h'];
        l = l+1;
    end
    xlabel('Distance from inlet (m)')
    if units == 2
        ylabel('VOC concentration (\mug/m^3)')
    else
        ylabel('VOC concentration (ppb)')
    end
    legend(lgnd_i)
end

% Biofilm profiles
if plot_bm_profiles && export_bm_profiles
    profiles.t = h5read(filepath,'/results/profiles/t');
    profiles.x = h5read(filepath,'/results/profiles/x');
    profiles.uo = h5read(filepath,'/results/profiles/uo');
    profiles.ui = h5read(filepath,'/results/profiles/ui');

    % Oxygen
    figure;
    l = 1;
    for i = 1:length(profiles.t)
        plot(profiles.x,profiles.uo(:,i))
        if i == 1
            hold all
        end
        lgnd_o{l} = ['t = ' num2str(profiles.t(i)/3600,'%.2f') ' h'];
        l = l+1;
    end
    xlabel('Distance biofilm surface (\mum)')
    ylabel('Oxygen concentration in the biofilm (ppm)')
    legend(lgnd_o)
    
    % VOC
    figure;
    l = 1;
    for i = 1:length(profiles.t)
        plot(profiles.x,profiles.ui(:,i)*cf)
        if i == 1
            hold all
        end
        lgnd_i{l} = ['t = ' num2str(profiles.t(i)/3600,'%.2f') ' h'];
        l = l+1;
    end
    xlabel('Distance from biofilm surface (\mum)')
    if units == 2
        ylabel('VOC concentration in the biofilm (\mug/m^3)')
    else
        ylabel('VOC concentration in the biofilm (ppb)')
    end
    legend(lgnd_i)
end
