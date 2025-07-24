% Start with a clean environment
% clear all; close all; clc;

addpath('/home/bqrosen/matlab/fieldtrip-20210825/') % Ensure this path is correct
ft_defaults;

%% ------------------------------------------------------------------------
% 1) Specify subject and file paths
% ------------------------------------------------------------------------
subj      = 'CA124';
bids_root = '/space/seh10/2/halgdev/projects/pchordiya/Cogigate/BIDS MEEG Data/COG_MEEG_EXP1_BIDS_SAMPLE'; % Ensure this path is correct

fif_file  = fullfile(bids_root, sprintf('sub-%s/ses-1/meg/sub-%s_ses-1_task-dur_run-01_meg.fif', subj, subj));
gz_mri    = fullfile(bids_root, sprintf('sub-%s/ses-1/anat/sub-%s_ses-1_T1w.nii.gz', subj, subj));
json_file = fullfile(bids_root, sprintf('sub-%s/ses-1/anat/sub-%s_ses-1_T1w.json', subj, subj));

%% ------------------------------------------------------------------------
% 1.5) Load Anatomical Landmark Coordinates (MRI Fiducials) from JSON
% These are in millimeters in the T1w space (BIDS standard).
% ------------------------------------------------------------------------
if exist(json_file, 'file')
    json_text = fileread(json_file);
    anat_landmarks_data = jsondecode(json_text);
    mri_fid_coords_raw = anat_landmarks_data.AnatomicalLandmarkCoordinates;

    target_fid_mri.nas = double(mri_fid_coords_raw.NAS(:)'); % Ensure double
    target_fid_mri.lpa = double(mri_fid_coords_raw.LPA(:)');
    target_fid_mri.rpa = double(mri_fid_coords_raw.RPA(:)');

    fprintf('Loaded MRI fiducials from JSON (NAS, LPA, RPA in mm, T1w space):\n');
    disp(['  NAS_mri: ', mat2str(target_fid_mri.nas,3)]);
    disp(['  LPA_mri: ', mat2str(target_fid_mri.lpa,3)]);
    disp(['  RPA_mri: ', mat2str(target_fid_mri.rpa,3)]);
else
    error('Anatomical landmarks JSON file not found: %s. Cannot proceed.', json_file);
end

%% ------------------------------------------------------------------------
% 2) Decompress and read the MRI
% ------------------------------------------------------------------------
fprintf('Decompressing and reading MRI...\n');
if exist(gz_mri(1:end-3), 'file') % Check if .nii already exists
    fprintf('  Found existing decompressed MRI: %s\n', gz_mri(1:end-3));
    mri_file = gz_mri(1:end-3);
else
    gunzip(gz_mri);
    mri_file = gz_mri(1:end-3);
    fprintf('  Decompressed to: %s\n', mri_file);
end
mri = ft_read_mri(mri_file);
% ft_read_mri usually sets mri.unit to 'mm' if it assumes mm based on header.
% The warning "assuming that the units are "mm"" suggests this happened.
if ~isfield(mri, 'unit') || ~strcmp(mri.unit, 'mm')
    fprintf('Warning: mri.unit is not "mm" (is "%s"). Setting to "mm" based on ft_read_mri typical behavior.\n', mri.unit);
    mri.unit = 'mm';
end
mri.coordsys = 'neuromag'; % Assert coordinate system for FieldTrip interpretation
mri.anatomy = double(mri.anatomy); % Ensure double precision
mri.transform = double(mri.transform); % Ensure double precision
fprintf('MRI read. Unit: %s, Coordsys: %s. Anatomy/Transform set to double.\n', mri.unit, mri.coordsys);



fprintf('Reading full sensor information (sens_all)...\n');
sens_all_raw = ft_read_sens(fif_file, 'readbids', 'no');
[sens_all, ~] = ensure_mm_robust(sens_all_raw, 'sens_all (electrodes)');
if ~isfield(sens_all, 'coordsys') || isempty(sens_all.coordsys), sens_all.coordsys = 'neuromag'; end
fprintf('  Sample sens_all.elecpos(1,:) after unit processing: %s (Unit: %s)\n', mat2str(sens_all.elecpos(1,:),3), sens_all.unit);

fprintf('Extracting EEG sensor information (eeg_sens)...\n');
eeg_idx_logical = [];
try % Robust EEG channel selection
    selected_output = ft_chantype(sens_all, 'eeg');
    if iscellstr(selected_output), [eeg_idx_logical, ~] = ismember(sens_all.label, selected_output);
    elseif isnumeric(selected_output) && ~islogical(selected_output), eeg_idx_logical = false(size(sens_all.label)); eeg_idx_logical(selected_output) = true;
    elseif islogical(selected_output), eeg_idx_logical = selected_output;
    else, error('ft_chantype unexpected output'); end
    if ~any(eeg_idx_logical), error('No EEG channels from ft_chantype'); end
catch ME_chantype
    fprintf('Warning: ft_chantype for EEG indices failed (%s). Using ft_channelselection.\n', ME_chantype.message);
    eeg_channel_labels_selected = ft_channelselection('EEG', sens_all.label);
    if isempty(eeg_channel_labels_selected), error('No EEG channels found via ft_channelselection.'); end
    [eeg_idx_logical, ~] = ismember(sens_all.label, eeg_channel_labels_selected);
    if ~any(eeg_idx_logical), error('No EEG channels identified.'); end
end

eeg_sens = [];
eeg_sens.label = sens_all.label(eeg_idx_logical);
if isfield(sens_all, 'elecpos') && ~isempty(sens_all.elecpos) && size(sens_all.elecpos,1) == numel(sens_all.label)
    eeg_sens.elecpos = sens_all.elecpos(eeg_idx_logical, :);
elseif isfield(sens_all, 'chanpos') && ~isempty(sens_all.chanpos) && size(sens_all.chanpos,1) == numel(sens_all.label)
    fprintf('Warning: sens_all.elecpos not found/suitable, using sens_all.chanpos for eeg_sens.elecpos\n');
    eeg_sens.elecpos = sens_all.chanpos(eeg_idx_logical, :);
else, error('Cannot get elecpos for eeg_sens.'); end
eeg_sens.chanpos = eeg_sens.elecpos; % For EEG, elecpos and chanpos are typically the same
eeg_sens.unit = sens_all.unit;       % Should be 'mm' from sens_all processing
eeg_sens.coordsys = sens_all.coordsys; % Should be 'neuromag'
if isfield(sens_all, 'chantype') && numel(sens_all.chantype) == numel(sens_all.label), eeg_sens.chantype = sens_all.chantype(eeg_idx_logical); else eeg_sens.chantype = repmat({'eeg'}, sum(eeg_idx_logical),1); end
fprintf('Successfully created eeg_sens with %d EEG channels. Unit: %s.\n', numel(eeg_sens.label), eeg_sens.unit);
fprintf('  Sample eeg_sens.elecpos(1,:): %s\n', mat2str(eeg_sens.elecpos(1,:),3));


fprintf('Reading headshape and EEG fiducials (headshape_dig)...\n');
headshape_dig_raw = ft_read_headshape(fif_file); % Your output showed this was read as 'cm'
[headshape_dig, ~] = ensure_mm_robust(headshape_dig_raw, 'headshape_dig (EEG fiducials)');
if ~isfield(headshape_dig, 'coordsys') || isempty(headshape_dig.coordsys), headshape_dig.coordsys = 'neuromag'; end
if ~isfield(headshape_dig, 'fid') || isempty(headshape_dig.fid.label)
    error('No fiducials found in headshape_dig data from FIF file. Check digitization.');
end
fprintf('  Sample headshape_dig.fid.pos(1,:) after unit processing: %s (Unit: %s)\n', mat2str(headshape_dig.fid.pos(1,:),3), headshape_dig.unit);

fprintf('Extracting template_fid_head (EEG fiducials)...\n');
template_fid_head = [];
fid_labels_to_find = {'nasion', 'nas', 'nz'; 'lpa', 'left', 'left posterior auricular'; 'rpa', 'right', 'right posterior auricular'};
fid_field_names = {'nas', 'lpa', 'rpa'};
for i = 1:size(fid_labels_to_find, 1)
    found = false;
    for j = 1:size(fid_labels_to_find, 2)
        idx = find(strcmpi(headshape_dig.fid.label, fid_labels_to_find{i,j}));
        if ~isempty(idx)
            template_fid_head.(fid_field_names{i}) = double(headshape_dig.fid.pos(idx(1), :)); % Ensure double
            found = true;
            break;
        end
    end
    if ~found, error('Could not find %s fiducial in headshape data. Available: %s', fid_field_names{i}, strjoin(headshape_dig.fid.label)); end
end
fprintf('Loaded EEG/Head fiducials (template_fid_head, should be mm):\n');
disp(['  NAS_eeg: ', mat2str(template_fid_head.nas,3)]);
disp(['  LPA_eeg: ', mat2str(template_fid_head.lpa,3)]);
disp(['  RPA_eeg: ', mat2str(template_fid_head.rpa,3)]);

% ... (all previous code for loading data and converting units - THIS IS ALL GOOD NOW) ...

%% ------------------------------------------------------------------------
% 4) Perform Fiducial-Based Alignment (Focus on INTERACTIVE METHOD for debugging)
% ------------------------------------------------------------------------
fprintf('MRI, Target (JSON) Fiducials, and Template (EEG) Fiducials are all in mm and double type.\n');
% (Your previous fprintf statements for displaying fiducial values are good here)

original_mri_transform = mri.transform; % Store for comparison

fprintf('Attempting INTERACTIVE realignment...\n');
cfg_interactive = [];
cfg_interactive.method   = 'interactive';
cfg_interactive.coordsys = mri.coordsys; % e.g., 'neuromag'

fprintf('\n>>> INSTRUCTIONS FOR INTERACTIVE REALIGNMENT GUI <<<\n');
fprintf('1. The FieldTrip realignment GUI will open.\n');
fprintf('2. In the GUI, you need to specify the fiducial points for BOTH the MRI and the HEAD (EEG).\n');
fprintf('3. For "MRI Fiducials" (this is your TARGET):\n');
fprintf('   NAS: %s\n', mat2str(target_fid_mri.nas,3));
fprintf('   LPA: %s\n', mat2str(target_fid_mri.lpa,3));
fprintf('   RPA: %s\n', mat2str(target_fid_mri.rpa,3));
fprintf('4. For "Head Fiducials" (this is your TEMPLATE, from EEG digitization):\n');
fprintf('   NAS: %s\n', mat2str(template_fid_head.nas,3));
fprintf('   LPA: %s\n', mat2str(template_fid_head.lpa,3));
fprintf('   RPA: %s\n', mat2str(template_fid_head.rpa,3));
fprintf('5. After entering the fiducials, use the GUI options to compute the alignment based on these fiducials.\n');
fprintf('6. You can then make further manual adjustments if needed.\n');
fprintf('7. When satisfied, close the GUI to apply the transformation.\n');
fprintf('------------------------------------------------------------------\n');

mri_aligned = []; % Initialize
try
    mri_aligned = ft_volumerealign(cfg_interactive, mri); % Pass original mri
    if isequal(mri_aligned.transform, original_mri_transform)
        fprintf('WARNING: Interactive ft_volumerealign did NOT change MRI.transform. Ensure you clicked "Apply" or that alignment was possible.\n');
        % Ask user if they want to proceed with unaligned or stop
        proceed_unaligned = input('Proceed with UNALIGNED MRI? (yes/no): ', 's');
        if ~strcmpi(proceed_unaligned, 'yes')
            error('User opted to stop due to failed interactive alignment.');
        end
        mri_aligned = mri; % Fallback to unaligned
    else
        fprintf('Interactive ft_volumerealign successfully updated MRI.transform.\n');
    end
catch ME_interactive_realign
    fprintf('ERROR during interactive ft_volumerealign: %s\n', ME_interactive_realign.message);
    mri_aligned = mri; % Fallback to unaligned MRI
    fprintf('CRITICAL: Using UNALIGNED MRI due to interactive alignment error.\n');
end

% ... (Rest of your script: Segmentation, Mesh Prep, Plotting using mri_aligned) ...
%% ------------------------------------------------------------------------
% 5) Segment the scalp surface FROM THE ALIGNED MRI and prepare a mesh
% ------------------------------------------------------------------------
fprintf('Segmenting scalp from ALIGNED MRI and preparing mesh...\n');
cfg_segment = [];
cfg_segment.output = {'scalp'};
cfg_segment.scalpsmooth = 5; % Added from your log, seems like a default you had
cfg_segment.scalpthreshold = 0.1; % Added from your log
segmented_scalp  = ft_volumesegment(cfg_segment, mri_aligned);

cfg_mesh = [];
cfg_mesh.method = 'projectmesh';
cfg_mesh.tissue = 'scalp';
cfg_mesh.numvertices = 10000;
scalp_mesh = ft_prepare_mesh(cfg_mesh, segmented_scalp);
if ~isfield(scalp_mesh, 'coordsys') || isempty(scalp_mesh.coordsys) && isfield(mri_aligned, 'coordsys')
    scalp_mesh.coordsys = mri_aligned.coordsys;
end
fprintf('Scalp mesh prepared. Unit: %s, Coordsys: %s\n', scalp_mesh.unit, scalp_mesh.coordsys); % scalp_mesh inherits mri_aligned.unit

%% ------------------------------------------------------------------------
% 6) Plot the scalp mesh and overlay EEG electrodes WITH LABELS
% ------------------------------------------------------------------------
fprintf('Plotting aligned data with electrode labels...\n');
figure('Color','w','Renderer','opengl', 'Name', sprintf('EEG Coregistration: %s (Interactive Alignment)', subj));

% Plot the scalp mesh
ft_plot_mesh(scalp_mesh, ...
    'facecolor', [0.9 0.85 0.8], ...
    'edgecolor', 'none', ...
    'facealpha', 0.9);
hold on;

% Plot the EEG sensors with labels
ft_plot_sens(eeg_sens, ...
    'elec', true, ...               % Plot electrodes
    'elecshape', 'sphere', ...      % Use spheres for electrodes
    'elecsize', 6, ...              % Size of the spheres (adjust as needed)
    'facecolor', 'b', ...           % Color of the electrode spheres
    'label', 'label', ...            % Display channel names from eeg_sens.label
    'fontsize', 8, ...         % Font size for the labels (adjust for readability)
    'labelcolor', [220 220 220]);       % Color for the labels (e.g., dark blue, adjust as needed)
                                    % Default is often black if not specified.

fprintf('--- Final check for plotting ---\n');
fprintf('eeg_sens unit: %s, coordsys: %s, sample elecpos(1,:): %s\n', eeg_sens.unit, eeg_sens.coordsys, mat2str(eeg_sens.elecpos(1,:),3) );
fprintf('scalp_mesh unit: %s, coordsys: %s, sample vertex(1,:): %s\n', scalp_mesh.unit, scalp_mesh.coordsys, mat2str(scalp_mesh.pos(1,:),3) );
fprintf('Median abs eeg_sens.elecpos: %.2f\n', median(abs(eeg_sens.elecpos(:))));
fprintf('Median abs scalp_mesh.pos: %.2f\n', median(abs(scalp_mesh.pos(:))));

% Plot EEG fiducials (template_fid_head, which are now mm)
fid_colors = {'g','g','g'}; % NAS, LPA, RPA
fid_fields = fieldnames(template_fid_head);
if isstruct(template_fid_head) && ~isempty(fid_fields) % Check if template_fid_head is a populated struct
    for i = 1:length(fid_fields)
        if isfield(template_fid_head, fid_fields{i}) % Check if specific fiducial exists
            fid_pos = template_fid_head.(fid_fields{i});
            plot3(fid_pos(1), fid_pos(2), fid_pos(3), 'o', 'MarkerFaceColor', fid_colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
            text(fid_pos(1)+5, fid_pos(2), fid_pos(3), upper(fid_fields{i}), 'Color', 'k', 'FontSize',10, 'FontWeight', 'bold');
        end
    end
end

camlight; lighting gouraud;
view([-90 20]); % Adjust view as needed
axis equal; axis off;
title(sprintf('EEG Electrodes on %s''s Scalp (Interactive Alignment)', subj), 'FontSize', 14);

fprintf('Script finished. Check plot for electrode labels.\n');





%% ------------------------------------------------------------------------
% 3) Read sensor information (sens_all) and Headshape/EEG Fiducials (headshape_dig)
%    AND ENSURE THEY ARE IN MILLIMETERS
% ------------------------------------------------------------------------

% --- Helper function for unit checking and conversion ---
    function [data_out, original_unit_best_guess] = ensure_mm_robust(data_in, name_str)
        data_out = data_in;
        original_unit_best_guess = 'unknown_initial'; % To track what we thought it was

        pos_data = []; % Get primary position data for magnitude check
        if isfield(data_out, 'elecpos') && ~isempty(data_out.elecpos), pos_data = data_out.elecpos;
        elseif isfield(data_out, 'chanpos') && ~isempty(data_out.chanpos), pos_data = data_out.chanpos;
        elseif isfield(data_out, 'pos') && ~isempty(data_out.pos), pos_data = data_out.pos;
        elseif isfield(data_out, 'pnt') && ~isempty(data_out.pnt), pos_data = data_out.pnt;
        elseif isfield(data_out, 'fid') && isfield(data_out.fid, 'pos') && ~isempty(data_out.fid.pos), pos_data = data_out.fid.pos;
        end

        reported_unit = '';
        if isfield(data_out, 'unit') && ~isempty(data_out.unit)
            reported_unit = lower(data_out.unit);
            fprintf('%s: Reported original unit: "%s".\n', name_str, data_out.unit);
        else
            fprintf('%s: Original unit field not set or empty.\n', name_str);
        end

        if strcmp(reported_unit, 'mm')
            fprintf('%s: Reported unit is already mm. Assuming values are correct.\n', name_str);
            original_unit_best_guess = 'mm';
        elseif strcmp(reported_unit, 'cm')
            fprintf('%s: Reported unit is cm. Converting to mm.\n', name_str);
            data_out = ft_convert_units(data_out, 'mm');
            original_unit_best_guess = 'cm';
        elseif strcmp(reported_unit, 'm')
            fprintf('%s: Reported unit is m. Converting to mm.\n', name_str);
            data_out = ft_convert_units(data_out, 'mm');
            original_unit_best_guess = 'm';
        else % Unit is empty, 'unknown', or something else -> use heuristic
            if ~isempty(pos_data)
                median_abs_val = median(abs(pos_data(:)));
                fprintf('%s: Unit was "%s". Median abs position value for heuristic: %.3f\n', name_str, reported_unit, median_abs_val);
                if median_abs_val > 1e-5 && median_abs_val < 2.0  % Likely meters
                    fprintf('%s: Magnitude (%.3f) suggests METERS. Setting unit to "m" for conversion.\n', name_str, median_abs_val);
                    data_out.unit = 'm'; % Set unit to allow ft_convert_units to recognize original scale
                    data_out = ft_convert_units(data_out, 'mm');
                    original_unit_best_guess = 'm (guessed)';
                elseif median_abs_val >= 2.0 && median_abs_val < 30.0 % Likely centimeters
                    fprintf('%s: Magnitude (%.3f) suggests CENTIMETERS. Setting unit to "cm" for conversion.\n', name_str, median_abs_val);
                    data_out.unit = 'cm';
                    data_out = ft_convert_units(data_out, 'mm');
                    original_unit_best_guess = 'cm (guessed)';
                else % Assume already mm or too ambiguous if very large or very small (e.g. zero)
                    fprintf('%s: Magnitude (%.3f) suggests already MILLIMETERS or is ambiguous/zero. Assuming mm.\n', name_str, median_abs_val);
                    data_out.unit = 'mm'; % Ensure unit field is set
                    original_unit_best_guess = 'mm (guessed/assumed)';
                end
            else
                fprintf('%s: No position data to guess units. Assuming mm.\n', name_str);
                data_out.unit = 'mm';
                original_unit_best_guess = 'mm (defaulted)';
            end
        end
        % Final check to ensure .unit field is 'mm'
        if ~isfield(data_out, 'unit') || ~strcmp(data_out.unit, 'mm')
            data_out.unit = 'mm';
        end
        fprintf('%s: Final unit after processing: %s.\n', name_str, data_out.unit);
    end
% --- End of helper function ---