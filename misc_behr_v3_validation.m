classdef misc_behr_v3_validation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties(Constant = true)
        behr_v2_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_v2-1C';
        gcas_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignRaw/DISCOVER-AQ_TX/B200/GCAS-SAO';
    end
    

    
    methods(Static = true)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property-like Methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function value = my_dir()
            value = fileparts(mfilename('fullpath'));
        end
        
        function value = validation_root_dir()
            value = fullfile(misc_behr_v3_validation.my_dir, 'Workspaces', 'Validation');
        end
        
        function value = profile_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', 'profile-structs.mat');
        end
        
        function value = gcas_comp_file()
            value = fullfile(misc_behr_update_plots.validation_root_dir, 'VCD-Comparison', 'gcas-structs.mat');
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generative validation methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function insitu_comparison = generate_vcd_comparison_structs(do_save)
            E = JLLErrors;
            
            if ~exist('do_save', 'var')
                do_save = ask_yn('Save the comparison files generated?');
            elseif ~isscalar(do_save) || ~islogical(do_save)
                E.badinput('DO_SAVE must be a scalar logical');
            end
            
            versions = {'v2','v3'};
            regions = {'us'};
            prof_modes = {'monthly', 'daily'};
            behr_field = 'BEHRColumnAmountNO2Trop'; % can switch to VisOnly later
            insitu_comparison = struct();
            for v=1:numel(versions)
                if strcmpi(versions{v}, 'v3')
                    n_regions = numel(regions);
                    n_profs = numel(prof_modes);
                else
                    n_regions = 1;
                    n_profs = 1;
                end
                
                version_comp = struct();
                
                for r=1:n_regions
                    for p=1:n_profs
                        if strcmpi(versions{v}, 'v3')
                            behr_dir = behr_paths.BEHRMatSubdir(regions{r}, prof_modes{p});
                            behr_prefix = regexp(behr_filename(today, prof_modes{p}, regions{r}), 'OMI_BEHR.*(?=\d\d\d\d\d\d\d\d)', 'match', 'once');
                        else
                            behr_dir = misc_behr_v3_validation.behr_v2_dir;
                            behr_prefix = 'OMI_BEHR_v2-1C_';
                        end
                        % First we do the four DISCOVER campaigns, which are nice
                        % because they are geared towards satellite validation with
                        % spirals clearly marked in the data.
                        spiral_campaigns = {'discover_md', 'discover_ca', 'discover_tx', 'discover_co'};
                        other_campaigns = {'seac4rs', 'dc3', 'arctas_carb', 'intex_b', 'soas'};
                        % For each campaign that doesn't identify the
                        % profiles in the data, we need to specify which
                        % of the precreated range files that identify
                        % profiles by the UTC ranges to use.
                        range_files = {'SEAC4RS_Profile_Ranges_Std.mat',...
                                       'DC3_Profile_Ranges_Std.mat',...
                                       'ARCTAS-CA Altitude Ranges Exclusive 3.mat',...
                                       'INTEXB_Profile_UTC_Ranges.mat',...
                                       'SOAS_Profile_UTC_Ranges.mat'};
                        time_windows = [1.5, 3]; % how far away from satellite overpass (in hours) the profile is allowed to be
                        
                        behr_comp = struct();
                        all_campaigns = [spiral_campaigns, other_campaigns];
                        campaign_use_ranges = [repmat({''}, size(spiral_campaigns)), repmat({'ranges'}, size(other_campaigns))];
                        campaign_range_files = [repmat({false}, size(spiral_campaigns)), range_files];
                        % Specify whether to use the LIF or
                        % NCAR/chemiluminesnce/non-LIF NO2 measurement.
                        no2_fields = {'lif','lif','lif','lif','lif','lif','lif','lif','cl'};
                        for a=1:numel(all_campaigns)
                            campaign = struct();
                            for b=1:numel(time_windows)
                                [prof_times, time_fieldname] = misc_behr_v3_validation.format_profile_time_range(time_windows(b));
                                comparison_params = {'campaign', all_campaigns{a},...
                                    'profile_input', campaign_use_ranges{a},...
                                    'ask_range', campaign_range_files{a},...
                                    'behr_dir', behr_dir,...
                                    'behr_prefix', behr_prefix,...
                                    'starttime', prof_times{1},...
                                    'endtime', prof_times{2},...
                                    'minheight', 0,...
                                    'minagl', 0.5,...
                                    'cloudtype', 'omi',...
                                    'cloudfrac', 0.2,...
                                    'rowanomaly', 'XTrackFlags',...
                                    'behrfield', behr_field,...
                                    'no2field', no2_fields{a},...
                                    'surf_pres_choice', 'yes',...
                                    'debug', 1};
                                try
                                    [sub.lon, sub.lat, sub.omi_no2, sub.behr_no2, sub.air_no2, sub.dbinfo, sub.dates] = Run_Spiral_Verification(comparison_params{:});
                                catch err
                                    if strcmp(err.identifier, 'Run_Spiral_Verification:run_failure')
                                        % The daily product will not have
                                        % data for some days which causes
                                        % an error in
                                        % Run_Spiral_Verification. 
                                        sub = [];
                                        fprintf('No BEHR %s data available for region = %s, prof_mode = %s\n', versions{v}, regions{r}, prof_modes{p});
                                    else
                                        rethrow(err);
                                    end
                                end
                                
                                campaign.(time_fieldname) = sub;
                            end
                            behr_comp.(all_campaigns{a}) = campaign;
                        end
                        
                        fn = sprintf('%s_%s', regions{r}, prof_modes{p});
                        version_comp.(fn) = behr_comp;
                    end
                end
                insitu_comparison.(versions{v}) = version_comp;
            end
            
            if do_save
                save(misc_behr_v3_validation.gcas_comp_file, '-struct', 'insitu_comparison');
            end
        end
        
        function [times, time_fieldname] = format_profile_time_range(win)
            base_time = datenum('2000-01-01 13:30:00');
            start_time = base_time - win/24;
            end_time = base_time + win/24;
            
            times = {datestr(start_time, 'HH:MM'), datestr(end_time, 'HH:MM')};
            time_fieldname = sprintf('t%s_%s', datestr(start_time, 'HHMM'), datestr(end_time, 'HHMM'));
        end
        
        function gcas_comparison = generate_gcas_comparison_struct(do_save)
            E = JLLErrors;
            
            if ~exist('do_save', 'var')
                do_save = ask_yn('Save the comparison files generated?');
            elseif ~isscalar(do_save) || ~islogical(do_save)
                E.badinput('DO_SAVE must be a scalar logical');
            end
            
            versions = {'v2','v3'};
            regions = {'us'};
            prof_modes = {'monthly', 'daily'};
            sat_fields = {'BEHRColumnAmountNO2Trop', 'ColumnAmountNO2Trop'}; % can add VisOnly later if desired
            all_campaigns = {'discover_tx'};
            time_windows = [1.5, 3];
            gcas_comparison = struct();
            for v=1:numel(versions)
                if strcmpi(versions{v}, 'v3')
                    n_regions = numel(regions);
                    n_profs = numel(prof_modes);
                else
                    n_regions = 1;
                    n_profs = 1;
                end
                
                version_comp = struct();
                
                for r=1:n_regions
                    for p=1:n_profs
                        if strcmpi(versions{v}, 'v3')
                            behr_dir = behr_paths.BEHRMatSubdir(regions{r}, prof_modes{p});
                        else
                            behr_dir = misc_behr_v3_validation.behr_v2_dir;
                        end
                        
                        behr_comp = struct();
                        
                        for a=1:numel(all_campaigns)
                            campaign = struct();
                            for b=1:numel(time_windows)
                                [~, time_fn] = misc_behr_v3_validation.format_profile_time_range(time_windows(b));
                                % Construct the options list
                                opts_list = {'campaign', all_campaigns{a},...
                                             'behr_dir', behr_dir,...
                                             'vectorize', false,...
                                             'cloud_prod', 'omi',...
                                             'cloud_frac_max', 0.2,...
                                             'row_anomaly', 'XTrackFlags',...
                                             'sat_fields', sat_fields,...
                                             'time_window', time_windows(b)};
                                try
                                    campaign.(time_fn) = run_gcas_verification(opts_list{:});
                                catch err
                                    if strcmp(err.identifier, 'load_behr_file_for_gcas:file_not_found')
                                        campaign.(time_fn) = [];
                                    else
                                        rethrow(err)
                                    end
                                end
                            end
                            
                            behr_comp.(all_campaigns{a}) = campaign;
                        end
                        
                        fn = sprintf('%s_%s', regions{r}, prof_modes{p});
                        version_comp.(fn) = behr_comp;
                    end
                end
                gcas_comparison.(versions{v}) = version_comp;
            end
            
            if do_save
                save(misc_behr_v3_validation.profile_comp_file, '-struct', 'gcas_comparison');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions %
        %%%%%%%%%%%%%%%%%%%%%%
        
        function plot_one_vcd_comparison(x_var, y_var, prof_mode, campaigns, time_range)
            % Plots one v2-v3 comparison against in situ aircraft data
            % derived columns. prof_mode may be 'monthly', 'daily' or
            % 'both'. campaign and time_range must match the fields in the
            % misc_behr_v3_validation.profile_comp_file. If omitted these
            % inputs will be asked interactively
            E = JLLErrors;
            
            comp_struct = load(misc_behr_v3_validation.profile_comp_file);
            
            allowed_vars = {'air_no2', 'omi_no2', 'behr_no2'};
            if ~exist('x_var', 'var')
                x_var = ask_multichoice('Which variable to plot on the x-axis?', allowed_vars, 'list', true);
            elseif ~ischar(x_var) || ~ismember(x_var, allowed_vars)
                E.badinput('X_VAR must be one of the strings: %s', strjoin(allowed_vars, ', '));
            end
            
            if ~exist('y_var', 'var')
                y_var = ask_multichoice('Which variable to plot on the y-axis?', allowed_vars, 'list', true);
            elseif ~ischar(y_var) || ~ismember(y_var, allowed_vars)
                E.badinput('Y_VAR must be one of the strings: %s', strjoin(allowed_vars, ', '));
            end
            
            allowed_prof_modes = {'monthly', 'daily', 'both'};
            if strcmpi(x_var, 'behr_no2') || strcmpi(y_var, 'behr_no2')
                if ~exist('prof_mode', 'var')
                    prof_mode = ask_multichoice('Which profile mode to use for v3 data?', allowed_prof_modes, 'list', true);
                elseif ~ismember(prof_mode, allowed_prof_modes)
                    E.badinput('PROF_MODE must be one of %s', strjoin(prof_mode, ', '));
                end
            else
                % There's no difference in the SP data between our monthly
                % and daily profile product
                prof_mode = 'monthly';
            end
            
            allowed_campaigns = fieldnames(comp_struct.v2.us_monthly);
            if ~exist('campaigns', 'var')
                campaigns = ask_multiselect('Which campaign(s) to plot?', allowed_campaigns);
            else
                if ischar(campaigns)
                    campaigns = {campaigns};
                elseif ~iscellstr(campaigns)
                    E.badinput('CAMPAIGN must be a string or cell array of strings')
                end
                if any(~ismember(campaigns, allowed_campaigns))
                    E.badinput('All strings in CAMPAIGN must be one of: %s', strjoin(allowed_campaigns, ', '));
                end
            end
            
            allowed_times = fieldnames(comp_struct.v2.us_monthly.discover_md);
            if ~exist('time_range', 'var')
                time_range = ask_multichoice('Which time range to use?', allowed_times, 'list', true);
            elseif ~ischar(time_range) || ~ismember(time_range, allowed_times)
                E.badinput('TIME_RANGE must be one of the strings: %s', strjoin(allowed_times, ', '));
            end
            
            if strcmpi(x_var, 'behr_no2') || strcmpi(y_var, 'behr_no2')
                v2_string = 'v2.1C';
                v3M_string = 'v3.0A (M)';
                v3D_string = 'v3.0A (D)';
            else
                v2_string = 'v2.1';
                v3M_string = 'v3.0';
                v3D_string = 'v3.0';
            end
            
            
            switch lower(prof_mode)
                case 'monthly'
                    data_structs = {comp_struct.v2.us_monthly, comp_struct.v3.us_monthly};
                    legend_strings = {v2_string, v3M_string};
                case 'daily'
                    data_structs = {comp_struct.v2.us_monthly, comp_struct.v3.us_daily};
                    legend_strings = {v2_string, v3D_string};
                case 'both'
                    data_structs = {comp_struct.v2.us_monthly, comp_struct.v3.us_monthly, comp_struct.v3.us_daily};
                    legend_strings = {v2_string, v3M_string, v3D_string};
                otherwise
                    E.notimplemented('prof_mode == %s', prof_mode);
            end
            
            
            
            
            data_lines = gobjects(size(data_structs));
            fit_lines = gobjects(size(data_structs));
            fit_legends = cell(size(data_structs));
            line_fmts = struct('color', {'k','b','r', [0, 0.5, 0]}, 'marker', {'s','x','^','o'});
            
            labels = struct('air_no2', 'Aircraft NO_2 VCD (molec. cm^{-2})',...
                'omi_no2', 'NASA NO_2 VCD (molec. cm^{-2})',...
                'behr_no2', 'BEHR NO_2 VCD (molec. cm^{-2})');
            
            figure;
            limits = [Inf, -Inf];
            for a=1:numel(data_structs)
                x_no2 = [];
                y_no2 = [];
                for b=1:numel(campaigns)
                    x_no2 = veccat(x_no2, data_structs{a}.(campaigns{b}).(time_range).(x_var));
                    y_no2 = veccat(y_no2, data_structs{a}.(campaigns{b}).(time_range).(y_var));
                    limits(1) = min([limits(1), min(x_no2), min(y_no2)]);
                    limits(2) = max([limits(2), max(x_no2), max(y_no2)]);
                end
                
                data_lines(a) = line(x_no2(:), y_no2(:), 'linestyle', 'none', 'color', line_fmts(a).color, 'marker', line_fmts(a).marker);
                [fit_x, fit_y, fit_legends{a}] = calc_fit_line(x_no2(:), y_no2(:), 'regression', 'RMA', 'xcoord', [-2e16, 2e16]);
                fit_lines(a) = line(fit_x, fit_y, 'color', line_fmts(a).color, 'linestyle', '--');
            end
            
            xylims(calc_plot_limits(limits, 'pow10'));
            
            xlabel(labels.(x_var));
            ylabel(labels.(y_var));
            
            all_lines = gobjects(numel(data_lines)*2, 1);
            all_lines(1:2:end) = data_lines;
            all_lines(2:2:end) = fit_lines;
            
            all_legends = cell(1, numel(data_structs)*2);
            all_legends(1:2:end) = legend_strings;
            all_legends(2:2:end) = fit_legends;
            
            legend(all_lines, all_legends);
            title_campaigns = strrep(campaigns, '_', '\_');
            title(strjoin(upper(title_campaigns), ', '));
        end
    end
    
end

