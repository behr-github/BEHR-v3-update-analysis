classdef misc_behr_v3_validation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties(Constant = true)
        behr_v2_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_v2-1C';
        wrf_v2_dir = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles';
        gcas_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignRaw/DISCOVER-AQ_TX/B200/GCAS-SAO';
        
        plot_colors = struct('aircraft', struct('raw', [0.5 0.5 0.5], 'avg', 'k'),... % black and grey for the aircraft data
                'v2', struct('raw', 'g', 'avg', [0 0.5 0]),... % greens for version 2
                'monthly', struct('raw', [1 0.75 0], 'avg', 'r'),... % red and orange for the monthly profiles
                'daily', struct('raw', 'c', 'avg', 'b')); % blue and cyan for the daily profiles
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
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', 'gcas-structs.mat');
        end
        
        function value = gcas_vec_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', 'gcas-vec-structs.mat');
        end
        
        function value = scd_comp_dir()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'SCD-Wind-Comparisons');
            if ~exist(value, 'dir')
                mkdir(value)
            end
        end
        
        function value = scd_comp_file(is_partial)
            if ~exist('is_partial', 'var') || ~is_partial
                partial_str = '';
            else
                partial_str = '_partial';
            end
            value = fullfile(misc_behr_v3_validation.scd_comp_dir, sprintf('scd_daily_%s%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'), partial_str));
        end
        
        function value = wrf_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'Profile-Comparison', 'wrf-match-structs.mat');
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
        
        function [gcas_comparison, gcas_comparison_vec] = generate_gcas_comparison_struct(do_save)
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
            gcas_comparison_vec = struct();
            for v=1:numel(versions)
                if strcmpi(versions{v}, 'v3')
                    n_regions = numel(regions);
                    n_profs = numel(prof_modes);
                else
                    n_regions = 1;
                    n_profs = 1;
                end
                
                version_comp = struct();
                version_comp_vec = struct();
                
                for r=1:n_regions
                    for p=1:n_profs
                        if strcmpi(versions{v}, 'v3')
                            behr_dir = behr_paths.BEHRMatSubdir(regions{r}, prof_modes{p});
                        else
                            behr_dir = misc_behr_v3_validation.behr_v2_dir;
                        end
                        
                        behr_comp = struct();
                        behr_comp_vec = struct();
                        
                        for a=1:numel(all_campaigns)
                            campaign = struct();
                            campaign_vec = struct();
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
                                    campaign_vec.(time_fn) = vectorize_gcas_matches(campaign.(time_fn));
                                catch err
                                    if strcmp(err.identifier, 'load_behr_file_for_gcas:file_not_found')
                                        campaign.(time_fn) = [];
                                        campaign_vec.(time_fn) = [];
                                    else
                                        rethrow(err)
                                    end
                                end
                            end
                            
                            behr_comp.(all_campaigns{a}) = campaign;
                            behr_comp_vec.(all_campaigns{a}) = campaign_vec;
                        end
                        
                        fn = sprintf('%s_%s', regions{r}, prof_modes{p});
                        version_comp.(fn) = behr_comp;
                        version_comp_vec.(fn) = behr_comp_vec;
                    end
                end
                gcas_comparison.(versions{v}) = version_comp;
                gcas_comparison_vec.(versions{v}) = version_comp_vec;
            end
            
            if do_save
                save(misc_behr_v3_validation.gcas_comp_file, '-struct', 'gcas_comparison');
                save(misc_behr_v3_validation.gcas_vec_comp_file, '-struct', 'gcas_comparison_vec');
            end
        end
        
        function generate_profile_comparison_struct(do_save)
            % This will generate a structure containing WRF data matched to
            % individual DISCOVER profiles, as well as information about
            % the OMI overpass time vs. the profile time
            % 
            % The hierarchy of the structure will be:
            %   profile type (daily, monthly, or version 2)
            %       campaign
            %           pXXXXXX (XXXXXX is the profile number)
            %               profdate (average profile UTC datetime)
            %               nearest OMI overpass (in space) time
            %               match
            
            if ~exist('do_save', 'var')
                do_save = ask_yn('Save the resulting match structures?');
            end
            
            profile_types = {'daily', 'monthly', 'v2'};
            campaigns = {'discover_md', 'discover_ca', 'discover_tx', 'discover_co'};
            %campaigns = {'discover_ca'};
            
            % Loop through each campaign, load each merge file, find each
            % profile, extract the necessary fields to match up the
            % aircraft and WRF data, and pass it to match_wrf2aircraft to
            % do the actual matching
            profile_comp_struct = make_empty_struct_from_cell(profile_types);
            
            for a=1:numel(campaigns)
                % List the merge files
                [merge_names, ~, merge_dir] = merge_field_names(campaigns{a});
                merges = dirff(fullfile(merge_dir, '*.mat'));
                all_merge_wrf_dirs = make_empty_struct_from_cell(profile_types, {{}});
                for b=1:numel(merges)
                    M = load(merges(b).name);
                    % Make raw structures for each profile that is listed
                    % in the file
                    Raws = misc_behr_v3_validation.make_raw_struct(M.Merge, merge_names, campaigns{a});
                    profile_fns = fieldnames(Raws);
                    
                    % Calculate the OMI overpass time from the v3 BEHR
                    % files, getting the closest overpass in space. We can
                    % use this to decide if a profile is important for the
                    % retrieval
                    OmiTimes = misc_behr_v3_validation.calc_nearest_omi_times(M.Merge.metadata.date, Raws);
                    
                    for p=1:numel(profile_types)
                        % Get the right WRF directory for the given profile
                        % type
                        if strcmpi(profile_types{p}, 'v2')
                            wrf_dir = misc_behr_v3_validation.wrf_v2_dir;
                        else
                            try
                                wrf_dir = find_wrf_path('us', profile_types{p}, M.Merge.metadata.date);
                            catch err
                                if strcmp(err.identifier, 'find_wrf_path:dir_does_not_exist')
                                    continue
                                else
                                    rethrow(err);
                                end
                            end
                        end
                        
                        if ~ismember(wrf_dir, all_merge_wrf_dirs.(profile_types{p}))
                            all_merge_wrf_dirs.(profile_types{p}){end+1} = wrf_dir;
                        end
                        
                        for f=1:numel(profile_fns)
                            Match = match_wrf2aircraft(Raws.(profile_fns{f}), wrf_dir, profile_types{p});
                            xlon = Match.wrf.xlon;
                            xlat = Match.wrf.xlat;
                            % These will be the same in every Match
                            % structure, including them greatly increases
                            % the size of the final file unnecessarily
                            Match.wrf = rmfield(Match.wrf, 'xlon');
                            Match.wrf = rmfield(Match.wrf, 'xlat');
                            profile_comp_struct.(profile_types{p}).(campaigns{a}).(profile_fns{f}).match = Match;
                            profile_comp_struct.(profile_types{p}).(campaigns{a}).(profile_fns{f}).prof_date = mean(Raws.(profile_fns{f}).dvec);
                            profile_comp_struct.(profile_types{p}).(campaigns{a}).(profile_fns{f}).omi_time = OmiTimes.(profile_fns{f});
                            profile_comp_struct.wrf_xlon = xlon;
                            profile_comp_struct.wrf_xlat = xlat;
                        end
                    end
                end
                
                % Match the WRF output to the entire P3 flight
                for p=1:numel(profile_types)
                    if ~isempty(all_merge_wrf_dirs.(profile_types{p}))
                        Match = match_wrf2campaigns.(campaigns{a})(profile_types{p}, all_merge_wrf_dirs.(profile_types{p}));
                        profile_comp_struct.(profile_types{p}).(campaigns{a}).All.match = Match;
                    end
                end
            end
            
            if do_save
                save(misc_behr_v3_validation.wrf_comp_file, '-struct', 'profile_comp_struct')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        % Utility functions %
        %%%%%%%%%%%%%%%%%%%%%
        
        function OmiTimes = calc_nearest_omi_times(this_date, Raws)
            Data = load_behr_file(this_date, 'monthly', 'us');
            % For each orbit, get the center longitude/latitude line. Also
            % get the average time
            center_lon = cell(size(Data));
            center_lat = cell(size(Data));
            avg_time = cell(size(Data));
            for a=1:numel(Data)
                n = round(size(Data(a).Longitude,2)/2);
                center_lon{a} = Data(a).Longitude(:,n);
                center_lat{a} = Data(a).Latitude(:,n);
                avg_time{a} = omi_time_conv(mean(Data(a).Time));
            end
            
            OmiTimes = make_empty_struct_from_cell(fieldnames(Raws));
            % For each profile, calculate the average lon/lat and figure
            % out which overpass comes the closest
            raw_fns = fieldnames(Raws);
            for f=1:numel(raw_fns)
                shortest_dist = Inf;
                omi_time = NaN;
                prof_avg_lon = nanmean(Raws.(raw_fns{a}).lon);
                prof_avg_lat = nanmean(Raws.(raw_fns{a}).lat);
                for a=1:numel(center_lon)
                    min_dist = min(sqrt( (prof_avg_lon - center_lon{a}).^2 + (prof_avg_lat - center_lat{a}).^2 ));
                    if min_dist < shortest_dist
                        shortest_dist = min_dist;
                        omi_time = avg_time{a};
                    end
                end
                OmiTimes.(raw_fns{f}) = omi_time;
            end
        end
        
        function Raws = make_raw_struct(Merge, merge_names, campaign)
            % The Raw structure for input into match_wrf2aircraft requires
            % the fields lon, lat, pres, dvec, and campaign. pres must be
            % pressure in hPa. dvec must be a vector that gives the date
            % and UTC time as a Matlab datenumber.
            
            % First get all the profile numbers
            profnums = remove_merge_fills(Merge, merge_names.profile_numbers);
            u_profnums = unique(profnums(profnums > 0));
            prof_fns = sprintfmulti('p%d', u_profnums);
            
            substruct = struct('lon', [], 'lat', [], 'pres', [], 'dvec', [], 'campaign', campaign);
            Raws = make_empty_struct_from_cell(prof_fns, substruct);
            
            % Read the other necessary variables. Read in UTC with NO2 to
            % ensure it has the same fill values.
            [AllRaw.no2, utc, ~, AllRaw.lon, AllRaw.lat] = remove_merge_fills(Merge, merge_names.no2_lif, 'unit', 'ppm',...
                'lon', merge_names.longitude, 'lat', merge_names.latitude);
            AllRaw.pres = remove_merge_fills(Merge, merge_names.pressure, 'unit', 'hPa');
            
            % UTC is given in seconds after midnight.
            if any(utc < 0)
                warning('Fill values in UTC time vector')
            end
            AllRaw.dvec = datenum(Merge.metadata.date) + utc ./ (24*60*60);
            
            % Now assign the proper subsets to the individual profile
            % fields
            data_fns = fieldnames(AllRaw);
            for a=1:numel(u_profnums)
                pp = profnums == u_profnums(a);
                for b=1:numel(data_fns)
                    Raws.(prof_fns{a}).(data_fns{b}) = AllRaw.(data_fns{b})(pp);
                end
            end
        end
        
        function [no2, xlon, xlat] = load_wrf_no2(prof_type, wrf_date, west_east_inds, north_south_inds, bottom_top_ind)
            if strcmpi(prof_type, 'v2')
                wrf_file_name = sprintf('m%02d_NO2_profile.mat', wrf_date);
                W = load(fullfile(misc_behr_v3_validation.wrf_v2_dir, wrf_file_name));
                no2 = permute(W.PROFILE.NO2_profile, [3 2 1]); % no2 in these files is arrange bottom_top, south_north, west_east; it needs to be the other way here
                xlon = W.PROFILE.Longitude'; % likewise lat and lon need transposed to have west_east in the first dimension
                xlat = W.PROFILE.Latitude';
            elseif strcmpi(prof_type, 'monthly')
                wrf_file_name = sprintf('WRF_BEHR_monthly_%02d.nc', wrf_date);
                [no2, xlon, xlat] = read_wrf_vars(find_wrf_path('us',prof_type,wrf_date), wrf_file_name, {'no2', 'XLONG', 'XLAT'});
            elseif strcmpi(prof_type, 'daily')
                wrf_file_name = sprintf('wrfout_d01_%s', datestr(wrf_date, 'yyyy-mm-dd_HH-MM-SS'));
                [no2, xlon, xlat] = read_wrf_vars(find_wrf_path('us',prof_type,wrf_date), wrf_file_name, {'no2', 'XLONG', 'XLAT'});
            else
                E.notimplemented('Loading WRF files for mode %s', prof_type);
            end
            
            no2 = no2(west_east_inds, north_south_inds, bottom_top_ind);
            xlon = xlon(west_east_inds, north_south_inds);
            xlat = xlat(west_east_inds, north_south_inds);
        end
        
        function locs = read_locs_file()
            locs_file = fullfile(misc_behr_v3_validation.my_dir, 'Workspaces', 'trend_locations.nc');
            ni = ncinfo(locs_file);
            ncvarnames = {ni.Variables.Name};
            struct_tmp_cell = cell(1, 2*numel(ncvarnames));
            for a=1:numel(ncvarnames)
                struct_tmp_cell{(a-1)*2 + 1} = ncvarnames{a};
                val = ncread(locs_file, ncvarnames{a});
                if ischar(val)
                    val = cellstr(val);
                else
                    val = num2cell(val);
                end
                
                struct_tmp_cell{a*2} = val;
            end
            
            % This makes it into a structure where each index is a location
            locs = struct(struct_tmp_cell{:});
        end
        
        function [xx, yy] = find_loc_indices(loc, lon, lat, radius)
            % LOC must be a scalar element of the locations structure, LON
            % and LAT must be 2D arrays of longitude and latitude
            % coordinates for an NO2 average or similar 2D field. RADIUS
            % must be a scalar number of grid cells in each direction to
            % get. If omitted, defaults to 0.
            E = JLLErrors;
    
            if ~exist('radius', 'var')
                radius = 0;
            end 
    
            sz = size(lon);
            
            r = sqrt((lon - loc.Longitude).^2 + (lat - loc.Latitude).^2);
            [~, i_min] = min(r(:));
            [xx, yy] = ind2sub(size(lon), i_min); 
    
            xx = (xx - radius):(xx + radius);
            xx = xx(xx > 0 & xx <= sz(1));
            yy = (yy - radius):(yy + radius);
            yy = yy(yy > 0 & yy <= sz(2));
        end 

        
        %%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions %
        %%%%%%%%%%%%%%%%%%%%%%
        
        function plot_one_insitu_comparison(varargin)
            allowed_vars = {'air_no2', 'omi_no2', 'behr_no2'};
            
            labels = struct('air_no2', 'Aircraft NO_2 VCD (molec. cm^{-2})',...
                'omi_no2', 'NASA NO_2 VCD (molec. cm^{-2})',...
                'behr_no2', 'BEHR NO_2 VCD (molec. cm^{-2})');
            
            misc_behr_v3_validation.plot_one_vcd_comparison(misc_behr_v3_validation.profile_comp_file, allowed_vars, labels, varargin{:});
        end
        
        function plot_one_gcas_comparison(varargin)
            allowed_vars = {'GCAS_NO2vcd', 'ColumnAmountNO2Trop', 'BEHRColumnAmountNO2Trop'};
            
            labels = struct('GCAS_NO2vcd', 'GCAS NO_2 VCD (molec. cm^{-2})',...
                'ColumnAmountNO2Trop', 'NASA NO_2 VCD (molec. cm^{-2})',...
                'BEHRColumnAmountNO2Trop', 'BEHR NO_2 VCD (molec. cm^{-2})');
            
            misc_behr_v3_validation.plot_one_vcd_comparison(misc_behr_v3_validation.gcas_vec_comp_file, allowed_vars, labels, varargin{:});
        end
        
    
        function plot_one_vcd_comparison(comp_file, allowed_vars, labels, x_var, y_var, prof_mode, campaigns, time_range)
            % This should be called from another method in this class that
            % provides the right comp_file, allowed_vars, and labels.
            %
            % Plots one v2-v3 comparison against in situ aircraft data
            % derived columns. prof_mode may be 'monthly', 'daily' or
            % 'both'. campaign and time_range must match the fields in the
            % misc_behr_v3_validation.profile_comp_file. If omitted these
            % inputs will be asked interactively
            E = JLLErrors;
            
            comp_struct = load(comp_file);
            
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
            if regcmpi(x_var, 'behr') || regcmpi(y_var, 'behr')
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
            
            allowed_times = fieldnames(comp_struct.v2.us_monthly.(allowed_campaigns{1}));
            if ~exist('time_range', 'var')
                time_range = ask_multichoice('Which time range to use?', allowed_times, 'list', true);
            elseif ~ischar(time_range) || ~ismember(time_range, allowed_times)
                E.badinput('TIME_RANGE must be one of the strings: %s', strjoin(allowed_times, ', '));
            end
            
            if regcmpi(x_var, 'behr') || regcmpi(y_var, 'behr')
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
        
        function plot_one_wrf_comparison(prof_types, campaign, prof_number, uncert_type)
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            allowed_prof_nums = fieldnames(wrfcomp.monthly.(campaign));
            if ~exist('prof_number', 'var')
                prof_number = ask_multichoice('Which profile?', allowed_prof_nums, 'list', true);
            else
                if isnumeric(prof_number) && isscalar(prof_number)
                    prof_number = sprintf('p%d', prof_number);
                elseif ~ischar(prof_number)
                    E.badinput('PROF_NUMBER must be a string or a scalar number')
                end
                
                if ~ismember(prof_number, allowed_prof_nums)
                    E.badinput('PROF_NUMBER must be one of these strings or the numeric component of it: %s', strjoin(allowed_prof_nums));
                end
            end
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_types', 'var')
                prof_types = ask_multiselect('Which profile types to include?', allowed_prof_types);
            elseif ~ismember(prof_types, allowed_prof_types)
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            allowed_uncert_types = {'raw','std', 'both', 'none'};
            if ~exist('uncert_type', 'var')
                uncert_type = ask_multichoice('How to show the variability in the profiles?', allowed_uncert_types, 'default', 'std');
            elseif ~ismember(uncert_type, allowed_uncert_types)
                E.badinput('UNCERT_TYPE must be one of: %s', strjoin(allowed_uncert_types, ', '));
            end
            
            % Do the actual plotting. First get each profile type's raw
            % data, and bin it to give the average shape. Convert from ppm
            % to ppt (hence the 1e6)
            for a=1:numel(prof_types)
                if a == 1
                    % The aircraft data should be the same for all profile
                    % types
                    no2.aircraft = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.data.no2*1e6;
                    pres.aircraft = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.data.pres;
                    [binned_no2.aircraft, binned_pres.aircraft, binned_no2_std.aircraft] = bin_omisp_pressure(pres.aircraft, no2.aircraft, 'mean');
                end
                
                no2.(prof_types{a}) = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.wrf.no2*1e6;
                pres.(prof_types{a}) = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.wrf.pres;
                [binned_no2.(prof_types{a}), binned_pres.(prof_types{a}), binned_no2_std.(prof_types{a})] = bin_omisp_pressure(pres.(prof_types{a}), no2.(prof_types{a}), 'mean');
            end
            
            % We'll always use the same colors for the various profile
            % types. Use lighter versions of the colors for the raw data
            plt_cols = misc_behr_v3_validation.plot_colors;
            
            figure;
            fns = fieldnames(no2);
            l = gobjects(numel(fns),1);
            legstr = cell(1, numel(fns));
            for a=1:numel(fns)
                if any(strcmpi(uncert_type, {'raw', 'both'}))
                    line(no2.(fns{a}), pres.(fns{a}), 'linestyle', 'none', 'marker', '.', 'color', plt_cols.(fns{a}).raw);
                end
                if any(strcmpi(uncert_type, {'std', 'both'}))
                    scatter_errorbars(binned_no2.(fns{a}), binned_pres.(fns{a}), binned_no2_std.(fns{a}), 'color', plt_cols.(fns{a}).avg, 'direction', 'x');
                end
                
                l(a) = line(binned_no2.(fns{a}), binned_pres.(fns{a}), 'color', plt_cols.(fns{a}).avg, 'linewidth', 2);
                legstr{a} = capitalize_words(fns{a});
            end
            
            legend(l, legstr);
            set(gca,'ydir','reverse','fontsize',14);
            xlabel('[NO_2] (pptv)');
            ylabel('Pressure (hPa)');
        end
        
        function plot_aircraft_on_wrf(prof_type, campaign, prof_number, wrf_level)
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            allowed_prof_nums = fieldnames(wrfcomp.monthly.(campaign));
            if ~exist('prof_number', 'var')
                prof_number = ask_multichoice('Which profile?', allowed_prof_nums, 'list', true);
            else
                if isnumeric(prof_number) && isscalar(prof_number)
                    prof_number = sprintf('p%d', prof_number);
                elseif ~ischar(prof_number)
                    E.badinput('PROF_NUMBER must be a string or a scalar number')
                end
                
                if ~ismember(prof_number, allowed_prof_nums)
                    E.badinput('PROF_NUMBER must be one of these strings or the numeric component of it: %s', strjoin(allowed_prof_nums));
                end
            end
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_type', 'var')
                prof_type = ask_multichoice('Which profile types to include?', allowed_prof_types);
            elseif ~ismember(prof_type, allowed_prof_types)
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            if ~exist('wrf_level', 'var')
                wrf_level = ask_number('Which WRF level to plot (0-29)? 0 will plot the lowest level that matched the aircraft.', 'testfxn', @(x) isscalar(x) && x >= 0 && x <= 29, 'testmsg', 'Must be between 0 and 29');
            end
            
            % Okay, this one's a little different. We need to look at the
            % date used for the WRF file (which is just a month number for
            % v2 and monthly) and load that WRF file, then extract the
            % subset of NO2 data to plot with PCOLOR() behind the aircraft
            % data to plot with SCATTER() on top of it.
            match = wrfcomp.(prof_type).(campaign).(prof_number).match;
            buffer = 3;
            we_inds = (min(match.indicies.west_east)-buffer):(max(match.indicies.west_east)+buffer);
            sn_inds = (min(match.indicies.south_north)-buffer):(max(match.indicies.south_north)+buffer);
            if wrf_level == 0
                wrf_level = min(match.indicies.bottom_top);
            end
            
            % Convert from ppm to ppt
            conv = 1e6;
            
            [wrf_no2, wrf_lon, wrf_lat] = misc_behr_v3_validation.load_wrf_no2(prof_type, match.wrf.time, we_inds, sn_inds, wrf_level);
            wrf_no2 = wrf_no2 * conv;
            
            % To make pcolor match up with the actual positioning of the
            % WRF grid cells better, we'll use the bottom left corner,
            % because pcolor would plot wrf_no2(i,j) with corners
            % lon(i,j):lon(i+1,j+1) and likewise for lat.
            [wrf_loncorn, wrf_latcorn] = wrf_grid_corners(wrf_lon, wrf_lat);
            wrf_loncorn = squeeze(wrf_loncorn(1,:,:));
            wrf_latcorn = squeeze(wrf_latcorn(1,:,:));
            
            % Do 10-second averaging to reduce the number of points to plot
            n_sec_per_avg = 10;
            air_lon = avg_n_elements(match.data.lon, n_sec_per_avg, 'op', 'nanmean');
            air_lat = avg_n_elements(match.data.lat, n_sec_per_avg, 'op', 'nanmean');
            air_no2 = avg_n_elements(match.data.no2, n_sec_per_avg, 'op', 'nanmean');
            air_pres = avg_n_elements(match.data.pres, n_sec_per_avg, 'op', 'nanmean');
            
            air_no2 = air_no2 * conv;
            
            % Sort everything so that the highest altitudes (lowest
            % pressure) are plotted first, so that the lower altitudes are
            % plotted on top
            [air_pres, xx_sort] = sort(air_pres);
            air_lon = air_lon(xx_sort);
            air_lat = air_lat(xx_sort);
            air_no2 = air_no2(xx_sort);
            
            % Scale pressure to make the scatter a reasonable size.
            size_from_pres = scale_to_range(air_pres, [50 200]);   
            
            figure;
            pcolor(wrf_loncorn, wrf_latcorn, wrf_no2);
            hold on;
            scatter(air_lon, air_lat, size_from_pres, air_no2, 'filled', 'markeredgecolor','k','linewidth',0.5);
            colorbar;
            set(gca,'fontsize',14);
            title(sprintf('%s: Profile %s vs. %s WRF', upper(strrep(campaign, '_', '-')), strrep(prof_number, 'p', ''), prof_type));
        end
        
        function plot_one_site_comparison(prof_types, campaign, site_number, uncert_type)
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            prof_fns = fieldnames(wrfcomp.monthly.(campaign));
            tmp = strrep(prof_fns, 'p', '');
            site_numbers = floor(str2double(tmp)/1000); % DISCOVER profile numbers are either snnn or sssnnn where s or sss is the site number and nnn the profile number at that location
            allowed_site_numbers_strs = cellfun(@num2str, num2cell(unique(site_numbers)), 'uniformoutput', false);
            if ~exist('site_number', 'var')
                site_number = ask_multichoice('Which profile?', allowed_site_numbers_strs, 'list', true);
                site_number = str2double(site_number);
            else
                if ~isnumeric(site_number) || ~isscalar(site_number)
                    E.badinput('SITE_NUMBER must be a scalar number')
                end
                
                if ~ismember(site_number, site_numbers)
                    E.badinput('SITE_NUMBER must be one of: %s', strjoin(allowed_site_numbers_strs));
                end
            end
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_types', 'var')
                prof_types = ask_multiselect('Which profile types to include?', allowed_prof_types);
            elseif ~ismember(prof_types, allowed_prof_types)
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            allowed_uncert_types = {'raw','std', 'both', 'none'};
            if ~exist('uncert_type', 'var')
                uncert_type = ask_multichoice('How to show the variability in the profiles?', allowed_uncert_types, 'default', 'std');
            elseif ~ismember(uncert_type, allowed_uncert_types)
                E.badinput('UNCERT_TYPE must be one of: %s', strjoin(allowed_uncert_types, ', '));
            end
            
            %%% PLOTTING %%%
            conv = 1e6; % convert NO2 from ppm to ppt
            
            % First make a list of all profiles for a given site then
            % concatenate all of the matched data. We will then bin that
            % and plot it.
            xx = site_numbers == site_number;
            prof_fns = prof_fns(xx);
            for a=1:numel(prof_types)
                if a == 1
                    no2.aircraft = [];
                    pres.aircraft = [];
                end
                no2.(prof_types{a}) = [];
                pres.(prof_types{a}) = [];
                
                for b=1:numel(prof_fns)
                    if a==1
                        no2.aircraft = veccat(no2.aircraft, wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.data.no2*conv, 'column');
                        pres.aircraft = veccat(pres.aircraft, wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.data.pres, 'column');
                    end
                    no2.(prof_types{a}) = veccat(no2.(prof_types{a}), wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.wrf.no2*conv, 'column');
                    pres.(prof_types{a}) = veccat(pres.(prof_types{a}), wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.wrf.pres, 'column');
                end
                
                % Bin
                if a==1
                    [binned_no2.aircraft, binned_pres.aircraft, binned_no2_std.aircraft] = bin_omisp_pressure(pres.aircraft, no2.aircraft, 'mean');
                end
                [binned_no2.(prof_types{a}), binned_pres.(prof_types{a}), binned_no2_std.(prof_types{a})] = bin_omisp_pressure(pres.(prof_types{a}), no2.(prof_types{a}), 'mean');
            end
            
            plt_cols = misc_behr_v3_validation.plot_colors;
            
            fns = fieldnames(no2);
            l = gobjects(numel(fns),1);
            legstr = cell(1, numel(fns));
            figure;
            for a=1:numel(fns)
                l(a) = line(binned_no2.(fns{a}), binned_pres.(fns{a}), 'color', plt_cols.(fns{a}).avg, 'linewidth', 2);
                legstr{a} = capitalize_words(fns{a});
            end
            set(gca,'ydir','reverse','fontsize',14);
            legend(l, legstr);
            title(sprintf('%s: Site %d', upper(strrep(campaign,'_','-')), site_number));
        end
        
        function plot_all_sites_comparison(prof_types, campaign, uncert_type)
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            prof_fns = fieldnames(wrfcomp.monthly.(campaign));
            tmp = strrep(prof_fns, 'p', '');
            site_numbers = unique(floor(str2double(tmp)/1000)); % DISCOVER profile numbers are either snnn or sssnnn where s or sss is the site number and nnn the profile number at that location
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_type', 'var')
                prof_types = ask_multiselect('Which profile types to include?', allowed_prof_types);
            elseif ~ismember(prof_types, allowed_prof_types)
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            allowed_uncert_types = {'raw','std', 'both', 'none'};
            if ~exist('uncert_type', 'var')
                uncert_type = ask_multichoice('How to show the variability in the profiles?', allowed_uncert_types, 'default', 'std');
            elseif ~ismember(uncert_type, allowed_uncert_types)
                E.badinput('UNCERT_TYPE must be one of: %s', strjoin(allowed_uncert_types, ', '));
            end
            
            for a=1:numel(site_numbers)
                misc_behr_v3_validation.plot_one_site_comparison(prof_types, campaign, site_numbers(a), uncert_type);
            end
        end
        
        function plot_scd_vs_wrf_columns(location_name, start_date, end_date)
            % Plots OMI tropospheric SCDs, WRF monthly, and WRF daily
            % columns for a given date range to show whether WRF is
            % capturing the wind direction correctly. If daily BEHR files
            % do not exist for a date in the range given, that day will
            % just be skipped.
            locs = misc_behr_v3_validation.read_locs_file();
            loc_names = {locs.ShortName};
            if ~exist('location_name', 'var')
                loc_ind = ask_multichoice('Choose a location', loc_names, 'list', true, 'index', true);
            elseif ~ismember(location_name, loc_names)
                E.badinput('LOCATION_NAME must be one of: %s', strjoin(loc_names, ', '));
            else
                loc_ind = strcmpi(location_name, loc_names);
            end
            
            plot_loc = locs(loc_ind);
            
            if ~exist('start_date', 'var')
                start_date = datenum(ask_date('Enter the start date'));
            else
                start_date = validate_date(start_date);
            end
            
            if ~exist('end_date', 'var')
                end_date = datenum(ask_date('Enter the end date'));
            else
                end_date = validate_date(end_date);
            end
            
            dvec = start_date:end_date;
            % Struct that holds the name of monthly file read in and its
            % column data
            monthly_file = '';
            
            for d=1:numel(dvec)
                % Load the daily BEHR file, if it exists
                try
                    Data = load_behr_file(dvec(d), 'daily', 'us');
                catch err
                    if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
                        fprintf('No BEHR daily profile file available for %s\n', datestr(dvec(d)));
                        continue
                    else
                        rethrow(err)
                    end
                end
                
                wrf_int_mode = 'box';
                
                % Load a new monthly file, if needed
                wrf_monthly_file_tmp = fullfile(find_wrf_path('us','monthly',dvec(d)), sprintf('WRF_BEHR_monthly_%02d.nc', month(dvec(d))));
                if ~strcmp(monthly_file, wrf_monthly_file_tmp)
                    monthly_file = wrf_monthly_file_tmp;
                    monthly_wrf_no2_vcds = compute_wrf_trop_columns(monthly_file, wrf_int_mode, 200);
                    monthly_wrf_lon = ncread(monthly_file, 'XLONG');
                    monthly_wrf_lat = ncread(monthly_file, 'XLAT');
                end
                
                % Get the WRF file name, but just retrieve the date b/c
                % files produced on the cluster will have different paths
                % and use the subset files, which aren't stored locally.
                [~, daily_wrf_file_tmp] = fileparts(Data(1).BEHRWRFFile);
                wrf_date = date_from_wrf_filenames(daily_wrf_file_tmp);
                daily_file = fullfile(find_wrf_path('us','daily',wrf_date), sprintf('wrfout_d01_%s', datestr(wrf_date, 'yyyy-mm-dd_HH-MM-SS')));
                daily_wrf_no2_vcds = compute_wrf_trop_columns(daily_file, wrf_int_mode, 200);
                daily_wrf_lon = ncread(daily_file, 'XLONG');
                daily_wrf_lat = ncread(daily_file, 'XLAT');
                
                for a=1:numel(Data)
                    % Loop over every swath. If the location is in the box
                    % defined by the four corner pixel centers of the swath
                    % plot it
                    corner_x = [Data(a).Longitude(1,1), Data(a).Longitude(1,end), Data(a).Longitude(end,end), Data(a).Longitude(end,1)];
                    corner_y = [Data(a).Latitude(1,1), Data(a).Latitude(1,end), Data(a).Latitude(end,end), Data(a).Latitude(end,1)];
                    if ~inpolygon(plot_loc.Longitude, plot_loc.Latitude, corner_x, corner_y);
                        continue
                    end
                    
                    badpix = mod(Data(a).BEHRQualityFlags, 2) ~= 0;
                    Data(a).BEHRColumnAmountNO2Trop(badpix) = NaN;
                    
                    plot_radius = 10; % number of grid cells or pixels around the location to plot
                    [xx_sat, yy_sat] = misc_behr_v3_validation.find_loc_indices(plot_loc, Data(a).Longitude, Data(a).Latitude, plot_radius);
                    [xx_monthly, yy_monthly] = misc_behr_v3_validation.find_loc_indices(plot_loc, monthly_wrf_lon, monthly_wrf_lat, plot_radius);
                    [xx_daily, yy_daily] = misc_behr_v3_validation.find_loc_indices(plot_loc, daily_wrf_lon, daily_wrf_lat, plot_radius);
                    
                    % Make the plot, 3 side-by-side figures
                    fig=figure;
                    fig.Position(3) = fig.Position(3)*2;
                    for p=1:3
                        subplot(1,3,p);
                        if p==1
                            % Convert VCD back to SCD
                            pcolor(squeeze(Data(a).FoV75CornerLongitude(1,xx_sat,yy_sat)), squeeze(Data(a).FoV75CornerLatitude(1,xx_sat,yy_sat)), Data(a).BEHRColumnAmountNO2Trop(xx_sat,yy_sat) .* Data(a).BEHRAMFTrop(xx_sat, yy_sat));
                            title('OMI Trop. SCD');
                        elseif p==2
                            pcolor(monthly_wrf_lon(xx_monthly, yy_monthly), monthly_wrf_lat(xx_monthly, yy_monthly), monthly_wrf_no2_vcds(xx_monthly, yy_monthly));
                            title('Monthly WRF');
                        elseif p==3
                            pcolor(daily_wrf_lon(xx_monthly, yy_monthly), daily_wrf_lat(xx_monthly, yy_monthly), daily_wrf_no2_vcds(xx_daily, yy_daily));
                            title('Daily WRF');
                        else
                            E.notimplemented('>3 plots')
                        end
                        
                        %caxis([0 1e16]);
                        colorbar;
                        line(plot_loc.Longitude, plot_loc.Latitude, 'linestyle', 'none', 'marker', 'p', 'color', 'w', 'markersize',16,'linewidth',2);
                        set(gca,'fontsize',14);
                    end
                end
            end
        end
        
        
        function results = record_scd_comparison()
            nsites = 20;
            ndays = 20;
            dvec = datenum('2012-01-01'):datenum('2012-12-31');
            locs = misc_behr_v3_validation.read_locs_file;
            % Remove the rural sites
            locs = locs(1:70);
            
            results = struct('loc_name', '', 'date', '', 'user_value', [], 'user_note', '');
            results = repmat(results, nsites*ndays, 1);
            
            eval_opts = {'Daily - good agreement with SCD', 'Daily - bad agreement with SCD', 'Similar to monthly - good agreement with SCD', 'Similar to monthly - bad agreement with SCD', 'Not enough data'};
            ned_ind = 5; % index for "Not enough data", if this is the response, we don't want to store the result.
            for a=1:nsites
                % Make a copy so that we can remove days as we test them
                isite = randi(numel(locs), 1);
                site_dvec = dvec;
                b = 1;
                while b <= ndays && ~isempty(site_dvec)
                    idate = randi(numel(site_dvec), 1);
                    misc_behr_v3_validation.plot_scd_vs_wrf_columns(locs(isite).ShortName, dvec(idate), dvec(idate));
                    try
                        user_ans = ask_multichoice('Evaluate this day', eval_opts, 'list', true, 'index', true);
                    catch err
                        if strcmp(err.identifier, 'ask_multichoice:user_cancel')
                            % If we cancel the work, save the results
                            % completed so far.
                            save(misc_behr_v3_validation.scd_comp_file(true), 'results');
                            return
                        else
                            rethrow(err)
                        end
                    end
                    % Whether or not we use this day, we don't want to
                    % repeat it, and we want to close the figures
                    close all
                    this_date = site_dvec(idate);
                    site_dvec(idate) = [];
                    if user_ans == ned_ind
                        % If the user responded "Not enough data", then
                        % don't store any result.
                        continue
                    end
                    
                    iresult = sub2ind([ndays, nsites], b, a);
                    results(iresult).loc_name = locs(isite).ShortName;
                    results(iresult).date = datestr(this_date);
                    results(iresult).user_value = user_ans;
                    results(iresult).user_note = input('Enter a note if you wish: ', 's');
                    b = b + 1;
                    fprintf('%d of %d completed...\n', iresult, nsites*ndays);
                end
                
                locs(isite) = [];
                
            end
            
            save(misc_behr_v3_validation.scd_comp_file, 'results');
        end
    end
    
end

