classdef misc_behr_update_plots
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        start_date = '2012-01-01';
        end_date = '2012-12-31';
        behr_modis_quality_best_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/6b-MODISQualityBest';
        behr_modis_quality_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/6-MODISQuality';
        behr_final_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/5-PSM';
        behr_final_cvm_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/5b-CVM';
        behr_nasa_brdf_vis_profs_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/4-NewProfs';
        behr_nasa_brdf_vis_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/3-NewVis';
        behr_nasa_brdf_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/2-BRDF';
        behr_nasa_only_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/1-NASAv3';
        behr_v2_1C_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/0-v2.1C';
        
        figs_root_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Figures';
    end
    
    methods(Static = true)
        function make_behr_modis_quality_best(do_overwrite)
            % Produces the final version, but with MODIS albedo quality
            % restricted to 2 or better.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            % Verify that the most up-to-date PSM code is active.
            G = GitChecker;
            G.addReqCommits(behr_paths.psm_dir, '7bd02b9');
            G.addReqCommits(behr_paths.python_interface, 'a217fcd');
            G.Strict = true;
            G.checkState();
            
            save_dir = misc_behr_update_plots.behr_modis_best_quality_dir;
            misc_behr_update_plots.make_behr_with_parameters('bbdc057', 'afd69ac', save_dir, do_overwrite, true);
        end
        
        function make_behr_modis_quality(do_overwrite)
            % Produces the final version, but with MODIS albedo quality
            % restricted to 2 or better.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            % Verify that the most up-to-date PSM code is active.
            G = GitChecker;
            G.addReqCommits(behr_paths.psm_dir, '7bd02b9');
            G.addReqCommits(behr_paths.python_interface, 'a217fcd');
            G.Strict = true;
            G.checkState();
            
            save_dir = misc_behr_update_plots.behr_modis_quality_dir;
            misc_behr_update_plots.make_behr_with_parameters('351da88', 'afd69ac', save_dir, do_overwrite, true);
        end
        
        function make_behr_final(do_overwrite)
            % Produces the final version, that includes:
            %   1) NASA SP v3 base
            %   2) BRDF land albedo
            %   3) COARTS sea albedo
            %   4) The new visible-only AMF
            %   5) PSM gridding
            %   6) New profiles (monthly and daily)
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            % Verify that the most up-to-date PSM code is active.
            G = GitChecker;
            G.addReqCommits(behr_paths.psm_dir, 'a861905');
            G.addReqCommits(fullfile(behr_repo_dir, '..', 'Python Interface'), 'a217fcd');
            G.Strict = true;
            G.checkState();
            
            save_dir = misc_behr_update_plots.behr_final_dir;
            misc_behr_update_plots.make_behr_with_parameters('ddb101b', save_dir, do_overwrite, true);
        end
        
        function make_behr_final_only_cvm(do_overwrite)
            G = GitChecker;
            G.addReqCommits(behr_paths.psm_dir, '49f77fc');
            G.checkState();
            
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            sp_dir = fullfile(misc_behr_update_plots.behr_final_dir, 'SP_Files');
            save_dir = misc_behr_update_plots.behr_final_cvm_dir;
            misc_behr_update_plots.make_behr_with_parameters('d71805e', save_dir, do_overwrite, true, sp_dir);
        end
        
        function make_behr_brdf_nasa_newvis_profs(do_overwrite)
            % Produces the version of BEHR that has:
            %   1) the NASA SP v3 base,
            %   2) the BRDF albedo,
            %   3) the COARTS sea albedo
            %   4) the new visible-only AMF
            %   5) the new profiles
            % but not the PSM gridding.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir;
            misc_behr_update_plots.make_behr_with_parameters('11177fd', save_dir, do_overwrite, true);
        end
        
        function make_behr_brdf_nasa_newvis(do_overwrite)
            % Produces the version of BEHR that has:
            %   1) the NASA SP v3 base,
            %   2) the BRDF albedo,
            %   3) the COARTS sea albedo
            %   4) the new visible-only AMF
            % but not the new profiles or the PSM gridding.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_brdf_vis_dir;
            misc_behr_update_plots.make_behr_with_parameters('d64a80f1', save_dir, do_overwrite, false);
        end
        
        function make_behr_brdf_and_nasa(do_overwrite)
            % Produces the version of BEHR that has:
            %   1) the NASA SP v3 base,
            %   2) the BRDF albedo
            %   3) the COARTS sea albedo
            % but not the new visible-only AMF, new profiles or the PSM
            % gridding.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_brdf_dir;
            misc_behr_update_plots.make_behr_with_parameters('2018ae4', save_dir, do_overwrite, false);
        end
        
        function make_behr_only_new_nasa(do_overwrite)
            % Produces the version of BEHR that has only the NASA SP v3
            % base. Everything else should be the same as in version 2,
            % except that we now use the NASA SP v3 as the base.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_only_dir;
            misc_behr_update_plots.make_behr_with_parameters('d889411', save_dir, do_overwrite, false);
        end
        
        
        function average_no2_vcds()
           
            do_overwrite = ask_yn('Overwrite existing average files?');
 
            % Structure array that describes how each incremental change
            % should be averaged. 3 fields: 'dir' is the root directory of
            % the outputs (which itself must contain the "MonthlyProfs" and
            % "DailyProfs" directories), 'use_new_avg', a boolean that, if
            % true will use psm_time_average instead of
            % no2_column_map_2014, and 'data_fields', a cell array of data
            % fields to average for that directory.
            average_info = cat(1,...
                struct('dir', misc_behr_update_plots.behr_v2_1C_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','MODISAlbedo'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_only_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','MODISAlbedo'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'MODISAlbedo'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_final_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_final_cvm_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','MODISAlbedo'}}),...
                struct('dir', misc_behr_update_plots.behr_modis_quality_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','MODISAlbedo'}})...
            );
        
            G_new_avg = GitChecker;
            G_new_avg.Strict = false;
            G_old_avg = GitChecker;
            G_old_avg.Strict = false;
            
            G_new_avg.addReqCommits(behr_analysis_repo_dir,'fffce59');
            G_old_avg.addCommitRange(behr_analysis_repo_dir,'f4dd1e4','f4dd1e4');
            
            do_new_avg = G_new_avg.checkState();
            do_old_avg = G_old_avg.checkState();
            
            if ~do_new_avg && ~do_old_avg
                E.callError('git_status', 'Not in the proper commit range to do either the old or new averaging');
            end
            
            for a=1:numel(average_info)
                for b=1:numel(average_info(a).data_fields)
                    if do_new_avg && average_info(a).use_new_avg
                        misc_behr_update_plots.average_behr_old(average_info(a).dir, average_info(a).data_fields{b}, do_overwrite, true);
                    elseif do_old_avg
                        misc_behr_update_plots.average_behr_old(average_info(a).dir, average_info(a).data_fields{b}, do_overwrite, false);
                    end
                end
            end
            
        end
        
        
        function plot_incr_figures()
            % This will need to go into each workspace directory, find the
            % average files, figure out which two files to compare, and
            % plot the difference.
            
            % These must be in the order they are to be compared
            incr_dirs = {misc_behr_update_plots.behr_v2_1C_dir,...
                         misc_behr_update_plots.behr_nasa_only_dir,...
                         misc_behr_update_plots.behr_nasa_brdf_dir,...
                         misc_behr_update_plots.behr_nasa_brdf_vis_dir,...
                         misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir,...
                         misc_behr_update_plots.behr_final_dir,...
                         misc_behr_update_plots.behr_final_cvm_dir,...
                         misc_behr_update_plots.behr_modis_quality_dir};
            
            % Use the "outside.mat" file to clean up the plots. Must be the
            % same length as incr_dirs
            use_outside = [true, true, true, true, false, false, false, false];
            
            outside = load(fullfile(behr_analysis_repo_dir, 'Utils', 'outside.mat'));
            outside = ~logical(outside.outside);
                     
            for a=1:numel(incr_dirs)
            %for a=numel(incr_dirs)
                if a < numel(incr_dirs)
                    base_ind = a;
                    new_ind = a+1;
                else
                    % We'll use the last index to make the
                    % "overall" change figures
                    base_ind = 1;
                    new_ind = a;
                end
                Fbase = dir(fullfile(incr_dirs{base_ind}, '*.mat'));
                base_names = {Fbase.name};
                Fnew = dir(fullfile(incr_dirs{new_ind}, '*.mat'));
                new_names = {Fnew.name};
                for b=1:numel(base_names)
                    if ismember(base_names{b}, new_names)
                        Dbase = load(fullfile(incr_dirs{base_ind}, base_names{b}));
                        Dnew = load(fullfile(incr_dirs{new_ind}, base_names{b}));
                        
                        if isscalar(Dbase.no2_vcds) || isscalar(Dnew.no2_vcds)
                            % psm_time_average returns a scalar NaN if, for
                            % whatever reason, no valid files are found.
                            continue
                        end
                        
                        if use_outside(base_ind)
                            Dbase.no2_vcds(outside) = nan;
                        end
                        if use_outside(new_ind)
                            Dnew.no2_vcds(outside) = nan;
                        end
                        
                        if ~isequal(size(Dnew.no2_vcds), size(Dbase.no2_vcds))
                            % We choose to interpolate to the new
                            % coordinates because the PSM grid is one
                            % smaller in both dimensions than the old
                            % method (due to the interpolation in
                            % no2_column_map). This way we don't have to
                            % worry about extrapolation, since the PSM grid
                            % should be entirely inside the old grid.
                            old_no2 = interp2(Dbase.lon_grid, Dbase.lat_grid, Dbase.no2_vcds, Dnew.lon_grid, Dnew.lat_grid);
                        else
                            old_no2 = Dbase.no2_vcds;
                        end
                        
                        abs_del = Dnew.no2_vcds - old_no2;
                        rel_del = reldiff(Dnew.no2_vcds, old_no2)*100;
                        
                        [title_str, quantity_name, unit_name] = misc_behr_update_plots.get_title_from_file_name(base_names{b}, incr_dirs{base_ind}, incr_dirs{new_ind});
                        save_dir = fullfile(misc_behr_update_plots.figs_root_dir, 'IncrementTests');
                        
                        misc_behr_update_plots.plot_change_avg(Dnew.lon_grid, Dnew.lat_grid, abs_del, false, title_str, quantity_name, unit_name, save_dir);
                        misc_behr_update_plots.plot_change_avg(Dnew.lon_grid, Dnew.lat_grid, rel_del, true, title_str, quantity_name, unit_name, save_dir);
                    end
                end
            end
        end
        
        
        function plot_day_v_month_figures()
            % This will, if available, plot the average difference between
            % the daily and monthly profile retrievals within a single
            % increment.
            
            % Only the increments with the new profiles have daily and
            % monthly profile retrievals
            incr_dirs = {misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir,...
                         misc_behr_update_plots.behr_final_dir};
                     
            for a=1:numel(incr_dirs)
                F = dir(fullfile(incr_dirs{a}, '*Daily.mat'));
                for b=1:numel(F)
                    monthly_file = strrep(F(b).name, 'Daily', 'Monthly');
                    if ~exist(fullfile(incr_dirs{a}, monthly_file), 'file')
                        continue
                    end
                    
                    Ddaily = load(fullfile(incr_dirs{a}, F(b).name));
                    Dmonthly = load(fullfile(incr_dirs{a}, monthly_file));
                    
                    abs_del = Ddaily.no2_vcds - Dmonthly.no2_vcds;
                    rel_del = reldiff(Ddaily.no2_vcds, Dmonthly.no2_vcds)*100;
                    
                    [~, quantity_name, unit_name, time_period] = misc_behr_update_plots.get_title_from_file_name(monthly_file, incr_dirs{a}, incr_dirs{a});
                    title_str = sprintf('Daily vs monthly profiles (%s, %s)', quantity_name, time_period); 
                    save_dir = fullfile(misc_behr_update_plots.figs_root_dir, 'DailyVsMonthly');
                    
                    misc_behr_update_plots.plot_change_avg(Ddaily.lon_grid, Ddaily.lat_grid, abs_del, false, title_str, quantity_name, unit_name, save_dir);
                    misc_behr_update_plots.plot_change_avg(Ddaily.lon_grid, Ddaily.lat_grid, rel_del, true, title_str, quantity_name, unit_name, save_dir);
                end
            end
        end
        
        
        function date_orbit = plot_psm_weights(start_date, end_date, use_daily)
            E = JLLErrors;
            
            if ~exist('start_date','var')
                start_date = datenum(ask_date('Enter the first date to examine weights for'));
            else
                start_date = validate_date(start_date);
            end
            if ~exist('end_date', 'var')
                end_date = datenum(ask_date('Enter the last date to examine weights for'));
            else
                end_date = validate_date(end_date);
            end
            
            if end_date < start_date
                E.badinput('END_DATE must be >= START_DATE');
            end
            
            use_daily = false; % for now, until we have daily profile PSM gridded data to examine
            
            if use_daily
                E.notimplemented('Use daily profile PSM product');
            else
                load_path = fullfile(misc_behr_update_plots.behr_final_dir, 'MonthlyProfs');
            end
            
            datevec = start_date:end_date;
            regex_pattern = 'OMI_BEHR_v\\d-\\d[A-Z](rev\\d)?_%s\\.mat';
            F = dir(fullfile(load_path, '*.mat'));
            file_names = {F.name};
            
            weights = nan(500,1200,numel(datevec)*4);
            no2_columns = nan(500,1200,numel(datevec)*4);
            date_orbit = cell(1, numel(datevec)*4);
            
            ind = 1;
            
            for d=1:numel(datevec)
                fprintf('Loading %d of %d\n', d, numel(datevec));
                
                file_pattern = sprintf(regex_pattern, datestr(datevec(d), 'yyyymmdd'));
                xx = ~cellfun('isempty', regexp(file_names, file_pattern)); % we want the error to happen if more than one file matches the pattern too
                if sum(xx) ~= 1
                    E.filenotfound(file_pattern);
                end
                
                D = load(fullfile(load_path, file_names{xx}));
                for a=1:numel(D.OMI)
                    if ind == 1
                        lon = D.OMI(a).Longitude;
                        lat = D.OMI(a).Latitude;
                    end
                    no2_columns(:,:,ind) = D.OMI(a).BEHRColumnAmountNO2Trop;
                    weights(:,:,ind) = D.OMI(a).BEHRColumnAmountNO2TropWeights;
                    date_orbit{ind} = sprintf('%s-Orbit%d', datestr(datevec(d),'yyyy-mm-dd'),a);
                    ind = ind+1;
                end
            end
            
            no2_columns(:,:,ind:end) = [];
            weights(:,:,ind:end) = [];
            
            [slon, slat] = state_outlines('not', 'ak', 'hi');
            plot_slice_gui(no2_columns, lon, lat, slon, slat, 'NO2');
            plot_slice_gui(weights, lon, lat, slon, slat, 'Weights');
        end
    end
    
    methods(Static = true, Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper functions for the production of the actual incr. changes %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function make_behr_with_parameters(req_commit, req_utils_commit, root_save_dir, do_overwrite, has_daily_profs, alt_sp_dir)
            % Runs BEHR, saving files in ROOT_SAVE_DIR and overwriting
            % according to DO_OVERWRITE. This runs BEHR assuming the
            % BEHR_main and read_omno2_v_aug2012 functions take all the
            % directories as parameters. ALT_SP_DIR may be omitted, if so,
            % this assumes that it needs to produce the SP files and place
            % them in fullfile(ROOT_SAVE_DIR, 'SP_Files'). If not omitted,
            % then it will assume it does not need to produce the SP files
            % and that they can be found in ALT_SP_DIR.
            E = JLLErrors;
            my_dir = fileparts(mfilename('fullpath'));
            behr_dir = behr_paths.behr_core;
            behr_utils = behr_paths.behr_utils;
            
            G = GitChecker;
            G.addCommitRange(behr_dir, req_commit, req_commit);
            G.addCommitRange(behr_utils, req_utils_commit, req_utils_commit);
            G.Strict = true;
            
            G.checkState();
            
            if ~exist(root_save_dir, 'dir')
                E.dir_dne(root_save_dir);
            end
            
            sp_save_dir = fullfile(root_save_dir, 'SP_Files');
            if ~exist(sp_save_dir, 'dir')
                mkdir(sp_save_dir);
            end
            
            if ~exist('alt_sp_dir', 'var')
                do_sp = true;
                alt_sp_dir = sp_save_dir;
            else
                if ~ischar(alt_sp_dir)
                    E.badinput('"alt_sp_dir" must be a string');
                elseif ~exist(alt_sp_dir, 'dir')
                    E.dir_dne('SP load directory "%s" does not exist', alt_sp_dir);
                end
                do_sp = false;
            end
            
            if do_sp
                read_omno2_v_aug2012('start', misc_behr_update_plots.start_date,...
                    'end', misc_behr_update_plots.end_date,...
                    'sp_mat_dir', sp_save_dir, 'overwrite', do_overwrite);
            else
                sp_files = dir(fullfile(alt_sp_dir,'*.mat'));
                sp_files = {sp_files.name};
                datevec = datenum(misc_behr_update_plots.start_date):datenum(misc_behr_update_plots.end_date);
                n_missing_files = 0;
                for d=1:numel(datevec)
                    this_pattern = sprintf('OMI_SP_v\\d-\\d[A-Z](rev\\d)?_%s.mat', datestr(datevec(d), 'yyyymmdd'));
                    n_missing_files = (sum(~cellfun('isempty', regexp(sp_files, this_pattern, 'once'))) ~= 1) + n_missing_files;
                end
                
                if n_missing_files > 0
                    E.filenotfound('%d OMI_SP .mat files missing from %s', n_missing_files, alt_sp_dir);
                end
            end
            
            % For the monthly profiles, we rely on the default profile
            % choice to be monthly. This avoids issues with older
            % increments that do not have the 'profile_mode' parameter.
            monthly_save_dir = fullfile(root_save_dir, 'MonthlyProfs');
            if ~exist(monthly_save_dir, 'dir')
                mkdir(monthly_save_dir);
            end
            BEHR_main('start', misc_behr_update_plots.start_date,...
                'end', misc_behr_update_plots.end_date,...
                'sp_mat_dir', alt_sp_dir,...
                'behr_mat_dir', monthly_save_dir,...
                'overwrite', do_overwrite);
            
            if has_daily_profs
                daily_save_dir = fullfile(root_save_dir, 'DailyProfs');
                if ~exist(daily_save_dir, 'dir')
                    mkdir(daily_save_dir);
                end
                BEHR_main('start', misc_behr_update_plots.start_date,...
                    'end', misc_behr_update_plots.end_date,...
                    'sp_mat_dir', alt_sp_dir,...
                    'behr_mat_dir', daily_save_dir,...
                    'profile_mode', 'daily',...
                    'overwrite', do_overwrite);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper functions for the incremental change averages %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function average_behr_old(data_dir, data_field, overwrite, use_new_avg)
            djf_start = {'2012-01-01', '2012-12-01'};
            djf_end = {'2012-02-29', '2012-12-31'};
            jja_start = '2012-06-01';
            jja_end = '2012-08-31';
            
            lon_lim = [-125 -65];
            lat_lim = [25 50];
            
            monthly_dir = fullfile(data_dir, 'MonthlyProfs');
            daily_dir = fullfile(data_dir, 'DailyProfs');
            save_name = fullfile(data_dir, sprintf('%s-DJF-Monthly.mat', data_field));
            
            if ~overwrite && exist(save_name, 'file')
                fprintf('%s exists\n', save_name);
            else
                if use_new_avg
                    [no2_vcds, lon_grid, lat_grid] = psm_time_average(djf_start, djf_end, lon_lim, lat_lim,...
                        'avgfield', data_field, 'behr_dir', monthly_dir);
                    count_grid = [];
                    avg_config = [];
                else
                    [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(djf_start, djf_end, lon_lim, lat_lim,...
                        'mapfield', data_field, 'behrdir', monthly_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                        'makefig', false);
                end
                fprintf('Saving as %s\n', save_name);
                save(save_name, 'no2_vcds', 'lon_grid', 'lat_grid', 'count_grid', 'avg_config');
            end
            
            save_name = fullfile(data_dir, sprintf('%s-JJA-Monthly.mat', data_field));
            if ~overwrite && exist(save_name, 'file')
                fprintf('%s exists\n', save_name);
            else
                if use_new_avg
                    [no2_vcds, lon_grid, lat_grid] = psm_time_average(jja_start, jja_end, lon_lim, lat_lim,...
                        'avgfield', data_field, 'behr_dir', monthly_dir);
                    count_grid = [];
                    avg_config = [];
                else
                    [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(jja_start, jja_end, lon_lim, lat_lim,...
                        'mapfield', data_field, 'behrdir', monthly_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                        'makefig', false);
                end
                
                fprintf('Saving as %s\n', save_name);
                save(save_name, 'no2_vcds', 'lon_grid', 'lat_grid', 'count_grid', 'avg_config');
            end
            
            % Only the new profiles and the version 3 final will have data
            % in the DailyProfs directory
            if exist(fullfile(daily_dir, 'OMI_BEHR_v2-1C_20120101.mat'), 'file')
                save_name = fullfile(data_dir, sprintf('%s-DJF-Daily.mat', data_field));
                if ~overwrite && exist(save_name, 'file')
                    fprintf('%s exists\n', save_name);
                else
                    if use_new_avg
                        [no2_vcds, lon_grid, lat_grid] = psm_time_average(djf_start, djf_end, lon_lim, lat_lim,...
                            'avgfield', data_field, 'behr_dir', daily_dir);
                        count_grid = [];
                        avg_config = [];
                    else
                        [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(djf_start, djf_end, lon_lim, lat_lim,...
                        'mapfield', data_field, 'behrdir', daily_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                        'makefig', false);
                    end
                
                    fprintf('Saving as %s\n', save_name);
                    save(save_name, 'no2_vcds', 'lon_grid', 'lat_grid', 'count_grid', 'avg_config');
                end
                
                save_name = fullfile(data_dir, sprintf('%s-JJA-Daily.mat', data_field));
                if ~overwrite && exist(save_name, 'file')
                    fprintf('%s exists\n', save_name);
                else
                    if use_new_avg
                        [no2_vcds, lon_grid, lat_grid] = psm_time_average(jja_start, jja_end, lon_lim, lat_lim,...
                            'avgfield', data_field, 'behr_dir', daily_dir);
                        count_grid = [];
                        avg_config = [];
                    else
                        [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(jja_start, jja_end, lon_lim, lat_lim,...
                        'mapfield', data_field, 'behrdir', daily_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                        'makefig', false);
                    end
                
                    fprintf('Saving as %s\n', save_name);
                    save(save_name, 'no2_vcds', 'lon_grid', 'lat_grid', 'count_grid', 'avg_config');
                end
            end
        end
        
        function [title_str, quantity_str, unit_str, time_period] = get_title_from_file_name(filename, base_dir, new_dir)
            % Parses one of the average .mat files names to figure out what
            % the title string should be and the colorbar limits.
            E = JLLErrors;
            
            [~, filename] = fileparts(filename); % remove any path or extension
            [~, base_dir] = fileparts(base_dir);
            [~, new_dir] = fileparts(new_dir);
            
            name_parts = strsplit(filename, '-');
            
            % Quantity - shorten it quite a bit so the title isn't insanely
            % long
            switch name_parts{1}
                case 'BEHRColumnAmountNO2Trop'
                    quantity_str = 'NO_2 VCD';
                    unit_str = 'molec. cm^{-2}';
                case 'BEHRColumnAmountNO2TropVisOnly'
                    quantity_str = 'Vis NO_2 VCD';
                    unit_str = 'molec. cm^{-2}';
                case 'MODISAlbedo'
                    quantity_str = 'Surf Refl';
                    unit_str = 'unitless';
                otherwise
                    E.notimplemented('No quantity title for data field "%s"', name_parts{1});
            end
            
            % Time period (DJF or JJA), should always be the second part of
            % the file name
            time_period = name_parts{2};
            
            % Monthly or daily profiles used
            profs_used = sprintf('%s Profs', name_parts{3});
            
            title_str = sprintf('%s vs %s - %s (%s, %s)', quantity_str, new_dir, base_dir, time_period, profs_used);
        end
        
        function plot_change_avg(lon_grid, lat_grid, val_diff, is_rel_diff, title_str, quantity_name, unit_name, save_dir)
            E = JLLErrors;
            
            if is_rel_diff
                save_name = sprintf('Relative difference %s', title_str);
                cb_label = sprintf(tex_in_printf('%%\Delta %s (%s)'), quantity_name, unit_name);
                cb_max_lims = [-100 100];
            else
                save_name = sprintf('Absolute difference %s', title_str);
                cb_label = sprintf(tex_in_printf('\Delta %s'), quantity_name);
                cb_max_lims = [-Inf Inf];
            end
            
            clims = calc_plot_limits(val_diff, 5, 'diff', 'pow10', 'max', cb_max_lims);
            
            fig=figure; pcolor(lon_grid, lat_grid, val_diff);
            shading flat; 
            cb=colorbar; 
            caxis(clims); 
            colormap(blue_red_cmap);
            state_outlines('k', 'not', 'ak', 'hi');
            title(title_str);
            
            
            cb.Label.String = cb_label;
            cb.FontSize = 14;
            
            if exist(misc_behr_update_plots.figs_root_dir, 'dir') && ~exist(save_dir, 'dir')
                mkdir(save_dir)
            elseif ~exist(misc_behr_update_plots.figs_root_dir, 'dir')
                E.dir_dne('The misc_behr_update_plots figures root directory "%s" does not exist', misc_behr_update_plots.figs_root_dir);
            end
            
            savefig(fig, fullfile(save_dir, save_name));
            saveas(fig, fullfile(save_dir, save_name), 'png');
            close(fig);
        end
    end
    
end

