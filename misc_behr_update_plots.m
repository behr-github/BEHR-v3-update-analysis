classdef misc_behr_update_plots
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        start_date = '2012-01-01';
        end_date = '2012-12-31';
        behr_final_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/6-Final';
        behr_nasa_brdf_vis_profs_wrftemp_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/5b-WRFTemp';
        behr_nasa_brdf_vis_profs_tempfix_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/5a-TempFix';
        behr_nasa_brdf_vis_profs_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/4-NewProfs';
        behr_nasa_brdf_vis_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/3-NewVis';
        behr_nasa_brdfD_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/2b-BRDF-D';
        behr_nasa_brdf_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/2-BRDF';
        behr_nasa_only_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/1-NASAv3';
        behr_v2_1C_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/0-v2.1C';
        
        bc_test_root = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/BC-Tests';
        
        wrf_path_v2 = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles';
        
        figs_root_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Figures';
        modis_truecol_dir = '/Volumes/share-sat/SAT/MODIS/MYD09/Winter2012/';
    end
    
    methods(Static = true)
        function make_behr_final(do_overwrite)
            % Produces the final version, with:
            %   * NASA v3 data
            %   * BRDF surf
            %   * New visible AMFs
            %   * New NO2 profiles
            %   * WRF temperature profiles
            %   * CVM gridding
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            % Verify that the most up-to-date PSM code is active.
            G = GitChecker;
            G.addReqCommits(behr_paths.psm_dir, '58ef9fb');
            G.addReqCommits(behr_paths.python_interface, 'a217fcd');
            G.Strict = true;
            G.checkState();
            
            save_dir = misc_behr_update_plots.behr_final_dir;
            misc_behr_update_plots.make_behr_with_parameters('6ddf77d', '298e5d9', save_dir, do_overwrite, true);
        end
        
        function make_behr_wrf_temp(do_overwrite)
            % Produces the final version, but with the longitude definition
            % in the temperature profiles fixed.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_brdf_vis_profs_wrftemp_dir;
            misc_behr_update_plots.make_behr_with_parameters('453c2d8', '88be28c', save_dir, do_overwrite, true);
        end
        
        function make_behr_temp_fix(do_overwrite)
            % Produces the final version, but with the longitude definition
            % in the temperature profiles fixed.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_brdf_vis_profs_tempfix_dir;
            misc_behr_update_plots.make_behr_with_parameters('1268384', '2ebbb06', save_dir, do_overwrite, true);
        end
        
        function make_behr_brdf_nasa_newvis_profs(do_overwrite)
            % Produces the version of BEHR that has:
            %   1) the NASA SP v3 base,
            %   2) the BRDF albedo,
            %   3) the COARTS sea albedo
            %   4) the new visible-only AMF
            %   5) the new profiles
            % but not the new CVM gridding.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir;
            misc_behr_update_plots.make_behr_with_parameters('e8c95e8', '9ec690e', save_dir, do_overwrite, true);
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
            misc_behr_update_plots.make_behr_with_parameters('e40a3fd', '2ac981a', save_dir, do_overwrite, false);
        end
        
        function make_behr_brdfD_and_nasa(do_overwrite)
            % Produces the version of BEHR that has:
            %   1) the NASA SP v3 base,
            %   2) the BRDF albedo using MCD43D{07,08,09,31} (instead of MCD43C1)
            %   3) the COARTS sea albedo
            % but not the new visible-only AMF, new profiles or the PSM
            % gridding.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_brdfD_dir;
            misc_behr_update_plots.make_behr_with_parameters('2298472', '2c692b5', save_dir, do_overwrite, false);
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
            misc_behr_update_plots.make_behr_with_parameters('4362069', 'c592651', save_dir, do_overwrite, false);
        end
        
        function make_behr_only_new_nasa(do_overwrite)
            % Produces the version of BEHR that has only the NASA SP v3
            % base. Everything else should be the same as in version 2,
            % except that we now use the NASA SP v3 as the base.
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            save_dir = misc_behr_update_plots.behr_nasa_only_dir;
            misc_behr_update_plots.make_behr_with_parameters('cfbad01', '8236170', save_dir, do_overwrite, false);
        end
        
        function make_bc_test_runs(do_overwrite)
            if ~exist('do_overwrite','var')
                do_overwrite = ask_yn('Overwrite existing files?');
            end
            
            moz_bc_dir = '/Volumes/share-wrf1/Tests/GeosBC-Test/MozBCs';
            geos_bc_dir = '/Volumes/share-wrf1/Tests/GeosBC-Test/GeosBCs';
            
            moz_save_dir = fullfile(misc_behr_update_plots.bc_test_root, 'MozBCs');
            if ~exist(moz_save_dir, 'dir')
                mkdir(moz_save_dir);
            end
            geos_save_dir = fullfile(misc_behr_update_plots.bc_test_root, 'GeosBCs');
            if ~exist(geos_save_dir, 'dir')
                mkdir(geos_save_dir);
            end
            
            bc_start_date = '2012-04-02';
            bc_end_date = '2012-04-08';
            
            % The SP files shouldn't change after the BRDF implementation
            sp_dir = fullfile(misc_behr_update_plots.behr_nasa_brdfD_dir, 'SP_Files');
            
            BEHR_main('start', bc_start_date,...
                'end', bc_end_date,...
                'sp_mat_dir', sp_dir,...
                'behr_mat_dir', moz_save_dir,...
                'no2_profile_path', moz_bc_dir,...
                'profile_mode', 'daily',...
                'overwrite', do_overwrite);
            
            BEHR_main('start', bc_start_date,...
                'end', bc_end_date,...
                'sp_mat_dir', sp_dir,...
                'behr_mat_dir', geos_save_dir,...
                'no2_profile_path', geos_bc_dir,...
                'profile_mode', 'daily',...
                'overwrite', do_overwrite);
        end
        
        function compare_bc_profiles
        end
        
        function compare_bc_amfs()
            % Assume that the files available are the same in both
            % directories
            files = dir(fullfile(misc_behr_update_plots.bc_test_root, 'MozBCs'));
            test_file = ask_multichoice('Choose file to compare', {files.name}, 'list', true);
            filter_clouds = ask_yn('Filter out cld. frac. > 0.2?');
            moz = load(fullfile(misc_behr_update_plots.bc_test_root, 'MozBCs', test_file), 'Data');
            geos = load(fullfile(misc_behr_update_plots.bc_test_root, 'GeosBCs', test_file), 'Data');
            
            file_date = datenum(regexp(test_file, '\d\d\d\d\d\d\d\d', 'match', 'once'), 'yyyymmdd');
            file_date = datestr(file_date);
            for a=1:numel(moz.Data)
                lon = moz.Data(a).Longitude;
                lat = moz.Data(a).Latitude;
                rdel = reldiff(moz.Data(a).BEHRAMFTrop, geos.Data(a).BEHRAMFTrop)*100;
                % Clouds should be the same in both, they are not affected
                % by the NO2 profile
                if filter_clouds
                    high_clds = moz.Data(a).CloudFraction > 0.2;
                    rdel(high_clds) = nan;
                end
                figure; 
                pcolor(lon, lat, rdel);
                shading flat;
                colorbar;
                state_outlines('k','not','ak','hi');
                title(sprintf('%s - Orbit %d', file_date, a));
            end
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
                struct('dir', misc_behr_update_plots.behr_nasa_brdfD_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'MODISAlbedo'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','BEHRAMFTrop','BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_tempfix_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','BEHRAMFTrop','BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_wrftemp_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','BEHRAMFTrop','BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_final_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'}})...
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
            
            E = JLLErrors;
            
            % These must be in the order they are to be compared
            incr_dirs = {misc_behr_update_plots.behr_v2_1C_dir,...
                         misc_behr_update_plots.behr_nasa_only_dir,...
                         misc_behr_update_plots.behr_nasa_brdf_dir,...
                         misc_behr_update_plots.behr_nasa_brdf_vis_dir,...
                         misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir,...
                         misc_behr_update_plots.behr_final_dir,...
                         misc_behr_update_plots.behr_final_cvm_dir,...
                         misc_behr_update_plots.behr_modis_quality_dir,...
                         misc_behr_update_plots.behr_modis_quality_best_dir};
            
            % Use the "outside.mat" file to clean up the plots. Must be the
            % same length as incr_dirs
            use_outside = [true, true, true, true, false, false, false, false, false];
            
            if numel(use_outside) ~= numel(incr_dirs)
                E.sizeMismatch('incr_dirs', 'use_outside');
            end
            
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
                        
                        [title_str, quantity_name, unit_name] = misc_behr_update_plots.get_title_from_file_name(base_names{b}, incr_dirs{base_ind}, incr_dirs{new_ind});
                        save_dir = fullfile(misc_behr_update_plots.figs_root_dir, 'IncrementTests');
                        
                        misc_behr_update_plots.plot_change_avg(Dbase.lon_grid, Dbase.lat_grid, Dbase.no2_vcds, Dnew.lon_grid, Dnew.lat_grid, Dnew.no2_vcds, false, title_str, quantity_name, unit_name, save_dir);
                        misc_behr_update_plots.plot_change_avg(Dbase.lon_grid, Dbase.lat_grid, Dbase.no2_vcds, Dnew.lon_grid, Dnew.lat_grid, Dnew.no2_vcds, true, title_str, quantity_name, unit_name, save_dir);
                    end
                end
            end
        end
        
        
        function plot_single_incr_diff()
            % Get all of the directories that would contain average files
            fns = fieldnames(misc_behr_update_plots);
            data_dirs = {};
            for a=1:numel(fns)
                the_dir = misc_behr_update_plots.(fns{a});
                if ~isempty(strfind(the_dir, '/IncrementTests/')) && exist(the_dir, 'dir')
                    data_dirs{end+1} = the_dir;
                end
            end
            
            base_dir = ask_multichoice('Select the base directory for the difference', data_dirs, 'list', true);
            new_opts = veccat(data_dirs, {'No difference'});
            new_dir = ask_multichoice('Select the new directory for the difference', new_opts, 'list', true);
            
            diff_bool = ~strcmpi(new_dir, 'No difference');
                
            
            % Now find all the data fields that the two directories have in
            % common
            Fbase = dir(fullfile(base_dir, '*.mat'));
            if diff_bool
                Fnew = dir(fullfile(new_dir, '*.mat'));
                xx_base = find_common_elements({Fbase.name}, {Fnew.name}, 'nodup');
            else
                xx_base = true(size(Fbase));
            end
            avg_files = {Fbase(xx_base).name};
            
            comparison_file = ask_multichoice('Select the file to compare:', avg_files, 'list', true);
            if diff_bool
                if ask_yn('Do a relative difference? (absolute if no)')
                    diff_type = 'rel';
                else
                    diff_type = 'abs';
                end
            else
                diff_type = 'none';
            end
            
            Dbase = load(fullfile(base_dir, comparison_file));
            if diff_bool
                Dnew = load(fullfile(new_dir, comparison_file));
                [title_str, quantity_name, unit_name] = misc_behr_update_plots.get_title_from_file_name(comparison_file, base_dir, new_dir);
            else
                Dnew = Dbase;
                Dbase.no2_vcds = zeros(size(Dbase.no2_vcds));
                [title_str, quantity_name, unit_name] = misc_behr_update_plots.get_title_from_file_name(comparison_file, base_dir);
            end
            
            
            misc_behr_update_plots.plot_change_avg(Dbase.lon_grid, Dbase.lat_grid, Dbase.no2_vcds, Dnew.lon_grid, Dnew.lat_grid, Dnew.no2_vcds, diff_type, title_str, quantity_name, unit_name, '');
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
        
        
        function plot_profile_differences(date_in, daily_or_monthly, utc_hour)
            % Plot the difference between the version 2 and version 3
            % profiles.
            
            if ~exist('date_in', 'var')
                date_in = ask_date('What date should we plot?');
            end
            if ~exist('daily_or_monthly', 'var')
                daily_or_monthly = ask_multichoice('Use daily or monthly new profiles?', {'daily', 'monthly'}, 'list', true);
            end
            if ~exist('utc_hour', 'var')
                if strcmpi(daily_or_monthly, 'daily')
                    utc_hour = ask_number('Enter the UTC hour to plot (0-23)', 'testfxn', @(x) isscalar(x) && x >= 0 && x <= 23, 'testmsg', 'Value must be between 0 and 23');
                else
                    utc_hour = 0;
                end
            end
            
            new_wrf_prof = find_wrf_path(daily_or_monthly, date_in);
            old_wrf_prof = misc_behr_update_plots.wrf_path_v2;
            compare_old_new_profiles(date_in, new_wrf_prof, old_wrf_prof, daily_or_monthly, utc_hour);
        end
        
        
        function modis_avg_quality(start_date, end_date)
            % Averages the quality flags in the MODIS MCD43C1 data between
            % the dates specified and plots it. If no inputs given, asks
            % for the dates interactively
            E = JLLErrors;
            
            if ~exist('start_date','var')
                start_date = ask_date('Enter the first date to consider for MODIS data');
            end
            if ~exist('end_date', 'var')
                end_date = ask_date('Enter the last date to consider for MODIS data');
            end
            
            lonlim = [-125 -65];
            latlim = [25 50];
            
            dvec = datenum(start_date):datenum(end_date);
            [band3_lons, band3_lats] = modis_cmg_latlon(0.05, [-180, 180], [-90, 90], 'grid');
            yy = band3_lats(:,1) >= min(latlim) & band3_lats(:,1) <= max(latlim);
            xx = band3_lons(1,:) >= min(lonlim) & band3_lons(1,:) <= max(lonlim);
            
            band3_lons = band3_lons(yy,xx);
            band3_lats = band3_lats(yy,xx);
            Mzero = zeros(size(band3_lons));
            modis_quantities = struct('BRDF_Quality', Mzero, 'Percent_Inputs', Mzero, 'Percent_Snow', Mzero, 'BRDF_Albedo_Uncertainty', Mzero);
            modis_counts = make_empty_struct_from_cell(fieldnames(modis_quantities), Mzero);
            fns = fieldnames(modis_quantities);
            for a=1:numel(dvec)
                fprintf('Loading file from %s\n', datestr(dvec(a)));
                % Figure out the MODIS file
                year_str = datestr(dvec(a), 'yyyy');
                modis_file_pat = sprintf('MCD43C1.A%s%03d.006*.hdf', year_str, modis_date_to_day(dvec(a)));
                F = dir(fullfile(behr_paths.mcd43c1_dir, year_str, modis_file_pat));
                if numel(F) > 1
                    E.toomanyfiles(modis_file_pat);
                elseif numel(F) < 1
                    E.filenotfound(modis_file_pat);
                else
                    mcd43_info = hdfinfo(fullfile(behr_paths.mcd43c1_dir, year_str, F(1).name));
                end
                for b=1:numel(fns)
                    wstate = warning('off','all'); % suppress warnings about missing scale factors, it's ok
                    val = hdfreadmodis(mcd43_info.Filename, hdfdsetname(mcd43_info,1,1,fns{b}));
                    warning(wstate);
                    val = val(yy,xx);
                    nans = isnan(val);
                    val(nans) = 0;
                    modis_quantities.(fns{b}) = modis_quantities.(fns{b}) + val;
                    modis_counts.(fns{b}) = modis_counts.(fns{b}) + double(~nans);
                end
            end
            
            for b=1:numel(fns)
                modis_quantities.(fns{b}) = modis_quantities.(fns{b}) ./ modis_counts.(fns{b});
                
                figure;
                pcolor(band3_lons, band3_lats, modis_quantities.(fns{b}));
                shading flat;
                colorbar;
                colormap(jet);
                state_outlines('k','not','ak','hi');
                title(sprintf('Avg. %s from %s to %s', fns{b}, start_date, end_date));
            end
        end
        
        function modis_picture(date_in)
            if ~exist('date_in', 'var')
                date_in = ask_date('Enter the date to plot the truecolor image for');
            end
            
            modis_pattern = sprintf('MYD09CMG.A%04d%03d*.hdf', year(date_in), modis_date_to_day(date_in));
            F = dirff(fullfile(misc_behr_update_plots.modis_truecol_dir, modis_pattern));
            if numel(F) ~= 1
                fprintf('Cannot find file for %s\n', datestr(date_in));
                return
            end
            
            modis_cmg_truecolor_image(F(1).name, [-125 -65], [25 50]);
            title(datestr(date_in));
        end
        
        function check_wrf_temp_profs(date_in, daily_or_monthly)
            % Plots temperature matched by rProfile_WRF to native pixels
            % with state outline so that I can inspect the temperature
            % profiles and make sure they match the geography correctly.
            E = JLLErrors;
            
            % Must have the commit that adds temperature profiles matching
            % to rProfile_WRF
            G = GitChecker;
            G.addReqCommits(behr_paths.behr_utils, '2a79b98');
            G.checkState();
            
            if ~exist('date_in', 'var')
                date_in = ask_date('Which date to use for the temperature profiles?');
            else
                date_in = validate_date(date_in);
            end
            
            dm_list = {'daily', 'monthly'};
            if ~exist('daily_or_monthly', 'var')
                daily_or_monthly = ask_multichoice('Use daily or monthly temperature profiles?', dm_list, 'list', true);
            elseif ~ismember(daily_or_monthly, dm_list)
                E.badinput('DAILY_OR_MONTHLY must be one of %s', strjoin(dm_list, ', '));
            end
            
            % We'll use the native data from the "final" directory
            behr_pattern = sprintf('OMI_SP_*_%s.mat', datestr(date_in, 'yyyymmdd'));
            F = dirff(fullfile(misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir, 'SP_Files', behr_pattern));
            if numel(F) ~= 1
                E.filenotfound(behr_pattern);
            end
            
            [slon, slat] = state_outlines('not', 'ak', 'hi');
            
            D = load(F(1).name);
            Data = D.Data;
            for a=1:numel(D.Data)
                lon = Data(a).Longitude;
                lat = Data(a).Latitude;
                loncorns = Data(a).FoV75CornerLongitude;
                latcorns = Data(a).FoV75CornerLatitude;
                time = Data(a).Time;
                surfPres = Data(a).GLOBETerpres;
                surfPres(surfPres>=1013)=1013;
                pressure = behr_pres_levels();
                [~, temperature] = rProfile_WRF(date_in, daily_or_monthly, loncorns, latcorns, time, surfPres, pressure);
                
                tstr = sprintf('%s - Swath %d', datestr(date_in), a);
                
                plot_slice_gui(permute(temperature, [2 3 1]), lon, lat, slon, slat, tstr);
                
                % Also plot surface pressure to help check where things are
                % NaNs
                figure; pcolor(lon, lat, surfPres);
                shading flat;
                colorbar;
                state_outlines('k','not','ak','hi');
                title(tstr);
            end
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
                    [no2_vcds, lon_grid, lat_grid] = behr_time_average(djf_start, djf_end, lon_lim, lat_lim,...
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
                    [no2_vcds, lon_grid, lat_grid] = behr_time_average(jja_start, jja_end, lon_lim, lat_lim,...
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
                        [no2_vcds, lon_grid, lat_grid] = behr_time_average(djf_start, djf_end, lon_lim, lat_lim,...
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
                        [no2_vcds, lon_grid, lat_grid] = behr_time_average(jja_start, jja_end, lon_lim, lat_lim,...
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
            
            is_diff = exist('new_dir', 'var');
            
            [~, filename] = fileparts(filename); % remove any path or extension
            [~, base_dir] = fileparts(base_dir);
            if is_diff
                [~, new_dir] = fileparts(new_dir);
            end
            
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
                case 'BEHRAMFTrop'
                    quantity_str = 'AMF';
                    unit_str = 'unitless';
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
            
            if is_diff
                title_str = sprintf('%s: %s vs %s (%s, %s)', quantity_str, new_dir, base_dir, time_period, profs_used);
            else
                title_str = sprintf('%s: %s (%s, %s)', quantity_str, base_dir, time_period, profs_used);
            end
        end
        
        function plot_change_avg(old_lon_grid, old_lat_grid, old_val, new_lon_grid, new_lat_grid, new_val, diff_type, title_str, quantity_name, unit_name, save_dir)
            % pass an empty string as save_dir to just plot the figure
            % without saving or closing it
            E = JLLErrors;
            
            if ~isequal(size(old_val), size(new_val))
                % Very simplistic check to see which is smaller, assumes
                % that the smaller one's coordinates are all inside the
                % larger one.
                if numel(old_val) < numel(new_val)
                    new_val = interp2(new_lon_grid, new_lat_grid, new_val, old_lon_grid, old_lat_grid);
                    lon_grid = old_lon_grid;
                    lat_grid = old_lat_grid;
                else
                    old_val = interp2(old_lon_grid, old_lat_grid, old_val, new_lon_grid, new_lat_grid);
                    lon_grid = new_lon_grid;
                    lat_grid = new_lat_grid;
                end
            else
                lon_grid = new_lon_grid;
                lat_grid = new_lat_grid;
            end
            
            if strcmpi(diff_type, 'rel')
                val_diff = reldiff(new_val, old_val)*100;
                save_name = sprintf('Relative difference %s', title_str);
                cb_label = sprintf(tex_in_printf('%%\Delta %s (%s)'), quantity_name, unit_name);
                clims = calc_plot_limits(val_diff, 5, 'diff', 'pow10', 'max', [-100 100]);
                cmap = blue_red_cmap;
            elseif strcmpi(diff_type, 'abs')
                val_diff = new_val - old_val;
                save_name = sprintf('Absolute difference %s', title_str);
                cb_label = sprintf(tex_in_printf('\Delta %s'), quantity_name);
                clims = calc_plot_limits(val_diff, 5, 'diff', 'pow10', 'max', [-Inf Inf]);
                cmap = blue_red_cmap;
            elseif strcmpi(diff_type, 'none')
                val_diff = new_val - old_val;
                save_name = title_str;
                cb_label = quantity_name;
                clims = calc_plot_limits(val_diff, 5, 'zero', 'pow10', 'max', [-Inf Inf]);
                cmap = parula;
            end
            
            fig=figure; pcolor(lon_grid, lat_grid, val_diff);
            shading flat; 
            cb=colorbar; 
            caxis(clims); 
            colormap(cmap);
            state_outlines('k', 'not', 'ak', 'hi');
            title(title_str);
            
            
            cb.Label.String = cb_label;
            cb.FontSize = 14;
            
            if ~isempty(save_dir)
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
    
end

