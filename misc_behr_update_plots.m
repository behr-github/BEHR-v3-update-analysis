classdef misc_behr_update_plots
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        start_date = '2012-01-01';
        end_date = '2012-12-31';
        behr_v3B_var_trop_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/7b-VarTrop';
        behr_v3B_daily_fix_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/7a-DailyFix';
        behr_final_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/6-Final';
        behr_nasa_brdf_vis_profs_wrftemp_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/5b-WRFTemp';
        behr_nasa_brdf_vis_profs_tempfix_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/5a-TempFix';
        behr_nasa_brdf_vis_profs_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/4-NewProfs';
        behr_nasa_brdf_vis_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/3-NewVis';
        behr_nasa_brdfD_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/2b-BRDF-D';
        behr_nasa_brdf_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/2-BRDF';
        behr_modisv6_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/1b-MODISv6';
        behr_nasa_only_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/1-NASAv3';
        behr_v2_1C_dir = '/Volumes/share-sat/SAT/BEHR/IncrementTests/0-v2.1C';
        
        bc_test_root = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/BC-Tests';
        surf_refl_root = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/Surf-Refl';
        wrf_root = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/WRF';
        
        wrf_path_v2 = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles';
        
        figs_root_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Figures';
        modis_truecol_dir = '/Volumes/share-sat/SAT/MODIS/MYD09/Winter2012/';
        modis_landcover_file = '/Volumes/share-sat/SAT/MODIS/MCD12C1/MCD12C1.A2012001.051.2013178154403.hdf';
    end
    
    methods(Static = true)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property like methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function d = my_dir()
            d = fileparts(mfilename('fullpath'));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Generation methods %
        %%%%%%%%%%%%%%%%%%%%%%
        
        function make_behr_v3B_var_tropopause(do_overwrite)
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            G = GitChecker;
            G.addReqCommits(behr_paths.psm_dir, '223e1e1');
            G.addReqCommits(behr_paths.python_interface, 'dcd79fc');
            G.addReqCommits(behr_paths.utils, '51a0869'); % new floodfill()
            G.addReqCommits(behr_paths.wrf_utils, '9272d08'); % need upgraded find_wrf_tropopause
            G.Strict = true;
            G.checkState();
            
            save_dir = misc_behr_update_plots.behr_v3B_var_trop_dir;
            misc_behr_update_plots.make_behr_with_parameters('8ed138f', 'eb1f53e', save_dir, do_overwrite, true);
        end
        

        function make_behr_v3B_daily_fix(do_overwrite)
            if ~exist('do_overwrite', 'var')
                do_overwrite = false;
            end
            
            G = GitChecker;
            G.addReqCommits(behr_paths.psm_dir, '223e1e1');
            G.addReqCommits(behr_paths.python_interface, 'dcd79fc');
            G.addReqCommits(behr_paths.utils, '6c401a8');
            G.addReqCommits(behr_paths.wrf_utils, 'd97776f');
            G.Strict = true;
            G.checkState();
            
            save_dir = misc_behr_update_plots.behr_v3B_daily_fix_dir;
            misc_behr_update_plots.make_behr_with_parameters('ea54fe3', 'a426124', save_dir, do_overwrite, true);
        end
        
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
        
        function average_geometry()
            dvec = datenum('2012-06-01'):datenum('2012-08-31');
            USGrid = GlobeGrid(0.05, 'domain', 'us');
            blank_mat = nan(size(USGrid.GridLon));
            Geometry = struct('SolarZenithAngle', blank_mat,...
                'ViewingZenithAngle', blank_mat,...
                'RelativeAzimuthAngle', blank_mat);
            areaweight = blank_mat;
            
            fns = fieldnames(Geometry);
            
            for d=1:numel(dvec)
                fprintf('Loading %s\n', datestr(dvec(d)));
                Data = load_behr_file(dvec(d), 'monthly', 'us');
                for a=1:numel(Data)
                    for f=1:numel(fns)
                        [this_val, this_weight] = cvm_generic_wrapper(Data(a).FoV75CornerLongitude, Data(a).FoV75CornerLatitude, Data(a).(fns{f}), USGrid);
                        if f==1
                            areaweight = nansum2(cat(3, areaweight, this_weight),3);
                        end
                        Geometry.(fns{f}) = nansum2(cat(3, Geometry.(fns{f}), this_val .* this_weight),3);
                    end
                end
            end
            
            for f=1:numel(fns)
                Geometry.(fns{f}) = Geometry.(fns{f}) ./ areaweight;
            end
            
            save(fullfile(misc_behr_update_plots.surf_refl_root, 'SummerGeometry.mat'), 'Geometry')
        end
        
        function average_brdf_coeffs()
            dvec = datenum('2012-06-01'):datenum('2012-08-31');
            for d=1:numel(dvec)
                fprintf('Reading %s\n', datestr(dvec(d)));
                b3data = read_modis_albedo(behr_paths.mcd43d_dir, dvec(d), [-125 -65], [25 50], 'DEBUG_LEVEL', 0);
                if d==1
                    f_iso = nan(size(b3data.iso));
                    f_iso_count = zeros(size(f_iso));
                    f_vol = nan(size(b3data.vol));
                    f_vol_count = zeros(size(f_vol));
                    f_geo = nan(size(b3data.geo));
                    f_geo_count = zeros(size(f_geo));
                    lon = b3data.lons;
                    lat = b3data.lats;
                end
                
                % Typically in BEHR I use albedo with quality better than
                % (<=) 3. So here we do not average coefficients with lower
                % (>) quality than 3. The lower the quality flag, the
                % better the quality.
                b3data.iso(b3data.quality > 3) = nan;
                b3data.vol(b3data.quality > 3) = nan;
                b3data.geo(b3data.quality > 3) = nan;
                
                f_iso = nansum2(cat(3, f_iso, b3data.iso),3);
                f_iso_count = f_iso_count + ~isnan(b3data.iso);
                f_vol = nansum2(cat(3, f_vol, b3data.vol),3);
                f_vol_count = f_vol_count + ~isnan(b3data.vol);
                f_geo = nansum2(cat(3, f_geo, b3data.geo),3);
                f_geo_count = f_geo_count + ~isnan(b3data.geo);
            end
            
            f_iso = f_iso ./ f_iso_count;
            f_vol = f_vol ./ f_vol_count;
            f_geo = f_geo ./ f_geo_count;
            
            fprintf('Saving...\n')
            save(fullfile(misc_behr_update_plots.surf_refl_root, 'SummerMODISCoefficients.mat'), '-v7.3', 'f_iso', 'f_vol', 'f_geo', 'lon', 'lat');
        end
        
        %%%%%%%%%%%%%%%%%%%%
        % Analysis methods %
        %%%%%%%%%%%%%%%%%%%%
        
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
                struct('dir', misc_behr_update_plots.behr_modisv6_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','MODISAlbedo'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdfD_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'MODISAlbedo'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','BEHRAMFTrop','BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_tempfix_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','BEHRAMFTrop','BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_wrftemp_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly','BEHRAMFTrop','BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_final_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'BEHRAMFTrop', 'BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_v3B_daily_fix_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'BEHRAMFTrop', 'BEHRAMFTropVisOnly'}}),...
                struct('dir', misc_behr_update_plots.behr_v3B_var_trop_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'BEHRAMFTrop', 'BEHRAMFTropVisOnly'}})...
                );
            
            % Removed: struct('dir', misc_behr_update_plots.behr_nasa_brdfD_dir, 'use_new_avg', false, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'MODISAlbedo'}}),...
            % since I skipped the MCD43C1 run this time.
            for a=1:numel(average_info)
                for b=1:numel(average_info(a).data_fields)
                    if average_info(a).use_new_avg
                        misc_behr_update_plots.average_behr(average_info(a).dir, average_info(a).data_fields{b}, do_overwrite, true);
                    else
                        misc_behr_update_plots.average_behr(average_info(a).dir, average_info(a).data_fields{b}, do_overwrite, false);
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
        
        
        function fig = plot_single_incr_diff(base_dir, new_dir, base_comparison_file, new_comparison_file, diff_type)
            E = JLLErrors;
            
            % Get all of the directories that would contain average files
            fns = fieldnames(misc_behr_update_plots);
            data_dirs = {};
            for a=1:numel(fns)
                the_dir = misc_behr_update_plots.(fns{a});
                if ~isempty(strfind(the_dir, '/IncrementTests/')) && exist(the_dir, 'dir')
                    data_dirs{end+1} = the_dir;
                end
            end
            
            % Check if base and new directories were given. If not, ask the
            % user which ones to use
            if ~exist('base_dir', 'var')
                base_dir = ask_multichoice('Select the base directory for the difference', data_dirs, 'list', true);
            elseif ~ismember(base_dir, data_dirs)
                E.badinput('BASE_DIR is not one of the recognized data directories');
            end
            
            new_opts = veccat(data_dirs, {'No difference'});
            if ~exist('new_dir', 'var')
                new_dir = ask_multichoice('Select the new directory for the difference', new_opts, 'list', true);
            elseif ~ismember(new_dir, new_opts)
                E.badinput('BASE_DIR is not one of the recognized data directories');
            end
            
            diff_bool = ~strcmpi(new_dir, 'No difference');
            
            
            % Now find all the data fields that the two directories have in
            % common
            Fbase = dir(fullfile(base_dir, '*.mat'));
            if diff_bool
                Fnew = dir(fullfile(new_dir, '*.mat'));
                base_avg_files = misc_behr_update_plots.make_common_files_list({Fbase.name}, {Fnew.name});
            end
            
            % Get user input of the file to compare (data field and season)
            % as well as the difference type. If doing a difference, allow
            % the user to compare daily and monthly files.
            if ~exist('base_comparison_file', 'var')
                base_comparison_file = ask_multichoice('Select the base file to compare against:', base_avg_files, 'list', true);
            elseif ~ismember(base_comparison_file, base_avg_files)
                E.badinput('BASE_COMPARISON_FILE is not one of the allowed files: %s', strjoin(base_avg_files, ', '));
            end
            
            if diff_bool
                % Only allow the user to compare files for the same
                % quantity and time period, but allow monthly/daily
                % difference
                new_avg_files = misc_behr_update_plots.make_common_files_list({Fnew.name}, {base_comparison_file});
                
                if ~exist('new_comparison_file', 'var')
                    new_comparison_file = ask_multichoice('Select the new file to compare:', new_avg_files, 'list', true);
                elseif ~ismember(new_comparison_file, new_avg_files)
                    E.badinput('NEW_COMPARISON_FILE is not one of the allowed files: %s', strjoin(new_avg_files, ', '));
                end
                
                allowed_diff_types = {'rel', 'abs'};
                if ~exist('diff_type', 'var')
                    if ask_yn('Do a relative difference? (absolute if no)')
                        diff_type = 'rel';
                    else
                        diff_type = 'abs';
                    end
                elseif ~ismember(diff_type, allowed_diff_types)
                    E.badinput('DIFF_TYPE must be one of: %s', strjoin(allowed_diff_types, ', '));
                end
            else
                diff_type = 'none';
            end
            
            Dbase = load(fullfile(base_dir, base_comparison_file));
            Dbase = misc_behr_update_plots.apply_outside(Dbase, base_dir);
            if diff_bool
                Dnew = load(fullfile(new_dir, new_comparison_file));
                Dnew = misc_behr_update_plots.apply_outside(Dnew, new_dir);
                [title_str, quantity_name, unit_name] = misc_behr_update_plots.get_title_from_file_name(base_comparison_file, base_dir, new_comparison_file, new_dir);
            else
                Dnew = Dbase;
                Dbase.no2_vcds = zeros(size(Dbase.no2_vcds));
                [title_str, quantity_name, unit_name] = misc_behr_update_plots.get_title_from_file_name(base_comparison_file, base_dir);
            end
            
            
            fig = misc_behr_update_plots.plot_change_avg(Dbase.lon_grid, Dbase.lat_grid, Dbase.no2_vcds, Dnew.lon_grid, Dnew.lat_grid, Dnew.no2_vcds, diff_type, title_str, quantity_name, unit_name, '');
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
        
        function plot_sza_range_over_year()
            dvec = datenum('2012-01-01'):datenum('2012-12-31');
            min_szas = nan(size(dvec));
            max_szas = nan(size(dvec));
            for a=1:numel(dvec)
                fprintf('%s\n', datestr(dvec(a)));
                Data = load_behr_file(dvec(a), 'monthly', 'us');
                szas = cat_sat_data(Data, 'SolarZenithAngle');
                min_szas(a) = min(szas(:));
                max_szas(a) = max(szas(:));
            end
            
            l = gobjects(2,1);
            figure;
            l(1) = line(dvec, min_szas, 'color', 'b', 'linewidth', 2);
            l(2) = line(dvec, max_szas, 'color', 'r', 'linewidth', 2);
            
            legend(l, {'Minimum SZA','Maximum SZA'});
        end
        
        function [prof_fig, prof_med_fig, freq_fig, fig_mean_bin, fig_med_bin, fig_map] = plot_no2_vs_cloudfrac(fxn_start_date, fxn_end_date, as_shape_factor, region, DEBUG_LEVEL)
            if ~exist('fxn_start_date', 'var')
                fxn_start_date = ask_date('Give the starting date');
            else
                validate_date(fxn_start_date);
            end
            if ~exist('fxn_end_date', 'var')
                fxn_end_date = ask_date('Give the ending date');
            else
                validate_date(fxn_end_date);
            end
            if ~exist('as_shape_factor', 'var')
                as_shape_factor = ask_yn('Calculate shape factor?');
            elseif ~isscalar(as_shape_factor) || ~islogical(as_shape_factor)
                E.badinput('AS_SHAPE_FACTOR must be a scalar logical');
            end
            
            allowed_regions = {'SE','NW'};
            if ~exist('region', 'var')
                region = ask_multichoice('Which region to restrict to?', allowed_regions, 'list', true);
            elseif ~ischar(region) || ~ismember(region, allowed_regions)
                E.badinput('REGION must be one of the strings: %s', strjoin(allowed_regions, ', '));
            end
            
            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 2;
            end
            
            lims.SE.lon = [-100 -75];
            lims.SE.lat = [25 37.5];
            lims.SE.color = 'b'; % color just used when plotting the regions
            lims.NW.lon = [-125 -100];
            lims.NW.lat = [37.5 50];
            lims.NW.color = 'r';
            
            try
                lonlim = lims.(region).lon;
                latlim = lims.(region).lat;
            catch err
                if strcmpi(err.identifier, 'MATLAB:nonExistentField')
                    E.notimplemented('Region %s has no lat/lon limits defined', region);
                else
                    rethrow(err)
                end
            end
            
            % Make a quick plot of the regions
            fig_map = figure;
            coasts = load('coast');
            state_outlines('k','not','ak','hi');
            line(coasts.long, coasts.lat, 'color', 'k');
            def_regions = fieldnames(lims);
            for a=1:numel(def_regions)
                x_all = lims.(def_regions{a}).lon([1 1 2 2 1]);
                y_all = lims.(def_regions{a}).lat([1 2 2 1 1]);
                line(x_all, y_all, 'color', lims.(def_regions{a}).color, 'linewidth', 2);
                text(mean(x_all), mean(y_all), def_regions{a}, 'FontSize', 18, 'BackgroundColor', [0.8 0.8 0.8]);
            end
            
            xlim([-130 -65]); ylim([22.5 52.5]);
            set(gca,'fontsize',14);
            
            omi_cloudfrac = [];
            all_pres_levs = [];
            behr_no2_daily = [];
            behr_no2_monthly = [];
            
            dvec = datenum(fxn_start_date):datenum(fxn_end_date);
            for d=1:numel(dvec)
                if DEBUG_LEVEL > 1
                    fprintf('Loading %s\n', datestr(dvec(d)));
                elseif DEBUG_LEVEL > 0
                    % Using this means we cannot have any other print
                    % statements in this loop
                    text_progress_bar(d, numel(dvec));
                end
                
                D = load(fullfile(misc_behr_update_plots.behr_final_dir, 'DailyProfs', behr_filename(dvec(d), 'daily', 'us')),'Data');
                DataDaily = D.Data;
                
                M = load(fullfile(misc_behr_update_plots.behr_final_dir, 'MonthlyProfs', behr_filename(dvec(d), 'monthly', 'us')), 'Data');
                DataMonthly = M.Data;
                
                for a=1:numel(DataDaily)
                    % We still need to reject row anomaly pixels because they
                    % affect the OMI cloud fraction. Cloud fraction and row
                    % anomaly should be the same between daily and monthly
                    % files.
                    xx = DataDaily(a).XTrackQualityFlags == 0 & DataDaily(a).Longitude > lonlim(1) & DataDaily(a).Longitude < lonlim(2) & DataDaily(a).Latitude > latlim(1) & DataDaily(a).Latitude < latlim(2);
                    if sum(xx(:)) == 0
                        continue
                    end
                    omi_cloudfrac = cat(1, omi_cloudfrac, DataDaily(a).CloudFraction(xx));
                    
                    % Get daily NO2
                    no2_profs = DataDaily(a).BEHRNO2apriori(:,xx);
                    pres_levs = DataDaily(a).BEHRPressureLevels(:,xx);
                    if as_shape_factor
                        vcds = integrate_published_behr_apriori(DataDaily(a));
                        vcds = vcds(xx);
                        for b=1:size(no2_profs,2)
                            no2_profs(:,b) = no2_profs(:,b) ./ vcds(b) * 1e14; %scale by 1e14 just to make the numbers not so crazy small
                        end
                    end
                    behr_no2_daily = cat(2, behr_no2_daily, no2_profs);
                    all_pres_levs = cat(2, all_pres_levs, pres_levs);
                    
                    % Now monthly
                    no2_profs_m = DataMonthly(a).BEHRNO2apriori(:,xx);
                    if as_shape_factor
                        vcds_m = integrate_published_behr_apriori(DataMonthly(a));
                        vcds_m = vcds_m(xx);
                        for b=1:size(no2_profs_m,2)
                            no2_profs_m(:,b) = no2_profs_m(:,b) ./ vcds_m(b) * 1e14; %scale by 1e14 just to make the numbers not so crazy small
                        end
                    end
                    behr_no2_monthly = cat(2, behr_no2_monthly, no2_profs_m);
                end
            end
            
            cc20 = omi_cloudfrac <= 0.2;
            cc_all = true(size(omi_cloudfrac));
            
            % Average profiles with standard deviations
            behr_no2_monthly_interp = nan(numel(behr_pres_levels()), size(behr_no2_monthly,2));
            behr_no2_daily_interp = nan(numel(behr_pres_levels()), size(behr_no2_daily,2));
            for a=1:size(behr_no2_monthly,2)
                % The whole log-log space thing shouldn't be necessary
                % because we're interpolating to pressure values that
                % already exist, this just removes the extra pressure
                % levels
                behr_no2_monthly_interp(:,a) = naninterp1(all_pres_levs(:,a), behr_no2_monthly(:,a), behr_pres_levels());
                behr_no2_daily_interp(:,a) = naninterp1(all_pres_levs(:,a), behr_no2_daily(:,a), behr_pres_levels());
            end
            
            if ~as_shape_factor
                x_str = '%s [NO_2] (mixing ratio)';
            else
                x_str = '%s NO_2 shape factor';
            end
            
            avg_no2_monthly = nanmean(behr_no2_monthly_interp(:,cc20),2);
            avg_no2_monthly_all = nanmean(behr_no2_monthly_interp(:,cc_all),2);
            avg_no2_daily = nanmean(behr_no2_daily_interp(:,cc20),2);
            avg_no2_daily_all = nanmean(behr_no2_daily_interp(:,cc_all),2);
            
            l=gobjects(4,1);
            prof_fig = figure;
            l(1)=line(avg_no2_monthly, behr_pres_levels(), 'color', 'b', 'linewidth', 2, 'marker', '^');
            scatter_errorbars(avg_no2_monthly, behr_pres_levels(), nanstd(behr_no2_monthly_interp(:,cc20),0,2), 'color', 'b','direction','x');
            l(2)=line(avg_no2_daily, behr_pres_levels(), 'color', 'r', 'linewidth', 2, 'marker', 'o');
            scatter_errorbars(avg_no2_daily, behr_pres_levels(), nanstd(behr_no2_daily_interp(:,cc20),0,2), 'color', 'r','direction','x');
            l(3)=line(avg_no2_monthly_all, behr_pres_levels(), 'color', [1 0 1], 'linewidth', 2, 'marker', 'v');
            scatter_errorbars(avg_no2_monthly_all, behr_pres_levels(), nanstd(behr_no2_monthly_interp(:,cc_all),0,2), 'color', [1 0 1], 'direction','x');
            l(4)=line(avg_no2_daily_all, behr_pres_levels(), 'color', [1 0.5 0], 'linewidth', 2, 'marker', 's');
            scatter_errorbars(avg_no2_daily_all, behr_pres_levels(), nanstd(behr_no2_daily_interp(:,cc_all),0,2), 'color', [1 0.5 0], 'direction','x');
            legend(l, {'Monthly CF < 20%', 'Daily CF < 20%', 'Monthly All CF', 'Daily All CF'});
            set(gca,'fontsize',14,'ydir','reverse');
            xlabel(sprintf(x_str,'Mean'));
            ylabel('Pressure (hPa)');
            title(region);
            
            med_no2_monthly = nanmedian(behr_no2_monthly_interp(:,cc20),2);
            quart_no2_monthly = quantile(behr_no2_monthly_interp(:,cc20), [0.25, 0.75],2);
            med_no2_daily = nanmedian(behr_no2_daily_interp(:,cc20),2);
            quart_no2_daily = quantile(behr_no2_daily_interp(:,cc20),[0.25,0.75],2);
            med_no2_monthly_all = nanmedian(behr_no2_monthly_interp(:,cc_all),2);
            quart_no2_monthly_all = quantile(behr_no2_monthly_interp(:,cc_all), [0.25, 0.75],2);
            med_no2_daily_all = nanmedian(behr_no2_daily_interp(:,cc_all),2);
            quart_no2_daily_all = quantile(behr_no2_daily_interp(:,cc_all),[0.25,0.75],2);
            
            l=gobjects(2,1);
            prof_med_fig = figure;
            l(1) = line(med_no2_monthly, behr_pres_levels(), 'color','b','linewidth',2,'marker','^');
            scatter_errorbars(med_no2_monthly, behr_pres_levels()-2, quart_no2_monthly(:,1), quart_no2_monthly(:,2), 'color','b','direction','x');
            l(2) = line(med_no2_daily, behr_pres_levels(), 'color', 'r', 'linewidth',2,'marker','o');
            scatter_errorbars(med_no2_daily, behr_pres_levels()+2, quart_no2_daily(:,1), quart_no2_daily(:,2), 'color', 'r', 'direction', 'x');
            l(3) = line(med_no2_monthly_all, behr_pres_levels(), 'color',[1 0 1],'linewidth',2,'marker','^');
            scatter_errorbars(med_no2_monthly_all, behr_pres_levels()-2, quart_no2_monthly_all(:,1), quart_no2_monthly_all(:,2), 'color',[1 0 1],'direction','x');
            l(4) = line(med_no2_daily_all, behr_pres_levels(), 'color', [1 0.5 0], 'linewidth',2,'marker','o');
            scatter_errorbars(med_no2_daily_all, behr_pres_levels()+2, quart_no2_daily_all(:,1), quart_no2_daily_all(:,2), 'color', [1 0.5 0], 'direction', 'x');
            legend(l, {'Monthly CF < 20%', 'Daily CF < 20%', 'Monthly All CF', 'Daily All CF'});
            set(gca,'fontsize',14,'ydir','reverse');
            xlabel(sprintf(x_str,'Median'));
            ylabel('Pressure (hPa)');
            title(region);
            
            % Frequency distribution plots of concentration or shape factor
            % for pixels with cloud fraction < 0.2
            
            pp = all_pres_levs > 400;
            behr_no2_monthly(pp) = nan;
            ut_no2_monthly = squeeze(nanmean(behr_no2_monthly,1));
            behr_no2_daily(pp) = nan;
            ut_no2_daily = squeeze(nanmean(behr_no2_daily,1));
            
            
            bin_edges = 0:1e-12:1e-9;
            if ~as_shape_factor
                x_str = '[NO_2] (mixing ratio, above 400 hPa)';
            else
                x_str = 'NO_2 shape factor (above 400 hPa)';
            end
            counts_monthly = histc(ut_no2_monthly(cc20), bin_edges);
            counts_monthly_all = histc(ut_no2_monthly(cc_all), bin_edges);
            counts_daily = histc(ut_no2_daily(cc20), bin_edges);
            counts_daily_all = histc(ut_no2_daily(cc_all), bin_edges);
            
            freq_fig = figure;
            line(bin_edges, counts_monthly, 'color', 'b', 'linewidth', 2);
            line(bin_edges, counts_daily, 'color', 'r', 'linewidth', 2);
            line(bin_edges, counts_monthly_all, 'color', [1 0 1], 'linewidth', 2);
            line(bin_edges, counts_daily_all, 'color', [1 0.5 0], 'linewidth', 2);
            legend('Monthly CF < 20%', 'Daily CF < 20%', 'Monthly All CF', 'Daily All CF');
            xlabel(x_str);
            ylabel('Counts');
            title(region)
            set(gca,'fontsize',14);
            
            % Mean and median UT NO2 binned by cloud fraction
            [binned_ut_no2_monthly, bin_centers] = bin_data(omi_cloudfrac, ut_no2_monthly, 0:0.05:1);
            binned_ut_no2_daily = bin_data(omi_cloudfrac, ut_no2_daily, 0:0.05:1);
            
            mean_ut_no2_monthly = cellfun(@nanmean,binned_ut_no2_monthly);
            std_ut_no2_monthly = cellfun(@nanstd, binned_ut_no2_monthly);
            mean_ut_no2_daily = cellfun(@nanmean, binned_ut_no2_daily);
            std_ut_no2_daily = cellfun(@nanstd, binned_ut_no2_daily);
            
            l = gobjects(2,1);
            fig_mean_bin = figure;
            offset = 0.01;
            l(1)=plot(bin_centers-offset, mean_ut_no2_monthly, 'bo');
            scatter_errorbars(bin_centers-offset, mean_ut_no2_monthly, std_ut_no2_monthly, 'color','b');
            hold on
            l(2)=plot(bin_centers+offset, mean_ut_no2_daily, 'r^');
            scatter_errorbars(bin_centers+offset, mean_ut_no2_daily, std_ut_no2_daily, 'color', 'r');
            xlabel('OMI Cloud Fration'); ylabel('Mean UT NO_2'); title(region)
            legend(l, {'Monthly','Daily'});
            
            med_ut_no2_monthly = cellfun(@nanmedian, binned_ut_no2_monthly);
            quart_ut_no2_monthly = cellfun(@(x) quantile(x, [0.25 0.75]), binned_ut_no2_monthly, 'uniformoutput', false);
            quart_ut_no2_monthly = cat(1, quart_ut_no2_monthly{:});
            med_ut_no2_daily = cellfun(@nanmedian, binned_ut_no2_daily);
            quart_ut_no2_daily = cellfun(@(x) quantile(x, [0.25 0.75]), binned_ut_no2_daily, 'uniformoutput', false);
            quart_ut_no2_daily = cat(1, quart_ut_no2_daily{:});
            
            l = gobjects(2,1);
            fig_med_bin = figure;
            l(1)=plot(bin_centers-offset, med_ut_no2_monthly, 'bo');
            scatter_errorbars(bin_centers-offset, med_ut_no2_monthly, quart_ut_no2_monthly(:,1), quart_ut_no2_monthly(:,2), 'color', 'b');
            hold on
            l(2)=plot(bin_centers+offset, med_ut_no2_daily, 'r^');
            scatter_errorbars(bin_centers+offset, med_ut_no2_daily, quart_ut_no2_daily(:,1), quart_ut_no2_daily(:,2), 'color', 'r');
            xlabel('OMI Cloud Fraction'); ylabel('Median UT NO_2'); title(region);
            legend(l, {'Monthly','Daily'});
        end
        
        function fig = plot_ground_cover()
            hdfi = hdfinfo(misc_behr_update_plots.modis_landcover_file);
            lc_type1 = double(hdfread(hdfi.Filename, hdfdsetname(hdfi,1,1,1)));
            [lon,lat,xx,yy] = modis_cmg_latlon(0.05,[-125 -65],[25 50],'grid');
            lc_type1 = lc_type1(yy,xx);
            
            % Get all but the "unclassified" land type
            lc_ticks = [hdfi.Vgroup.Vgroup(1).SDS(1).Attributes(5:end-1).Value];
            lc_labels = {hdfi.Vgroup.Vgroup(1).SDS(1).Attributes(5:end-1).Name};
            
            fig=figure;
            map_ax = pcolor(lon,lat,lc_type1);
            shading flat;
            colormap(jet);
            state_outlines('k');
            cb=colorbar;
            cb.Ticks = lc_ticks;
            cb.TickLabels = lc_labels;
            
            set(gca,'fontsize',16);
            cb.FontSize = 10;
            fig.Position(3) = fig.Position(3)*1.5;
            fig.Children(2).Position(3) = 0.63;
        end
        
        function [fig_s3, fig_rad_v_geo, fig_clear, fig_cloudy, fig_tot_v_vis] = plot_vis_scatter(plot_gui)
            if ~exist('plot_gui', 'var')
                plot_gui = true;
            end
            
            %vis_start_date = '2012-07-01';
            vis_start_date = '2012-06-01';
            %vis_end_date = '2012-07-14';
            vis_end_date = '2012-08-31';
            fields = {'BEHRColumnAmountNO2TropVisOnly', 'BEHRColumnAmountNO2Trop','CloudRadianceFraction','CloudFraction','CloudPressure','BEHRNO2apriori','XTrackQualityFlags'};
            [no2_vis_old, no2_tot_old, cld_rad_frac_old, cld_frac_old, cld_pres_old, apriori_old, xtrack_old, pixel_info_old] = cat_sat_data(fullfile(misc_behr_update_plots.behr_nasa_brdfD_dir, 'MonthlyProfs'), fields, 'startdate', vis_start_date, 'enddate', vis_end_date);
            [no2_vis_new, no2_tot_new, ~, ~, ~, ~, xtrack_new] = cat_sat_data(fullfile(misc_behr_update_plots.behr_nasa_brdf_vis_dir, 'MonthlyProfs'), fields, 'startdate', vis_start_date, 'enddate', vis_end_date);
            
            row_anom = xtrack_old > 0 | xtrack_new > 0;
            
            no2_vis_old(row_anom) = [];
            no2_tot_old(row_anom) = [];
            cld_rad_frac_old(row_anom) = [];
            cld_frac_old(row_anom) = [];
            cld_pres_old(row_anom) = [];
            apriori_old = behr_apriori_surface(apriori_old);
            apriori_old(row_anom) = [];
            no2_vis_new(row_anom) = [];
            no2_tot_new(row_anom) = [];
            pixel_info_old(row_anom) = [];
            
            rdel = reldiff(no2_vis_new, no2_vis_old)*100;
            
            %{
            figure;
            scatter(no2_vis_old, no2_vis_new, [], rdel*100);
            xlabel('Old visible VCD');
            ylabel('New visible VCD');
            cb=colorbar;
            cb.Label.String = 'Percent change';
            set(gca,'fontsize',16);
            xylims([0 3e16]);
            
            figure;
            line(no2_tot_old, no2_vis_old, 'marker', 'o', 'color', 'k', 'linestyle', 'none');
            xlabel('Old total VCD');
            ylabel('Old visible VCD');
            set(gca,'fontsize',16);
            
            figure;
            line(no2_tot_new, no2_vis_new, 'marker', 'o', 'color', 'k', 'linestyle', 'none');
            xlabel('New total VCD');
            ylabel('New visible VCD');
            set(gca,'fontsize',16);
            
            figure;
            scatter3(cld_rad_frac_old, cld_pres_old, rdel, [], apriori_old*1e9);
            xlabel('Cloud radiance fraction');
            ylabel('Cloud pressure');
            zlabel('Percent change in vis NO2');
            cb = colorbar;
            cb.Label.String = 'Surface apriori [NO_2] (ppb)';
            set(gca,'fontsize',16);
            
            figure;
            scatter(rdel, cld_rad_frac_old, [], cld_frac_old);
            xlabel('Percent change in vis NO2');
            ylabel('Cloud radiance fraction');
            cb=colorbar;
            cb.Label.String = 'Geo. cloud frac.';
            set(gca,'fontsize',16);
            
            figure;
            plot(rdel, cld_pres_old, 'ko');
            xlabel('Percent change in vis NO2');
            ylabel('Cloud pressure');
            set(gca,'fontsize',16);
            
            figure;
            plot(rdel, apriori_old, 'ko');
            xlabel('Percent change in vis NO2');
            ylabel('Surface apriori [NO_2]');
            set(gca,'fontsize',16);
            %}
            
            % This is to answer why the difference in old and new visible
            % columns gets so large
            xstr = 'Cloud radiance fraction';
            ystr = 'Cloud pressure';
            zstr = 'Percent change in vis NO_2';
            cstr = 'Surface apriori [NO_2] (ppb)';
            if plot_gui
                scatter3_slice_gui(cld_rad_frac_old, cld_pres_old, rdel, [], apriori_old*1e9,...
                    'xlabel',xstr,'ylabel',ystr,'zlabel',zstr,'clabel',cstr);
            else
                fig_s3=figure;
                scatter3(cld_rad_frac_old, cld_pres_old, rdel, [], apriori_old*1e9);
                cb = colorbar;
                xlabel(xstr);
                ylabel(ystr);
                zlabel(zstr);
                cb.Label.String = cstr;
                set(gca,'fontsize',16);
            end
            
            % This should show that only when CF < 1 do we get really large
            % changes, otherwise, the old and new vis VCDs should be
            % similar.
            fig_rad_v_geo = figure; 
            scatter(cld_rad_frac_old, cld_frac_old, [], rdel);
            line([0 1], [0 1], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
            xlabel('Cloud radiance fraction');
            ylabel('Cloud fraction');
            cb = colorbar;
            cb.Label.String = 'Percent change in vis NO_2';
            caxis([0 100]);
            set(gca,'fontsize',16);
            
            xx_clear = cld_frac_old == 0 & cld_rad_frac_old ==0;
            fig_clear = figure;
            plot(no2_vis_old(xx_clear), no2_vis_new(xx_clear), 'ko');
            xlabel('v2.1C V_{vis}');
            ylabel('v3.0A V_{vis}');
            title('Wholly clear conditions');
            set(gca,'fontsize',16);
            
            xx_cloudy = cld_frac_old >= 1.0 & cld_rad_frac_old >= 1.0;
            fig_cloudy = figure;
            plot(no2_vis_old(xx_cloudy), no2_vis_new(xx_cloudy), 'ko');
            xlabel('v2.1C V_{vis}');
            ylabel('v3.0A V_{vis}');
            title('Wholly cloudy conditions');
            set(gca,'fontsize',16);
            
            % This is to answer why the new visible column is sometimes 0
            % when the total is not.
            xx = no2_vis_new == 0 & no2_tot_new ~= 0;
            fig_tot_v_vis = figure;
            scatter(cld_frac_old(xx), cld_pres_old(xx), [], no2_vis_new(xx));
            cb = colorbar;
            cb.Label.String = 'New vis NO2';
            xlabel('Cloud fraction');
            ylabel('Cloud pressure');
            set(gca,'fontsize',16);
        end
        
        function [fig_rel, fig_abs] = plot_bsa_to_brdf_indiv_pixel_changes(regenerate_from_behr_files)
            
            if ~exist('regenerate_from_behr_files', 'var')
                regenerate_from_behr_files = true;
            end
            
            concatenated_data_file = fullfile(misc_behr_update_plots.surf_refl_root, 'IndividualPixels_BSA_BRDF.mat');
            
            if regenerate_from_behr_files    
                dates.jan_feb = {'startdate', '2012-01-01', 'enddate', '2012-02-29'};
                dates.dec = {'startdate', '2012-12-01', 'enddate', '2012-12-31'};
                dates.jun_jul_aug = {'startdate', '2012-06-01', 'enddate', '2012-08-31'};
                
                data.bsa_v5 = struct('path', fullfile(misc_behr_update_plots.behr_nasa_only_dir, 'SP_Files'));
                data.bsa_v6 = struct('path', fullfile(misc_behr_update_plots.behr_modisv6_dir, 'SP_Files'));
                data.brdf_v6 = struct('path', fullfile(misc_behr_update_plots.behr_nasa_brdfD_dir, 'SP_Files'));
                
                dates_fns = fieldnames(dates);
                data_fns = fieldnames(data);
                
                for a=1:numel(dates_fns)
                    fprintf('Loading ocean flags for %s\n', upper(dates_fns{a}));
                    ocean_flag.(dates_fns{a}) = cat_sat_data(data.brdf_v6.path, 'AlbedoOceanFlag', dates.(dates_fns{a}){:});
                    for b=1:numel(data_fns)
                        fprintf('Loading %s for %s\n', data_fns{b}, upper(dates_fns{a}));
                        data.(data_fns{b}).(dates_fns{a}) = cat_sat_data(data.(data_fns{b}).path, 'MODISAlbedo', dates.(dates_fns{a}){:});
                    end
                end
                
                SRData.bsa_v5_jja = data.bsa_v5.jun_jul_aug;
                SRData.bsa_v6_jja = data.bsa_v6.jun_jul_aug;
                SRData.brdf_v6_jja = data.brdf_v6.jun_jul_aug;
                SRData.ocean_flag_jja = ocean_flag.jun_jul_aug;
                
                SRData.bsa_v5_djf = cat(1, data.bsa_v5.jan_feb, data.bsa_v5.dec);
                SRData.bsa_v6_djf = cat(1, data.bsa_v6.jan_feb, data.bsa_v6.dec);
                SRData.brdf_v6_djf = cat(1, data.brdf_v6.jan_feb, data.brdf_v6.dec);
                SRData.ocean_flag_djf = cat(1, ocean_flag.jan_feb, ocean_flag.dec);
                
                save(concatenated_data_file, '-v7.3', '-struct', 'SRData');
            else
                SRData = load(concatenated_data_file); 
            end
            
            v5v6_adel_djf = SRData.bsa_v6_djf(~SRData.ocean_flag_djf) - SRData.bsa_v5_djf(~SRData.ocean_flag_djf);
            v5v6_rdel_djf = reldiff(SRData.bsa_v6_djf(~SRData.ocean_flag_djf), SRData.bsa_v5_djf(~SRData.ocean_flag_djf), true)*100;
            
            v5v6_adel_jja = SRData.bsa_v6_jja(~SRData.ocean_flag_jja) - SRData.bsa_v5_jja(~SRData.ocean_flag_jja);
            v5v6_rdel_jja = reldiff(SRData.bsa_v6_jja(~SRData.ocean_flag_jja), SRData.bsa_v5_jja(~SRData.ocean_flag_jja), true)*100;
            
            bsa_brdf_adel_djf = SRData.brdf_v6_djf(~SRData.ocean_flag_djf) - SRData.bsa_v6_djf(~SRData.ocean_flag_djf);
            bsa_brdf_rdel_djf = reldiff(SRData.brdf_v6_djf(~SRData.ocean_flag_djf), SRData.bsa_v6_djf(~SRData.ocean_flag_djf), true)*100;
            
            bsa_brdf_adel_jja = SRData.brdf_v6_jja(~SRData.ocean_flag_jja) - SRData.bsa_v6_jja(~SRData.ocean_flag_jja);
            bsa_brdf_rdel_jja = reldiff(SRData.brdf_v6_jja(~SRData.ocean_flag_jja), SRData.bsa_v6_jja(~SRData.ocean_flag_jja), true)*100;
            
            fig_rel = figure;
            % Plot without outliers because they span a huge range (that's
            % what " 'symbol', '' " does).
            boxplot(padcat(2, v5v6_rdel_jja, v5v6_rdel_djf, bsa_brdf_rdel_jja, bsa_brdf_rdel_djf), 'symbol','');
            set(gca,'XTickLabel',{'BSA v6 - BSA v5, JJA', 'BSA v6 - BSA v5, DJF', 'BRF v6 - BSA v6, JJA', 'BRF v6 - BSA v6, DJF'},...
                'fontsize', 16,...
                'XTickLabelRotation',30,...
                'YGrid','on');
            ylabel('%\Delta surface reflectance');
            ylim([-80 80]);
            
            % The absolute differences can keep the outliers because they
            % won't be quite so ridiculous
            fig_abs = figure;
            boxplot(padcat(2, v5v6_adel_jja, v5v6_adel_djf, bsa_brdf_adel_jja, bsa_brdf_adel_djf));
            set(gca,'XTickLabel',{'BSA v6 - BSA v5, JJA', 'BSA v6 - BSA v5, DJF', 'BRDF v6 - BSA v6, JJA', 'BRDF v6 - BSA v6, DJF'},...
                'fontsize', 16,...
                'XTickLabelRotation',30,...
                'YGrid','on');
            ylabel('%\Delta surface reflectance');
        end
        
        function plot_dc3_avg_profiles()
            Opts = struct('match_file', 'WRF-DC3-Matched-Data.mat', 'match_dir', '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/WRF',...
                'quantity', 'NO2', 'lats_filter', 'all', 'strat_filter', false, 'fresh_filter', false);
            profs = misc_wrf_chem_comp_plots('make-hybrid',Opts);
            figure;
            plot(profs.dc3_prof*1e3, profs.dc3_pres);
            hold on
            plot(profs.wrf_prof*1e3, profs.wrf_pres);
            set(gca,'fontsize',16,'ydir','reverse');
            xlabel('[NO_2] (ppbv)');
            ylabel('Pressure (hPa)');
            legend('DC3','WRF');
        end
        
        function plot_sl7_rmse(redo_rmse_calcs)
            if ~exist('redo_rmse_calcs','var')
                redo_rmse_calcs = true;
            end
            
            sl7_rmse_file = fullfile(misc_behr_update_plots.wrf_root,'sacramento_sl7_rmse.mat');
            
            if redo_rmse_calcs
                pbase = '/Volumes/share-wrf1/FocusedRuns/Sacramento-SL6-check/DateWrongCheck';
                pnew = '/Volumes/share-wrf1/FocusedRuns/Sacramento-BEARPEX-FullR2SMH-12km/Output';
                [W.rmse, W.wrf_dates] = wrf_rmse(pbase,pnew,{'no','no2','o3'},'normalization','none');
                save(sl7_rmse_file, '-struct', 'W');
            else
                W = load(sl7_rmse_file);
            end
            
            figure;
            hold on
            plot(W.wrf_dates, W.rmse.no*1e6);
            plot(W.wrf_dates, W.rmse.no2*1e6);
            plot(W.wrf_dates, W.rmse.o3*1e6);
            xlim([min(W.wrf_dates), max(W.wrf_dates)]);
            legend('NO','NO_2','O_3');
            datetick('x', 'mm/dd', 'keepticks','keeplimits');
            ylabel('RMSE (pptv)');
            set(gca,'fontsize',16,'yscale','log','ygrid','on')
            
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
        
        function average_behr(data_dir, data_field, overwrite, use_new_avg)
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
                    % The new files may have different name patterns than
                    % the old ones (e.g. final has the pattern
                    % OMI_BEHR-PROFS_REGION_v3-0A_yyyymmdd)
                    [no2_vcds, lon_grid, lat_grid] = behr_time_average(djf_start, djf_end, lon_lim, lat_lim,...
                        'avgfield', data_field, 'behr_dir', monthly_dir, 'filepattern', 'OMI_BEHR*.mat');
                    count_grid = [];
                    avg_config = [];
                else
                    % Since I've moved no2_column_map_2014 into the private
                    % subdirectory of the analysis repo, this class can
                    % access it, but the rest of the Matlab functions
                    % cannot. I also put the dependencies that changed in
                    % there two, suitably renamed.
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
                        'avgfield', data_field, 'behr_dir', monthly_dir, 'filepattern', 'OMI_BEHR*.mat');
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
            F = dir(fullfile(daily_dir, '*.mat'));
            if ~isempty(F)
                save_name = fullfile(data_dir, sprintf('%s-DJF-Daily.mat', data_field));
                if ~overwrite && exist(save_name, 'file')
                    fprintf('%s exists\n', save_name);
                else
                    if use_new_avg
                        [no2_vcds, lon_grid, lat_grid] = behr_time_average(djf_start, djf_end, lon_lim, lat_lim,...
                            'avgfield', data_field, 'behr_dir', daily_dir, 'filepattern', 'OMI_BEHR*.mat');
                        count_grid = [];
                        avg_config = [];
                    else
                        file_stem = regexp(F(1).name, 'OMI_BEHR.+(?=\d\d\d\d\d\d\d\d)', 'match', 'once');
                        [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(djf_start, djf_end, lon_lim, lat_lim,...
                            'mapfield', data_field, 'behrdir', daily_dir, 'fileprefix', file_stem,...
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
                            'avgfield', data_field, 'behr_dir', daily_dir, 'filepattern', 'OMI_BEHR*.mat');
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
        
        function [title_str, quantity_str, unit_str, time_period] = get_title_from_file_name(base_filename, base_dir, new_filename, new_dir)
            % Parses one of the average .mat files names to figure out what
            % the title string should be and the colorbar limits.
            E = JLLErrors;
            
            is_diff = nargin > 2;
            
            [~, base_filename] = fileparts(base_filename); % remove any path or extension
            [~, base_dir] = fileparts(base_dir);
            if is_diff
                [~, new_filename] = fileparts(new_filename);
                [~, new_dir] = fileparts(new_dir);
            end
            
            base_name_parts = strsplit(base_filename, '-');
            new_name_parts = strsplit(new_filename, '-');
            % Quantity - shorten it quite a bit so the title isn't insanely
            % long
            switch base_name_parts{1}
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
                    E.notimplemented('No quantity title for data field "%s"', base_name_parts{1});
            end
            
            % Time period (DJF or JJA), should always be the second part of
            % the file name
            time_period = base_name_parts{2};
            
            % Monthly or daily profiles used
            if ~is_diff || strcmp(base_name_parts{3}, new_name_parts{3});
                profs_used = sprintf('%s Profs', base_name_parts{3});
            else
                profs_used = sprintf('%s vs %s Profs', new_name_parts{3}, base_name_parts{3});
            end
            
            if is_diff
                title_str = sprintf('%s: %s vs %s (%s, %s)', quantity_str, new_dir, base_dir, time_period, profs_used);
            else
                title_str = sprintf('%s: %s (%s, %s)', quantity_str, base_dir, time_period, profs_used);
            end
        end
        
        function fig = plot_change_avg(old_lon_grid, old_lat_grid, old_val, new_lon_grid, new_lat_grid, new_val, diff_type, title_str, quantity_name, unit_name, save_dir)
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
            set(gca,'fontsize',14,'xtick',-120:10:-70,'ytick',25:5:50); % ensure consistent ticks
            
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Other utility functions %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function full_dir = fullpath_incr_dir(dir_in)
            % If it does not start with "/" or "./", then consider it just
            % the end folder of the increment directory and prepend the
            % usual increment directory front path. Otherwise, it's either
            % a full absolute or relative path and should be left as-is
            if regcmp(dir_in, '^\.?/')
                full_dir = dir_in;
                return
            end
            
            % Assume that the "final" directory will split into the proper
            % increment path stem and the final subdir
            incr_dir_stem = fileparts(misc_behr_update_plots.behr_final_dir);
            full_dir = fullfile(incr_dir_stem, dir_in);
        end
        
        function Data = apply_outside(Data, data_dir)
            use_outside_dirs = {misc_behr_update_plots.behr_v2_1C_dir, misc_behr_update_plots.behr_nasa_only_dir,...
                misc_behr_update_plots.behr_nasa_brdf_dir, misc_behr_update_plots.behr_nasa_brdfD_dir,...
                misc_behr_update_plots.behr_nasa_brdf_vis_dir, misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir,...
                misc_behr_update_plots.behr_nasa_brdf_vis_profs_tempfix_dir, misc_behr_update_plots.behr_nasa_brdf_vis_profs_wrftemp_dir};
            
            if ismember(data_dir, use_outside_dirs)
                outside = load(fullfile(behr_analysis_repo_dir, 'Utils', 'outside.mat'));
                outside = ~logical(outside.outside);
                Data.no2_vcds(outside) = nan;
            end
        end
        
        function [file_list] = make_common_files_list(list_of_files, list_to_match)
            % Compare lists of old and new average files. Find common data
            % fields and time periods, but don't match for monthly/daily
            % profiles
            new_files_tmp = regexprep(list_to_match, '(Daily|Monthly)', '');
            old_files_tmp = regexprep(list_of_files, '(Daily|Monthly)', '');
            xx_base = find_common_elements(old_files_tmp, new_files_tmp, 'nodup');
            file_list = list_of_files(xx_base);
        end
    end
    
end
