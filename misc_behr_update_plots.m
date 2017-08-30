classdef misc_behr_update_plots
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        start_date = '2012-01-01';
        end_date = '2012-12-31';
        behr_final_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/5-PSM';
        behr_nasa_brdf_vis_profs_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/4-NewProfs';
        behr_nasa_brdf_vis_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/3-NewVis';
        behr_nasa_brdf_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/2-BRDF';
        behr_nasa_only_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/1-NASAv3';
        behr_v2_1C_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/0-v2.1C';
        
        figs_root_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Figures';
    end
    
    methods(Static = true)
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
                struct('dir', misc_behr_update_plots.behr_final_dir, 'use_new_avg', true, 'data_fields', {{'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'}})...
            );
        
            for a=1:numel(average_info)
                for b=1:numel(average_info(a).data_fields)
                    if average_info(a).use_new_avg
                        % pass for now
                    else
                        misc_behr_update_plots.average_behr_old(average_info(a).dir, average_info(a).data_fields{b}, do_overwrite);
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
                         misc_behr_update_plots.behr_final_dir};
                     
            for a=1:numel(incr_dirs)-1
                Fbase = dir(fullfile(incr_dirs{a}, '*.mat'));
                base_names = {Fbase.name};
                Fnew = dir(fullfile(incr_dirs{a+1}, '*.mat'));
                new_names = {Fnew.name};
                for b=1:numel(base_names)
                    if ismember(base_names{b}, new_names)
                        Dbase = load(fullfile(incr_dirs{a}, base_names{b}));
                        Dnew = load(fullfile(incr_dirs{a+1}, base_names{b}));
                        
                        abs_del = Dnew.no2_vcds - Dbase.no2_vcds;
                        rel_del = reldiff(Dnew.no2_vcds, Dbase.no2_vcds)*100;
                        
                        misc_behr_update_plots.plot_change_avg(Dbase.lon_grid, Dbase.lat_grid, abs_del, false, incr_dirs{a}, incr_dirs{a+1}, base_names{b});
                        misc_behr_update_plots.plot_change_avg(Dbase.lon_grid, Dbase.lat_grid, rel_del, true, incr_dirs{a}, incr_dirs{a+1}, base_names{b});
                    end
                end
            end
        end
    end
    
    methods(Static = true, Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper functions for the production of the actual incr. changes %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function make_behr_with_parameters(req_commit, root_save_dir, do_overwrite, has_daily_profs)
            % Runs BEHR, saving files in ROOT_SAVE_DIR and overwriting
            % according to DO_OVERWRITE. This runs BEHR assuming the
            % BEHR_main and read_omno2_v_aug2012 functions take all the
            % directories as parameters.
            E = JLLErrors;
            my_dir = fileparts(mfilename('fullpath'));
            behr_dir = fullfile(my_dir, '..', 'BEHR');
            
            G = GitChecker;
            G.addCommitRange(behr_dir, req_commit, req_commit);
            G.Strict = true;
            
            G.checkState();
            
            if ~exist(root_save_dir, 'dir')
                E.dir_dne(root_save_dir);
            end
            
            sp_save_dir = fullfile(root_save_dir, 'SP_Files');
            if ~exist(sp_save_dir, 'dir')
                mkdir(sp_save_dir);
            end
            
            read_omno2_v_aug2012('start', misc_behr_update_plots.start_date,...
                'end', misc_behr_update_plots.end_date,...
                'sp_mat_dir', sp_save_dir, 'overwrite', do_overwrite);
            
            % For the monthly profiles, we rely on the default profile
            % choice to be monthly. This avoids issues with older
            % increments that do not have the 'profile_mode' parameter.
            monthly_save_dir = fullfile(root_save_dir, 'MonthlyProfs');
            if ~exist(monthly_save_dir, 'dir')
                mkdir(monthly_save_dir);
            end
            BEHR_main('start', misc_behr_update_plots.start_date,...
                'end', misc_behr_update_plots.end_date,...
                'sp_mat_dir', sp_save_dir,...
                'behr_mat_dir', monthly_save_dir,...
                'overwrite', do_overwrite);
            
            if has_daily_profs
                daily_save_dir = fullfile(root_save_dir, 'DailyProfs');
                if ~exist(daily_save_dir, 'dir')
                    mkdir(daily_save_dir);
                end
                BEHR_main('start', misc_behr_update_plots.start_date,...
                    'end', misc_behr_update_plots.end_date,...
                    'sp_mat_dir', sp_save_dir,...
                    'behr_mat_dir', daily_save_dir,...
                    'profile_mode', 'daily',...
                    'overwrite', do_overwrite);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper functions for the incremental change averages %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function average_behr_old(data_dir, data_field, overwrite)
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
                [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(djf_start, djf_end, lon_lim, lat_lim,...
                    'mapfield', data_field, 'behrdir', monthly_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                    'makefig', false);
                fprintf('Saving as %s\n', save_name);
                save(save_name, 'no2_vcds', 'lon_grid', 'lat_grid', 'count_grid', 'avg_config');
            end
            
            save_name = fullfile(data_dir, sprintf('%s-JJA-Monthly.mat', data_field));
            if ~overwrite && exist(save_name, 'file')
                fprintf('%s exists\n', save_name);
            else
                [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(jja_start, jja_end, lon_lim, lat_lim,...
                    'mapfield', data_field, 'behrdir', monthly_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                    'makefig', false);
                
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
                    [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(djf_start, djf_end, lon_lim, lat_lim,...
                    'mapfield', data_field, 'behrdir', daily_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                    'makefig', false);
                
                    fprintf('Saving as %s\n', save_name);
                    save(save_name, 'no2_vcds', 'lon_grid', 'lat_grid', 'count_grid', 'avg_config');
                end
                
                save_name = fullfile(data_dir, sprintf('%s-JJA-Daily.mat', data_field));
                if ~overwrite && exist(save_name, 'file')
                    fprintf('%s exists\n', save_name);
                else
                    [~, no2_vcds, lon_grid, lat_grid, count_grid, avg_config] = no2_column_map_2014(jja_start, jja_end, lon_lim, lat_lim,...
                    'mapfield', data_field, 'behrdir', daily_dir, 'fileprefix', 'OMI_BEHR_v2-1C_',...
                    'makefig', false);
                
                    fprintf('Saving as %s\n', save_name);
                    save(save_name, 'no2_vcds', 'lon_grid', 'lat_grid', 'count_grid', 'avg_config');
                end
            end
        end
        
        function [title_str, quantity_str, unit_str] = get_title_from_file_name(filename, base_dir, new_dir)
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
        
        function plot_change_avg(lon_grid, lat_grid, val_diff, is_rel_diff, base_dir, new_dir, file_name)
            E = JLLErrors;
            
            [title_str, quantity_name, unit_name] = misc_behr_update_plots.get_title_from_file_name(file_name, base_dir, new_dir);
            
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
            
            save_dir = fullfile(misc_behr_update_plots.figs_root_dir, 'IncrementTests');
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

