classdef misc_behr_update_plots
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        start_date = '2012-01-01';
        end_date = '2012-01-01';
        behr_final_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/5-PSM';
        behr_nasa_brdf_vis_profs_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/4-NewProfs';
        behr_nasa_brdf_vis_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/3-NewVis';
        behr_nasa_brdf_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/2-BRDF';
        behr_nasa_only_dir = '/Users/Josh/Documents/MATLAB/BEHR-v3-analysis/Workspaces/IncrementTests/1-NASAv3';
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
            G.addReqCommits(fullfile(behr_repo_dir, '..', 'Python Interface'), 'b6bc6e4');
            G.Strict = true;
            G.checkState();
            
            save_dir = misc_behr_update_plots.behr_final_dir;
            misc_behr_update_plots.make_behr_with_parameters('631c643', save_dir, do_overwrite, true);
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
    end
    
    methods(Static = true, Access = private)
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
    end
    
end

