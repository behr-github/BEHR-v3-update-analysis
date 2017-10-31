classdef misc_behr_v3_validation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties(Constant = true)
        behr_v2_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_v2-1C';
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
        
        function value = spiral_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', 'spiral-structs.mat');
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main validation methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                        other_campaigns = {'seac4rs', 'dc3', 'arctas-carb', 'intex-b', 'soas'};
                        time_windows = [1.5, 3]; % how far away from satellite overpass (in hours) the profile is allowed to be
                        
                        behr_comp = struct();
                        all_campaigns = [spiral_campaigns, other_campaigns];
                        campaign_use_ranges = [repmat({''}, size(spiral_campaigns)), repmat({'ranges'}, size(other_campaigns))];
                        for a=1:numel(spiral_campaigns)
                            campaign = struct();
                            for b=1:numel(time_windows)
                                [prof_times, time_fieldname] = misc_behr_v3_validation.format_profile_time_range(time_windows(b));
                                comparison_params = {'campaign', all_campaigns{a},...
                                    'profile_input', campaign_use_ranges{a},...
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
                            behr_comp.(spiral_campaigns{a}) = campaign;
                        end
                        
                        fn = sprintf('%s_%s', regions{r}, prof_modes{p});
                        version_comp.(fn) = behr_comp;
                    end
                end
                insitu_comparison.(versions{v}) = version_comp;
            end
            
            if do_save
                save(misc_behr_v3_validation.spiral_comp_file, '-struct', 'insitu_comparison');
            end
        end
        
        function [times, time_fieldname] = format_profile_time_range(win)
            base_time = datenum('2000-01-01 13:30:00');
            start_time = base_time - win/24;
            end_time = base_time + win/24;
            
            times = {datestr(start_time, 'HH:MM'), datestr(end_time, 'HH:MM')};
            time_fieldname = sprintf('t%s_%s', datestr(start_time, 'HHMM'), datestr(end_time, 'HHMM'));
        end
    end
    
end

