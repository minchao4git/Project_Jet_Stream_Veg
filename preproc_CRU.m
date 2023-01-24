function [] = preproc_CRU()
    
%     vars_in={'frs', 't2m','rad'}; % frost days
    vars_in={'frs'}; % frost days
    vars_conv={'' , ''   ,''}; 
    
    % Local climate data from ERA5
    global data_src;
    global domain_def;
    global chosen_prd;
    global DATA_CRU_05rs_out;
    
    % Input variables
    % frs
    
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'CRU', 'RS_05');
    
    % period for complete annual data
    data_sy=1901; % data start year (from Jan.)
    data_ey=2019; % data start year (from Jan.)
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);

        
    % check if the array dataset has created
    if size(DATA_CRU_05rs_out,2) <= 1
        % Create 
        % Monthly data          nlon         nlat   nweek/nmonth  nyear   ndatasource
        DATA_CRU_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1);
    end

    for v=1:size(vars_in,2)
        
        dir_in=sprintf('%s/CRU/CRU-TS404/mon/orig/%s',data_src, vars_in{v});
        
        for y=yr_s:yr_e
            
            % ---- Read data ---
            ndecade=ceil((data_ey-data_sy+1)/10);
            for decade_id=1:ndecade

                data_decade_sy=data_sy+(decade_id-1)*10;
                % cater the name of last decade file, which uses the last
                % data year for file name
                if decade_id==ndecade
                    data_decade_ey=data_ey;
                else
                    data_decade_ey=(data_sy-1)+decade_id*10;
                end

                % Search the corresponding decade file for input, skip the
                % file not in the target period
                if ~(y>=data_decade_sy && y<=data_decade_ey)
                    continue;
                end
                
                fname=sprintf('cru_ts4.04.%d.%d.%s.dat.nc',data_decade_sy,data_decade_ey,vars_in{v});

                fprintf(sprintf('--> Processing year %d in file : %s\n',y, fname));

                fullname=sprintf('%s/%s',dir_in, fname);
                nc_var = ncgeodataset(fullname);

                % ---- Vegetation indices  ---
                % vars (month, lat, lon)
                dtmp=squeeze(double(nc_var.data(vars_in{v})));

                % Unit conversion
                dtmp=eval(['dtmp' vars_conv{v}]);
                
                yi=y-yr_s+1;
                dyi=y-data_decade_sy; % starting from 0
                
                for m=1:12 % to permute, how to loop by month here
%                     DATA_CRU_05rs_out(:,:, 1:12, yi, v)=squeeze(dtmp((dyi*12+1):(dyi+1)*12, y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
                    DATA_CRU_05rs_out(:,:, m, yi, v)=squeeze(dtmp(dyi*12+m, y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
                end
            end
        end
    end
end

% figure;
% imagescn(squeeze(dtmp(1,:,:)));