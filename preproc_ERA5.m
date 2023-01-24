function [data_out] = preproc_ERA5(domain_def)
    
    vars_in={       't2m', 'ssrd', 'pr'                  ,   'z'   ,     'sst'};
    vars_conv={'- 273.15',     '', '*days_of_mon(m)*1000','/9.8'   , '- 273.15'};
%     vars_in={       't2m', 'ssrd', 'pr'                  ,   'u'};
%     vars_conv={'- 273.15',     '', '*days_of_mon(m)*1000',   ''};
    % z: specifically referring to geopotential height at 500 hPa, original data is
    % geopotential height, need to /9.8 to convert it to the height
    % u: specifically referring to zonal wind speed at 300 hPa,
    
    %             1  2  3  4  5  6  7  8  9 10 11,12
    days_of_mon=[31,28,31,30,31,30,31,31,30,31,30,31];
    
    % Local climate data from ERA5
    global data_src;
    global chosen_prd;
    
    % Input variables
    % t2m
    
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'ERA5', 'RS_05');
    
    % period for complete annual data
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    ny=yr_e-yr_s+1;
        
    % check if the array dataset has created
%     if size(DATA_ERA5_05rs_out,2) <= 1
        % Create 
        % Monthly data:              nlon         nlat   nweek/nmonth  nyear   nvariable
        data_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1, size(vars_in,2));
%     end

%     for v=1:1
    for v=1:size(vars_in,2)
        for y=yr_s:yr_e
            % ---- Read data ---
            yi=y-yr_s+1;

            switch vars_in{v}
                case 'z'
                    dir_in=sprintf('%s/ECMWF/ERA5/raw/atmosphere/mon/pressure_levels/remap/0.5deg/500/%s/',data_src,vars_in{v});
                    fname=sprintf('%s_ECMWF-ERA5_rean_pressure_levels_500_mon_%d-%d_0.5deg_remapycon.nc',vars_in{v},y,y);
                case 'sst'
                    dir_in=sprintf('%s/ECMWF/ERA5/raw/surface/mon/single_level/remap/0.5deg/%s',data_src,vars_in{v});
                    fname=sprintf('%s_ECMWF-ERA5_rean_single_level_mon_%d-%d_0.5deg_remapycon.nc',vars_in{v},y,y);
                otherwise
                    dir_in=sprintf('%s/ECMWF/ERA5-Land/remap/0.5deg/%s',data_src,vars_in{v});
                    fname=sprintf('%s_ECMWF-ERA5-Land_rean_mon_%d-%d_0.5deg_remapycon.nc',vars_in{v},y,y);
            end

            fprintf(sprintf('--> Processing file : %s\n',fname));

            fullname=sprintf('%s/%s',dir_in, fname);
            nc_var = ncgeodataset(fullname);

            % ---- Climate variables  ---
            % vars (month, lat, lon)
            dtmp=squeeze(double(nc_var.data(vars_in{v})));
            [s1 s2 s3]=size(dtmp);
            dtmp_uconv=nan(s2,s3);
            
            for m=1:12
                % Unit conversion
                dtmp_uconv=squeeze(eval(['dtmp(m,:,:)' vars_conv{v}]));
                
                data_out(:,:, m, yi, v)=squeeze(dtmp_uconv(y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
            end

            clearvars dtmp dir_in fname;
        end
    end
    
    % Temporal code for u at 300hpa, to be deleted
    fname='/Data/obs/ECMWF/ERA5/raw/atmosphere/mon/pressure_levels/remap/0.5deg/300/u/u_300hPa_global_1979_2020_lon-180_lon180.nc';
    nc_var = ncgeodataset(fname);
    
	fprintf(sprintf('--> Processing U300 : %s\n',fname));
    dtmp=squeeze(double(nc_var.data('u')));
    [nm nlat nlon]=size(dtmp);
    v=v+1;
    for y=yr_s:yr_e
        for m=1:12
            mi=(y-1979)*12+m;
            data_out(:,:, m, y-yr_s+1, v)=squeeze(dtmp(mi,y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))'; % Note the rotate sign here
        end
    end
    
	fprintf(sprintf('--> Processing zonal-mean of U300 : %s\n',fname));
    % ======== THIS IS BASED ON MONTHLY DATA ========
    if 1==2
        
        % Calculate the zonal wind speed of the defined domain
        dtmp_lat=squeeze(double(nc_var.data('lat'))); % start from negative values
        dtmp_lon=squeeze(double(nc_var.data('lon'))); % start from negative values

        % lati and loni of the defined JLI domain
        lat_def=[15 75];lon_def=[-60 40];
        dtmp_zm=nan(nm, nlat, nlon);

        lati_rng=nan(2,1);loni_rng=nan(2,1);
        for n=1:2
            % get lat index
            for i=1:length(dtmp_lat)
                if dtmp_lat(i)>=lat_def(1,n)
                    lati_rng(n,1)=i;
                    break;
                end
            end

            % get lon index
            for i=1:length(dtmp_lon)
                if dtmp_lon(i)>=lon_def(1,n)
                    loni_rng(n,1)=i;
                    break;
                end
            end
        end

        for m=1:nm
            for i=1:nlat
                for j=1:nlon
                    if i >= lati_rng(1,1) && i <= lati_rng(2,1) && j >= loni_rng(1,1) && j <= loni_rng(2,1)
                        dtmp_zm(m,i,j)=squeeze(nanmean(dtmp(m,i,loni_rng(1,1):loni_rng(2,1)),3));
                    end
                end
            end
        end
        
        v=v+1;
        for y=yr_s:yr_e
            for m=1:12
                mi=(y-1979)*12+m;
                data_out(:,:, m, y-yr_s+1, v)=squeeze(dtmp_zm(mi,y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))'; % Note the rotate sign here
            end
        end
    end

    % ======== THIS IS BASED ON MONTH-MEAN OF DAILY DATA ========
    fname='./data/u_zon_mon_mean_1979_2020.nc'; % monthly-mean data
    nc_var = ncgeodataset(fname);
	dtmp_lat=squeeze(double(nc_var.data('lat'))); % 
	dtmp_u=squeeze(double(nc_var.data('u'))); %
    
    % NOTE: + 0.5 here to align the number of lat
    ind_dmn=dtmp_lat >= (domain_def(3)+0.5) & dtmp_lat <= domain_def(4);

    % get first and last lat index
    xs=min(find(ind_dmn==1));
    xe=max(find(ind_dmn==1));
    dtmp_u_array=permute(repmat(dtmp_u(:,xs:xe),[1 1 dmn_lon_n_g]), [1 3 2]);
    
    v=v+1;
    for y=yr_s:yr_e
        for m=1:12
            mi=(y-1979)*12+m;
            data_out(:,:, m, y-yr_s+1, v)=squeeze(dtmp_u_array(mi,:,:)); % Note the rotate sign here
        end
    end
end