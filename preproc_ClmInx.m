function [] = preproc_ClmInx(use_gs,dsid, clmidx_name)
    % Input Jet Stream Index/Zonal Jet Index and Jet Speed Index
    % Preprocessing teleconnection climate indices
    
    % --- Input ---
    global chosen_prd nyr;
    
    ni=2; % number of the indics
    
    data_type='daily';
    if strcmp(data_type,'daily')
        ndays_m=[31 28 31 30 31 30 31 31 30 31 30 31];
        eday_m=cumsum(ndays_m);sday_m=eday_m-ndays_m+1;

        filename=sprintf('%s/data/era5_jet_indices_North_Atlantic_1979_2020.nc', '.');
        nc_var = ncgeodataset(sprintf('%s',filename));
        clmidx_tmp=nc_var.data(clmidx_name)';
        
        clmidx_tmp(2,:)=nc_var.data('jsi')'; % Additionally input Jet Speed Index

        % Convert daily data to monthly
        clmidx_rs=reshape(clmidx_tmp, [2 365 42]);
        for n=1:ni
            for m=1:12
                clmidx_rs_m(m,:,n)=nanmean(clmidx_rs(n,sday_m(m):eday_m(m), :),2);
            end
        end
    else
        if strcmp(data_type,'monthly')
            filename=sprintf('%s/data/era5_jet_indices_mmu_North_Atlantic_1979_2020.nc', '.');
            nc_var = ncgeodataset(sprintf('%s',filename));
            clmidx_tmp=nc_var.data(clmidx_name)';

            clmidx_rs_m=reshape(clmidx_tmp,[12 length(clmidx_tmp)/12]);
        end
    end

    % Deseasonalize
    clmidx_m_ds=nan(12*42, ni);
    for n=1:ni
        clmidx_m=clmidx_rs_m(:,:,n);
        clmidx_m=clmidx_m(:); % convert to ts
        a=repmat(squeeze(nanmean(clmidx_rs_m(:,:,n),2)),[1 42]);
        clmidx_m_clm=a(:); % convert to ts climatology
        clmidx_m_ds(:,n)=clmidx_m-clmidx_m_clm; % ds as ts
    end
    
    global DATA_VI_05rs_out;
    global nyr_gs;
    global sos_map lgs_map;
    
    % --- Output ---
    global DATA_CLMINX_out DATA_CLMINX_out_gs DATA_CLMINX_out_gs1;
    
    % Climate indices data with lag
    %                   n lags x n year x 12 variable (year, month, indx1, indx2,
    %                   indx3, indx4, indx5, indx6, indx7, indx8, indx9, indx10)
    nlag=6;
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    data_sy=1979;
    
	[s1 s2 s3 s4 s5]=size(DATA_VI_05rs_out);
    
    % ---- Get growing season data based on the defined growing season
	%                   lon/lat,lon/lat,month*year, nvar
    DATA_CLMINX_out=nan(s1,s2,s3*s4,ni);
    DATA_CLMINX_out_gs=nan(s1,s2,s3,nyr_gs,nlag,ni);
    
    for n=1:ni
        dtmp=squeeze(clmidx_m_ds(((yr_s-data_sy)*12+1):((yr_e-data_sy+1)*12),n));
        DATA_CLMINX_out(:,:,1:s3*s4,n)=permute(repmat(dtmp,[1 s1 s2]),[2 3 1]);
    end
    
    if use_gs

        % Growing season definition based on calculated SOS
        for v=1:ni
            
            fprintf(sprintf('--> Processing climate index : %d\n',v));
            
            for i=1:s1
                for j=1:s2
                    % Growing season definition based on astronomic month starting from December
                    % gssn_s=12;

                    gssn_s=sos_map(i,j,dsid);
                    sosmi=2; % sos month index in array
                    
                    for lag=0:(nlag-1)
                        gssn_s1=gssn_s-(sosmi-1)-lag;
                        
                        if gssn_s1<1
                            gssn_s1=gssn_s1+12;
                        end

                        if ~isnan(gssn_s)
                            
                            for y=1:nyr_gs
                                m_s=(y-1)*12+gssn_s1;
                                m_e=m_s+11;
                                
                                DATA_CLMINX_out_gs(i,j,1:12,y,lag+1,v)=squeeze(DATA_CLMINX_out(i,j,m_s:m_e,v));
                            end % year
                        end
                    end % lag
                    
                end % lon/lat
            end % lon/lat
        end % variables/indics
        
        % season 1 growing season dataset
        DATA_CLMINX_out_gs1=nan(s1,s2,12,nyr_gs,nlag, ni);
        for v=1:ni
            for i=1:s1
                for j=1:s2

                    m_lgs=lgs_map(i,j,dsid);
                    if ~isnan(m_lgs)

                        if m_lgs < 12
                            gs_prd=(1:m_lgs)+(sosmi-1);
                        else
                            gs_prd=1:12;
                        end

                        % mean of actual growing season (season 1)
                        DATA_CLMINX_out_gs1(i,j,gs_prd,:,:,v)=DATA_CLMINX_out_gs(i,j,gs_prd,:,:,v); 
                    end
                end
            end
        end
        
        ny=nyr_gs;
        dtmp=DATA_CLMINX_out_gs1;
    else
        ny=nyr;
        dtmp=DATA_CLMINX_out;
    end
end

