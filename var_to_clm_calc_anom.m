function []=var_to_clm_calc_anom(dsid)
    % ------------------------------------------------------------------
    % Correlation coefficient between variable of interests and teleconnection climate indices
    % ------------------------------------------------------------------
    
    % --- Input ---
    global DATA_VI_05rs_dtds_out_gs1 DATA_CLMINX_out_gs1;
    global DATA_CLM_05rs_dtds_out_gs1;
    global lgs_map;
    
    % --- Output ---
	global vh_clm_anom vl_clm_anom;
    
    dtmp(:,:,:,:,1)=DATA_VI_05rs_dtds_out_gs1(:,:,:,:,dsid); % VI anomalies
    dtmp(:,:,:,:,2)=DATA_CLMINX_out_gs1(:,:,:,:,1,1); % Climate Index anomalies, only lag-0 of the first climate index

    dtmp(:,:,:,:,3)=DATA_CLM_05rs_dtds_out_gs1(:,:,:,:,4); % now using zonal mean of 300hPa u anomalies
%     dtmp(:,:,:,:,3)=DATA_CLMINX_out_gs1(:,:,:,:,1,2); % Jet Speed Index, only lag-0
    
    [s1 s2 dummpy nyr dummpy]=size(dtmp);
    
    % Definition of seasons (winter is actually starting from November, 1st month of the growing season)
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[],[] [] [] [] []};
    nm=size(m_rng,2); % number of months/seasons
    nclmv=2; % now only include climate index and wind speed/jet speed index
    nexi=4; % around 10% of the sample size, which is 32 values for 32-years growing season 
            % (one mean value for a subperiod of a year)
    
    vh_clm_anom=nan(s1,s2,nm,nclmv);
    vl_clm_anom=nan(s1,s2,nm,nclmv);
    
	for i=1:s1
        for j=1:s2
            for clmv=1:nclmv
                for m=1:nm

                    m_lgs=lgs_map(i,j,dsid);
                    [m_rng{13} m_rng{14} m_rng{17} m_rng{18}]=get_ssn12(m_lgs); % also get the end of growing season
                    m_rng{15}=[m_rng{13} m_rng{14}];  % all growing season months
                    m_rng{16}=[2 3];  % first two growing season

                    % get the top and lowest VI anomalies
                    v_mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m},:,1),3));
                    if sum(isnan(v_mon_ts))>=10
                        continue;
                    end
                    
                    [dummy,maxI] = maxk(v_mon_ts,nexi);
                    [dummy,minI] = mink(v_mon_ts,nexi);

                    % get the corresponding climate composites
                    vh_clm_anom(i,j,m,clmv)=nanmean(nanmean(dtmp(i,j,m_rng{m},maxI,clmv+1),3),4);
                    vl_clm_anom(i,j,m,clmv)=nanmean(nanmean(dtmp(i,j,m_rng{m},minI,clmv+1),3),4);

                end % month
            end %clmv
        end % j
	end % i
end
