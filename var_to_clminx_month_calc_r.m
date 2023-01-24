function []=var_to_clminx_month_calc_r(dsid)
    % ------------------------------------------------------------------
    % Correlation coefficient between variable of interests and teleconnection climate indices
    % ------------------------------------------------------------------
    
    % --- Input ---
    global DATA_VI_05rs_dtds_out_gs1 DATA_CLMINX_out_gs1;
    global lgs_map;
    global DATA_CLM_05rs_dtds_out_gs1;
    
    % --- Output ---
	global v_clminx_r v_clminx_p;
    
%     dtmp=DATA_VI_05rs_dtds_out_gs1; % Growing season months
    dtmp=DATA_CLM_05rs_dtds_out_gs1;
    dtmp(:,:,:,:,4)=DATA_VI_05rs_dtds_out_gs1(:,:,:,:,dsid);

    [s1 s2 s3 nyr s5]=size(dtmp);
    % Definition of seasons (winter is actually starting from November, 1st month of the growing season)
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[],[] [] [] [] []};
    nm=size(m_rng,2); % number of months/seasons
    nv=1; % Only for Jet Stream Index
    nclmv=4; % including EVI2
    nlag=6;
    v_clminx_r=nan(s1,s2,nm,nlag, nclmv);
    v_clminx_p=nan(s1,s2,nm,nlag, nclmv);
    for clmv=1:4
        for v=1:nv % not use now

            fprintf(sprintf('=== > Variable : %d \n',clmv));

            for lag=0:(nlag-1)
                for m=1:nm
                    clminx_mon_ts=nan;

                    fprintf(sprintf(' * Calculating lag %d for month %d ...\n',lag, m));

                    for i=1:s1
                       for j=1:s2

                           m_lgs=lgs_map(i,j,dsid);
                           [m_rng{13} m_rng{14} m_rng{17} m_rng{18}]=get_ssn12(m_lgs); % also get the end and the peak of growing season
                           m_rng{15}=[m_rng{13} m_rng{14}];  % all growing season months
                           m_rng{16}=[2 3];  % first two growing season

                           v_mon_ts=squeeze(nanmean(dtmp(i,j,m_rng{m},:,clmv),3));
                           % For the growing season month, and seasons
                           clminx_mon_ts=squeeze(nanmean(DATA_CLMINX_out_gs1(i,j,m_rng{m},:,lag+1,v),3));

                           if (nansum(abs(v_mon_ts))>0) && sum(isnan(clminx_mon_ts))<nyr/2 % skip the ocean and nan values

                              [r, p]=corrcoef(v_mon_ts,clminx_mon_ts);
                              v_clminx_r(i,j,m,lag+1,clmv)=r(1,2);
                              v_clminx_p(i,j,m,lag+1,clmv)=p(1,2);

                           end

                       end % j
                    end % i
                end % month
            end % lag
        end % v
    end %clmv
    
end
