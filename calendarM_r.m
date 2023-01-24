function calendarM_r(clmidx_name)
    dsid=2;
    global DATA_VI_05rs_sm_out DATA_SMRoot_05rs_out DATA_ERA5_05rs_out;
    global chosen_prd nyr;
    global lons lats ax_show;
    global clmidx_name;

    data_type='daily';
    if strcmp(data_type,'daily')
        ndays_m=[31 28 31 30 31 30 31 31 30 31 30 31];
        eday_m=cumsum(ndays_m);sday_m=eday_m-ndays_m+1;
        % Preprocessing teleconnection climate indices

        filename=sprintf('%s/data/era5_jet_indices_North_Atlantic_1979_2020.nc', '.');
        nc_var = ncgeodataset(sprintf('%s',filename));
        jli_tmp=nc_var.data(clmidx_name)';

        % Convert daily data to monthly
        jli_rs=reshape(jli_tmp, [365 42]);
        for m=1:12
            jli_rs_m(m,:)=nanmean(jli_rs(sday_m(m):eday_m(m), :));
        end
    else
        if strcmp(data_type,'monthly')
%             filename=sprintf('%s/data/era5_jet_indices_mm_North_Atlantic_1979_2020_corr.nc', '.');
            filename=sprintf('%s/data/era5_jet_indices_mmu_North_Atlantic_1979_2020.nc', '.');
            nc_var = ncgeodataset(sprintf('%s',filename));
            jli_tmp=nc_var.data(clmidx_name)';

            jli_rs_m=reshape(jli_tmp,[12 length(jli_tmp)/12]);
        end
    end

    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    data_sy=1979;
    dtmp=permute(repmat(jli_rs_m,[1 1 92 77]),[3 4 1 2]);
    DATA_CLMINX_rs_out=dtmp(:,:,:,(yr_s-data_sy+1):(yr_e-data_sy+1));
    DATA_SMRoot_05rs_out1=reshape(DATA_SMRoot_05rs_out,[92 77 12 33]);

    % Detrend and Deseasonalize
    DATA_CLMINX_rs_dtds_out=data_dsdt(DATA_CLMINX_rs_out);
    DATA_ERA5_05rs_dtds_out=data_dsdt(DATA_ERA5_05rs_out);
    DATA_SMRoot_05rs_dtds_out=data_dsdt(DATA_SMRoot_05rs_out1);

    DATA_05rs_dtds_out=nan(92,77,12,33,4);
    DATA_05rs_dtds_out(:,:,:,:,1:2)=DATA_ERA5_05rs_dtds_out(:,:,:,:,1:2);
    DATA_05rs_dtds_out(:,:,:,:,3)=DATA_SMRoot_05rs_dtds_out;
    DATA_05rs_dtds_out(:,:,:,:,4)=data_dsdt(DATA_VI_05rs_sm_out(:,:,:,:,dsid,4));

    %% Calculate the correlation
    [s1 s2 s3 s4 s5]=size(DATA_CLMINX_rs_dtds_out);
    v_clminx_r_cln_m=nan(92,77,16,6);
    v_clminx_p_cln_m=nan(92,77,16,6);
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[12 1 2],[3 4 5] [6 7 8] [9 10 11]};
    for v=1:4
        fprintf(sprintf('--> Processing variable %d\n',v));
        for i=1:s1
            for j=1:s2
                for m=1:16

                   clminx_mon_ts=squeeze(nanmean(DATA_CLMINX_rs_dtds_out(i,j,m_rng{m},:),3));

                   var_mon_ts=squeeze(nanmean(DATA_05rs_dtds_out(i,j,m_rng{m},:,v),3));

                   if (nansum(abs(var_mon_ts))>0) && sum(isnan(clminx_mon_ts))<nyr/2 % skip the ocean and nan values

                      [r, p]=corrcoef(var_mon_ts,clminx_mon_ts);
                      if abs(r(1,2))>0
                          v_clminx_r_cln_m(i,j,m,v)=r(1,2);
                          v_clminx_p_cln_m(i,j,m,v)=p(1,2);
                      end

                   end
                end
            end
        end
    end

    %% Plotting by seasons
    y_pos=0.40;
    x_pos=-0.1;
    c_thr=0.600001;
    
    figure('color','white','Position',[781   136   716   833]);
    gap=[0.02 0.02]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(4,4,gap,marg_h,marg_w);
    yscal=0.96; xscal=0.015;
    cb_xR=0.8;
    cb_yR=0.4;
    cb_wR=4;
    ax_show.frame=0;
    ncb=25;
    
    var_nam={'Temp.','Rad.','SM','EVI2'};
    ssn_nam={'DJF','MAM','JJA','SON'};
    labels={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)', 'n)','o)','p)'};
    
    [s1 s2 s3 s4]=size(v_clminx_p_cln_m);
    lons_sig=repmat(lons,[1 1 s3 s4]);
    lats_sig=repmat(lats,[1 1 s3 s4]);
    lons_sig(isnan(v_clminx_p_cln_m))=nan;
    lats_sig(isnan(v_clminx_p_cln_m))=nan;
    lons_sig(v_clminx_p_cln_m>=0.05)=nan;
    lats_sig(v_clminx_p_cln_m>=0.05)=nan;
    intv=2;
    
	n=0;
    for v=1:4
        f=0;
        for m=[2 3 4 1]
            n=n+1;
            f=f+1;
            axes(ha((v-1)*4+f)); 
            
            bg_show=squeeze(flipud(rot90(v_clminx_r_cln_m(:,:,m+12,v))));
            bg_show(isnan(bg_show))=0;
            geoplot(ax_show,bg_show);
            crng=[-c_thr c_thr];
            colormap(gca, flipud(cbrewer('div','RdBu',ncb)));
            caxis([crng(1) crng(2)]);
            
            % Mark signficant points
            a=squeeze(lats_sig(1:intv:end,1:intv:end,m+12,v));
            b=squeeze(lons_sig(1:intv:end,1:intv:end,m+12,v));
            plotm(b(:), a(:), 'Color',[50 50 50]/255,'LineStyle','none', 'Marker', '+', 'MarkerSize',3, 'LineWidth',1);
            
            if v==1
                title(sprintf('%s',ssn_nam{m}));
            end

            if f==1
                a=get(gca);
                gca_w = (a.XLim(2)-a.XLim(1));
                gca_h = (a.YLim(2)-a.YLim(1));
                text(a.XLim(1)+gca_w*x_pos, a.YLim(1)+gca_h*y_pos,sprintf('%s~%s(r)',var_nam{v},upper(clmidx_name)),'FontSize',11,'FontName','Arial','Rotation',90);
            end

            % labels for sub-plots
            yscal=0.92; xscal=0.03;
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{n},'FontSize',11,'FontName','Arial');
            
            if v==4 && m==f
                cb3=colorbar('Southoutside');
                resizeCB(cb3, cb_xR, cb_yR, cb_wR, 0.6, '',19,10,8);
            end
            
        end
    end
    %%
end

