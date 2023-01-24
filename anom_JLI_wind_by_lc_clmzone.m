function anom_JLI_wind_by_lc_clmzone(m)

    % ===== R between Veg. and JLI grouped by landclass and climate zone ===== 
    global lc_dom_grp all_subgrp;
    global clm_cl_x;
    global vh_clm_anom vl_clm_anom;
    global clmidx_name;

    [s1 s2]=size(all_subgrp);

    % Calculation
    data_calc=nan(s1,s2,3);
    clm_cl_x=nan(7,8);

    % number of gridpoint
    clm_cl_cnt=nan(7,8,5);

    % group in bar chart
    clm_cl_rpt=nan(7,8,5);
    clm_cl_stderr=nan(8,3,5);
    clm_lc_p=nan(8,3,5);

    nds=1; % to store r and p values

    jli_u_scalar=2;
    data_calc(:,:,1,1)=squeeze(vh_clm_anom(:,:,m,1)); % high VI, JLI
    data_calc(:,:,2,1)=squeeze(vh_clm_anom(:,:,m,2))*jli_u_scalar; % high VI, u300
    data_calc(:,:,3,1)=squeeze(vl_clm_anom(:,:,m,1)); % low  VI, JLI
    data_calc(:,:,4,1)=squeeze(vl_clm_anom(:,:,m,2))*jli_u_scalar; % low  VI, u300
    
    [s1 s2 nv s4]=size(data_calc);
    
    ylim_rng=[-4.5 2.5];
    for clm=1:7
        for lc=1:8

            % spatial points
            amask=all_subgrp; amask(amask~=clm)=nan;amask(amask==clm)=1;
            bmask=lc_dom_grp; bmask(bmask~=lc)=nan;bmask(bmask==lc)=1;
            dtmp=data_calc.*repmat(amask,[1 1 nv nds]).*repmat(bmask,[1 1 nv nds]);

            for v=1:nv  % variable presented in the calculation

                % The mean values across space
                rtmp=squeeze(dtmp(:,:,v,1)); % anomalies

                nanidc=~isnan(rtmp) & (rtmp~=0);
                rtmp=rtmp(nanidc);
                clm_cl_rpt(clm,lc,v,1)=nanmean(rtmp);
                
                % number of gridpoints
                clm_cl_cnt(clm,lc,v)=length(rtmp(:));
                clm_cl_stderr(clm,lc,v)=std(rtmp(:),0,1,'omitnan');

                % One-sample t-test to check the significant difference
                % from zero
                [h_r,p_r,ci_r,stats_r] = ttest(rtmp);
                clm_cl_rpt(clm,lc,v,2)=h_r;
            end
        end
    end

    lc_names={'EF','DF','MF','WS', 'GRS','WL', 'CRP', 'SIB'};
	clmz_name={'Boreal, Scan.+Fin.','Boreal, East. Eur.', 'Temp., UK+France.','Temp., Central Eur.','Temp., East Eur.','Medit., Iberia','Central-East Medit.' };
    labels={'a)','b)','c)','d)','e)','f)','g)'};
    gsm_str={'GS(1)','GS(2)','GS(3)','GS(4)','GS(5)','GS(6)','GS(7)','GS(8)','GS(9)','GS(10)','GS(11)','GS(12)','EGS','LGS','Entire GS'};
    nv_show=1:4;
    figure('Position',[353  343  1097  546],'Color','w');

    gap_h=0.07; gap_w=0.08;
    gap=[gap_h gap_w]; marg_h=[0.08 0.1]; marg_w=[0.05 0.05];
    ha = tight_subplot(3,3,gap,marg_h,marg_w);
    
    clm_grp_grp={[1:2],[3:5],[6:7]};
    tile_num={[1 4],[2 5 8],[3 6]};
    % --- climate group ---
    for gi=1:3
        ti=0;
        for clm_grp=clm_grp_grp{gi}
            ti=ti+1;

            % sorted by JLI anomalies of high VI of land class
            [B sort_i(:,clm_grp)]=sort(clm_cl_rpt(clm_grp,:,1,1),'descend','MissingPlacement','last');
            xtmp=sort_i(:,clm_grp);

            % --- climatezone A ---
            axes(ha(tile_num{gi}(ti)));

            n=0;
            x=nan; % array to store land class index
            vi_sort=1; % based on JLI anomalies
            for i=1:length(xtmp)
              if ~isnan(clm_cl_rpt(clm_grp,xtmp(i),vi_sort,1)) && clm_cl_cnt(clm_grp,xtmp(i),vi_sort) >= 10 % skim based on number of valid points
                 n=n+1;
                 x(n)=xtmp(i);
                 clm_cl_x(clm_grp,n)=xtmp(i); % output for decomp_VegClm_relation_plot
              end
            end

            hold on
            width=0.8;
            xoff=0;
            yoff1=0.1;
            yoff2=2;
            yoff3=0.28;

            % Show relevant varaible in the barchat column
            b1=bar((1:length(x))-xoff, squeeze(clm_cl_rpt(clm_grp,x,[1 2 3 4],1)),width);
            b1(1).FaceColor = [0 85 0]/255;
            b1(2).FaceColor = [0 85 255]/255;
            b1(3).FaceColor = [170 0 127]/255;
            b1(4).FaceColor = [255 85 0]/255;

            sig_ind=squeeze(clm_cl_rpt(clm_grp,x,[1 2 3 4],2));
            sig_ind(sig_ind==0)=nan;

            % All bar in the group
            for bi=1:nv
                x_bar=[b1(1,bi).XEndPoints]';
                y_bar=[b1(1,bi).YEndPoints]';
                y_bar_sig=(~isnan(y_bar)).*yoff2.*sig_ind(:,bi);
                plot(x_bar, y_bar_sig, '*', 'color','k');

            end

            % The second and the third bar in the group
            x_bar_2_4=[b1(1,2).XEndPoints; b1(1,3).XEndPoints; b1(1,4).XEndPoints]';

            % Draw a line for 2nd and the 3rd bar if they are significantly
            colororder([0 0 0;0.5 0.5 0.5]);
            
            set(gca, 'XTick', 1:length(x));
            set(gca, 'XTickLabel',{lc_names{x}});
            xtickangle(0);

            plot([-100 100],[0 0],'-k');
            
            hold off
            box on;
            ylabel(sprintf('%s(-)',upper(clmidx_name)));
            
            sgtitle(sprintf('%s %s and wind anomalies by climate zone and land class\n',gsm_str{m},upper(clmidx_name)),'Fontsize',12, 'FontWeight','bold');
            title(clmz_name{clm_grp},'Fontsize',9, 'FontWeight','normal');

            if clm_grp==2 || clm_grp==5 || clm_grp==7 
                xlabel('Landclass');
            end
            xlim([0.2 6.8]);
            if clm_grp==1
                ylim(ylim_rng);
            elseif clm_grp==2
                ylim(ylim_rng);
            else
                ylim(ylim_rng);
            end

            yyaxis right;
            ylabel('ZSI (m/s)');
            ylim(ylim_rng/jli_u_scalar);
            
            if clm_grp==3
                legend(b1,{'High VI, JLI', 'High VI, ZSI', 'Low VI, JLI', 'Low VI, ZSI'});
                legend boxoff;
            end
            
            % labels for sub-plots
            yscal=0.92; xscal=0.03;
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{clm_grp},'FontSize',11,'FontName','Arial');

        end
    end
end