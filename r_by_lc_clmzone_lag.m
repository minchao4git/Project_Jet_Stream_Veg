function r_by_lc_clmzone()

    % ===== R between Veg. and JLI grouped by landclass and climate zone ===== 
    global lc_dom_grp all_subgrp v_clminx_r v_clminx_p;
    global clm_cl_x;
    global clmidx_name;

    [s1 s2]=size(all_subgrp);

    % Calculation
    data_calc=nan(s1,s2,3);
    clm_cl_x=nan(7,8);

    % number of gridpoint
    clm_cl_cnt=nan(7,8,5);

    % group in bar chart
    clm_cl_rpt=nan(7,8,3,3);
    clm_cl_stderr=nan(8,3,5);
    clm_lc_p=nan(8,3,5);

    nds=2; % to store r and p values
    v_in=4; % 4 is EVI2
    nv=3;
    lags=1:6;
    nlag=length(lags);
    data_calc(:,:,1,1:nlag,1)=squeeze(v_clminx_r(:,:,13,lags,v_in)); % 1st half
    data_calc(:,:,2,1:nlag,1)=squeeze(v_clminx_r(:,:,14,lags,v_in)); % 2nd half
    data_calc(:,:,3,1:nlag,1)=squeeze(v_clminx_r(:,:,15,lags,v_in)); % 15 for entire GS
    
    data_calc(:,:,1,1:nlag,2)=squeeze(v_clminx_p(:,:,13,lags,v_in)); % 1st half
    data_calc(:,:,2,1:nlag,2)=squeeze(v_clminx_p(:,:,14,lags,v_in)); % 2nd half
    data_calc(:,:,3,1:nlag,2)=squeeze(v_clminx_p(:,:,15,lags,v_in)); % 15 for entire GS

    for clm=1:7
        for lc=1:8

            % spatial points
            amask=all_subgrp; amask(amask~=clm)=nan;amask(amask==clm)=1;
            bmask=lc_dom_grp; bmask(bmask~=lc)=nan;bmask(bmask==lc)=1;
            dtmp=data_calc.*repmat(amask,[1 1 nv nlag nds]).*repmat(bmask,[1 1 nv nlag nds]);

            rlags=nan(1,nlag);
            for lag=1:nlag
                rtmp=dtmp(:,:,3,lag,1); % get the maximum based on the entire GS values

                nanidc=~isnan(rtmp) & (rtmp~=0);
                rtmp=rtmp(nanidc);
                rlags(lag)=nanmean(rtmp);
            end
            [dummy maxi]=max(abs(rlags),[],'omitnan');
            
            % GS periods
            for v=1:3  % variable presented in the calculation

                % mean across data sources
                btmp=squeeze(nanmean(dtmp(:,:,v,maxi,1),4));  % ds=1, r values
                nanidc=~isnan(btmp) & (btmp~=0);
                btmp=btmp(nanidc);
                
                clm_cl_rpt(clm,lc,v,1)=nanmean(btmp); % instead of lag-0 value, store the maximum values across lags 
               
                % number of gridpoints
                clm_cl_cnt(clm,lc,v)=length(btmp);
                % Standard deviation across space
                clm_cl_stderr(clm,lc,v)=std(btmp,0,1,'omitnan');

                % Confidence interval across space
                clm_lc_p(clm,lc,v)= clm_cl_stderr(clm,lc,v); % Here I use 1 SD as error bar instead
                
                % One-sample t-test to check the significant difference
                % from zero
                [h_r,p_r,ci_r,stats_r] = ttest(btmp);
                clm_cl_rpt(clm,lc,v,2)=h_r;
                
                clm_cl_rpt(clm,lc,v,3)=maxi; % store the lag with maximum abs values
                
            end
            
            % Two-sample t-test to check EGS-LGS difference
            a=dtmp(:,:,1,1);a=a(~isnan(a)); % EGS
            b=dtmp(:,:,2,1);b=b(~isnan(b)); % LGS
            clm_cl_sig(clm,lc)=ttest2(a(:),b(:));
        end
    end

    % ===== plotting ===== 
    color_lc_grp=[[0 129 27];
                      [124 255 100]; 
                      [142 204 51];
                      [214 236 163]; 
                      [244 181 120]; 
                      [0 104 150]; 
                      [255 236 131]; 
                      [170 255 255]
        ]/255;

    lc_names={'EF','DF','MF','WS', 'GRS','WL', 'CRP', 'SIB'};
    clmz_name={'Boreal, Scan.+Fin.','Boreal, East. Eur.', 'Temp., UK+France.','Temp., Central Eur.','Temp., East Eur.','Medit., Iberia','Central-East Medit.' };
    labels={'a)','b)','c)','d)','e)','f)','g)'};
    nv_show=1:4;
    figure('Position',[461   377   984   546],'Color','w');

    gap_h=0.07; gap_w=0.05;
    gap=[gap_h gap_w]; marg_h=[0.08 0.1]; marg_w=[0.05 0.05];
    ha = tight_subplot(3,3,gap,marg_h,marg_w);

    clm_grp_grp={[1:2],[3:5],[6:7]};
    tile_num={[1 4],[2 5 8],[3 6]};
    % --- climate group ---
    for gi=1:3
        ti=0;
        for clm_grp=clm_grp_grp{gi}
            ti=ti+1;

            % sorted by annual values of land class
            [B sort_i(:,clm_grp)]=sort(clm_cl_rpt(clm_grp,:,3,1),'descend','MissingPlacement','last');
            xtmp=sort_i(:,clm_grp);

            % --- climatezone A ---
%             nexttile(tile_num{gi}(ti));
            axes(ha(tile_num{gi}(ti)));

            n=0;
            x=nan; % array to store land class index
            for i=1:length(xtmp)
              if ~isnan(clm_cl_rpt(clm_grp,xtmp(i),3,1)) && clm_cl_cnt(clm_grp,xtmp(i),3) >= 10 % skim based on number of valid points
                 n=n+1;
                 x(n)=xtmp(i);
                 clm_cl_x(clm_grp,n)=xtmp(i); % output for decomp_VegClm_relation_plot
              end
            end

            hold on
            width=0.8;
            xoff=0;
            yoff1=0.1;
            yoff2=0.35;
            yoff3=0.28;

            % Show relevant varaible in the barchat column
            b1=bar((1:length(x))-xoff, squeeze(clm_cl_rpt(clm_grp,x,[3 1 2],1)),width);
            b1(2).FaceColor = [22 170 116]/255;
            b1(3).FaceColor = [255 232 99]/255;

            sig_ind=squeeze(clm_cl_rpt(clm_grp,x,[3 1 2],2));
            sig_ind(sig_ind==0)=nan;

            % All bar in the group
            for bi=1:3
                x_bar=[b1(1,bi).XEndPoints]';
                y_bar=[b1(1,bi).YEndPoints]';
                y_bar_sig=(~isnan(y_bar)).*yoff2.*sig_ind(:,bi);
                plot(x_bar, y_bar_sig, '*', 'color','k');

            end

            % The second and the third bar in the group
            x_bar_2_3=[b1(1,2).XEndPoints; b1(1,3).XEndPoints]';

            % Draw a line for 2nd and the 3rd bar if they are significantly
            % different
            for i=1:1:length(x) % sequential i is the index in the selected lc array
                if clm_cl_sig(clm_grp,x(i))==1
                    plot(x_bar_2_3(i,:), [yoff2 yoff2]-0.05,'-r','LineWidth',2);
                    plot(x_bar_2_3(i,:), [yoff2 yoff2]-0.05,'.r','MarkerSize',15);
                end
                
                text(x_bar(i,1)-0.6,-0.23,sprintf('lag-%d',clm_cl_rpt(clm_grp,x(i),[3],3)-1)); % convert to the real lag value

            end

            set(gca, 'XTick', 1:length(x));
            set(gca, 'XTickLabel',{lc_names{x}});
            xtickangle(0);

            plot([-100 100],[0 0],'-k');

            hold off
            box on;
            ylabel('r');

            sgtitle(sprintf('GS-based Max. EVI~%s(r) by climate zone and land class\n\n',upper(clmidx_name)),'Fontsize',12, 'FontWeight','bold');
            title(clmz_name{clm_grp},'Fontsize',9, 'FontWeight','normal');

            if clm_grp==2 || clm_grp==5 || clm_grp==7 
                xlabel('Landclass');
            end
            xlim([0.2 6.8]);
            if clm_grp==1
                ylim([-0.3 0.4]);
            elseif clm_grp==2
                ylim([-0.3 0.4]);
            else
                ylim([-0.3 0.4]);
            end

            if clm_grp==3
                legend(b1,{'Entire GS', 'Early GS', 'Late GS'});
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