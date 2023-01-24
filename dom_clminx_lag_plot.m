function dom_clminx_lag_plot(m_show3, mlabels)
%%
    global v_clminx_r v_clminx_p;
	global lons lats ax_show;
    global clmidx_name;

    % Pick the climate indices with highest correlation coefficient
    % 6 lags x 10 climate indices
    [s1 s2 nmn nlag nv]=size(v_clminx_r);

    % ===== Pre-calculation =====
    clearvars dtmp d_mx d_mxi dtmp_rs;
    
    % Mask out the insignificant point
    clmi_chosen=[4];    % EVI2
    nclmi=length(clmi_chosen);
    dtmp=v_clminx_r(:,:,:,:,clmi_chosen);
    dtmp_p=v_clminx_p(:,:,:,:,clmi_chosen);
    % significant points
%     dtmp(v_clminx_p(:,:,:,:,clmi_chosen)>0.05)=nan;

    % Get the position of the maximum abs(r) values
    dtmp_rs = reshape(dtmp,[s1 s2 nmn nlag*nclmi]);
    dtmp_p_rs = reshape(dtmp_p,[s1 s2 nmn nlag*nclmi]);
    [d_mx d_mxi]=nanmax(abs(dtmp_rs),[],4);
    
    % Note: matlab return 1 for the indices for the all-nan array, we needs to mask
    % them out
    d_mxi(isnan(d_mx))=0;
    
    d_mxi=squeeze(d_mxi);
    map_mxclm_r=nan(s1,s2,nmn);
    map_mxclm_p=nan(s1,s2,nmn);
    for i=1:s1
        for j=1:s2
            for m=1:nmn
                mxi=d_mxi(i,j,m);
                if mxi > 0
                    map_mxclm_r(i,j,m)=dtmp_rs(i,j,m,mxi); % original r values based on the position
                    map_mxclm_p(i,j,m)=dtmp_p_rs(i,j,m,mxi); % original p values based on the position
                end
            end
        end
    end

    mycolormap=[
        255 255 255 ;
        173	255	47	;   % POL
        202	255	112	;
        188	238	104	;
        162	205	90	;
        110	139	61	;
        85	107	47	;
        ];
    
    % Significant points
    lons_sig=repmat(lons,[1 1 nmn]);
    lats_sig=repmat(lats,[1 1 nmn]);
    lons_sig(isnan(map_mxclm_p))=nan;
    lats_sig(isnan(map_mxclm_p))=nan;
    lons_sig(map_mxclm_p>=0.05)=nan;
    lats_sig(map_mxclm_p>=0.05)=nan;
%%

	% ===== Dominant map =====
%     labels1={'a)','b)','c)'};
%     labels2={'d)','e)','f)'};
    labels1={'a)','b)'};
    labels2={'c)','d)'};

    figure('color','w', 'Position',[ 773   242   571   570]);
    yscal=0.9; xscal=0.05;
    gap_h=0.02; gap_w=0.005;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.15];
    ha = tight_subplot(2,length(m_show3),gap,marg_h,marg_w);
    ax_show.frame=0;

    x_off=0.1;
    y_off=-0.05;
    h_ratio=1.3;
    
    n=0;
    for m=m_show3 % IMPARTANT TO START FROM 2 here! because sosmi=2!!
        n=n+1;
        axes(ha(n));

        bg_show=squeeze(map_mxclm_r(:,:,m))';
        bg_show(isnan(bg_show))=0;
        geoplot(ax_show,bg_show);
            
        colormap(gca, flipud(cbrewer('div','RdBu',21)));
        caxis([-0.8 0.80001]);
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
        title(mlabels{n});
        
        % Mark significant points
        a=squeeze(lats_sig(1:2:end,1:2:end,m));
        b=squeeze(lons_sig(1:2:end,1:2:end,m));
        plotm(b(:), a(:), 'Color',[50 50 50]/255,'LineStyle','none', 'Marker', '+', 'MarkerSize',4, 'LineWidth',1);
        
        % Label
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels1{n},'FontSize',11,'FontName','Arial');
        
        if n==length(m_show3)
            cbh2 = colorbar('eastoutside');
            cb_x=cbh2.Position(1); cb_y=cbh2.Position(2);
            cb_w=cbh2.Position(3); cb_h=cbh2.Position(4);
            set(cbh2, 'AxisLocationMode','manual');
            set(cbh2, 'Position',[cb_x+x_off cb_y+y_off cb_w*0.8 cb_h*h_ratio]);
        end
        
        if n==1
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)-gca_w*0.1, a.YLim(1)+gca_h*0.25,sprintf('Max. r (VI~%s)',upper(clmidx_name)),'FontSize',12,'FontName','Arial','Rotation',90);
        end
    end
    
    n=0;
    for m=m_show3 % IMPARTANT TO START FROM 2 here! because sosmi=2!!
        n=n+1;
        axes(ha(length(m_show3)+n));
        
        bg_show=squeeze((squeeze(d_mxi(:,:,m))'+0.5));
        bg_show(isnan(bg_show))=0;
        geoplot(ax_show,bg_show);
        
        colormap(gca, mycolormap(1:(nclmi*6)+1,:)/255);
        caxis([0 nclmi*6+1]); % IMPORTANT TO INCLUDE 0 HERE, which corresponds to the white color
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
        
        % Label
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels2{n},'FontSize',11,'FontName','Arial');
        
        if n==length(m_show3)
            cbh1 = colorbar('eastoutside');
            cbh1.Ticks = [0.5 1.5:1:6.5];
            cbh1.TickLabels={'no data','lag-0','lag-1','lag-2','lag-3','lag-4','lag-5'};
            cb_x=cbh1.Position(1); cb_y=cbh1.Position(2);
            cb_w=cbh1.Position(3); cb_h=cbh1.Position(4);
            set(cbh1, 'AxisLocationMode','manual');
            set(cbh1, 'Position',[cb_x+x_off cb_y+y_off cb_w*0.8 cb_h*h_ratio]);
        end
        
        if n==1
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)-gca_w*0.1, a.YLim(1)+gca_h*0.15,'Dominant lag','FontSize',12,'FontName','Arial','Rotation',90);
        end
    end
    
%     export_fig 'dom_clminx_lag_plot_EGS_LGS_EntireGS.png' -png -r600;
%     export_fig 'dom_clminx_lag_plot_SOS_plusmonths.png' -png -r600;
    export_fig(sprintf('dom_clminx_lag_plot_%s.png', clmidx_name), '-png', '-r600');

end

