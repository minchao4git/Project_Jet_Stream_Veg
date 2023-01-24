function JLI_CLIM_Veg_SI(m)
    global v_clminx_r v_clminx_p;
    global lons lats ax_show;
    %% === For local climatic drivers ===
    y_pos=0.40;
    x_pos=-0.05;
    c_thr=0.600001;

    figure('color','w','Position',[935   132   575   890]);
    gap=[0.02 0.02]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(4,3,gap,marg_h,marg_w);
    yscal=0.96; xscal=0.015;
    cb_xR=0.55;
    cb_yR=0.4;
    cb_wR=3.5;
    ax_show.frame=0;
    ncb=25;
    vnames={'Temp.', 'Rad.', 'SM', 'EVI2'};
    labels={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
    [s1 s2 s3 s4 s5]=size(v_clminx_p);
    % Significant points
    lons_sig=repmat(lons,[1 1 s3 s4 s5]);
    lats_sig=repmat(lats,[1 1 s3 s4 s5]);
    lons_sig(isnan(v_clminx_p))=nan;
    lats_sig(isnan(v_clminx_p))=nan;
    lons_sig(v_clminx_p>=0.05)=nan;
    lats_sig(v_clminx_p>=0.05)=nan;

    for lag=1:3
        for v=1:4 % % Temp., Radiation, Soil moisture, EVI2
            axes(ha((v-1)*3+lag)); 
            bg_show=squeeze(flipud(rot90(squeeze(v_clminx_r(:,:,m,lag,v)))));
            bg_show(isnan(bg_show))=0;
            geoplot(ax_show,bg_show);
            crng=[-c_thr c_thr];
            colormap(gca, flipud(cbrewer('div','RdBu',ncb)));
            caxis([crng(1) crng(2)]);

            % Mark significant points
            if v==4
                intv=1;
            else
                intv=2;
            end
            a=squeeze(lats_sig(1:intv:end,1:intv:end,m,lag,v));
            b=squeeze(lons_sig(1:intv:end,1:intv:end,m,lag,v));
%             plotm(b(:), a(:), 'Color',[80 80 80]/255,'LineStyle','none', 'Marker', '+', 'MarkerSize',3);
            plotm(b(:), a(:), 'Color',[50 50 50]/255,'LineStyle','none', 'Marker', '+', 'MarkerSize',4, 'LineWidth',0.9);

            if v==1
                title(sprintf('lag-%d',lag-1));
            end

            if lag==1
                a=get(gca);
                gca_w = (a.XLim(2)-a.XLim(1));
                gca_h = (a.YLim(2)-a.YLim(1));
                text(a.XLim(1)+gca_w*x_pos, a.YLim(1)+gca_h*y_pos,vnames{v},'FontSize',11,'FontName','Arial','Rotation',90);
            end

            if lag==2 && v==4
                cb3=colorbar('Southoutside');
                resizeCB(cb3, cb_xR, cb_yR, cb_wR, 0.6, '',19,10,8);
            end

            % labels for sub-plots
            yscal=0.92; xscal=0.03;
            a=get(gca);
            gca_w = (a.XLim(2)-a.XLim(1));
            gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{(v-1)*3+lag},'FontSize',11,'FontName','Arial');

        end
    end
end