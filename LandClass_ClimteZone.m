function LandClass_ClimteZone()
    global DATA_LC_out DATA_KG_CLMZ_05rs_out;
    global lc_dom_grp clmzone_grp all_subgrp;
    global ax_show;
    %% Climate zone and landclass grouping
    % ====> IGBP Land cover
    % > The calculation is based on the regridded dataset from MODIS IGBP land
    % cover type
    % > Two appraoches are used to achieve the dominant landcover type, and
    % both yield similar results
    % > The landcover type based on 0.5 grid and 0.05 grid are quite similar,
    % and also consistent with landclass from IGBP: 
    % https://e4ftl01.cr.usgs.gov/MOTA/MCD12C1.006/2008.01.01/BROWSE.MCD12C1.A2008001.006.2018053184623.1.jpg

    % ----------------------------
    %  Eight groups of land class 
    % ----------------------------
    % Skip igbp 1: water boday

    % Method I: Freqency approach: get the stable land cover
    [s1 s2 s3 s4]=size(DATA_LC_out);
    nyr_lc=s3;
    for i=1:s1
        for j=1:s2
            for y=1:nyr_lc
                [dump lc_dom_ts(i,j,y)]=max(DATA_LC_out(i,j,y,:));
            end
        end
    end

    lc_dom=mode(lc_dom_ts,3); % get the most freqent value
    lc_dom_stable=nan(s1,s2);
    for i=1:s1
        for j=1:s2

            if sum((find(lc_dom_ts(i,j,:)==lc_dom(i,j))>0))==nyr_lc
                lc_dom_stable(i,j)=lc_dom(i,j);
            end
        end
    end
    % 
    % % Method II: Mean approach
    % for i=1:s1
    %     for j=1:s2
    %         [dump lc_dom(i,j)]=max(nanmean(DATA_LC_out(i,j,:,:),3));
    %     end
    % end

    % grouping
    lc_dom_grp=nan(s1,s2);
    for i=1:s1
        for j=1:s2
            switch lc_dom(i,j)
                % evergreen
                case {2 3}
                    lc_grp=1;
                % decidous
                case {4 5}
                    lc_grp=2;
                % mixed forest
                case {6}
                    lc_grp=3;
                % closed shrubland, open shrubland and woody savannas
                case {7,8,9}
                    lc_grp=4;
                % grasslands, savannas
                case {10, 11}
                    lc_grp=5;
                % permanent wetland
                case {12}
                    lc_grp=6;
                % crop land, natural vegetation mosaic (we skip 14: Urban here)
                case {13,15}
                    lc_grp=7;
                % snow and ice, barren sparsely veg.
                case {16,17}
                    lc_grp=8;
                otherwise
                    lc_grp=nan;
            end

            lc_dom_grp(i,j)=lc_grp;
        end
    end

    color_lc_grp=[[0 129 27]; 
                  [124 255 100]; 
                  [142 204 51]; 
                  [214 236 163]; 
                  [244 181 120]; 
                  [0 104 150]; 
                  [255 236 131]; 
                  [170 255 255]
    ]/255;
    %%
    % -----------------------------------
    %     Three groups of climate zones 
    % -----------------------------------
    % - also devided into several sub-groups for each climate zones.
    % Boreal: 110: Scandinavia + Finland, 120 Eastern Europe, 130 The Alps
    % Temperate: 210: UK and France, 220, Continental Europe, 230, Eastern Europe
    % Mediterranean: 310: Iberia, 320: Central-Eastern Med
    
     % grouping
    clmzone_grp=nan(s1,s2);
    for i=1:s1
        for j=1:s2
            switch DATA_KG_CLMZ_05rs_out(i,j)
                % Dfb, Dfc, ET (snow, fully humid, warm and cool summer, polar tundra)
                case {42, 43, 62}
                    cl_grp=1;
                % Cfa, Cfb, Cfc (warm temperate, fully humid, hot-warm-cool summer)
                case {31, 32, 33}
                    cl_grp=2;
                % Csa, Csb, Csc (warm temperate, summer dry, hot-warm-cool summer)
                % BSk (arid steppe cold arid)
                case {34,   35,  36, 26}
                    cl_grp=3;
                otherwise
                    cl_grp=nan;
            end

            clmzone_grp(i,j)=cl_grp;
        end
    end
    % You can skip the group with number of grid point < 10
    color_clm_grp=[[226 68 211]
                  [0 121 24];
                  [255 254 118];
    ]/255;

%% Sub climate group
   
    clmzone_subgrp=nan(s1,s2);
    % === Boreal ===
	a=clmzone_grp; a(a~=1)=nan;
    b1=nan(s1,s2);b2=nan(s1,s2);
    
    % For subgroup 110: Scandinavia
	tmp=nan(s1,s2);
    tmp(35:66,47:77)=1;
    tmp(67:82,52:77)=1;
    tmp(82:92,54:77)=1;
    b1((tmp+a)>0)=1;
    
    % For subgroup 120: Eastern Europe
	tmp=nan(s1,s2);
    tmp(59:92,15:77)=1;
    b1_tmp=b1; b1_tmp(isnan(b1))=0; b1_tmp(~isnan(b1))=nan;
    b2=tmp + a  + b1_tmp;
    b2(~isnan(b2))= 2;

    % The rest is the Alps
    b1_tmp=b1; b1_tmp(isnan(b1))=0; b1_tmp(~isnan(b1))=nan;
    b2_tmp=b2; b2_tmp(isnan(b2))=0; b2_tmp(~isnan(b2))=nan;
    b12=b1_tmp + b2_tmp;
    b3=a+b12;
    b3(~isnan(b3))= nan;  % Leave out the Alps
    
    % === Temperate ===
	a=clmzone_grp; a(a~=2)=nan;
    
    t1=nan(s1,s2);t11=nan(s1,s2);t12=nan(s1,s2);t2=nan(s1,s2);t3=nan(s1,s2);
    % For subgroup 210: UK, France and Coastal Norway
	tmp=nan(s1,s2);
    tmp(33:53,48:66)=1;
    t11((tmp+a)>0)=3;
    
	tmp=nan(s1,s2);
    tmp(3:39,10:52)=1;
    t1_tmp=t11; t1_tmp(isnan(t11))=0; t1_tmp(~isnan(t11))=nan;
    t12=tmp + a  + t1_tmp;
    t12(~isnan(t12))= 3;
    
    t11_tmp=t11; t11_tmp(isnan(t11))=0; t11_tmp(~isnan(t11))=nan;
    t12_tmp=t12; t12_tmp(isnan(t12))=0; t12_tmp(~isnan(t12))=nan;
    t1=t12_tmp;
    t1(isnan(t1))= 3;
    % also include the UK islands
    t1(23:24, 53:54)=3;
    t1(13, 56)=3;
    t1(12:13, 57)=3;
    t1(38:39, 49)=3;
    t1(50:53,48:49)=nan; % Parts of Sweden should belong to central temperate, see it nan to reserve for the later group
	t1(t1==0)= nan;
    
    % For subgroup 220: Continental Europe
	tmp=nan(s1,s2);
    tmp(40:65,10:52)=1;
    
    t1_tmp=t1; t1_tmp(isnan(t1))=0; t1_tmp(~isnan(t1))=nan;
    t2=tmp + a  + t1_tmp;
    t2(~isnan(t2))= 4;
    
    % The rest is 230: Eastern Europe
    t1_tmp=t1; t1_tmp(isnan(t1))=0; t1_tmp(~isnan(t1))=nan;
    t2_tmp=t2; t2_tmp(isnan(t2))=0; t2_tmp(~isnan(t2))=nan;
    t12=t1_tmp + t2_tmp;
    t3=a+t12;
    t3(~isnan(t3))= 5;
    for i=1:s1
        for j=1:s2
            if t3(i,j)>0 && j>=50 % Leave out the Norwegian coast
                t3(i,j)=nan;
            end
        end
    end
    
    % === Mediterranian ===
	a=clmzone_grp; a(a~=3)=nan;
    m1=nan(s1,s2);m2=nan(s1,s2);
    
    % For subgroup 310: Iberia
	tmp=nan(s1,s2);
    tmp(5:30,4:22)=1;
    m1((tmp+a)>0)=6;
    m1(26:34,3:6)=nan;     % remove north Africa
    
    % For subgroup 320: 
	tmp=nan(s1,s2);
    tmp(31:92,1:27)=1;
    m2((tmp+a)>0)=7;
    m2(30:49,1:7)=nan;     % remove north Africa

    % combine all
    b1_tmp=b1; b1_tmp(isnan(b1))=0;
    b2_tmp=b2; b2_tmp(isnan(b2))=0;
    b3_tmp=b3; b3_tmp(isnan(b3))=0;
    
    t1_tmp=t1; t1_tmp(isnan(t1))=0;
    t2_tmp=t2; t2_tmp(isnan(t2))=0;
    t3_tmp=t3; t3_tmp(isnan(t3))=0;
    
    m1_tmp=m1; m1_tmp(isnan(m1))=0;
    m2_tmp=m2; m2_tmp(isnan(m2))=0;
    a=clmzone_grp; a(a>0)=0;
    
    all_subgrp= b1_tmp + b2_tmp + b3_tmp + t1_tmp + t2_tmp + t3_tmp + m1_tmp + m2_tmp;
    %%
    if 1==2
        figure('color','white','Position',[555   451   613   422]);
        gap_h=0.005; gap_w=0.005;
        gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
        ha = tight_subplot(1,2,gap,marg_h,marg_w);
        yscal=0.92; xscal=0.03;
        cb_xR=1.0;
        cb_yR=1.6;
        cb_wR=1.0;
        ax_show.frame=0;

        % --- Land class ---
        axes(ha(1));
    %     imagescn(rot90(lc_dom_grp));
        geoplot(ax_show, lc_dom_grp');

        colormap(gca,color_lc_grp);
        caxis([0.5 8.5]);
        cbh1=colorbar('southoutside');
        cbh1.Ticks = [1:8];
        cbh1.TickLabels={'EF','DF','MF','WS', 'GRS','WL', 'CRP', 'SIB'};
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
        title('Landclass','FontSize',11);
        resizeCB(cbh1, cb_xR, cb_yR, cb_wR, 0.6, '',0.32,30,9);

        % --- Climate zone ---
        axes(ha(2));
        bg_show=clmzone_grp;
        bg_show(isnan(bg_show))=3; % nan in North Africa, set it as Medit.

        geoplot(ax_show, bg_show');

        colormap(gca,color_clm_grp);
        caxis([0.5 3.5]);
        cbh1=colorbar('southoutside');
        cbh1.Ticks = [1:3];
        cbh1.TickLabels={'Boreal','Temperate', 'Medit.'};
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'b)','FontSize',11,'FontName','Arial');
        title('Climate zone','FontSize',11);
        resizeCB(cbh1, cb_xR, cb_yR, cb_wR, 0.6, '',0.32,30,9);

    %     export_fig 'landclass_and_climatezone.png' -png -r600;
    end
    
    % Plot climate subgroup
    color_clm_grp=[[226 68 211];
                   [170 0 255];
                  [0 121 24];
                  [0 85 127];
                  [0 85 255];
                  [255 170 0];
                  [255 254 118];
                  [180 180 180];
    ]/255;
    
    figure('color','white','Position',[922   558   712   346]);
    gap_h=0.005; gap_w=0.005;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(1,1,gap,marg_h,marg_w);
    yscal=0.92; xscal=0.03;
    cb_xR=1.0;
    cb_yR=1.6;
    cb_wR=1.0;
    ax_show.frame=0;
    
  % --- Climate zone ---
    axes(ha(1));
    bg_show=all_subgrp;
    bg_show(bg_show==0)=9; % nan in North Africa, set it as Medit., but not in our analysis

    geoplot(ax_show, bg_show');

    colormap(gca,color_clm_grp);
    caxis([0.5 8.5]);
    cbh1=colorbar('eastoutside');
    cbh1.Ticks = [1:8];
    cbh1.TickLabels={'Boreal, Scandinavia+Finland',...
                     'Boreal, Eastern Europe',...
                     'Temperate, UK and France', ...
                     'Temperate, Central Europe', ...
                     'Temperate, Eastern Europe', ...
                     'Medit., Iberia', ...
                     'Medit., Central-Eastern Medit.', ...
                     'Not for analysis', ...
                     };
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
    title('Climate zone','FontSize',11);
    resizeCB(cbh1, cb_xR, cb_yR, cb_wR, 0.6, '',0.32,30,9);
end
