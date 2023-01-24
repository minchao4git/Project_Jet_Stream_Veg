function [data_out]=data_dsdt(data_in)

        % Dimension in data_in: lon/lat lon/lat mm yy nv
        [s1 s2 s3 ny nv]=size(data_in);
        % dtmp_movm dimension: lon/lat lon/lat yy ds
        dtmp_movm=movmean(squeeze(nanmean(data_in,3)),5, 3); % 5-year running mean of annual mean

        % Detrend the monthly data with the detrended annual mean
        atmp=repmat(dtmp_movm,[1 1 1 1 12]);
        data_in=data_in-permute(atmp,[1 2 5 3 4]);
        
        % Deseasonalize with the climatology seasonal cycle
        for v = 1:nv
            % manually deseason by substracting climatology seasonal cycle
            % lon/lat lon/lat m
            climssn=squeeze(nanmean(data_in(:,:,:,:,v),4));
            data_in(:,:,:,:,v)=data_in(:,:,:,:,v)-repmat(climssn,[1 1 1 ny]);
        end
        data_out=data_in;
end