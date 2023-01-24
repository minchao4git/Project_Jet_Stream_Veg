% get months for 1st and 2nd half season
function [ssn1 ssn2 ssn_end gs_mid]=get_ssn12(m_lgs)

    if m_lgs ==12 % less likely
       sosmi=1;
    else
       sosmi=2; % normal case: sos month index in array
    end
    gs_mid=floor(m_lgs/2);

    % dynamical adjust the 1st and 2nd half of the growing
    % season
    if gs_mid > 0
        ssn1=[1:gs_mid]+(sosmi-1);
        ssn2=(max(ssn1)+1):(m_lgs+(sosmi-1));
        ssn_end=m_lgs+(sosmi-1);
    else
        ssn1=sosmi;
        ssn2=sosmi;
        ssn_end=sosmi;
        gs_mid=sosmi;
    end
end