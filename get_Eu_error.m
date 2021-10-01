function out = get_Eu_error(r_targ,r_resp,theta_targ,theta_resp)

for i=1:numel(r_targ)
    
    x_targ = r_targ{i}.*sind(theta_targ{i}) ;
    y_targ = r_targ{i}.*cosd(theta_targ{i}) ;
    
    x_resp = r_resp{i}.*sind(theta_resp{i}) ;
    y_resp = r_resp{i}.*cosd(theta_resp{i}) ;
    
    out{i} = sqrt((x_targ-x_resp).^2 + (y_targ-y_resp).^2);
   
end

end