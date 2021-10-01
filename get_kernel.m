function [kernel,ts] = get_kernel(Yresp,Amplitudes,Yphant)

ts = 0:.006:2;
t_sq=.024:.006:2; %basis function width


for tj=1:length(t_sq)-1
    
    sq_indx = find(ts >= t_sq(tj),1);
    sq = zeros(1,length(ts));
    if sq_indx+3 > numel(ts)
        sq(sq_indx-4:numel(ts)) = 1;
    else
        sq(sq_indx-4:sq_indx+3) = 1;
    end
    sq_funs(:,tj)=sq';
end



D_v=[];
for i=1:length(Amplitudes)
    v_trl = sq_funs*Amplitudes(i);
    D_v = [D_v;v_trl];
end

b1=regress(Yresp,D_v);
rsp1= sq_funs*b1;

b2=regress(Yphant,D_v);
rsp2= sq_funs*b2;

rsp1 = medfilt3(rsp1);
rsp2=medfilt3(rsp2);


kernel = rsp1- rsp2; 
end