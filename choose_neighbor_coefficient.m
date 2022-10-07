function [S_weight,S_number]=choose_neighbor_coefficient(Z,num_good,num_ordinary)

% num_good: 8 here, the number of good neighbors
% num_ordinary: 20 here, the number of ordinary neighbors
% the optimal parameter set varies with different datasets


Z = Z - diag(diag(Z));

z=sum(Z');
w = zeros(size(Z,1));
for i=1:size(Z,1)
        w(i,:)=Z(i,:)./z(i);
end
C=zeros(size(w,1));
S_number_temp = zeros(size(w,1),num_ordinary);
S_weight_temp = zeros(size(w,1),num_ordinary);
for i=1:size(w,1)
    C(i,:)=w(i,:)+w(:,i)';
    
    for j=1:num_ordinary 
        [p,q]=max(C(i,:));
        S_number_temp(i,j)=q;
        S_weight_temp(i,j)=p;
        C(i,q)=0; 
    end
    
end

        
C=zeros(size(w,1));
num=0;
S_weight = zeros(size(w,1));
S_number = zeros(size(w,1),num_good);
for i=1:size(w,1)
    
    C(i,:)=w(i,:)+w(:,i)';
    
    for j=1:num_good 
        [p,q]=max(C(i,:));
        
        
        for k=1:size(w,1)
           if  find(S_number_temp(q,:)==i)
               S_weight(i,q)=p;
               S_number(i,j)=q;
               C(i,q)=0; 
               break;
           else
               num=num+1;
               C(i,q)=0; 
               [p,q]=max(C(i,:));
           end
        end
    end
end

for k=1:size(w,1)
   idx=find(S_number(k,:)==0);
   S_number(k,idx)=k;
   S_weight(k,idx)=1;
end
