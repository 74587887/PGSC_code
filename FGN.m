function [res,Tim]=FGN(X,s,alpha,beta,gamma)
% s is the true class label.
niter=20;  % num of iteration
tic;
mv=size(X,1);
n=size(X{1},1);
Z=eye(n);
Zv=repmat(Z,[1,1,mv]);
c=length(unique(s));
wv=ones(mv,1)/mv;
%options = optimset( 'Algorithm','interior-point-convex','Display','off');
% for i=1:mv
%     f=X{i};
%     Zv(:,:,i)=(f*f'+alpha*eye(n))\(f*f');
% end
for ii=1:niter
    disp(['now this is the ',num2str(ii),' th iter']);    
    % disp(['ii = ',num2str(ii)]);
    Z(Z<0)=0;
    Z= (Z+Z')/2;
    Zold=Z;
    D = diag(sum(Z));
    L = D-Z;     
    
    [F, temp, ev]=eig1(L, c, 0);    

    
    for i=1:mv
        f=X{i};
        Zv(:,:,i)=(f*f'+alpha*eye(n)+beta*wv(i)*eye(n))\(beta*wv(i)*Z+f*f');
        T=Zv(:,:,i);
        T(T<0)=0;
        % fprintf('now T is obtained, we are choosing value')
        [T,~]=choose_neighbor_coefficient(T,8,20); 
        T=(T+T')/2;
        Zv(:,:,i)=T;
        wv(i)=1/2/norm(Zv(:,:,i)-Z,'fro');
        M(:,:,i)=wv(i)*Zv(:,:,i);
        
    end
    
    parfor ij=1:n
        all=distance(F,n,ij);      
        Z(:,ij)=(sum(M(:,ij,:),3)-gamma*all'/(4*beta))/sum(wv);
    end

    actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 1, 'display', 'off');
    
    Tim=toc;
    if ii>5 && ((norm(Z-Zold,'fro')/norm(Zold,'fro'))<1e-3)
        break
    end
end 
res = EvaluationMetrics(actual_ids,s);




end

function [all]=distance(F,n,ij);
  for ji=1:n
            all(ji)=(norm(F(ij,:)-F(ji,:)))^2; 
  end
end
