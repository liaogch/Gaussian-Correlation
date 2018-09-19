function [ Threshold,index,cor_var] = CollectorStrategy_DataCorrelation( Wsum,target_M,M,N,T,W,r,r_d )
%for the data collector



s = sum(W(target_M,:));

W(target_M,:) = W(target_M,:)*Wsum/s;
W(:,target_M)=W(target_M,:)';




[Threshold,index,cor_var] = Get_Threshold( M,N,T,W,r,r_d );





end

