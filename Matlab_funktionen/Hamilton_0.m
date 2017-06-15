function H = Hamilton_0(J,E )
    H=diag(E);
    H(1,length(H))=J;
    H(length(H),1)=J;
    for i=1:length(H)-1
        H(i+1,i)=J;
        H(i,i+1)=J;
    end

end

%
% function H=Hamilton(J,E,Zustande)
% [n,m]=size(Zustande);
%   H=zeros(m,m);
%         for i=1:m
%             Zustand=Zustande(:,i);
%             Ergebnisse=zeros(n,n*2);
%             for j=1:(n-1)
%               Ergebnisse(:,j)=CTC(j+1,j,Zustand);
%             end
%             Ergebnisse(:,n)=CTC(1,n,Zustand);
%             for j=(n+1):(2*n-1)
%               Ergebnisse(:,j)=CTC(j-n,j+1-n,Zustand);
%             end
%             Ergebnisse(:,n*2)=CTC(n,1,Zustand);
%             for k=1:n*2
%               for l=1:m
%               if (Ergebnisse(:,k)==Zustand);
%                 H(i,l)+=0;
%               else
%                 if (Ergebnisse(:,k)==Zustande(:,l))
%                     H(i,l)+=J;
%                 end
%               end
%               end
%             end
%         end
% end
