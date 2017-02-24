function x = trimbits(y,bit)



xi = typecast(y,'uint32');
%xb = bitget(xi,1:32)==true;

%xi = bitset(xi,1,rand>.5);
if bit>1
for b=1:bit-1
    %xi = bitset(xi,b,rand>.5);
    xi = bitset(xi,b,false);
end
end
xi = bitset(xi,bit,false);
%xi = bitset(xi,bit,rand>.5);
x = typecast(xi,'single');



% nup=0;
% for i=1:N
% 
%     %convert to sign-mag
%     xb = bitget(Xi(i,:),1:32)==true;
%     
%     
%     tmp = uint8(0);
%     for b=1:8
%         tmp = bitset(tmp,b,xb(23+b));
%     end
%         
%     if tmp<127
%         xb(31)=true;
%     else
%         xb(31)=false;
%     end
%     tmp = uint8(abs(127-int16(tmp)));
%     
%     
%     xb(24:30) = bitget(tmp,1:7)==true;
%     
%     
%     Xb(i,:) = xb;
%     
%     if i/N>nup
%         nup=nup+.01
%     end
%     
% end

