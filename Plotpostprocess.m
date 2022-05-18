clear Xpost Ypost Cpost Zpost
close all
clc
q=1;
k=1;%initialize count
figure1=figure(1);
for x=startX:scanResolution:endX
    for y=startY:scanResolution:endY
        kstart=k;
        for j=startData:endData
            k=k+1;
        end
        kend=k;
        [pks,locs]=findpeaks(envelope(C(kstart:min(kend,length(C)))),'MinPeakProminence',400); % tune last value for each new scan
        for l=1:length(pks)
            if locs(l)>0 && locs(l)<100000 %use this to isolate different depths
                Xpost(q)=x;
                Ypost(q)=y;
                Cpost(q)=pks(l);
                Zpost(q)=locs(l);
                q=q+1;
            end
        end
    end
    x
end
figure
scatter3(Xpost,Ypost,-Zpost*0.34301*0.001*timeIntervalNanoseconds*10,[],Cpost,'.')%includes calculation based on speed of sound through water
title('Depth Plot')
xlabel('mm')
ylabel('mm')
zlabel('.1 mm')
axis equal
a = colorbar;
a.Label.String = 'mV';