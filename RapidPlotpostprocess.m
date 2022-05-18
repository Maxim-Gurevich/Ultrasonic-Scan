clear Xpost Ypost Cpost Zpost
close all
clc
q=1;
k=1;%initialize count
for x=startX:scanResolutionX:endX
    kstart=k;
    for w=1:width(chA)
        for j=1:length(data)
            k=k+1;
        end
        kend=k;
        percentDone=k*100/(width(chA)*length(data)*((endX-startX)/scanResolutionX))
        [pks,locs]=findpeaks(envelope(C(kstart:min(kend,length(C))),rmsWindow,'rms'),'MinPeakProminence',minPeakProminence);
        for l=1:length(pks)
            if locs(l)>0 && locs(l)<100000
                Xpost(q)=x;
                Ypost(q)=startY-(startY-endY)/width(chA)*w;
                Cpost(q)=pks(l);
                Zpost(q)=locs(l);
                q=q+1;
            end
        end
    end
end
figure()
scatter3(Xpost,Ypost,-Zpost*0.34301*0.001*timeIntervalNanoseconds*10,[],Cpost,'.')
title('Depth Plot')
xlabel('mm')
ylabel('mm')
zlabel('.1 mm')
%axis equal
a = colorbar;
a.Label.String = 'mV'