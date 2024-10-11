load(['./Results/source/sub-CC120120/source_rec_results.mat'],'source_roi_data','labels');

T=2;
fs=100;
segleng = T*fs; %2 sec
segshift =  T*fs; %1 sec
epleng = T*fs; %2 sec
df=1/(segleng/fs);
desired_f=50;
maxfreqbins= desired_f/df  ;
eps=[30,150,size(source_roi_data,3)];

for B=0:1
%0 is Bispectrum and 1 for Bicoherence
max_b=zeros(1,size(source_roi_data,1));
all_b={};

for e=1:3
    nep=eps(e);
    %1,5 and whole [30,150,size(source_roi_data,3)]minuntes
for i = 1:size(source_roi_data,1)
    data=source_roi_data(i,:,1:nep);
    data=data(:,:)';
    [cs,csnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreqbins);

    if B==1
         cs=cs./csnr;
    end

    cs=abs(cs);
    b=reshape(cs(1,:,:),size(cs,2),[]);

    all_b{e,i}=b;

    if max(b,[],'all')>max_b(i)
        max_b(i)=max(b,[],'all');
    end
end
end

for e=1:3
    nep=eps(e);
for i = 1:size(source_roi_data,1)
    fig=figure('Visible','off');
    imagesc(all_b{e,i});
    colormap hot
    colorbar;
    name=labels{i};
    name=[name ' ' num2str(nep*2) ' S'];
    title(name);
    xticks = linspace(1, maxfreqbins, 6);  % 6 evenly spaced ticks (0, 10, 20, 30, 40, 50 Hz)
    xticklabels = 0:10:desired_f;  % Corresponding labels from 0 to 50 Hz
    % Apply the xticks and corresponding labels
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    set(gca, 'YTick', xticks, 'YTickLabel', xticklabels);  % Same for y-axis if needed
    % Add labels and title (optional)
    xlabel('Frequency (Hz)');
    ylabel('Frequency (Hz)');

    if B==1
        clim([0 1])
        saveas(fig,['./roi_pac/Bicoh/' name '.png'])
    else 
        clim([0 max_b(i)])
        saveas(fig,['./roi_pac/' name '.png'])
    end
    close all
end
end
end
