% read sequencing trace file and determine the methylation level

function alignedSeq = bisulfitePCRSeqAnalysis(SCFfile, PCRseq, PCRstart, PCRend, BSFstart, BSFend, imgfile, methylMotif)

if nargin < 8
    methylMotif = {'GC'};
end

[sample, probability] = scfread(SCFfile);
% peak data background parameter
bkg.fraction = 0.3; % this fraction at the bottom are cosidered as background
bkg.scale = 2; % can scale the bck level by this number.

BSFseq = probability.base';

% display('Go to the MegAlign file and find the aligned region')
% PCRstart = input('The start of the aligned original seq (>10bp): ', 's');
% BSFstart = input('The start of the aligned bisulfite seq (>10bp): ', 's');
% PCRend = input('The end of the aligned original seq (>10bp): ', 's');
% BSFend = input('The end of the aligned bisulfite seq (>10bp): ', 's');

k1 = strfind(PCRseq,PCRstart);
k2 = strfind(PCRseq,PCRend)+length(PCRend)-1;
% disp(k1)
% disp(k2)
PCRalign = PCRseq(k1:k2);

k3 = strfind(BSFseq,BSFstart);
k4 = strfind(BSFseq,BSFend)+length(BSFend);

[k1, k2, k3, k4]

BSFalign = BSFseq(k3:k4);
peakPos = probability.peak_index(k3:k4);
peakA = sample.A(peakPos);
peakC = sample.C(peakPos);
peakG = sample.G(peakPos);
peakT = sample.T(peakPos);

% get rid of the background noise in the peak data
bkg_length = int32((k4-k3)*bkg.fraction);
tempA = sort(peakA); bkgA = mean(tempA(1:bkg_length))*bkg.scale;
tempC = sort(peakC); bkgC = mean(tempC(1:bkg_length))*bkg.scale;
tempG = sort(peakG); bkgG = mean(tempG(1:bkg_length))*bkg.scale;
tempT = sort(peakT); bkgT = mean(tempT(1:bkg_length))*bkg.scale;

peakA = max(peakA-bkgA,0); % if the number is less than 0, set it to 0.
peakC = max(peakC-bkgC,0);
peakG = max(peakG-bkgG,0);
peakT = max(peakT-bkgT,0);

% alignment data
group = {PCRalign,BSFalign};
% alignedSeq = showAlignmentMulti(group);

% find all the "C"s in the methylMotif
motifC = [];
for i=1:length(methylMotif)
    motif = methylMotif{i};  
    Cpos = strfind(motif,'C');
    motifC = [motifC, strfind(PCRalign,motif)+Cpos-1]; % find all the occurence of the motif in the original sequence    
end

motifC = unique(motifC);
totalPeak = double(peakA(motifC)+peakC(motifC)+peakG(motifC)+peakT(motifC));
probA_motifC = double(peakA(motifC))./totalPeak;
probC_motifC = double(peakC(motifC))./totalPeak;
probG_motifC = double(peakG(motifC))./totalPeak;
probT_motifC = double(peakT(motifC))./totalPeak;
[(motifC+k1)',probT_motifC]

disp('plotting');
fig=figure(); % 'visible', 'off');
subplot(2,1,1);
H=bar(motifC+k1,[probT_motifC,probA_motifC,probG_motifC,probC_motifC],'stacked');
ch = get(gca, 'children');
set(ch(1), 'facecolor', 'r');
set(ch(2), 'facecolor', 'y');
set(ch(3), 'facecolor', 'g');
set(ch(4), 'facecolor', 'b');

AX=legend(H, {'T', 'A', 'G', 'C'}, 'Location', 'Best', 'FontSize', 8);
% P=findobj(gca,'type','patch');
% myC= [0 0 1
%       0 1 0
%       1 1 0
%       1 0 0 ];
% for n= 1 : length(P) 
%     set(P(n),'FaceColor',myC(n,:));
%     set(P(n),'EdgeColor',myC(n,:));
% end
%[double(probA(motifC)),double(probC(motifC)),double(probG(motifC)),double(probT(motifC))]

% find all the other "C"s (not in the methylMotif)
temp = strfind(PCRalign,'C');
non_motifC = setdiff(temp,motifC); 

totalPeak_non = double(peakA(non_motifC)+peakC(non_motifC)+peakG(non_motifC)+peakT(non_motifC));
probA_non_motifC = double(peakA(non_motifC))./totalPeak_non;
probC_non_motifC = double(peakC(non_motifC))./totalPeak_non;
probG_non_motifC = double(peakG(non_motifC))./totalPeak_non;
probT_non_motifC = double(peakT(non_motifC))./totalPeak_non;

subplot(2,1,2);
H=bar(non_motifC+k1,[probT_non_motifC,probA_non_motifC,probG_non_motifC,probC_non_motifC],'stacked');
ch = get(gca, 'children');
set(ch(1), 'facecolor', 'r');
set(ch(2), 'facecolor', 'y');
set(ch(3), 'facecolor', 'g');
set(ch(4), 'facecolor', 'b');
AX=legend(H, {'T', 'A', 'G', 'C'}, 'Location', 'Best', 'FontSize', 8);

% P=findobj(gca,'type','patch');
% for n= 1 : length(P) 
%     set(P(n),'FaceColor',myC(n,:));
%     set(P(n),'EdgeColor',myC(n,:));
% end
saveas(fig, imgfile)

% [(non_motifC+k1)',probT_non_motifC]




