
function seqFinal = showAlignmentMulti(seq)

%lineWidth = 100;

for i=1:length(seq)-1
    seq1 = seq{i};
    seq2 = seq{i+1};
    for j=1:length(seq1)
       if strcmp(seq1(j),seq2(j))
           seqAlig(i,j)='|';
       else
           seqAlig(i,j)='X';
       end
    end
end

seqFinal = '';

for i=1:length(seq)
    temp = seq{i};
    length(temp)
    seqFinal(i*2-1,:)=temp;
end

for i=1:length(seq)-1
    seqFinal(i*2,:)=seqAlig(i,:);
end

    
