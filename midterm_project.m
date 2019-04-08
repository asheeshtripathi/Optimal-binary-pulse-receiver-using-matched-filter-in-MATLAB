clear;
clc;
receiverData = load('midtermData.mat');
receivedSignal = receiverData.signal;
a = 0:length(receiverData.signal)-1;
z = receiverData.p_sampled;
figure(1),plot(a,receiverData.signal);
figure(2),plot(z);
numSymbols = receiverData.numSymbols;
samplesPerSymbol = receiverData.samplesPerSymbol;
p_pulseEnergy = receiverData.p_pulseEnergy;
idx = 1;
op = zeros(1,248);
for i = 1:10:2471
 y = conv(receiverData.signal(i:i+9),z);
 k =   y(10);
 if k > 0
     op(idx) = 1;
     idx = idx +1;
 else
     op(idx) = 0;
     idx = idx +1;
 end
end
c = dec2bin(op);
d = reshape(c',1,numel(c));
ch = char(bin2dec(reshape(d,8,[]).')).';
disp(ch);
no_of_zeros = numel(op)-nnz(op);
probability_of_zero = (no_of_zeros/numel(op));
probability_of_one = nnz(op)/numel(op);
source_entropy = (probability_of_one)*log2(1/probability_of_one)+(probability_of_zero)*log2(1/probability_of_zero);
efficiency1 = source_entropy/(1*(log2(2)));
G1 = [1,1,1,0,1,0,1,0,1];
gL = length(G1);
m = zeros(1,gL-1);
M1=[op,m];
mL = length(M1);
count = 0;
while((mL-count) >= gL)
msg9 = M1(1:gL);
rem = mod((msg9 + G1),2)
M1(1:gL) = rem;
j=1;
shift = 0;
while(j<=gL)
if(rem(j)~=0)
break;
else shift = j;
j = j + 1;
end
end
count = count + shift;
M1(1:shift) = [];
end
j = 0;
value = 0;
chksuml = length(M1);
for j = 1:chksuml % convert binary to decimal
if(M1(j) == 1)
value = value + (2^(chksuml-j));
end
end
crc8 = dec2hex(value);
disp(crc8);
 no_of_00 = 0;
 no_of_01 = 0;
 no_of_10 = 0;
 no_of_11 = 0;
 idx2 = 1;
 
for e = 1:2:length(op)-1
    if op(e) == 0
        if op(e+1) == 0
            no_of_00 = no_of_00 + 1;
            alp(idx2) = 1;
            
        else
            no_of_01 = no_of_01 + 1;
            alp(idx2) = 2;
           
        end
    else
        if op(e+1) == 0
            no_of_10 =  no_of_10 + 1;
            alp(idx2) = 3;
           
        else
             no_of_11 =  no_of_11 + 1;
             alp(idx2) = 4;
            
        end
    end
    idx2 = idx2 + 1;
end
p_00 = no_of_00 / (length(op)/2);
p_01 = no_of_01 / (length(op)/2);
p_10 = no_of_10 / (length(op)/2);
p_11 = no_of_11 / (length(op)/2);
entropy1 = (p_00 * log2(1/p_00)) + (p_01 * log2(1/p_01)) + (p_10 * log2(1/p_10)) + (p_11 * log2(1/p_11));
sym = (1:4);
prob = [p_00 p_01 p_10 p_11];
[dict, avglen] = huffmandict(sym,prob);
dict{1,2}
dict{2,2}
dict{3,2} 
dict{4,2}
alp1 = alp';
comp = huffmanenco(alp1,dict);
efficiency2 = entropy1/(avglen*(log2(2)));



