% takes a matrix of intcal13 and speleothem data (cal age BP, error, 
% a14Cmeas, error) and returns dcp and error.

function [results] = dcpcalc(speleo)

format shortG

% intcal columns = cal age, 14C age, error, D14C, error
intcal = load('IntCal13.txt');

lambda = 1/8267;

rounded = zeros(1,1);

for i =1:1:size(speleo)
    if(speleo(i,1)<=13900)
        rounded = vertcat(rounded, round(speleo(i,1)/5)*5);
    elseif((speleo(i,1)>13900) && (speleo(i,1)<=25000))
            rounded = vertcat(rounded, round(speleo(i,1)/10)*10);
    elseif((speleo(i,1)>25000) && (speleo(i,1)<=50000))
                rounded = vertcat(rounded, round(speleo(i,1)/20)*20);
    end
end 

rounded = rounded(2:end);

list = zeros(1,1);
for i=1:size(rounded(:,1))
    list(i) = find(rounded(i,1)==intcal(:,1));
end

list=list';

values = zeros(1,5);

for i=1:size(list(:,1))
    values(i,:) = intcal(list(i),:);
end

a14Catminit = (((values(:,4)/1000)+1)*100);
a14Catminit_error = values(:,5)/10;
atminit = horzcat(values, a14Catminit, a14Catminit_error);

a14Cinit = speleo(:,3)./exp(-1*lambda*speleo(:,1));

a = ((speleo(:,4)./speleo(:,3)).^2);
b = ((lambda*speleo(:,2)).^2);
a14Cinit_error = a14Cinit.*sqrt(a+b);

dcp = (1-a14Cinit./a14Catminit).*100;
a1 = (a14Cinit_error./a14Cinit).^2;
b1 = (atminit(:,7)./atminit(:,6)).^2;
c1 = (a14Cinit./atminit(:,6))*100;
dcp_error = c1.*sqrt(a1+b1);

results = horzcat(dcp, dcp_error);
