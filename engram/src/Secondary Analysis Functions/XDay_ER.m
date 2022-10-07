function [Date1, Date2, ER_Correlation] = XDay_ER(ER1, ER2, firstDate, secondDate)

Date1 = firstDate;
Date2 = secondDate;


a = size(ER1,1);
b = size(ER2,1);

if a - b == 1
    ER2 = cat(1, ER2, 0);
elseif b - a == 1
    ER1 = cat(1, ER1, 0);
elseif abs(a-b) ~= 1 && a>b
    ER2 = cat(1, ER2, zeros(a-b, 1));
elseif abs(a-b) ~= 1 && b>a
    ER1 = cat(1, ER1, zeros(b-a, 1));
else
    disp('both event rate vectors have same size, good to go')
end

for n = 1:size(ER1,1)
    if isnan(ER1(n,1)) || isnan(ER2(n,1))
        ER1(n,1) = nan;
        ER2(n,1) = nan;
    elseif ER1(n,1) == 0 || ER2(n,1) == 0
        ER1(n,1) = nan;
        ER2(n,1) = nan;
    end
end
ER1(isnan(ER1)) = [];
ER2(isnan(ER2)) = [];

% calculate ER correlation
ER_Correlation = corr(ER1(:), ER2(:));

end