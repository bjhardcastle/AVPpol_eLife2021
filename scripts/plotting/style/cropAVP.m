%[alpha, v] = crop_borders(alpha, 0, 1, options.crop_amounts);
%A = A(v(1):v(2),v(3):v(4),:);
[alpha, vA, vB] = crop_borders(1, 0, 0,[nan, nan, nan, nan]);
if ~any(isnan(vB)) % positive padding
    B = repmat(uint8(zeros(1,1,size(A,3))),size(alpha));
    B(vB(1):vB(2), vB(3):vB(4), :) = A(vA(1):vA(2), vA(3):vA(4), :); % ADDED BY OH
    A = B;
else  % negative padding
    A = A(vA(1):vA(2), vA(3):vA(4), :);
end
