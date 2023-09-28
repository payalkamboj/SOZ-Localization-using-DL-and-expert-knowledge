function y = completeInside(B,A)
%% A is the bigger one
p = inpoly2(B',A');

if(sum(p) == size(B,2))
    y = 1;
else
    y = 0;
end
    