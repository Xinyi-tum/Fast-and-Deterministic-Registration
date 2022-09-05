function out_B = Divide(branch)

a = branch(1:3);
b = branch(4:6);
c = 0.5*(a+b);
m = [a,c,b];
for i = 1:8
    out(1,i) = m(1,bitget(i,1)+1);
    out(2,i) = m(2,bitget(i,2)+1);
    out(3,i) = m(3,bitget(i,3)+1);
    out(4,i) = m(1,bitget(i,1)+2);
    out(5,i) = m(2,bitget(i,2)+2);
    out(6,i) = m(3,bitget(i,3)+2);
end
out_B = out;
end