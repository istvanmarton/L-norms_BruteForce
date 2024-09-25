n=17;
r=100;
fileName = sprintf("mtx%d.txt", n);
A = randi(2*r+1,n,n)-r-1;
dlmwrite(fileName, A, ' ')
