clc;
A = Utility.generateRandomIntMatrix(50, 6000);
B = Utility.generateRandomIntMatrix(6000, 50);
C = Utility.generateRandomIntMatrix(50, 50);
nruns = 100;

tic;
for i = 1:nruns
    tr1 = trace(A * B * C);
end
elapsedTime = toc;
fprintf('Tr(ABC): %.4f seconds\n', elapsedTime);


tic;
for i = 1:nruns
    tr2 = trace(B * C * A);
end
elapsedTime = toc;
fprintf('Tr(BCA): %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    tr3 = trace(C * A * B);
end
elapsedTime = toc;
fprintf('Tr(CAB): %.4f seconds\n', elapsedTime);



tic;
for i = 1:nruns
    D = (A * B)';
    tr4 = sum(sum(D .* C));
end
elapsedTime = toc;
fprintf('Elapsed time (AB)''C: %.4f seconds\n', elapsedTime);


tic;
for i = 1:nruns
    D = A * B;
    tr5 = sum(sum(C' * D));
end
elapsedTime = toc;
fprintf('Elapsed time (AB)C'': %.4f seconds\n', elapsedTime);


tic;
for i = 1:nruns
    D = (B * C)';
    tr6 = sum(sum(D .* A));
end
elapsedTime = toc;
fprintf('Elapsed time (BC)''A: %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    D = B * C;
    tr7 = sum(sum(A' .* D));
end
elapsedTime = toc;
fprintf('Elapsed time (BC)A'': %.4f seconds\n', elapsedTime);


return;
tic;
for i = 1:nruns
    D = C * A;
    tr8 = D(:)' * B(:);
end
elapsedTime = toc;
fprintf('Elapsed time (CA)''B: %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    D = C * A;
    tr9 = B(:)' * D(:);
end
elapsedTime = toc;
fprintf('Elapsed time (CA)B'': %.4f seconds\n', elapsedTime);