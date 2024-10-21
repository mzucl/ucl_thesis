clc;
A = Utility.generateRandomIntMatrix(50, 6000);
B = Utility.generateRandomIntMatrix(6000, 50);
nruns = 10;

tic;
for i = 1:nruns
    tr1 = trace(A * B);
end
elapsedTime = toc;
fprintf('Elapsed time Tr(AB): %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    tr2 = trace(B * A);
end
elapsedTime = toc;
fprintf('Elapsed time Tr(BA): %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    D = A';
    tr3 = sum(D(:) .* B(:));
end
elapsedTime = toc;
fprintf('Elapsed time A''B: %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    D = B';
    tr4 = sum(D(:) .* A(:));
end
elapsedTime = toc;
fprintf('Elapsed time B''A: %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    tr5 = sum(sum(A' .* B));
end
elapsedTime = toc;
fprintf('Elapsed time A'' .* B: %.4f seconds\n', elapsedTime);

tic;
for i = 1:nruns
    tr6 = sum(sum(B' .* A));
end
elapsedTime = toc;
fprintf('Elapsed time B'' .* A: %.4f seconds\n', elapsedTime);