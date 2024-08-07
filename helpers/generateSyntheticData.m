function data = generateSyntheticData(numPoints, d, stdDevs)
    if length(stdDevs) ~= d
        error(['Error in ' mfilename ': The length of stdDevs must be equal to the number of dimensions d.']);
    end

    covMatrix = diag(stdDevs.^2);

    mu = ones(1, d);

    data = mvnrnd(mu, covMatrix, numPoints);
end
