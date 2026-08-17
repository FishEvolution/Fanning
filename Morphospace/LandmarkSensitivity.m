function LandmarkSensitivity(image, N)

nSamples = length(N);

for i = 1:nSamples
    P = autoLdmk(image,N(i));
    harm = efa1(N(i),P,12,100);
    if i > 1
        D(i) = norm(harm(:)-lastHarm(:));
    end
    lastHarm = harm;
end

figure
plot(N,D);