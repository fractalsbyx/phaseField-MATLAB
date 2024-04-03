close all;

models = ["gentry", "gradprod", "cubic"];
dt = 0.02;
cycles = 500;
t = dt*(1:cycles);
figure(1)
hold on
for i=1:length(models)
    mode = models(i);
    fileID = fopen("../data/IntegratedForce_" + mode + ".dat");
        force = fread(fileID, cycles, 'double')./dt;
    fclose(fileID);
    plot(t, force, 'LineWidth',3);
end
legend(models)