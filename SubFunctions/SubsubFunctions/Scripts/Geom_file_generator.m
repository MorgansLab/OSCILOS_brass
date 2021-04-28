% Script to sprout different discretisations of the geometry
% Clearly sprout is the wrong word for this
% I love sprouts

addpath('./SubFunctions/SubsubFunctions')
MaxFracRange = [0.1, 0.0005];
MinFracRange = [0.02, 0.0001];
StartNum = 41;
EndNum = 60;
Geoms = StartNum:EndNum;
NumGeoms = length(Geoms);

MaxFracs = logspace(log10(MaxFracRange(1)),log10(MaxFracRange(2)),NumGeoms);
MinFracs = logspace(log10(MinFracRange(1)),log10(MinFracRange(2)),NumGeoms);

infile = "C:/Users/amacl/OneDrive - Imperial College London/UROP 2020/Data/Norman 2013/Cornophone.csv";

for i = 1:NumGeoms
    outfile = sprintf("./Inputs/TWC/20TWC%d.txt",Geoms(i));
    Geom_file_utility(infile, outfile, MaxFracs(i), MinFracs(i));
end