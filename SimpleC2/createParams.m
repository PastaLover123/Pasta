function [ params ] = createParams()

params = struct();

%The dataset number, also serves to name the set and determine seed
params.setNum=9999;

params.baseSeed = params.setNum*1000;
%params.totRuns = 1;

params.is3D = true; %true for 3D, false for 2D
params.L = 64; %width of grid, ideally should be power of 2

%array of dimensions of the grid
if params.is3D
    params.griddim = [params.L params.L params.L];
else
    params.griddim = [params.L params.L];
end

%what initial psi is. Default is liquid.
params.specificinit = '';%'seed2D';


%step parameters. maxStep should be a multiple of outStep.
params.maxStep = 21500; params.outStep = 1000; params.printStep = params.outStep;

%Remember to set this to false for remote runs
params.dispOutput=true;


params.output_dir = 'C:\Users\justin\Downloads\kenpfc_nuclear_pasta\kenpfc_nuclear_pasta';

%Whether to output as hdf5 file. If true, each run gets its own h5 file. If
%false, all the runs are stored together as a cell array in one mat file.
%For small runs/small sizes, better to output as mat file, since it's setup
%in that case to store all data in memory until the run is done, thus
%minimizing disk access and file compression.
params.fullout_hdf5=false; %whether to output the full field at every outstep in hdf5
%chunksize for hdf5 data storage (needs a bit of black magic to choose these better...???)
if params.is3D
    params.chunksize = [params.L params.L 1 1];
else
    params.chunksize = [params.L params.L 1];
end
params.resizing_hdf5 = true; %if true, uses an 'unlimited' dimension in the time direction. Otherwise, sets the size of the file based on maxStep.

%Params for making a rough first-pass movie, won't look as good as a properly formatted movie made with finished data
params.makeMovie = false;
params.movieFramerate = 5;%frames per second

%'lattice' constants for the three ordered phases
%(careful: 3D and 2D in Nik's book seem incorrect...?)
params.a3D=2*pi*sqrt(2);
params.a2D=2*pi*2/sqrt(3);
params.a1D=2*pi;

%k-space vectors for 1-mode approximation
params.q1 = [-1/2; -sqrt(3)/2];
params.q2 = [1; 0];
params.q3 = [-1/2; +sqrt(3)/2];

%length and time steps
params.dx = params.a1D/10;
params.dt = 1/2;%1/10;

params.n0 = 0.5;% try these: 0.95;0.8;0.50;0.38;0.20;.05;


%effective temperature parameters, Del_B = Bl - Bx
params.Bx = 1.0;
params.Bl = params.Bx + 0.000;%0.155

%init noise
params.chi = 0.01;

%numerical stability parameter, implements a splitting term to ensure better stability
params.stabP = 1;

%k-space cutoffs, for testing effects of higher modes
params.kcut_active = false;
params.kcutoff = 1.1;%value of k, maybe should be something else... for now, recall k=1.0 is radius of the ring.

params.startStrain = 20000;
params.stopStrain = 21500;
params.strain_rate = -1e-4;

end

