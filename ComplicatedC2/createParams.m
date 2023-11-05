function [ params ] = createParams()

params = struct();

%The dataset number, also serves to name the set and determine seed
params.setNum=9999;

params.baseSeed = params.setNum*1000;
%params.totRuns = 1;

params.is3D = true; %true for 3D, false for 2D
params.L = 128; %width of grid, ideally should be power of 2, also dimensionless
%64, 128, 256, 512, 1024

%array of dimensions of the grid
if params.is3D
    params.griddim = [params.L params.L params.L];
else
    params.griddim = [params.L params.L];
end

%what initial psi is. Default is liquid.
params.specificinit = '';%'seed2D'; did found the params to be call anywhere, should try to put it to random thing to see if it will break


%step parameters. maxStep should be a multiple of outStep.
params.maxStep = 100000; params.outStep = 100; params.printStep = params.outStep; % caution don't forget to modify the max step to always finish in the steady state for the time step chosen

%Remember to set this to false for remote runs
params.dispOutput=true;



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
params.dx = params.a1D/10; % DIMENSIONLESS
%0.6, 0.3, 0.24, 0.12, 0.06
params.dt = 1/500; %, DIMENSIONLESS
% 1e-4 seem to work fine for small simulation, but if increase L, then it
% is breaking, must try bigger range to confirm good value of time step
params.n0 = 0.6;% try these: 0.95;0.8;0.50;0.38;0.20;.05; 
% assume n=(rhobar-rho)/rhobar (rhobar is the average density), so how to choose rho
%rhobar n0= rhobar -rho => rhobar = rho/(1-n0); what is rho??? How can I
%relate it to the density inside the pasta and what is the physical meaning
%of n0 OR rhobar is the liquid density at coexistence, if it is the case,
%what is that density of the liquid at coexistence, we have two possible
%liquid-ish, nuclear matter at fix density (know and perfect as it is) OR
%the liquid before the phase which is not 100% pourcente clear what it is
%OR the true liquid density at that phase (SHOULD GO CHECK IF THE DENSITY ALWAYS THE SAME IN THE LIQUID ????)
% ASSUME THE DENSITY BEFORE NUCLER PASTA TO BE THE CLASSICAL DENSITY STILL
% PRESENT IN THE DIFFERENR REGION

params.output_dir = 'G:\My Drive\McGill\Hackathon\noise';
params.name_k_file = "G:\My Drive\McGill\Hackathon\kenpfc_nuclear_pasta\kvalue";
params.name_effectivec2_file = "G:\My Drive\McGill\Hackathon\kenpfc_nuclear_pasta\effectiveC2";

% file:///C:/Users/Justin/Downloads/SSRN-id4568913.pdf  
%init noise
params.chi = 0.05;

%numerical stability parameter, implements a splitting term to ensure better stability
params.stabP = 1;

%k-space cutoffs, for testing effects of higher modes
params.kcut_active = false;
params.kcutoff = 1.1;%value of k, maybe should be something else... for now, recall k=1.0 is radius of the ring.

%Name prof copolymer
%An-Chang Shi, physics

% dynamics noise
params.sigma = 0.08;

end

