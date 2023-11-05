
tic;disp('[[[[[ Preparing for simulation run... ]]]]]');

%premadeParamID = -1;

params=createParams(); %even if we want to use previously made params, we need to use the params generator script for the location of data output
if exist('premadeParamID','var') && premadeParamID>0 %replace the params with a previously saved one, if needed
    filename_params = strcat('o_',num2str(premadeParamID),'_params.mat');
    fullfilename_params = fullfile(params.output_dir,filename_params);
    load(fullfilename_params,'params');
end
clear('premadeParamID');%to avoid accidents during testing

toc;tic;disp(['[[[[[ Run ID is ' num2str(params.setNum) ' ]]]]]']);

gpuDevice([]);
gpu=gpuDevice();
parallel.gpu.rng(params.baseSeed, 'Philox4x32-10');



%Prepare operators on GPU
toc;tic;disp('[[[[[ Preparing operators... ]]]]]');
dx = params.dx; dy = params.dx; dz = params.dx;
[ operators_h, operators_d ] = prepareOperators( params, gpu, dx, dy, dz );

%Prepare output files
toc;tic;disp('[[[[[ Preparing output... ]]]]]');
if params.fullout_hdf5
    filename_data = strcat('o_',num2str(params.setNum),'_data.h5');
    fullfilename_data = fullfile(params.output_dir,filename_data);
    if params.resizing_hdf5
        hdf5timedim = inf;
    else
        hdf5timedim = params.maxStep/params.outStep;
    end
    h5create(fullfilename_data,'/psi',cat(2,params.griddim,hdf5timedim),'ChunkSize',params.chunksize,'Deflate',9);
end
filename_params = strcat('o_',num2str(params.setNum),'_params.mat');
fullfilename_params = fullfile(params.output_dir,filename_params);
save(fullfilename_params,'params');

if params.makeMovie
    filename_video1=strcat('s',num2str(params.setNum),'_simulmov_psi.avi');
    fullfilename_video1 = fullfile(params.output_dir,filename_video1);
    v1 = VideoWriter(fullfilename_video1);
    v1.FrameRate = params.movieFramerate;
    open(v1)
    filename_video2=strcat('s',num2str(params.setNum),'_simulmov_psiMF.avi');
    fullfilename_video2 = fullfile(params.output_dir,filename_video2);
    v2 = VideoWriter(fullfilename_video2);
    v2.FrameRate = params.movieFramerate;
    open(v2)
end


%Prepare simulation fields
toc;tic;disp('[[[[[ Preparing simulation fields... ]]]]]');
fields_h = struct();
fields_h.psi   = zeros(params.griddim);
fields_h.psi_F = complex(zeros(params.griddim));
status = struct(); %contains non-field values that change during the simulation, probably doesn't require loading on GPU as long as it's not large data?
status.timestep = 0;
status.interfacePos = 0;
status.debug = 0;
status.equlibrate = false;

%Move host fields to gpu device
toc;tic;disp('[[[[[ Moving fields to GPU device... ]]]]]');
fields_d = struct();
fields_d.psi = gpuArray(fields_h.psi);wait(gpu);
fields_d.psi_F = gpuArray(fields_h.psi_F);wait(gpu);

%Initialize simulation fields
toc;tic;disp('[[[[[ Initializing simulation fields... ]]]]]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial state

fields_d.psi = params.chi*(rand(params.griddim,'gpuArray')*2-1) + params.n0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields_d.psi_F = fftn(fields_d.psi);
%fields_d.psi_SOD = fields_d.psi;
% fields_d.psi_F_SOD = fields_d.psi_F;

%return

if params.dispOutput
    if params.is3D
        imagesc(fields_d.psi(:,:,1));colorbar
    else
        imagesc(fields_d.psi);colorbar
    end
    pause(1.0);
end

%return

toc;tic;disp('[[[[[ Starting main simulation loop... ]]]]]');
counter = 0;
for step=1:params.maxStep
    
    status.timestep = step;

    if step>params.startStrain && step<=params.stopStrain
        dt = params.dt;
        dx = dx * (1.0 + params.strain_rate*dt);
        dy = dy/sqrt(1.0 + params.strain_rate*dt);
        dz = dz/sqrt(1.0 + params.strain_rate*dt);
        [ operators_h, operators_d ] = prepareOperators( params, gpu, dx, dy, dz );
    end

    [fields_d,status] = stepPFC( params, gpu, fields_d, operators_d, status );
    
    if mod(step,params.outStep)==0
        fields_h.psi = gather(fields_d.psi);wait(gpu);
        fields_h.psi_F = gather(fields_d.psi_F);wait(gpu);

        isovalue = 0.5*max(max(max(fields_d.psi)));

    clf;
    isosurface(fields_d.psi,isovalue)
    isocaps(fields_d.psi,isovalue)
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'ztick',[])
    set(gca,'zticklabel',[])
    
    xlabel('X (fm)', 'fontSize', 12)
    ylabel('Y (fm)', 'fontSize', 12)
    zlabel('Z (fm)', 'fontSize', 12)
    title('Density Field (dimensionless)', 'fontSize', 12)
    colorbar
    %set(get(colorbar,'title'), 'fontSize', 14)
    
    xlim([1 64])
    ylim([1 64])
    zlim([1 64])
    
    set(gca, 'linewidth', 1.5, 'fontsize', 16)
    %set(gcf, 'Position',  [100, 100, 800, 700])
    
    view(45, 45)
    axis equal
    name = strcat(".\plots\figure_",num2str(counter),".png");
    counter=counter+1;
    saveas(gcf, name, 'png')
        
        if nnz(isnan(fields_h.psi))>0
            error(['ERROR: OBTAINED NAN VALUES BY STEP ' num2str(step)])
        end
        
        if params.fullout_hdf5
            if params.is3D
                writestart=[1 1 1 step/params.outStep];
            else
                writestart=[1 1 step/params.outStep];
            end
            h5write(fullfilename_data,'/psi',fields_h.psi,writestart,cat(2,params.griddim,1));
        end
        
        if params.makeMovie
            if params.is3D
                g = exp(operators_h.kLap(:,:,1)/(2*(0.2)^2));
                data1 = fields_h.psi(:,:,1);
                data2 = real(ifftn(fftn(data1).*g));
            else
                g = exp(operators_h.kLap/(2*(0.2)^2));
                data1 = fields_h.psi;
                data2 = real(ifftn(fftn(data1).*g));
            end
            minv1 = min(data1(:));
            maxv1 = max(data1(:));
            frame=(data1-minv1)/(maxv1-minv1);
            writeVideo(v1,frame);
            minv2 = min(data2(:));
            maxv2 = max(data2(:));
            frame=(data2-minv2)/(maxv2-minv2);
            writeVideo(v2,frame);
        end
        
        if params.dispOutput
            if params.is3D
                %tmp = fields_h.psi_F(:,:,1);tmp(1,1) = 0;
                %imagesc(abs(tmp));colorbar
                imagesc(fields_h.psi(:,:,1));colorbar
            else
                imagesc(fields_h.psi);colorbar
            end
            
            pause(.1);
        end
        
    end
    
    
    
    if mod(step,params.printStep)==0
        timed = toc;
        
        etatime = (params.maxStep - step)/params.printStep*timed;
        disp(  ['### Step ' num2str(step) ' / ' num2str(params.maxStep) ' (' num2str(round(step/params.maxStep*100)) '%) ###']  );
        disp( ['### Time since last print out: ' num2str(timed) ' seconds. ### ETA: ' char(duration(0,0,etatime)) ' ###' ] )

        tic
    end
    
    
    
end

if params.makeMovie
    close(v1)
    close(v2)
end

disp('[[[[[ Done. ]]]]]');















