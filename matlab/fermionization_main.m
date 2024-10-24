%% Main file controlling the simulations for finding the ground state of a combined system of N interacting particles immersed into a BEC in 1D 
% Depending on the available hardware the outer loops should be parallelized as much as possible, e.g. by running separate instances for each combination of N = [2, 3] and opt = ["superfluid", "pinned"].
% Also by further replacing the loop over gMIX by gMIX = gMIXarr(${SLURM_ARRAY_TASK_ID}) in each instance when running the code as an array job on a cluster using slurm. 
function fermionization_main()
    % Number of immersed particles
    for N = [2 3]

        % Load simulation parameters and pack into struct
        run(sprintf('fermionization_N%d_parameters',N));
        simpara = v2struct(N,Ngrid,NBEC,gBEC,posmax,LTG,wall,dt,cutoff,steps,samples,xi,gMIXarr,interactions); 

        % Simulations are performed twice. Once starting in the superfluid phase and slowly ramping up the immersed intra-particle interaction strength from gTG = 0.
        % And then once the other way round, starting in the pinned phase and slowly decreasing gTG.
        for opt = ["superfluid" "pinned"]

            % Flip interactions array for ramping down gTG if starting in the pinned phase 
            simpara.opt = opt;
            if opt == "pinned"
                simpara.interactions = fliplr(interactions);
            end

            % Loop over interspecies interaction strength 
            for j = 1:length(gMIXarr)

                % Set parameter and filename 
                simpara.gMIX = gMIXarr(j);
                simpara.fname = sprintf('data/fermionization_N%d_%s_data/fermionization_N%d_gMIX%3.2f_gTGramp_%s.mat',N,opt,N,gMIXarr(j),opt);
                simpara.fname = fullfile(fileparts(pwd),simpara.fname); 

                try % Resume previous simulation or skip if already finished
                    load(simpara.fname,"data","simpara");
                    if length(data) == length(simpara.interactions) 
                            continue % Replace by 'exit' if parallelizing gMIX loop
                    else
                        load(simpara.fname,"wfi","psi_ini");
                    end
                    simpara.start = length(data) + 1;
                    % Main simulation loop ramping over gTG 
                    fermionization_ramp(simpara,wfi,psi_ini); 
                    
                catch ME % Otherwise start new simulation 
                    fprintf("%s\nStarting from the beginning...\n",ME.message);
                    simpara.start = 1;

                    % Setting position grid and initial homogeneous BEC wave function
                    [x,dx] = fftdef(posmax,Ngrid);
                    simpara.dx = dx;
                    wfi = sqrt(NBEC/(Ngrid*dx))*ones(Ngrid,1);
                    simpara.ETG_ini = 0.5*N*(pi/LTG)^2;

                    % Setting initial wave function similar to the expected ground state in each phase for speed-up
                    switch N
                        case 2
                            switch opt
                                case "pinned"
                                    psi_ini = exp(-0.5*((x-xi).^2 + ((x+xi)').^2)) + exp(-0.5*((x+xi).^2 + ((x-xi)').^2));
                                case "superfluid"
                                    psi_ini = exp(-0.5*(x.^2 + (x').^2));
                            end
                            psi_ini(abs(x)>LTG/2,:) = 0; psi_ini(:,abs(x)>LTG/2) = 0; 
                            psi_ini = psi_ini./sqrt(sum(sum(abs(psi_ini).^2))*dx*dx);
                        case 3
                            [xm,ym,zm] = meshgrid(x,x,x);
                            switch opt
                                case "pinned"
                                    psi_ini = exp(-0.5*((xm-xi).^2 + ym.^2 + (zm+xi).^2)) + exp(-0.5*(xm.^2 + (ym-xi).^2 + (zm+xi).^2)) + exp(-0.5*((xm+xi).^2 + (ym-xi).^2 + zm.^2)) + exp(-0.5*((xm-xi).^2 + zm.^2 + (ym+xi).^2)) + exp(-0.5*(xm.^2 + (zm-xi).^2 + (ym+xi).^2)) + exp(-0.5*((xm+xi).^2 + (zm-xi).^2 + ym.^2));
                                case "superfluid"
                                    psi_ini = exp(-0.5*(xm.^2 + ym.^2 + zm.^2));
                            end
                            psi_ini(abs(x)>LTG/2,:,:) = 0; psi_ini(:,abs(x)>LTG/2,:) = 0; psi_ini(:,:,abs(x)>LTG/2) = 0;
                            psi_ini = psi_ini./sqrt(sum(abs(psi_ini).^2,'all')*dx^3);
                    end

                    % Main simulation loop ramping over gTG 
                    fermionization_ramp(simpara,wfi,psi_ini);
                end

            end
        end
    end

    % Continue to data processing
    fermionization_analysis()
end