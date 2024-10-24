%% Process data from the simulation runs  
function fermionization_analysis()

    % Collecting energy per particle and coherence from simulation data
    for N = [2 3]

        % Load parameters
        run(sprintf('fermionization_N%d_parameters',N));

        for opt = ["superfluid" "pinned"]

            energy = NaN(length(gMIXarr),length(interactions));
            coherence = NaN(length(gMIXarr),length(interactions)); 

            for i = 1:length(gMIXarr)

                gMIX = gMIXarr(i);
                fname = sprintf('data/fermionization_N%d_%s_data/fermionization_N%d_gMIX%3.2f_gTGramp_%s.mat',N,opt,N,gMIXarr(i),opt);
                fname = fullfile(fileparts(pwd),fname); 

                try
                    load(fname,"data");
                    mubar = 0.5*gBEC*(NBEC + gMIX*N/gBEC)/posmax;

                    for j = 1:length(interactions)
                        if j > length(data)
                            energy(i,j) = NaN;
                            coherence(i,j) = NaN; 
                        else
                            energy(i,j) = (data(j).energy(end) - mubar^2*posmax/gBEC)/N; % Energy per particle, adjusted for background BEC energy
                            coherence(i,j) = data(j).occupationnr(1);
                        end

                    end
                catch
                    energy(i,:) = NaN;
                    coherence(i,:) = NaN; 
                end
            end
            
            eval(sprintf('energy_N%d_%s = energy;',N,opt));
            eval(sprintf('coherence_N%d_%s = coherence;',N,opt));
        end    
    end
    clearvars -except gMIXarr interactions energy_N* coherence_N* 

    % Determining phase diagrams 
    for N = [2 3]

        % Load parameters
        run(sprintf('fermionization_N%d_parameters',N));

        energy_num = NaN(length(gMIXarr),length(interactions));
        energy_ana = NaN(length(gMIXarr),length(interactions));
        coherence_num = NaN(length(gMIXarr),length(interactions));
        gTG_crit = []; 

        for opt = ["superfluid" "pinned"] 
            eval(sprintf('energy_%s = energy_N%d_%s;',opt,N,opt));
            eval(sprintf('coherence_%s = coherence_N%d_%s;',opt,N,opt));
        end

        for i = 1:length(gMIXarr)

            gMIX = gMIXarr(i);

            % Analytical energies for pinned and superfluid phase according to the model detailed in the paper 
            mubar = 0.5*gBEC*(NBEC + gMIX*N/gBEC)/posmax;
            a0 = 0.5*gMIX^2/gBEC; 
            eps = 6*a0^2/(5*mubar);
            ap = a0*(sqrt(1+2*eps)-1)/eps;

            % Pinned-state energy per particle, independent from gTG 
            Epin = gMIX^2*ap^3/(30*mubar*gBEC) + ap^2/6 - gMIX^2*ap/(6*gBEC);

            % Superfluid-state energy per particle, model only valid if resulting asf > 0
            asf = a0*(sqrt(1 + 2*N^2*eps*(1 - gBEC*interactions/gMIX^2*(1 - 1/N))) - 1)/(N*eps);
            energy_ana(i,:) = gMIX^2*N*asf.^3/(30*mubar*gBEC) + asf.^2/6 - N*gMIX^2*asf/(6*gBEC) + (N-1)*(interactions.*asf)/6;
            energy_ana(i,asf<0) = NaN;

            % We use the analytical value for Epin to determine the phase boundary in order to minimize finite-size effects. See paper for details. 
            ind = energy_superfluid(i,:) < Epin;
            energy_num(i,ind) = energy_superfluid(i,ind);
            energy_num(i,~ind) = energy_pinned(i,~ind);
            coherence_num(i,ind) = coherence_superfluid(i,ind);
            coherence_num(i,~ind) = coherence_pinned(i,~ind);

            % Intraspecies interaction strength at phase boundary 
            if isempty(find(ind,1,'last'))
                gTG_crit(i).ind = [];
                gTG_crit(i).val = NaN;
            else
                gTG_crit(i).ind = find(ind,1,'last') + 1;
                gTG_crit(i).val = interactions(find(ind,1,'last') + 1);
            end
        end
        
        eval(sprintf('energy_N%d_ana = energy_ana;',N));
        eval(sprintf('energy_N%d_num = energy_num;',N));
        eval(sprintf('coherence_N%d_num = coherence_num;',N));
        eval(sprintf('gTG_crit_N%d = gTG_crit;',N));
    end
    clearvars -except gMIXarr interactions energy_N* coherence_N* gTG_crit_N*

    % Calculating numerical and analytical density data
    for gMIX = 2
        for N = [2 3]

            % Load parameters
            run(sprintf('fermionization_N%d_parameters',N));
            x = fftdef(posmax,Ngrid);

            % Collect numerical data from simulations 
            for opt = ["superfluid" "pinned"] 
                density = NaN(Ngrid,length(interactions));
                fname = sprintf('data/fermionization_N%d_%s_data/fermionization_N%d_gMIX%3.2f_gTGramp_%s.mat',N,opt,N,gMIX,opt);
                fname = fullfile(fileparts(pwd),fname);
                try
                    load(fname,"data");
                    for i = 1:length(interactions)
                        density(:,i) = data(i).rho;
                    end
                catch
                end
                eval(sprintf('density_%s = density;',opt));
                eval(sprintf('density_N%d_gMIX%1.0f_%s = density;',N,gMIX,opt));
            end
            clear density fname i data

            % Compose density according to phase boundary 
            density_num = NaN(Ngrid,length(interactions));
            eval(sprintf('ind = gTG_crit_N%d(gMIXarr==gMIX).ind;',N));
            density_num(:,1:ind-1) = density_superfluid(:,1:ind-1); 
            density_num(:,ind:end) = density_pinned(:,ind:end); 

            eval(sprintf('rho_N%d_gMIX%1.0f_sf = density_num(:,1);',N,gMIX));
            eval(sprintf('rho_N%d_gMIX%1.0f_crit = density_num(:,ind-1);',N,gMIX));
            eval(sprintf('rho_N%d_gMIX%1.0f_pin = density_num(:,end);',N,gMIX));
            
            % Rescale density used in heatmap plots for better visibility 
            for i = 1:length(interactions)
                density_num(:,i) = density_num(:,i)/max(density_num(:,i));
            end
            eval(sprintf('density_N%d_gMIX%1.0f_num = density_num;',N,gMIX));

            % Calculate analytical approximations to the density 
            mubar = 0.5*gBEC*(NBEC + gMIX*N/gBEC)/posmax;
            a0 = 0.5*gMIX^2/gBEC; 
            eps = 6*a0^2/(5*mubar);

            % Pinned-state model 
            ap = a0*(sqrt(1+2*eps)-1)/eps;
            xi = x(find(islocalmax(density_num(:,end)),1,'last')); % align position of density peaks of the analytical model with numerical results for easier comparison 
            switch N
                case 2
                    rho_pin_ana = 0.5*ap*(1./cosh(ap*(x-xi)).^2 + 1./cosh(ap*(x+xi)).^2);
                case 3
                    rho_pin_ana = 0.5*ap*(1./cosh(ap*(x-xi)).^2 + 1./cosh(ap*x).^2 + 1./cosh(ap*(x+xi)).^2);
            end

            % Superfluid-state model for non-interacting case 
            gTG = 0;
            asf = a0*(sqrt(1 + 2*N^2*eps*(1 - gBEC*gTG/gMIX^2*(1 - 1/N))) - 1)/(N*eps);
            rho_sf_ana = 0.5*N*asf./cosh(asf*x).^2;
            
            eval(sprintf('rho_N%d_gMIX%1.0f_pin_ana = rho_pin_ana;',N,gMIX));
            eval(sprintf('rho_N%d_gMIX%1.0f_sf_ana = rho_sf_ana;',N,gMIX));
        end
    end

    clearvars -except gMIXarr interactions energy_N* coherence_N* gTG_crit_N* rho_N* density_N*

    % Save and proceed to figure creation 
    save(fullfile(fileparts(pwd),'data/fermionization_data_processed.mat'))
    fermionization_figures()
end