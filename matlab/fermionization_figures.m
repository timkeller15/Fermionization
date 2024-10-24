%% Create data files for the publication figures from processed data
function fermionization_figures()
    
    % Load data 
    load(fullfile(fileparts(pwd),'data/fermionization_data_processed.mat'))

    % Fig. 1 - Write data file for energy per particle as a function of intra-species interaction for gMIX = 2 and gMIX = 3
    for gMIX = [2 3]
        fname = sprintf('data/fermionization_gMIX%1.0f_gTGramp_energy.dat',gMIX);
        fname = fullfile(fileparts(pwd),fname); 
        header = ["gTG","E_N2_num","E_N2_ana","E_N2_sf","E_N3_num","E_N3_ana","E_N3_sf"];
        dataout = [interactions.' energy_N2_num(gMIXarr==gMIX,:).' energy_N2_ana(gMIXarr==gMIX,:).' energy_N2_superfluid(gMIXarr==gMIX,:).' energy_N3_num(gMIXarr==gMIX,:).' energy_N3_ana(gMIXarr==gMIX,:).' energy_N3_superfluid(gMIXarr==gMIX,:).'];
        dataout = [header; dataout];
        dataout(ismissing(dataout)) = 'NaN';
        writematrix(dataout,fname,'Delimiter','tab');
    end

    % Fig. 2 - Heatmap and line plots for immersed component density for gMIX = 2
    load(fullfile(fileparts(pwd),'data/smoothcoolwarm.mat')); % Colormap 
    for gMIX = 2
        for N = [2 3]

            % Load parameters
            run(sprintf('fermionization_N%d_parameters',N));
            x = fftdef(posmax,Ngrid);

            % Heatmap plots
            fname = sprintf('data/fermionization_N%d_gMIX%1.0f_density.png',N,gMIX);
            fname = fullfile(fileparts(pwd),fname);
            figure;
            set(gca,'position',[0 0 1 1]); hold all
            eval(sprintf('surf(interactions,x,density_N%d_gMIX%1.0f_num);',N,gMIX));
            % Restrict axis to relevant regimes 
            switch N
                case 2
                    axis([0 30 -3 3]);
                case 3
                    axis([0 18 -5 5]);
            end
            shading flat; view(2);
            colormap(smoothcoolwarm);
            clim([0 1]);
            saveas(gcf,fname);
            close(gcf);

            % Line Plots
            fname = sprintf('data/fermionization_N%d_gMIX%1.0f_gTGramp_densities.dat',N,gMIX);
            fname = fullfile(fileparts(pwd),fname);
            
            % Write data file for density line plots 
            header = ["x","rho_sf","rho_sf_ana","rho_pin","rho_pin_ana","rho_crit"];
            eval(sprintf('dataout = [x rho_N%d_gMIX%1.0f_sf rho_N%d_gMIX%1.0f_sf_ana rho_N%d_gMIX%1.0f_pin rho_N%d_gMIX%1.0f_pin_ana rho_N%d_gMIX%1.0f_crit];',N,gMIX,N,gMIX,N,gMIX,N,gMIX,N,gMIX));
            dataout = [header; dataout];
            dataout(ismissing(dataout)) = 'NaN';
            writematrix(dataout,fname,'Delimiter','tab');
        end
    end

    % Fig. 3 - Write data file for coherence as a function of intra-species interaction for gMIX = 2
    for gMIX = 2
        fname = sprintf('data/fermionization_gMIX%1.0f_gTGramp_coherence.dat',gMIX);
        fname = fullfile(fileparts(pwd),fname); 
        header = ["gTG","coherence_N2","coherence_N3"];
        dataout = [interactions.' coherence_N2_num(gMIXarr==gMIX,:).' coherence_N3_num(gMIXarr==gMIX,:).'];
        dataout = [header; dataout];
        dataout(ismissing(dataout)) = 'NaN';
        writematrix(dataout,fname,'Delimiter','tab');
    end

    % Fig. 4 - Create heatmap plot for phase diagram
    for N = [2 3]
        fname = sprintf('data/fermionization_N%d_phasediag.png',N);
        fname = fullfile(fileparts(pwd),fname); 
        figure;
        set(gca,'position',[0 0 1 1]); hold all;
        % We checked that also numerically the system behavior is symmetric with respect to the sign of gMIX and therefore only mirror the results for positive gMIX
        eval(sprintf('surf(interactions,gMIXarr,coherence_N%d_num);',N));
        hold all; 
        eval(sprintf('surf(interactions,-gMIXarr,coherence_N%d_num);',N));
        shading flat; view(2);
        colormap(bone(1024));
        % Restrict axis to relevant regimes 
        switch N
            case 2
                gmax = 100;
            case 3
                gmax = 40;
        end
        axis([0 gmax -5 5]); 
        clim([1/N 1]);
        saveas(gcf,fname); 
        close(gcf);
    end
end