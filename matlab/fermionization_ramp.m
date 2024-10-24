    %% Main simulation loop studying the fermionization of an interacting N-particle system immersed into a BEC by ramping over gTG 
    function fermionization_ramp(simpara,wfi,psi_ini) 

    % Unpack parameters
    v2struct(simpara);

    % Load results from previous simulation if available
    if start > 1
        load(fname,"data");
    end

    % Ramping up/down intraspecies interaction strength depending on regime
    for i = start:length(interactions)

        % Setting parameters 
        fprintf('gTG = %3.2f:\n',interactions(i));
        funcname = sprintf('fermionization_N%d_groundstate',N);
        simpara.gTG = interactions(i);

        % Finding ground state
        tic; d = feval(funcname,simpara,wfi,psi_ini); toc 
        d.gTG = interactions(i);

        % Setting final wave functions from current run as initial states for next one for speed-up
        wfi = d.wf;
        psi_ini = d.psi;
        simpara.ETG_ini = d.energy_TG(end);

        % Diagonalize reduced single-particle density matrix (RSPDM) for occupation numbers / coherence
        switch N
            case 2
                rspdm = (d.psi*d.psi').*dx;
            case 3
                rspdm = (reshape(d.psi,Ngrid,Ngrid*Ngrid)*reshape(d.psi,Ngrid,Ngrid*Ngrid)').*dx^2;
        end

        EV = eig(rspdm);
        occupationnr = flipud(EV)*dx;
        d.occupationnr = occupationnr;
        clear EV rspdm occupationnr

        % Immersed system's wave function is removed to reduce necessary storage space
        d.psi = [];
        data(i) = d;
        clear d

        save(fname);
    end

    % Clean up and if necessary flip back data and intraspecies interaction array for easier processing
    simpara = rmfield(simpara,{'gTG','ETG_ini','start'});

    if opt == "pinned"
        simpara.interactions = fliplr(interactions);
        data = fliplr(data);
    end

    save(fname,"simpara","data");
end