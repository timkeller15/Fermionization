%% Find the combined system ground state for N = 3 particles immersed into the BEC via the split-step method in imaginary time
function data = fermionization_N3_groundstate(simpara,wfi,psi_ini) 

    % Unpack parameters
    v2struct(simpara);

    % Construct position and momentum grids
    [x,dx,px,~] = fftdef(posmax,Ngrid);
    [pxm,pym,pzm] = meshgrid(px,px,px);

    % Set external potentials
    Vtrap = zeros(Ngrid,1);
    TG_trap = (gTG/dx)*(repmat(eye(Ngrid),1,1,Ngrid) + repmat(reshape(eye(Ngrid),Ngrid,1,Ngrid),1,Ngrid,1) + repmat(reshape(eye(Ngrid),1,Ngrid,Ngrid),Ngrid,1,1));  
    TG_trap(abs(x)>LTG/2,:,:) = wall;
    TG_trap(:,abs(x)>LTG/2,:) = wall;
    TG_trap(:,:,abs(x)>LTG/2) = wall;

    samplestep = floor(steps/samples);

    % Initial wave functions and immersed line density
    wf = wfi;
    psi = psi_ini;
    rho = N*sum(abs(psi).^2,[2 3])*dx^2; 

    % Allocate observable arrays 
    samples = samples + 2;
    energy = zeros(1,samples);
    energy_BEC = zeros(1,samples);
    energy_TG = zeros(1,samples);
    Natom = zeros(2,samples);
    time = zeros(1,samples);

    % Operators for imaginary time split-step algorithm 
    Ekin = exp(-0.5*dt*px.^2);
    Ekin_TG = exp(-0.5*dt*(pxm.^2 + pym.^2 + pzm.^2));

    % Calculate initial energy and check normalization 
    Ki = 0.5*conj(wfi).*ifft(px.^2.*fft(wfi));
    Vi = (Vtrap + gMIX*rho).*abs(wfi).^2 + 0.5*gBEC.*abs(wfi).^4;
    Ei = real(sum(Ki + Vi))*dx;

    energy_BEC(1) = Ei;
    energy_TG(1) = ETG_ini;
    energy(1) = energy_BEC(1) + energy_TG(1);
    Natom(1,1) = sum(abs(wf).^2)*dx;
    Natom(2,1) = sum(abs(psi).^2,'all')*dx^3;

    count = 1;

    for i=1:steps
        
        % Step 1 - Imaginary time split-step for BEC component
        V = Vtrap + gMIX*rho + gBEC*abs(wf).^2;
        wf = exp(-0.5*dt*V).*wf;

        % Step 2
        wf = ifft(Ekin.*fft(wf));

        % Step 3
        wf = exp(-0.5*dt*V).*wf;

        % Normalize
        wf = wf/sqrt(sum(conj(wf).*wf)*dx/NBEC);
 
         % Step 1 - Update potential and imaginary time split-step for immersed component 
        Vint = gMIX*(repmat(abs(wf).^2,1,Ngrid,Ngrid) + repmat(abs(wf.').^2,Ngrid,1,Ngrid) + repmat(reshape(abs(wf).^2,1,1,Ngrid),Ngrid,Ngrid,1));   
        V_TG = TG_trap + Vint;

        psi = exp(-0.5*dt*V_TG).*psi;

        % Step 2
        psi = ifftn(Ekin_TG.*fftn(psi));

        % Step 3
        psi = exp(-0.5*dt*V_TG).*psi;

        % Normalize
        psi = psi/sqrt(sum(conj(psi).*psi,'all')*dx^3);
        rho = N*sum(abs(psi).^2,[2 3])*dx^2; 

        % Calculate observables 
        if mod(i,samplestep) == 0
            count = count + 1;
            time(count) = dt*i;

            % Energy 
            K_BEC = 0.5*conj(wf).*ifft(px.^2.*fft(wf));
            V_BEC = (Vtrap + gMIX*rho).*abs(wf).^2 + 0.5*gBEC.*abs(wf).^4;
            K_TG = 0.5*conj(psi).*ifftn((pxm.^2 + pym.^2 + pzm.^2).*fftn(psi));
            V_TG = TG_trap.*abs(psi).^2; % Not counting VInt twice (only added to BEC energy)

            energy_BEC(count) = real(sum(K_BEC + V_BEC))*dx;
            energy_TG(count) = real(sum(K_TG + V_TG,'all'))*dx^3;
            energy(count) = energy_BEC(count) + energy_TG(count);

            % Normalization 
            Natom(1,count) = sum(abs(wf).^2)*dx;
            Natom(2,count) = sum(abs(psi).^2,'all')*dx^3;

            % Check for convergence and exit if reached 
            if count > 2
                diff = abs(energy(count) - energy(count-1));
                if diff < cutoff
                    time = time(1:count);
                    energy_BEC = energy_BEC(1:count);
                    energy_TG = energy_TG(1:count);
                    energy = energy(1:count);
                    Natom = Natom(:,1:count);
                    steps = i;
                    % Only keep real part of wave function to save storage space 
                    if sum(imag(psi),'all') < 1e-6
                        psi = real(psi);
                    end
                    data = v2struct(wf,psi,rho,energy,energy_BEC,energy_TG,time,Natom);
                    return
                end
            end
        end
    end

    % Calculate final energy and normalization 
    K_BEC = 0.5*conj(wf).*ifft(px.^2.*fft(wf));
    V_BEC = (Vtrap + gMIX*rho).*abs(wf).^2 + 0.5*gBEC.*abs(wf).^4;
    K_TG = 0.5*conj(psi).*ifftn((pxm.^2 + pym.^2 + pzm.^2).*fftn(psi));
    V_TG = TG_trap.*abs(psi).^2; % Not counting VInt twice (only added to BEC energy)

    energy_BEC(end) = real(sum(K_BEC + V_BEC))*dx;
    energy_TG(end) = real(sum(K_TG + V_TG,'all'))*dx^3;
    energy(end) = energy_BEC(end) + energy_TG(end);
    Natom(1,end) = sum(abs(wf).^2)*dx;
    Natom(2,end) = sum(abs(psi).^2,'all')*dx^3;
    time(end) = dt*i;

    % Only keep real part of wave function to save storage space 
    if sum(imag(psi),'all') < 1e-6
        psi = real(psi);
    end

    data = v2struct(wf,psi,rho,energy,energy_BEC,energy_TG,time,Natom);
end