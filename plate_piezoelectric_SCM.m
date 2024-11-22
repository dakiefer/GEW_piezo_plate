% Dispersion calculation in a piezoelectric plate.
% 
% _↑z_______________________________________ top surface
%  →x   --> k  guided wave    c,e,epsr,rho   plate
% __________________________________________ bottom surface
% 
% This is an extension to https://github.com/dakiefer/GEW_dispersion_script 
% 
% The computation is performed in terms of the mechanical displacements u and
% the electrostatic potential Phi. At angular frequency w, the governing
% equations are
% 
% div(T) + w^2 rho u = 0,    (balance of linear momentum)
% div(D) = 0,                (charge-free material)
% 
% with mass-density rho. The stress T and electric flux D are given by 
% 
% T = c:grad(u) + grad(Phi).e , 
% D = e:grad(u) - E.grad(Phi) . 
% 
% with c:     4th-order stiffness tensor
%      E:     2nd-order permittivity tensor
%      e:     3rd-order piezoelectric constitutive tensor 
% 
% Depends on DMSUITE's chebdif.m by Weideman and Reddy [2] (included): 
% https://mathworks.com/matlabcentral/fileexchange/29-dmsuite
% 
% see also:
% [1] D. A. Kiefer, "Elastodynamic quasi-guided waves for transit-time ultrasonic flow
% metering", ser. FAU Forschungen, Reihe B, Medizin, Naturwissenschaft, Technik, 
% vol. 42. Erlangen: FAU University Press, 2022, doi: 10.25593/978-3-96147-550-6
% [2] J. A. Weideman and S. C. Reddy, “A MATLAB Differentiation Matrix Suite,”
% ACM Trans. Math. Softw. 26(4), 465–519, Dec. 2000, doi: 10.1145/365723.365727.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
%        Clemens Grünsteidl, RECENDT, Linz, Austria

mat = jsondecode(fileread('lithium_niobate_128Ycut.json')); % load material data
h = 0.51e-3;     % thickness in m
N = 20;          % discretization: polynomial order of interpolants
BCtop = 'open';  % "open" | "shorted". Electrical boundary condition at top surface
BCbot = 'open';  % "open" | "shorted". Electrical boundary condition at bottom surface
theta = pi/2;    % propagation direction in the plane of the plate (from X)
k = linspace(1e-2, 15, 120).'/h; % wavenumbers on which to solve

% convert material parameters from Voigt notated matrices to tensors
c = voigt2tensor4(mat.C); % stiffness in Pa = V/m*A*s/m^2, size [3x3x3x3]
rho = mat.rho;            % mass density in kg/m^3 = A*s^3*V/m^5
e = voigt2tensor3(mat.e); % charge constants in A*s/m^2, size [3x3x3]
eps0 = 8.8541878188e-12;  % permittivity of free space in A*s/(V*m)
E = mat.epsr*eps0;        % permittivity in A*s/(V*m), size [3x3]

% rotate to propagation direction:
c = rotate(c,-theta);
e = rotate(e,-theta);
E = rotate(E,-theta);

% normalize parameters:
c0 = c(1,2,1,2); rho0 = rho; E0 = E(1,1); e0 = sqrt(E0*c0); fh0 = sqrt(c0/rho0); % choose
cn = c/c0; rhon = rho/rho0; En = E/E0; en = e/e0; eps0n = eps0/E0;  % normalize

%% relevant material matrices (projection to propagation coordinates): 
udof = [1 2 3]; % displacement degrees of freedom (for Lamb-polarization choose [1 3])
cxx = squeeze(cn(1,udof,udof,1));
cxz = squeeze(cn(1,udof,udof,3));
czx = squeeze(cn(3,udof,udof,1));
czz = squeeze(cn(3,udof,udof,3));
exx = squeeze(en(1,1,udof));
exz = squeeze(en(1,3,udof));
ezx = squeeze(en(3,1,udof));
ezz = squeeze(en(3,3,udof));
Exx = En(1,1);
Exz = En(1,3);
Ezx = En(3,1);
Ezz = En(3,3);

%% assemble mechanical and electrical tensors
A  = [cxx,    exx;      % ~ (ik)^2
      exx.'  -Exx];
B  = [czx,    exz;      % ~ (ik)^1
      ezx.'  -Ezx];
C  = [czz,    ezz;      % ~ (ik)^0
      ezz.'  -Ezz];
MM  = blkdiag(rhon*eye(size(cxx)), 0);       % ~ w^2 (mechanics only, singular)
V   = zeros(size(B)); V(end,end) = 1i*eps0n; % For imposing continuity of D across surfaces

%% discretization 
[~, Dy_dash] = chebdif(N, 2); % create differentiation matrices
D1 = -2*Dy_dash(:,:,1);       % differentiation on unit domain
D2 = 4*Dy_dash(:,:,2);        % second order derivative
Id = eye(size(D1));           % identity matrix for discretization

% discrete wave operators:
L2 = kron(A, Id); 
L1 = kron(B + B.', D1); 
L0 = kron(C, D2); 
M =  kron(MM, Id);

% discrete boundary operators for inhomogeneous Neumann BC (1: bottom, N: top):
B1 = kron(B, Id([1 N],:)) + kron(V, [-1; 1].*Id([1 N],:)); % different sign for top/bottom in V
B0 = kron(C, D1([1, N],:));

%% incorporate boundary conditions
% incorporate Neumann BCs (continuity of electric flux D across the surfaces of the plate):
Nudof = length(udof);
dofBC = [(0:Nudof)*N+1; (1:Nudof+1)*N]; % equation numbers where to impose BCs
dofPhiBot = dofBC(1,end); dofPhiTop = dofBC(2,end); % equation numbers for electrical BCs
L2(dofBC,:) = 0; L1(dofBC,:) = B1; L0(dofBC,:) = B0; M(dofBC,:) = 0;

% metalized surfaces (Dirichlet BCs for the electric potential):
nfix = []; % list of degrees of freedom (DOF) where to apply shorted conditions
if string(BCtop) == "shorted", nfix = [nfix, dofPhiTop]; end % add Phi-DOF at top to list
if string(BCbot) == "shorted", nfix = [nfix, dofPhiBot]; end % add Phi-DOF at bottom to list
L2(nfix,:) = []; L1(nfix,:) = []; L0(nfix,:) = []; M(nfix,:) = []; 
L2(:,nfix) = []; L1(:,nfix) = []; L0(:,nfix) = []; M(:,nfix) = [];

%% solve for real frequencies:
kh = k(:)*h; % wavenumber-thickness will be imposed to compute the corresponding frequency-thicknesses
whn = nan(length(kh), size(M,1));      % allocate array for normalized angular frequencies
Psi =   zeros([length(kh), size(M)]);  % allocate for eigenvectors
tic
for ii = 1:length(kh)
    [Psii, wh2i] = eig(kh(ii)^2*L2 - 1i*kh(ii)*L1 - L0, M, 'vector'); % solve
    whi = real(sqrt(wh2i));
    whi(abs(whi) == 0) = nan; whi(isinf(whi)) = nan; % remove "spurious" solutions
    [whi, ind] = sort(whi);  % sort ascending
    Psii = Psii(:,ind);      % sort like the frequencies 
    whn(ii,:) = whi;         % save frequencies/eigenvalues
    Psi(ii,:,:) = Psii;      % save eigenvectors
end
w = whn*fh0/h; % convert to rad/s
chron = toc; fprintf('nK: %d, nF: %d, elapsed time: %g, time per point: %g. ms\n', size(w, 1), size(w, 2), chron, chron/length(w(:))*1e3);

%% plot 
figure; hold on;
kk = k.*ones(size(w)); % expand to same size as w
plot(kk/1e3, w/2/pi/1e6, 'k-');
xlim([0, 25]), ylim([0, 30])
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz'),






%% helper functions
function c = voigt2tensor4(C)
    % VOIGT2TENSOR returns the 4th order stiffness tensor c corresponding to the Voigt
    % notated matrix C. This function is from https://github.com/dakiefer/GEWtool
    % 
    % see also: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
    % (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
    %
    % 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    validateattributes(C,{'numeric'},{'size',[6 6]});
    if ~issymmetric(C), error('Voigt matrix C must be symmetric.'); end
    c = zeros(3, 3, 3, 3);
    ij = [ {1},{1}; {2},{2}; {3},{3}; {2},{3}; {3},{1}; {1},{2}]; % voigt -> tensor indices
    ji = [ {1},{1}; {2},{2}; {3},{3}; {3},{2}; {1},{3}; {2},{1}]; % voigt -> tensor indices
    for I = 1:6
        for J = I:6 % start from I (symmetry in C) -> explicitly assign major symmetry below
            % minor symmetries: 
            c(ij{I,:}, ij{J,:}) = C(I, J); % cijkl
            c(ij{I,:}, ji{J,:}) = C(I, J); % cijlk
            c(ji{I,:}, ij{J,:}) = C(I, J); % cjikl
            c(ji{I,:}, ji{J,:}) = C(I, J); % cjilk
            % major symmetry:
            c(ij{J,:}, ij{I,:}) = C(I, J); % cklij
            c(ij{J,:}, ji{I,:}) = C(I, J); % cklji
            c(ji{J,:}, ij{I,:}) = C(I, J); % clkij
            c(ji{J,:}, ji{I,:}) = C(I, J); % clkji
        end
    end
end % function

function e = voigt2tensor3(E)
    % VOIGT2TENSOR3 - returns the 3rd order piezoelectric constitutive tensor e
    % corresponding to the Voigt notated matrix E. This function is from 
    % https://github.com/dakiefer/GEWtool
    % 
    % see also: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
    % (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
    %
    % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    validateattributes(E,{'numeric'},{'size',[3 6]});
    e = zeros(3, 3, 3);
    jk = [ {1},{1}; {2},{2}; {3},{3}; {2},{3}; {3},{1}; {1},{2}]; % voigt -> tensor indices
    kj = [ {1},{1}; {2},{2}; {3},{3}; {3},{2}; {1},{3}; {2},{1}]; % voigt -> tensor indices
    % assign e_i{jk} from E_i{J}, where J in {1..6} -> {jk} in {11,22,33,23,32,13,31,12,21} :
    for i = 1:3
        for J = 1:6
            e(i, jk{J,:}) = E(i, J);
            e(i, kj{J,:}) = E(i, J);
        end
    end
end % function

function B = rotate(A,theta)
    % rotate - rotate nth-order tensor A theta radian about ez (plate normal). 
    % This function is adapted from https://github.com/dakiefer/GEWtool
    %
    % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    R = [cos(theta), sin(theta), 0; 
        -sin(theta), cos(theta), 0; 
            0,         0,        1];
    sA = size(A);
    if iscolumn(A), sA = sA(1); end % correct weird matlab size
    n = length(sA); % order of tensor
    for dim = 1:n % transform each of the basis vectors of the tensor A one by one
        A = permute(A, [1:dim-1, n+1, dim:n]); % add one dimension before "dim"
        RR = shiftdim(R, -(dim-1)); % for contraction of R with the appropriate dimension of A
        A = squeeze(sum(RR.*A, dim+1)); % note: permute has shifted "dim" one to the right
    end
    B = A;
end
