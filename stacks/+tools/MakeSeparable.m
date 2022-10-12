function [Col,Row,Scalar,Model,Residual,SingValues,Constant] = MakeSeparable(MatrixIn,ShowGraphics,GuessSign)
% MakeSeparable finds the best separable approximation to a matrix
% 
% [Row,Col,Scalar,Model,Residual] = MakeSeparable(MatrixIn)
% uses singular value decomposition to obtain the best separable model of
% MatrixIn.  It returns the best row, the best column,the best scaling
% factor, the best model (Row*Col*bestscale) and the residual
% between the best model and MatrixIn
%
% [Row,Col,Scalar,Model,Residual,SingValues] = MakeSeparable() returns all
% the singular values. 
%
% [Row,Col,Scalar,Model,Residual,SingValues,Constant] = MakeSeparable()
% also finds the best optimal additive constant to remove from the matrix.
%
% [] = MakeSeparable(MatrixIn, ShowGraphics) lets you specify whether to
% show graphics (1) or not (DEFAULT: 0).
% 
% [] = MakeSeparable(MatrixIn, ShowGraphics, GuessSign) 
% if you set GuessSign to 1 it does it best to
% give you Col and Row that have the right signs...
%
%Originally called MakeBestSeparableModel
%
% Example:
% mycol = normrnd(0,1,200,1);
% myrow = exp(-(-10:10).^2/5)+1;
% data = mycol*myrow;
% [Col,Row,Scalar,Model,Residual] = MakeSeparable(data,1);
%   
% 2006-03 RAF wrote it.
% 2006-05 RAF added the eigenspectrum as an output
% 2007-06 MC commented out "axes equal"
% 2009-05 MC imposed removing the grand mean first
% 2009-07 SK made it return the singular values
% 2010-05 MC removed the grand mean subtraction (user does it if they want)
% 2013-12 MC added option to try to fix the sign
% 2015-05 MC added option to find the best additive constant

%% parse the input parameters

if nargin<3 
    GuessSign = [];
end

if nargin<2
    ShowGraphics = [];
end

if nargin<1
    error('MakeSeparable requires a matrix as an input.');
end

if isempty(ShowGraphics)
    ShowGraphics = 0;
end

if isempty(GuessSign)
    GuessSign = 0;
end

%% find the best constant

if nargout>=7
    Constant = fminsearch(@(alpha) SeparableError(MatrixIn,alpha),0);
else
    Constant = 0;
end

%%

[U,S,V] = svd(MatrixIn - Constant,'econ');
Row = U(:,1);
Col = V(:,1)';
Scalar = S(1,1);
Model = Row*Col*Scalar + Constant;
Residual = MatrixIn - Model;

SingValues = nan(1, size(S,1));
for i = 1:size(S,1)
    SingValues(i) = S(i,i);
end

if GuessSign
    nR = size(MatrixIn,1);
    Projections = zeros(nR,1);
    for iR = 1:nR
        Projections(iR) = MatrixIn(iR,:)'\Col';
    end
    RespSign = sign(mean(Projections));
    Row = Row * RespSign;
    Col = Col * RespSign;
end

%%

if ShowGraphics
    %1st principal Component figure
    lims(1,:) = [min(MatrixIn(:)) max(MatrixIn(:))];
    lims(2,:) = [min(Model(:)) max(Model(:))];
    lims(3,:) = [min(Residual(:)) max(Residual(:))];
    clims = [min(lims(:,1)) max(lims(:,2))];
    figure('Color',[1 1 1]);
    subplot(1,3,1);
    imagesc(MatrixIn,clims);
    title('Original Matrix');axis tight; % axis equal;
    subplot(1,3,2);
    imagesc(Model,clims);
    title('Best Separable Model');axis tight; % axis equal;
    subplot(1,3,3);
    imagesc(Residual,clims);
    title('Residual');axis tight; % axis equal;
    colormap bone;
    
    %Eigenspectrum figure
    PC = linspace(1,size(S,1),size(S,1));
    for il = 1:size(S,1)
        VarAccount(il) = S(il,il)./sum(S(:));
    end
    CumAccount = cumsum(VarAccount);
    figure('Color',[1 1 1]);
    stem(PC,VarAccount);hold on;plot(PC,CumAccount,'r-');
    title('Eigenspectrum');
    xlabel('Principal Component #');
    ylabel('% Variance accounted');
    legend('Individual','Cumulative');
end

return

%% Code to test it:

ShowGraphics = 1;
FieldMat = repmat(linspace(0,1,50),50,1);
MatrixIn = randn(50,50)+10.*sin(2.*pi.*FieldMat);
[Col,Row,Scalar,Model,Residual] = MakeSeparable(MatrixIn,ShowGraphics);

