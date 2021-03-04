resnormi = resnorm;

% get residual norm
resnorm = norm(dU(:),2)./(1e-16+norm(W(:)+WBG(:),2)) ...
        + norm(dW(:),2)./(1e-16+norm(W(:)+WBG(:),2)) ...
        + norm(dP(:),2)./(1e-16+norm(P(:),2)) ...
        + norm(df(:),2)./(1e-16+norm(f(:),2));

if it<=max(100,nup) || resnorm>resnorm0; resnorm0 = resnorm + 1e-32; end  % reset reference residual
 
% get total mass error
if bnchmrk
    err = norm((U(:)-U_mms(:)),2)./(1e-16+norm(U_mms(:),2)) ...
        + norm((W(:)-W_mms(:)),2)./(1e-16+norm(W_mms(:),2)) ...
        + norm((P(:)-P_mms(:)),2)./(1e-16+norm(P_mms(:),2));
else
    if step<=1 && it<=max(100,nup); fmass0 = sum(f(:)); end
    fmass = sum(f(:));
    err   = abs(fmass-fmass0)/fmass0;
end

% report iterations
if     it >=  0  && it <  10
    fprintf(1,'    ---  it =      %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,err);
elseif it >= 10  && it < 100
    fprintf(1,'    ---  it =     %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,err);
elseif it >= 100 && it < 1000
    fprintf(1,'    ---  it =    %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,err);
elseif it >= 1000 && it < 10000
    fprintf(1,'    ---  it =   %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,err);
elseif it >= 10000
    fprintf(1,'    ---  it =  %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,err);
end 

% plot convergence of outer iterations
if plot_cv
    figure(100); if it<=max(100,nup); clf; else; hold on; end
    plot(it,log10(resnorm),'r.','MarkerSize',15); box on; axis tight;
    
%     figure(101);
%     subplot(1,3,1); imagesc(-dWi); axis equal tight; box on; colorbar; colormap(ocean);
%     subplot(1,3,2); imagesc( dUi); axis equal tight; box on; colorbar; colormap(ocean);
%     subplot(1,3,3); imagesc( dPi); axis equal tight; box on; colorbar; colormap(ocean); 
    
    drawnow;
end