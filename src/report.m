resnormi = resnorm;

% get residual norm
resnorm = norm(dU(:),2)./(1e-16+norm(W(:)+WBG(:),2)) ...
        + norm(dW(:),2)./(1e-16+norm(W(:)+WBG(:),2)) ...
        + norm(dP(:),2)./(1e-16+norm(P(:),2));

if it<=max(100,nup) || resnorm>resnorm0; resnorm0 = resnorm + 1e-32; end  % reset reference residual

% report iterations
if     it >=  0  && it <  10
    fprintf(1,'    ---  it =      %d;   abs res = %4.4e;   rel res = %4.4e \n',it,resnorm,resnorm/resnorm0);
elseif it >= 10  && it < 100
    fprintf(1,'    ---  it =     %d;   abs res = %4.4e;   rel res = %4.4e \n',it,resnorm,resnorm/resnorm0);
elseif it >= 100 && it < 1000
    fprintf(1,'    ---  it =    %d;   abs res = %4.4e;   rel res = %4.4e \n',it,resnorm,resnorm/resnorm0);
elseif it >= 1000 && it < 10000
    fprintf(1,'    ---  it =   %d;   abs res = %4.4e;   rel res = %4.4e \n',it,resnorm,resnorm/resnorm0);
elseif it >= 10000
    fprintf(1,'    ---  it =  %d;   abs res = %4.4e;   rel res = %4.4e; \n',it,resnorm,resnorm/resnorm0);
end 

if any(isnan([U(:);W(:);P(:);f(:)]))
    error('!!! Solver diverged !!!')
end

% plot convergence of outer iterations
if plot_cv
    figure(100); if it<=max(100,nup); clf; else; hold on; end
    plot(it,log10(resnorm),'r.','MarkerSize',15); box on; axis tight;
    drawnow;
end