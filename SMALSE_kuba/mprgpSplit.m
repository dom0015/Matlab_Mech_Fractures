function [gf,gc,gr] = mprgpSplit(x,g,mprgpctx)
  bounddiff = x - mprgpctx.lb;
  activeset = (abs(bounddiff) <= mprgpctx.settol);
  freeset = ~activeset;
  gf = freeset.*g; % free
  gc = min(activeset.*g,0); % chopped
  gr = min(bounddiff/mprgpctx.abar,gf); % free reduced
end
