function [fuw,fuT,fwT,fuu,fww,fTT,f]=compute_spectra(un,wn,Tn,fs)
 Nwin=length(un);

 %------ Compute the co-spectra
 Nwind=floor(Nwin/16);
 window=hamming(Nwind);
 [fuu1,~]=cpsd(un,un,window);
 [fww1,~]=cpsd(wn,wn,window);
 [fTT1,~]=cpsd(Tn,Tn,window);
 [fuw1,~]=cpsd(un,wn,window);
 [fwT1,~]=cpsd(wn,Tn,window);
 [fuT1,fn]=cpsd(un,Tn,window);
 
 dfn=fn(2)-fn(1);
 f=(fn/max(fn))*fs*length(fn)/Nwind;
 df=f(2)-f(1);

 fuu=real(fuu1)*dfn/df;
 fww=real(fww1)*dfn/df;
 fTT=real(fTT1)*dfn/df;
 fuw=real(fuw1)*dfn/df;
 fwT=real(fwT1)*dfn/df;
 fuT=real(fuT1)*dfn/df;
 