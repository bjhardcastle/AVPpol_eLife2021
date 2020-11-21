function superPBTrial(objarray)

for oidx = 1:length(objarray)
   if ismember(4,objarray(oidx).Exps) && ...
           ( ~isfield('polOffBetweenTrials',objarray(oidx).MIP.pSet(4)) || objarray(oidx).MIP.pSet(4).polOffBetweenTrials == 0 )
      
       getFrames(objarray(oidx).MIP)
       trialidx = objarray(oidx).MIP.TrialStartFrame( objarray(oidx).MIP.TrialSeqNum == 4 );
       fm = [];
       fa = [];
      for tidx = 1:length(trialidx)
       t0 = objarray(oidx).MIP.TrialStartFrame(tidx);
       t1 = objarray(oidx).MIP.TrialEndFrame(tidx);
          fm(:,:,tidx) = max(objarray(oidx).MIP.Frames(:,:,t0:t1),[],3);
          fa(:,:,tidx) = mean(objarray(oidx).MIP.Frames(:,:,t0:t1),3);
          
      end
   end
   
   figure
   for n = 1:size(fm,3)
       ax(n) = subplot(size(fm,3)/6,6,n);
       imagesc(fm(:,:,n)  - mean(fa,3))
       ax(n).CLim = [0 max(fa(:))];
       axis off
       
   end
   tightfig(gcf)
end

end