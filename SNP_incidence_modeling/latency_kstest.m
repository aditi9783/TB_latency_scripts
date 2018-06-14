  % Creates a 10x10 Magic square 
  % To do two-tailed 2sample ks test for 2d data 
  % This file should be run from same folder that has the kstest_2s_2d.m file
  datapath = ''; #write your data path here
  truefname = fullfile(datapath,'nmut_monthwise_latencypairs.csv');
  truedata = csvread(truefname);
  models = ["Constant", "Linear", "Expo", "RapidExpo_2", "RapidExpo_10"];
  outf = fopen('pvalues_ks2s2d_monthwise_latencymodels.txt', 'w');
  for m = [1:1:length(models)]
    samp_pvals = []; % list of p-values from comparing true data to each sampled data from a given latency-model
    for i = [0:1:999] % there are 100 randomly sampled distributions 0-99
      %sampf = fullfile(datapath,models(m),[models(m),'_poisson_26samples_',int2str(i),'.txt']);  
      sampf = datapath+models(m)+"/"+models(m)+"_poisson_25samples_"+int2str(i)+".txt";  
      sampdata = csvread(sampf);
      [hval, pval, ksval] = kstest_2s_2d(truedata, sampdata);
      samp_pvals = [samp_pvals,pval];
    end
    fprintf(outf,'%s,', models(m));
    fprintf(outf,'%d,',samp_pvals(1:length(samp_pvals)-1));
    fprintf(outf,'%d\n',samp_pvals(end));
  end
  fclose(outf);
  exit
