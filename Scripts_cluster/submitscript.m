% script to submit job to cluster
sched = parcluster('CARL');
sched.AdditionalProperties.memory='10G';
sched.AdditionalProperties.runtime='04:00:00';
sched.AdditionalProperties.diskspace='50G';
job_prec = batch(sched, 'cluster', 'Pool', 1, 'AttachedFiles', {'Beelen_small_2PDE_7_stoch.txtbc', 'IQMstochsim2.m'});