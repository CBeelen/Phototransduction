% script to submit job to cluster
sched = parcluster('CARL');
sched.AdditionalProperties.memory='10G';
sched.AdditionalProperties.runtime='04:00:00';
sched.AdditionalProperties.diskspace='50G';
job_small = batch(sched, 'cluster_rand', 'Pool', 1, 'AttachedFiles', {'Beelen_small_2PDE_7_stoch.txtbc', 'Beelen_small_2PDE_7_stoch_prec.txtbc', 'IQMstochsim2.m'});