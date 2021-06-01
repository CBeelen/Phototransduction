% script to submit job to cluster
sched = parcluster('CARL');
sched.AdditionalProperties.memory='10G';
sched.AdditionalProperties.runtime='04:00:00';
sched.AdditionalProperties.diskspace='50G';
job_small_RKKO_1 = batch(sched, 'cluster_small_RKKO', 'Pool', 1, 'AttachedFiles', {'Beelen_small_2PDE_7_stoch_RKKO.txtbc', 'IQMstochsim2.m'});