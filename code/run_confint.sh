#!/bin/bash
#$ -N confint
#$ -o ../job_output/$JOB_NAME-$JOB_ID.log
#$ -j y

# epsa range
epsa_min=-3
epsa_gap=1
epsa_max=1
# epsy range
epsy_min=-8
epsy_gap=1
epsy_max=-8
# betau
betau_min=0
betau_gap=2
betau_max=4
# betav
betav_min=0
betav_gap=2
betav_max=0
# d
d_min=0
d_gap=2
d_max=0
# n
n_min=8
n_gap=1
n_max=8
# n
n_test_min=0
n_test_gap=2
n_test_max=0
# nsim
nsim=500
# beta0
beta0='1-3-5'
xmat='a'
max_iter=10000
jobid=$JOB_ID
scaleuvtest=1
power_beta=1
confint=1
confint_A=1
hulc=0
obs_ratio_min=1
obs_ratio_gap=.125
obs_ratio_max=1
batch_size=100

echo "job_id, nsim, epsa_min, epsa_gap, epsa_max, epsy_min, epsy_gap, epsy_max, beta0, betau_min, betau_gap, betau_max, betav_min, betav_gap, betav_max, d_min, d_gap, d_max, n_min, n_gap, n_max, n_test_min, n_test_gap, n_test_max, obs_ratio_min, obs_ratio_gap, obs_ratio_max, xmat, max_iter, scaleuvtest, power_beta, confint, hulc" >> "../output/${JOB_NAME}_${JOB_ID}.log"
echo "$jobid, $nsim, $epsa_min, $epsa_gap, $epsa_max, $epsy_min, $epsy_gap, $epsy_max, $beta0, $betau_min, $betau_gap, $betau_max, $betav_min, $betav_gap, $betav_max, $d_min, $d_gap, $d_max, $n_min, $n_gap, $n_max, $n_test_min, $n_test_gap, $n_test_max, $obs_ratio_min, $obs_ratio_gap, $obs_ratio_max, $xmat, $max_iter, $scaleuvtest, $power_beta, $confint, $hulc" >> "../output/${JOB_NAME}_${JOB_ID}.log"


for obs_ratio in $(seq $obs_ratio_min $obs_ratio_gap $obs_ratio_max)
do
for n_train in $(seq $n_min $n_gap $n_max)
do
	for n_test in $(seq $n_test_min $n_test_gap $n_test_max)
	do
		for d in $(seq $d_min $d_gap $d_max)
		do
			for epsa in $(seq $epsa_min $epsa_gap $epsa_max)
			do
				for epsy in $(seq $epsy_min $epsy_gap $epsy_max)
				do
					for betau in $(seq $betau_min $betau_gap $betau_max)
					do
						for betav in $(seq $betav_min $betav_gap $betav_max)
						do
						    taskname="${JOB_NAME}-$epsa-$epsy-$n_train"
						    echo $taskname
						   	qsub -q short.q -N $taskname -t 1-5 -o '../job_output/$JOB_NAME-$JOB_ID-$TASK_ID.log' -j y -b y "Rscript --vanilla confint.R --n_train $n_train --n_test $n_test --d=$d --beta0 $beta0 --betau $betau --betav $betav --epsa $epsa --epsy $epsy --nsim $nsim --ibatch=\${SGE_TASK_ID} --max_iter $max_iter --xmat $xmat --jobid $jobid --power_beta $power_beta --confint $confint --hulc $hulc --observed_ratio $obs_ratio --batch_size $batch_size"
						done
					done
				done
			done
		done
	done
done
done



