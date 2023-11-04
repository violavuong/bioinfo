# ---------------------------------------
# submitting a batch job with HTCondor
# ---------------------------------------

## do NOT work in the root
cp -r /data/BDP1_2022/condor/ . #copying the data from the condor dir
cd condor

ls #myexec.sh is the executable file 
./myexec.sh world.txt 0 #executing the file

condor_submit first_batch.job #condor submission
condor_q #checking the status of the job
condor_q -better-analyze  
condor_q -better-analyze  <job_ID>
condor_history <job_ID> 

ls -ltr #checking all files are outputted
cat condor.out
cat condor.error
cat condor.log
