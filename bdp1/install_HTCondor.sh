# ---------------------------
# installing HTCondor
# ---------------------------

# dependencies
yum install wget
wget https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
yum localinstall epel-release-latest-7.noarch.rpm
yum clean all

# repositories and packages
wget http://research.cs.wisc.edu/htcondor/yum/repo.d/htcondor-stable-rhel7.repo
cp htcondor-stable-rhel7.repo /etc/yum.repos.d/
yum install condor-all
wget http://htcondor.org/yum/RPM-GPG-KEY-HTCondor
rpm --import RPM-GPG-KEY-HTCondor

# basic configurations
vim /etc/condor/condor_config 

## ----------------------------------
## change the IP to your master IP
CONDOR_HOST = #private master IP
## add on the master
DAEMON_LIST = COLLECTOR, MASTER, NEGOTIATOR, STARTD, SCHEDD #master's jobs
## add these lines on the nodes
DAEMON_LIST = MASTER, STARTD
HOSTALLOW_READ = *
HOSTALLOW_WRITE = *
HOSTALLOW_ADMINISTRATOR = *	
## ---------------------------------

# enabling Condor
systemctl enable condor
systemctl start condor
systemctl status condor

condor_status # to check condor status
