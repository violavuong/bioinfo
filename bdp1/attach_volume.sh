# ----------------------
# in your linux shell 
# ----------------------

fdisk /dev/sdb #creating a new disk

## formatting the disk (with default settings)
## an interactive interface will open to create a partition in the disk, use the following commands
### p - checks if the disk was created with default settings; 
### n - creates the primary partition;
### w - creates a partition table. 

fdisk -l #list of the attached disk 

## creating a file system on top of the partition
mkfs.ext4 /dev/sdb1   #creating the ext4 filesystem
mkdir /data           #creating a mountpoint for the new filesystem
vim /etc/fstab        #editing the fstab file by adding: /dev/sdb1     /data  ext4 defaults 0 0

mount -a #mounting the file system
df -h #checking the mounted file system

