# -----------------------
# on the server/master
# -----------------------

## dependencies
yum install nfs-utils rpcbind
systemctl enable nfs-server
systemctl enable rpcbind
systemctl enable nfs-lock
systemctl enable nfs-idmap
systemctl start rpcbind
systemctl start nfs-server
systemctl start nfs-lock
systemctl start nfs-idmap
systemctl status nfs

vim /etc/exports #adding: /data  <client_internal_IP>(rw,sync,no_wdelay)
exportfs -r #checking for any errors (if nothing is printed, everything went well)
exportfs #printing the address

# ----------------------
# on the client/slave
# ----------------------

yum install nfs-utils
mkdir /data
mount -t nfs -o ro,nosuid <internal_IP>:/data /data
ll /data/ #checking the dir content as list
umount /data

vim /etc/fstab 

## adding the following lines: 
## UUID=1fc4211c-2271-43b7-92f0-3fbdbe1c2f2f /                       xfs     defaults        0 0
## UUID=DA08-90CB          /boot/efi               vfat    defaults,uid=0,gid=0,umask=0077,shortname=winnt 0 0
## <server_private_IP>:/data /data   nfs defaults        0 0

mount -a
df #printing all file system
