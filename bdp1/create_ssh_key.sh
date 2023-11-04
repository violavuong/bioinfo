# --------------------------
# in your linux shell
# always (!) use the public key (.pub) 
# --------------------------

## <key_filename>: the name of your SSH key
## <user>: the username of the VM (root cannot be used)

ssh-keygen -t rsa -f <key_filename> -C <user> -b 2048
