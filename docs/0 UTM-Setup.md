
# UTM Setup
The following instructions outline the process for setting up a pulsar environment using a Ubuntu image on a MacOS computer. These instructions are based on UTM Version 4.3.5 (87). 

## Image Setup
After installing UTM, follow the instructions in the [UTM docs](https://docs.getutm.app/guides/ubuntu/) for the setup of the Ubuntu image. For these instructions, the `UTM Server for ARM` image was selected (Ubuntu 22.04).

During the image installer, create and save a hostname and a username. The ones used for the following instructions are the following:
```
Host: csirolinux
User: csirouser
```

## Set up SSH access
On your local computer, use `ssh-keygen` to create a ssh key in the `~/.ssh` folder. The key generation process will require that a unique key name and optional passcode is entered. The private keyname for the following instructions is `csirolinux_ed25519`. The hostname is an IP addressed that is accessed by running `ip a` on the remote image.

```
Host csirolinux
	User csirouser
	Hostname X.X.X.X 
	PubKeyAuthentication yes
	IdentityFile ~/.ssh/csirolinux_ed25519
```

On the remote, add the public key to the file `~/.ssh/authorized_keys` file. The public can be transferred using the Shared folder configuration in UTM. This should allow you to ssh into the image and set up the remote development environment.

## Headless Configuration
Once the above is completed, the image can be started headlessly by going to the image settings, right clicking on the `Display` section and selecting `Remove`. This will prevent the display screen from popping up when the remote image is started. 

## VSCode Remote SSH
A remote development environment can be set up locally in VSCode by downloading the extension `Remote - SSH` and selecting the remote image. If the local `.ssh/config` is correctly, the `csirolinux` Host should come up automatically. This will allow you to easily configure and develop on the remoute Ubuntu image.  