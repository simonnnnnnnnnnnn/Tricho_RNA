# Tricho_RNA


**Objective:** Creation of an RNAseq pipeline with Nextflow for _Trichoplax adhaerens_

## Setup of Execution Environment

The Pipelines were executed on a virtual machine running Ubuntu 22.04, the VM itself ran on SimpleVMs Servers.

### Requirements

 - In order to run a Virtual Machine on the SimpleVM servers the project (and the resources) must be approved by de.NBI
 - next the user must create a private Key, this key is needed to log into the virtual machine(s) using ssh.
 - the provate key has to be downloaded, the file must have sufficiently restricted access, I made root the owner of this file using the following command:
```b̀ash
sudo chown root:root /path/to/my/private/key
```
 - the ssh command to access the virtual machine is displayed in [simplevm](https://simplevm.denbi.de/portal/webapp/#/profile) after configuring the virtual machine

### Configuring the Virtual Machine

The configuration of the virtual machine was entirely done on SimpleVMs Portal, I chose the "Flavor" de.NBI medium.
This Flavor came with the following specs:
 - 14VCPUs
 - 32 GB RAM
 - 50 GB Root Disk space

I chose Ubuntu 22.04 de.NBI as the Image (2025-05-06)

### Configuration of the VM

In order to use the virtual machine I installed some software and ran some other commands. All commands dealing with the configuration of the VM are displayed below (ordered chronologically)

```b̀ash
mkdir project
nano test.txt
sudo apt install git
sudo apt install docker
sudo apt install openjdk-21-jdk
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv /home/ubuntu/nextflow /bin
nextflow info
```

The commands to install nextflow follow [Nextflows Instructions](https://training.nextflow.io/latest/envsetup/02_local/#download-nextflow) for the local environment (2025-05-07)
