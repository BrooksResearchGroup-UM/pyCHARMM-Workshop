# Using Pre-built Environment with jupyter-lab on Gollum

If you have access to ***gollum***, you can use jupyter lab to run the tutorial examples to start a jupyter-server on ***gollum/satyr***.

**NOTE:** Since each jupyter server needs to bind to a dedicated port on gollum, you will need to choose a unique port number for your own jupyter server. It should be a dynamic port in the range from 49152 to 65535. (Some port numbers are already taken e.g. 65432)

To use a compute node for GPU computation with jupyter lab, you will need to login to your favorite compute node on gollum. It is always better to check whether the node is busy by checking the status with

```bash
# to check partition
sinfo -N -o '%10N %10T %18E %.8e %.8O %.15C' -p gpu    
# To check a single node
sinfo -N -o '%10N %10T %18E %.8e %.8O %.15C' -n gollum152 
```

For this example, I use ***gollum152***.

## On Gollum

```shell
# login to gollum 152 from head node
ssh gollum152 
intra_ip=`ifconfig | grep "inet 192.*" | awk '{print $2;}'`
module load pycharmm/0.4
jupyter lab --no-browser --ip $intra_ip --port xxxxx 
```

## On Local Machine

```shell
ssh -N -f -L localhost:xxxxx:gollum152:xxxxx gollum
```

Where `xxxxx` is the port number of your choice, and you should also substitute ***gollum152*** with the compute node of your choice (But of course you can also just use gollum152, and we can share the same compute node). You should be able to see an url like `http://127.0.0.1:xxxxx/lab` from the jupyter lab stdout. Copy and paste the url in a web browser, and you should be connected to jupyter lab on gollum with a pycharmm kernel (Do not use the url starting with 192). To check if you are on the correct node as well as the current GPU usage, you can type in a jupyter cell

```jupyter
%%bash
hostname
nvidia-smi
```
