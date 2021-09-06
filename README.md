# PYTHIA

[![Docker Pulls](https://img.shields.io/docker/pulls/dsimb/pythia.svg)](https://hub.docker.com/r/dsimb/pythia)

(C) Gabriel Cretin, Tatiana Galochkina, Alexandre de Brevern, Jean-Christophe Gelly  
https://doi.org/10.3390/ijms22168831

**Deep Learning Approach For Local Protein Conformation Prediction**


## Abstract

Protein Blocks (PBs) is a structural alphabet describing the local protein conformation with higher precision than classical secondary structures [1]. PBs correspond to 16 structural conformational states, which can be adopted by five consecutive amino acids. Encoding of complex protein structures (3D) in a PB sequence (1D) has already been successfully applied to protein structure alignment and protein structure prediction [2,3]. In the current study we developed a deep learning model for prediction of the protein local conformations in terms of PB directly from the amino acid sequence. Each amino acid is encoded by 58 physico-chemical properties [4] and a position-specific substitution matrix (PSSM) generated by PSI-BLAST. We performed a 10-fold cross-validation on a non-redundant dataset of 9638 protein chains from the Protein Databank. Prediction was performed using a deep residual inception‐inside‐inception neural architecture based on convolutional block attention modules. The developed model named PYTHIA (Predicting Any Conformation at High Accuracy) clearly outperforms the reference method for PB prediction LOCUSTRA [5]. The mean accuracy (Q16) equals 70% for PYTHIA and 61% for LOCUSTRA. Furthermore, PYTHIA outperforms LOCUSTRA on every PB class even for the smallest ones such as ‘g’ (MCC equal 0.209 for PYTHIA vs. 0.154 for LOCUSTRA) and ‘j’ (MCC equal 0.315 for PYTHIA vs. 0.223 for LOCUSTRA).

## Dependencies

### 1. Install Docker

**Linux**

Docker provides a convenient script to install itself (Linux: Ubuntu, Debian|Raspbian, CentOS|RHEL)
```term
$ curl -fsSL https://get.docker.com -o get-docker.sh
$ sudo sh get-docker.sh

# Add yourself to Docker group
$ sudo usermod -aG docker <user>
# This will reload your group assignments,
# avoiding the need to logout and log back in
$ su - $USER
```

**Mac**  

[Get Docker Desktop for Mac](https://docs.docker.com/docker-for-mac/install/)  

**Windows**  

[Get Docker Desktop for Windows](https://docs.docker.com/docker-for-windows/install/)  



### Download this repository  

```term
$ git clone https://github.com/DSIMB/PYTHIA
$ cd PYTHIA
```

### Database

HHblits requires a sequence database e.g. uniref30.  
If you don't have a sequence database follow the steps explained on HHblits [extensive user-guide](https://github.com/soedinglab/hh-suite/wiki#hh-suite-databases).  

**TL;DR**:  
1. Download the [latest release](http://wwwuser.gwdg.de/~compbiol/uniclust/current_release/) of database from the HHblits repository into an empty directory using command:  
`wget http://wwwuser.gwdg.de/~compbiol/uniclust/[date]]/UniRef[date]_hhsuite.tar.gz`

2. Extract the files using command:  
`tar xzvf UniRef[date]_hhsuite.tar.gz`


### Download the docker image  

You can download the latest build of Medusa docker image (recommended):  

```
$ docker pull dsimb/pythia
```

or build it yourself from the git repository:  

```
$ docker build -t dsimb/pythia .
```

### Run the Docker container  
  
**IMPORTANT**  
  
| Option | Adapt to your local paths (no trailing slash '/') |   | Do not modify |
|--------|---------------------------------------------------|---|---------------|
| -v     | /path/to/database                                 | : | /database     |
| -v     | $(pwd)                                            | : | /project      |


#### 1 - Set paths
```term
# Change paths and names accordingly
PATH_DATABASE=/path/to/database
# Path to PYTHIA git repository
PATH_PYTHIA=$(pwd)  
```

#### 2 - Run docker image  
```
$ docker run -it --rm \
    # Launch docker as user's id
    --user "$(id -u):$(id -g)" \  
    # Bind mount as read-only the database for HHblits (ex: /home/gabriel/UniRef)
    -v ${PATH_DATABASE}:/database:ro \  
    # Bind mount the project's folder
    -v ${PATH_PYTHIA}:/project \  
    # The name of the container we launch
    dsimb/pythia \  
    # Fasta file containing the target sequence (path relative to project)
    -i ./data/sequence.fasta \  
    # Name of the database prefix for HHblits (c.f --help for more details)
    # If database file is "uniclust30_2016_09_a3m.ffindex", set option to prefix "uniclust30_2016_09"
    -d uniclust30_2016_09 \
    # Use either the 'balanced' or 'global' accuracy model
    # Directory which will contain results (path relative to project)
    -o ./results  
```
Or, this one-liner, for convenience
```term
$ docker run -it --user "$(id -u):$(id -g)" -v ${PATH_DATABASE}:/database:ro -v ${PATH_PYTHIA}:/project dsimb/pythia -i ./data/sequence.fasta -d uniclust30_2016_09 -l balanced -o ./results
```

### 2. Install Conda environment

Alternatively to Docker, you can run PYTHIA inside a conda environment.  
To create the environment:  
```
$ conda env create -f environment.yaml
$ conda activate pythia
```  

Now you can launch PYTHIA:  
```
./launch_PYTHIA.sh -i ./data/T1027-D1.fasta -d /path/to/database/uniclust30_2016_09 -l balanced -o ./results -c 0 -m 0
```  

#### **ATTENTION**  

Please note that the option -d value depends on the way you launch PYTHIA:  
- With Docker: provide only the name of the database, ex: uniclust30_2016_09  
- With Conda: provide the absolute path to the database with the name, ex: /home/gabriel/databases/UNICLUST/uniclust30_2016_09  


### Example for running batch prediction  

If you have a multifasta file and you wish to run PYTHIA protein flexibility prediction for all of your sequences,
you can follow these few steps to split your multifasta file into single fasta files and run the docker container
for each of them individually.  

```term
$ ls
multi.fasta
$ cat multi.fasta
>SEQUENCE_1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
$ awk '/^>/{s=++d".fasta"} {print > s}' multi.fasta
1.fasta 2.fasta
$ for seq in ./*.fasta; do docker run -it --user "$(id -u):$(id -g)" -v ${PATH_DATABASE}:/database:ro -v ${PATH_PYTHIA}:/project dsimb/pythia -i ./$seq -d uniclust30_2016_09 -l balanced -o ./results -c 0 -o 0; done
```




### Help

```term
$ docker run dsimb/pythia
or
$ docker run dsimb/pythia --help
At least 4 parameters required.
Usage:
        -i    | --seq         (Required)     Path to input Fasta sequence file. The path is relative to the project folder.
        -l    | --model       (Required)     Use "balanced" or "global" accuracy model.
        -d    | --database    (Required)     Name of the database for HHBlits.
                                             Ex: UniRef30_2020_06
                                             If your path is '/path/to/Uniclust/UniRef30_2020_03_a3m.ffdata'
                                             please provide: -d UniRef30_2020_03
        -o    | --outdir      (Required)     Path to the output directory.
        -c    | --cpus        (Optionnal)    Number of CPUs to use (for HHblits). Default is 2. Set to 0 for all available memory.
        -m    | --memory      (Optionnal)    Maximum RAM to use in Gb (for HHblits). Default is 3. Set to 0 for all.
        -h    | --help        (Optionnal)    Brings up this help
```

### Example

```term
$ docker run -it --user "$(id -u):$(id -g)" -v /home/cretin/uniclust:/database:ro -v $(pwd):/project dsimb/pythia -i ./data/sequence.fasta -d uniclust30_2016_09 -l balanced -o ./results -c 0 -m 0
Run HHblits ... done
Run HHfilter ... done
Reformat ... done
Transform alignment to frequence matrix ... done
Create AAindex and one-hot encodings ... done
Translate encodings into vectors ... done
Merge vectors ... done
Run predictions ... done

Results can be found in ./results/PYTHIA.jy5J-20201208215144

$ ls ./results/PYTHIA.jy5J-2020120
logs  pythia_job_results.tar.gz  prediction
```

### Reference  

Gabriel Cretin, Tatiana Galochkina, Alexandre G. de Brevern, and Jean-Christophe Gelly. 2021. "PYTHIA: Deep Learning Approach for Local Protein Conformation Prediction" International Journal of Molecular Sciences 22, no. 16: 8831. https://doi.org/10.3390/ijms22168831

### Issues  

If you encounter any issue, do not hesitate to open an issue.  
