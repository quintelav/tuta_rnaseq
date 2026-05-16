 1  sudo apt update && sudo apt upgrade -y
    2  sudo apt install -y wget curl git build-essential
    3  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    4  bash Miniconda3-latest-Linux-x86_64.sh
    5  source ~/.bashrc
    6  conda --version
    7  conda 
    8  conda --version
    9  uname -a
   10  cd
   11  cd ~
   12  bash Miniconda3-latest-Linux-x86_64.sh
   13  source ~/.bashrc
   14  conda --version
   15  ls ~
   16  ~/miniconda3/bin/conda --version
   17  ~/miniconda3/bin/conda init bash
   18  source ~/.bashrc
   19  conda --version
   20  echo $PATH
   21  conda create -n rnaseq python=3.10
   22  conda install -c conda-forge -c bioconda fastqc
   23  conda create -n rnaseq python=3.10
   24  conda activate rnaseq
   25  conda config --add channels defaults
   26  conda config --add channels bioconda
   27  conda config --add channels conda-forge
   28  conda config --set channel_priority strict
   29  conda install fastqc
   30  conda install hisat2
   31  conda install samtools
   32  pip install pandas numpy matplotlib seaborn jupyter biopython
   33  mkdir -p ~/projects/tuta_rnaseq/{data_raw,data_processed,scripts,results,figures,metadata}
   34  cd ~/projects/tuta_rnaseq
   35  code .
   36  fastqc --version
   37  samtools --version
   38  conda install sra-tools
   39  prefetch
   40  fasterq-dump
   41  du -sh raw_data
   42  ps aux | grep -E 'prefetch|fasterq|gzip'
   43  pkill -9 prefetch
   44  pkill -9 fasterq-dump
   45  pkill -9 gzip
   46  du -sh raw_data
   47  ls -lh raw_data
   48  rm raw_data/SRR35985984*
   49  rm raw_data/SRR35985987_2.fastq.gz
   50  ls -lh raw_data
   51  pkill -9 prefetch
   52  pkill -9 fasterq-dump
   53  pkill -9 gzip
   54  ps aux | grep -E 'prefetch|fasterq|gzip'
   55  ls -lh raw_data
   56  ls -lh raw_data/SRR35985987_1.fastq.gz
   57  ls -lh raw_data/SRR35985988_1.fastq.gz
   58  ls -lh raw_data/SRR35985988_2.fastq.gz
   59  ls -lh raw_data
   60  rm raw_data/SRR35985987_1.fastq.gz
   61  ls -lh raw_data
   62  rm raw_data/SRR35985988_1.fastq.gz
   63  rm raw_data/SRR35985988_2.fastq.gz
   64  ls -lh raw_data
   65  grep SRR35985987 metadata/runinfo.csv
   66  rm raw_data/SRR35985987*
   67  ls -lh raw_data
   68  ps aux | grep -E 'prefetch|fasterq|gzip|salmon'
   69  ls scripts
   70  pkill -9 fasterq-dump
   71  ps aux | grep -E 'prefetch|fasterq|gzip|salmon'
   72  pkill -9 gzip
   73  ps aux | grep -E 'prefetch|fasterq|gzip|salmon'
   74  pkill -9 prefetch
   75  pkill -9 fasterq-dump
   76  pkill -9 gzip
   77  ps aux | grep -E 'prefetch|fasterq|gzip|salmon'
   78  ls -lh raw_data
   79  rm raw_data/SRR35985985_1.fastq
   80  rm raw_data/SRR35985985_1.fastq.gz
   81  rm raw_data/SRR35985985_2.fastq
   82  ls -lh raw_data
   83  ps aux | grep -E 'prefetch|fasterq|gzip|salmon'
   84  tree -L 2
   85  sudo apt install tree
   86  tree -L 2
   87  ls scripts
   88  pwd
   89  ls
   90  cd scripts
   91  ls
   92  cd ~/projects/tuta_rnaseq
   93  cd scripts
   94  nano README.txt
   95  history
   96  history > commands_history.txt
   97  ls
   98  touch scripts/download_sra.sh
   99  touch scripts/run_fastqc.sh
  100  touch scripts/run_fastp.sh
  101  touch scripts/run_salmon.sh
  102  ls scripts
  103  mkdir -p scripts
  104  touch scripts/download_sra.sh
  105  touch scripts/run_fastqc.sh
  106  touch scripts/run_fastp.sh
  107  touch scripts/run_salmon.sh
  108  ls scripts
  109  nano pipeline_notes.md
  110  history
