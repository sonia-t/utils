# Adventures in downloading data from basespace 

# only have to do this once 
# install bs (installed on Dragonet)
# see https://blog.basespace.illumina.com/2016/07/06/new-basespacecli-tools/
bs authenticate # you'll need a basespace account to login on the web

bs list runs  # this is useful ! 
# | 107307xx | 150213_M03428_0014_000000000-ACGxx | 150213XX  
bs list projects 

# try to cp this run using the runid 
bs cp -v  conf://default/Run/107307xx XX021115
# downloads a shit ton of other stuff, but no fastq
find  XX021115_bs_cp/ -name \*fastq\*  # nada
more XX021115_bs_cp/SampleSheet.csv  # it *IS* the right run

# try a third party solution: BaseSpaceRunDownloader
wget https://gist.githubusercontent.com/rlesca01/7ce2ca0c35c7ff97a215/raw/0eeaa8cc1b3eff00babf398a82a31f4b0946f5bb/BaseSpaceRunDownloader_v2a.py
AccessToken = 16fce877599b4bb69f43bb8a8dXXXX # get this from your config ~/.basespace/*cfg or thru website
# supposed to be able to run this with project name, not run id 
python BaseSpaceRunDownloader_v2a.py -r XX021115 -a  $AccessToken
# NameError: name 'ProjectDir' is not defined
# try again with project name 
python BaseSpaceRunDownloader_v2a.py -r 150213XX -a  $AccessToken  # this worked !
# see 150213XX-20295296/XX021115-Fx6-22302747/Data/Intensities/BaseCalls/XX021115-Fx6_S1_L001_R1_001.fastq.gz
# however, it gets ONLY the fastq 
find 150213XX-202952xx/ -name  \*csv # nada 


# BEST OPTION for flexibility: bs mount 
bs mount /shared/soniat/basespace_/
find  /shared/soniat/basespace_/Runs/150213XX/ -name \*fastq\*  # I can see them, should be easy to transfer. 
# confirm its the same run I tried and failed to download fastq from before 
diff /shared/soniat/basespace_/Runs/150213XX/Files/SampleSheet.csv XX021115_bs_cp/SampleSheet.csv
# now just cp or rsync the files from the mount to local 

# another third party option
git clone https://github.com/nh13/basespace-invaders
# ran into interrupted download issues ...WARNING:BaseSpacePy.model.MultipartFileTransfer:Task failed after too many retries

