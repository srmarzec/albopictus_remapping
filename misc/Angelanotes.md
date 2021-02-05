## Downloading trimmomatic into home directory
# enter directory you want to store trimmomatic in: /home/zz220
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

# you can test everything worked with this next comand giving you the version number of 0.39
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version

# downloading the SRA toolkit
Chose Ubuntu: wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
#could not get export path to work
