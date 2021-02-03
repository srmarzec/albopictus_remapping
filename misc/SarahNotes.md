# Sarah's Notes
Things Sarah has tried, wants to keep track of

## Downloading the right version of trimmomatic
We will be using version 0.39
```
# enter directory you want to store trimmomatic in
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

# you can test everything worked with this next comand giving you the version number of 0.39
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar -version
```
