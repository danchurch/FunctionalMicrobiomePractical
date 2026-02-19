
## for macs, you will need to log in directly with SSH from your terminal

## to log in with a terminal:
ssh -X ubuntu@129.70.51.6 -p 31993 ## this is for my VM. Your port will be different. 

## you may need to tell ssh  where you put your private ssh key:

ssh -X -i /home/daniel/.ssh/./officeComp_denbi26 ubuntu@129.70.51.6 -p 31993 
## this is for my VM. Your port will be different. 

## for windows users, you can try this, and/or you can use MobaXterm


## let's practice some commands and getting around in a linux filetree.

pwd 

ls 

cd  /

## to get back to home drive:

cd 

## most of our work will actually not be in the "home" directory. 
## Look at the storage volume:

cd /vol/funmic



### and a note about filepaths in linux

## relative file path
cd 

cd /vol

cd funMicStorage

cd ..

pwd

## absolute file path:

cd /vol/funMicStorage

## what's the difference between these?
## (relative and absolute file paths)

## get back to your home directory, and back into the someFiles directory

## the file path is:
/home/ubuntu/someFiles

## how will you get there?:


## and what is in there?:


## the copy command is safe:

cp A12-05.fastq  anotherFile.fastq





## the mv (move) command is not so safe (but useful!)

mv anotherFile.fastq thisIsTheSameFile.fastq




mv thisIsTheSameFile.fastq doNotOverwriteMePlease.txt



## the rm (remove) command is also dangerous!!

rm doNotOverwriteMePlease.txt


mkdir newDirectory



mkdir /vol/funMicStorage/zoop



rmdir 




echo 




cat 




rm 




sudo 




top 




man top
top --help
top -h


head 

tail 

less 

## special characters




## ~  ##



## .  ##



## .. ## 













## > ##

echo "this a file" 

echo "this a file" > file.txt



## |  ##

#one program | second program

bbmap.sh -h

less



## /  ##

## see filepaths

## \  ##

ls -ltrh

ls -l -t -r -h

ls -l \
   -t \
   -r \
   -h

## &  ##



top



top &





## = and $ ##

X="this is a variable"
echo $X




## *  ##

## go back to someFiles

ls 

ls A12*

ls doNot*

ls *Please.txt



## # ## 

## this is a note 







## loops ##

for FILENAME in *
do
  echo here is a file: $FILENAME
done



#### advanced stuff ####

### more loops ###

## we can use files as lists:

## go back to someFiles directory

for i in $(cat A12-05.fastq)
do
  echo "this is data?" $i
done

## or look at some numbers...

for i in {0..10}; do
  echo "this is a number?" $i
done

### find command ###

## did you lose a file? find it with find!

cd 

find . -name doNotOverwriteMePlease.txt


find /vol/funmic/datasets -name sip2026.csv



## find one of the files in the someFiles directory
## try using MobaXterm and put it on your local computer. 
## Make sure you can find it on both computers.

## if you are using a mac, you'll need to do your 
## transfer the old-fashioned way, from the terminal
## on your home computer (log off of de.NBI VM or 
## open a new terminal

## scp
getFile=/home/ubuntu/someFiles/A12-05.fastq
putDir=/home/daniel/Documents/teaching/funmic/scratchpad
scp -i /home/daniel/.ssh -P 31993 -r ubuntu@129.70.51.6:$getFile $putDir

## or rsync:

getFile=/home/ubuntu/someFiles/A12-05.fastq
putDir=/home/daniel/Documents/teaching/funmic/scratchpad
rsync -auv \
  --progress \
  -e "ssh -p 31993" \
  ubuntu@129.70.51.6:$getFile $putDir

 

