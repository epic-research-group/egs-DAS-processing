#! /bin/bash
# count sort and list the number of files in subdirectories inside the current directory
du -a | cut -d/ -f2 | sort | uniq -c | sort -nr > some_file.txt

# or this 
find . -type f | cut -d/ -f2 | sort | uniq -c > some_file.txt

# count the number of unique occurrences of a number in column 1 of the output file above
awk '{A[$1]++}END{for(i in A)print i,A[i]}' EGS_CASSM_fc.txt

# delete all file from a directory and its subdirectories of a certain extension
find . -name "*.ext" -type f -delete
# find and list all files of a certain extension from directory and subdirs
find . -name "*.ext" -type f
