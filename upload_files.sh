#!/bin/bash

# Load environment variables
source .env

# Check if the correct number of arguments is provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <file1> <file2> ... <fileN> <remote_destination>"
    exit 1
fi

# Get the remote destination from the last argument
remote_destination="${!#}"

# Loop through all the arguments except the last one
for ((i=1; i<=$#-1; i++)); do
    file="${!i}"
    # Execute the expect script to upload the file
    sshpass -p $PASSWORD_SERVER scp -r "$file" st76s_6@chome.metz.supelec.fr:/usr/users/st76s/st76s_6/SEM3D_ST7_project/"$remote_destination"
    echo "File $file uploaded to $remote_destination"
done
