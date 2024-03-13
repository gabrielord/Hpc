#!/bin/bash

# Load environment variables
source .env

# Check if the correct number of arguments is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <file1> <file2> ... <fileN> "
    exit 1
fi

# Loop through all the arguments except the last one
for ((i=1; i<=$#; i++)); do
    file="${!i}"
    # We remove the eventual initial . and / from the file name
    file=$(echo $file | sed 's/^\.\///')
    echo "Uploading file $file"
    # Execute the expect script to upload the file
    sshpass -p $PASSWORD_SERVER scp -r "$file" st76s_6@chome.metz.supelec.fr:/usr/users/st76s/st76s_6/SEM3D_ST7_project/"$file"
    echo "==> File $file uploaded successfully!"
done