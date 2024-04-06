#!/bin/bash

# Load environment variables
source .env

# Check if the correct number of arguments is provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <file1> <file2> ... <fileN> <local_destination>"
    exit 1
fi

# Get the remote destination from the last argument
local_destination="${!#}"

# Function to download sinfle file
download_file() {
    local file="$1"
    local local_destination="$2"
    echo "Downloading file $file to $local_destination"
    #sshpass -p "$PASSWORD_SERVER" scp -r st76s_6@chome.metz.supelec.fr:/usr/users/st76s/st76s_6/SEM3D_ST7_project/"$file"  "$local_destination"
    scp -r st76s_6@chome.metz.supelec.fr:/usr/users/st76s/st76s_6/SEM3D_ST7_project/"$file"  "$local_destination"
    echo "==> File $file uploaded to $local_destination"
}

# Loop through all the arguments except the last one
for ((i=1; i<=$#-1; i++)); do
    source="${!i}"
    
    download_file "$source" "$local_destination"
done
