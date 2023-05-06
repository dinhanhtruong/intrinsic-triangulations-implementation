#!/bin/bash

# Check that two arguments were provided
if [ $# -ne 3 ]; then
  echo "Usage: $0 <name of img dir> <output name> <delete img dir = 0 or 1>"
  exit 1
fi

dir_name=$1
out_name=$2
delete_imgs=$3

echo "Compiling images in $dir_name to video $out_name"

ffmpeg -framerate 30 -pattern_type glob -i "$dir_name/*.png" -c:v libx264 -pix_fmt yuv420p $out_name;

echo "Done creating $out_name"

if [ $delete_imgs == "1" ]; then
    echo "Deleting $dir_name"
    rm -rf $dir_name
fi

echo "Success!"

# Exit with a success status code
exit 0