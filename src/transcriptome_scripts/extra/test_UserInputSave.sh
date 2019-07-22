#!/bin/bash

echo "Please put in the base directory:"
read base
echo "Please put in the output folder name:"
read output

echo "Outputs saving to : " $base$output

if [ -d "$base$output" ]; then
    echo "Directory Already Exists, please use another name"
else
    echo "Directory Created"
    mkdir "$base$output"
fi

