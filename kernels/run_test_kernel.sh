#!/bin/bash

#This script was optimized with the assistance of ChatGPT by OpenAI.
echo `python --version`
echo "mpmath version: $(python -c 'import mpmath; print(mpmath.__version__)')"

START_TIME=$(date +%s)
./laplace_inverse.py 0.25 'pstrain' 100
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

echo "Execution time: $ELAPSED_TIME seconds"

# Define local file and URL

for comp in "00" "01" "11"; do
    LOCAL_FILE="./nu0.25_h${comp}.txt"
    URL="https://gitlab.com/uguca/uguca/-/raw/main/kernels/nu0.25_h${comp}.txt?ref_type=heads&inline=false"

    # Download the remote file temporarily
    TMP_FILE="./original_nu0.25_h${comp}.txt"
    curl -sSL "$URL" -o "$TMP_FILE"

    # Check if the download was successful
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download the file from $URL"
        exit 1
    fi

    # Compare the files
    diff "$LOCAL_FILE" "$TMP_FILE" > /dev/null
    if [ $? -eq 0 ]; then
        echo "${comp} Files match."
    else
        echo "${comp} Files differ."
    fi

    # Remove temporary file
    rm "$TMP_FILE"

done

