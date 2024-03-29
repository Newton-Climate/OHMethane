#!/bin/bash

# Directory where you want to save the downloaded files
download_dir=$(pwd)

# Ensure the download directory exists
mkdir -p "$download_dir"
cd "$download_dir"

# Read each URL from the urls.txt file
while IFS= read -r url; do
    wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies -O "$(basename "$url")" "$url"
done < "merra_urls.txt"

