import requests

def read_urls_from_file(file_path):
    """Read URLs from a given file."""
    with open(file_path, 'r') as file:
        urls = [line.strip() for line in file.readlines() if line.strip()]
    return urls

def estimate_download_size(urls):
    """Estimate the total download size of the URLs."""
    total_size = 0
    for url in urls:
        try:
            response = requests.head(url)
            if response.status_code == 200:
                total_size += int(response.headers.get('content-length', 0))
            else:
                print(f"Error accessing {url}: Status code {response.status_code}")
        except Exception as e:
            print(f"Exception occurred for {url}: {e}")
    return total_size

def download_files(urls):
    """Download files from the given URLs."""
    for url in urls:
        local_filename = url.split('/')[-1]
        try:
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print(f"Downloaded {local_filename}")
        except Exception as e:
            print(f"Failed to download {local_filename}: {e}")

def main(file_path):
    """Main function to estimate download size and download files."""
    urls = read_urls_from_file(file_path)
    total_size = estimate_download_size(urls)
    print(f"Total estimated download size: {total_size / 1024 / 1024:.2f} MB")
    
    confirm = input("Do you want to proceed with the download? (yes/no): ").lower()
    if confirm == 'yes':
        download_files(urls)
    else:
        print("Download cancelled.")

if __name__ == "__main__":
    file_path = input("Enter the path to your file: ")
    main(file_path)

