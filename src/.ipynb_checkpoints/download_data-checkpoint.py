import os
import re
import argparse
import requests
from bs4 import BeautifulSoup
from urllib.parse import unquote

def get_remote_entries(base_url):
    resp = requests.get(base_url)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, 'html.parser')
    entries = []
    for link in soup.find_all('a', href=True):
        href = link['href']
        if href.endswith('/'):
            continue
        name = unquote(href)
        entries.append((href, name))
    return entries

def download_matching(base_url, pattern, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    entries = get_remote_entries(base_url)
    regex = re.compile(pattern)
    matches = [(href, name) for href, name in entries if regex.match(name)]
    
    if not matches:
        return
    
    for href, name in matches:
        dest_path = os.path.join(dest_dir, name)
        if os.path.exists(dest_path):
            continue
        print(f"--- Downloading {name}")
        url = f"{base_url.rstrip('/')}/{href}"
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(dest_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

def main():
    parser = argparse.ArgumentParser(
    )
    parser.add_argument('--base_url', default='https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/')
    parser.add_argument('--pattern', default='^ELG_LOPnotqso_.*_clustering\\.(?:ran|dat)\\.fits$'
)
    parser.add_argument('--dest_dir', default='/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/raw/',)
    args = parser.parse_args()
    download_matching(args.base_url, args.pattern, args.dest_dir)

if __name__ == "__main__":
    main()