import shutil
# import urllib.request as request
# from contextlib import closing
import urllib3
import json


def download_url(url, filename):
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    with open(filename, 'wb') as f:
        f.write(response.data)


def get_url_json(url):
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    if data := response.data.decode('utf-8'):
        return json.loads(data)
