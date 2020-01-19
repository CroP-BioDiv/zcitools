import shutil
import urllib.request as request
from contextlib import closing


def download_url(url, filename):
    with closing(request.urlopen(url)) as r:
        with open(filename, 'wb') as f:
            shutil.copyfileobj(r, f)
