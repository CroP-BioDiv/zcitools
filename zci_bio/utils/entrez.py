import datetime
import time
from .import_methods import import_bio_entrez


class Entrez:
    def __init__(self, calls_per_sec=None):
        self.entrez = import_bio_entrez()
        self.last_call = None
        #
        if not calls_per_sec:
            calls_per_sec = 10 if self.entrez.email else 3
        assert calls_per_sec > 0, calls_per_sec
        self.between_calls = datetime.timedelta(milliseconds=round(1000 / calls_per_sec) + 1)

    def _sleep(self):
        # Wait
        if self.last_call:
            d = datetime.datetime.now() - self.last_call
            if d < self.between_calls:
                d = self.between_calls - d
                print(f'  Entrez: sleeping for {d.microseconds // 1000} ms')
                time.sleep(d.microseconds / 1000000)
        #
        self.last_call = datetime.datetime.now()

    def _read(self, handle):
        records = self.entrez.read(handle) if handle else None
        handle.close()
        return records

    #
    def efetch(self, filename, **kwargs):
        # Saves into filename
        self._sleep()
        with self.entrez.efetch(**kwargs) as handle:
            with open(filename, 'w') as out:
                out.write(handle.read())

    def esearch(self, **kwargs):
        self._sleep()
        return self._read(self.entrez.esearch(**kwargs))

    def esummary(self, **kwargs):
        self._sleep()
        return self._read(self.entrez.esummary(**kwargs))

    #
    def search_summary(self, db, retmax=None, **kwargs):
        search = self.esearch(db=db, usehistory='y', retmax=1, **kwargs)  # retmax=1: to have less network traffic
        if int(search['Count']):
            if retmax:
                return self.esummary(db=db, query_key=search['QueryKey'], WebEnv=search['WebEnv'], retmax=retmax)
            return self.esummary(db=db, query_key=search['QueryKey'], WebEnv=search['WebEnv'])

    def search_count(self, db, info=True, **kwargs):
        if info:
            print(f'Entrez search count ({db}): {", ".join(f"{k}={v}" for k, v in kwargs.items())}')
        search = self.esearch(db=db, usehistory='y', retmax=1, **kwargs)  # retmax=1: to have less network traffic
        return int(search['Count'])
