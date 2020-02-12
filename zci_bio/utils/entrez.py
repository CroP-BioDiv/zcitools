import datetime
import time
from .import_methods import import_bio_entrez


class Entrez:
    def __init__(self, calls_per_sec=3):
        assert calls_per_sec > 0, calls_per_sec
        self.entrez = import_bio_entrez()
        self.last_call = None
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
        records = self.entrez.read(handle)
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
