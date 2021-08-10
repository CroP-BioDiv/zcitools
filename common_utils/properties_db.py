import os
import json
from .exceptions import ZCItoolsValueError

_create_cmd = """
CREATE TABLE properties (
    key1 VARCHAR(50) NOT NULL,
    key2 VARCHAR(50) NOT NULL,
    data_type VARCHAR(5) NOT NULL,
    data TEXT,
    PRIMARY KEY (key1, key2));
"""


class PropertiesDB:
    def _json_dumps(x):
        return json.dumps(x, default=str)

    _to_db_methods = dict(
        NoneType=('none', lambda _: None),
        list=('json', _json_dumps),
        dict=('json', _json_dumps),
    )
    _from_db_methods = dict(
        int=int,
        float=float,
        none=lambda _: None,
        json=json.loads
    )

    def __init__(self, dbfile=None, directory=None):
        if dbfile:
            self._dbfile = dbfile
        else:
            _dir = directory or os.environ.get('ZCI_COMMON_DB')
            if not _dir:
                raise ZCItoolsValueError('Properties db filename is not specified!')
            self._dbfile = os.path.join(_dir, 'properties_db.sqlite')
        self._db = None

    def db(self):
        if not self._db:
            import sqlite3
            create_struct = not os.path.exists(self._dbfile)
            self._db = sqlite3.connect(self._dbfile)
            if create_struct:
                self._db.execute(_create_cmd)
                self._db.commit()
        return self._db

    def set_property(self, key1, key2, data):
        dt = type(data).__name__
        dt_m = self._to_db_methods.get(dt)
        if dt_m:
            data_type, m = dt_m
            data = m(data)
        else:
            data_type = dt
            data = str(data)
        db = self.db()
        db.execute(f"REPLACE INTO properties VALUES (?, ?, ?, ?);", (key1, key2, data_type, data))
        db.commit()

    def get_property(self, key1, key2):
        cursor = self.db().execute(f"SELECT data_type, data FROM properties WHERE key1 = ? AND key2 = ?", (key1, key2))
        rec = cursor.fetchone()
        if rec:
            data_type, data = rec
            m = self._from_db_methods.get(data_type)
            return m(data) if m else data

    def fetch_property(self, key1, key2, _callable, *args, **kwargs):
        if (vals := self.get_property(key1, key2)) is None:
            if (vals := _callable(*args, **kwargs)) is not None:
                self.set_property(key1, key2, vals)
        return vals

    #
    def not_stored_keys1(self, keys1, key2):
        k1 = ','.join(f"'{k}'" for k in keys1)
        sql = f"SELECT key1 FROM properties WHERE key1 IN ({k1}) AND key2 = '{key2}'"
        cursor = self.db().execute(sql)
        return set(keys1) - set(r[0] for r in cursor.fetchall())

    def get_properties_keys1(self, keys1, key2):
        k1 = ','.join(f"'{k}'" for k in keys1)
        sql = f"SELECT key1, data_type, data FROM properties WHERE key1 IN ({k1}) AND key2 = '{key2}'"
        cursor = self.db().execute(sql)
        return dict((k, (m(data) if (m := self._from_db_methods.get(data_type)) else data))
                    for k, data_type, data in cursor.fetchall())

    def fetch_properties_keys1(self, keys1, key2, _callable, *args, **kwargs):
        data = self.get_properties_keys1(keys1, key2)
        for k1 in sorted(set(keys1) - set(data.keys())):
            data[k1] = vals = _callable(k1, *args, **kwargs)
            if vals is not None:
                self.set_property(k1, key2, vals)
        return data

    #
    def get_properties_key2_like(self, key1, key2):
        cursor = self.db().execute(
            "SELECT key2, data_type, data FROM properties WHERE key1 = ? AND key2 LIKE ?", (key1, key2))
        return dict((k, (m(data) if (m := self._from_db_methods.get(data_type)) else data))
                    for k, data_type, data in cursor.fetchall())
