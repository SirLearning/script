import re
import sqlite3
import uuid


class Type:
    DTA = 'DTA'
    DTC = 'DTC'
    DTH = 'DTH'
    DTM = 'DTM'
    DTT = 'DTT'
    DTX = 'DTX'
    DHH = 'DHH'
    DXX = 'DXX'
    RLC = 'RLC'
    RLG = 'RLG'
    RLX = 'RLX'
    RIX = 'RIX'
    RSX = 'RSX'
    XXX = 'XXX'


def to_lower(x):
    return str(x).lower()


class GffNode:
    id_patt = '\tID=([\w.-]+);?'
    key_hangkes = [to_lower,]

    def set_node(self, node: str):
        self.node = node
        self.keys = []
        self._nodes = []
        return self

    def set_key_handel(cls, handle):
        cls.key_handles.append(handle)

    def do_key_handel(cls, keys: str):
        for handel in cls.key_handles:
            keys = list(map(handel, keys))
        return keys

    def _parser_key(self, item: str):
        r = re.findall(self.id_patt, item, re.I | re.S)
        return r or None

    def _set_keys(self, keys: str):
        if not keys:
            return
        keys = self.do_key_handel(keys)
        self.keys = set(keys)

    def _get_key(self, item: str):
        return self.do_key_handel(self._parser_key(item))
