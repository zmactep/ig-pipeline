import json
import copy
from common.enums import DomainType, Direction


class SplitterMID(object):
    def __init__(self, dt=DomainType.undefined, direct=Direction.undefined):
        self.d_type = dt
        self.direction = direct
        self.mid = ""
        self.primer = ""

    def dumpd(self):
        a = copy.copy(self)
        a.d_type = a.d_type.name
        a.direction = a.direction.name
        return a.__dict__

    def dumps(self):
        return json.dumps(self.dumpd())

    def loadd(self, d):
        self.__dict__ = d
        self.d_type = getattr(DomainType, self.d_type)
        self.direction = getattr(Direction, self.direction)

    def loads(self, s):
        self.loadd(json.loads(s))


class SplitterConfiguration(object):
    def __init__(self):
        self.mids = []

    def dumps(self):
        a = []
        for i in self.mids:
            a.append(i.dumpd())
        return json.dumps({"mids": a})

    def loads(self, s):
        a = json.loads(s)
        for i in a["mids"]:
            sm = SplitterMID()
            sm.loadd(i)
            self.mids.append(sm)


def generate_defaults():
    sc = SplitterConfiguration()
    for d_type in DomainType:
        if d_type == DomainType.undefined:
            continue
        for direction in Direction:
            if direction == Direction.undefined:
                continue
            mid = SplitterMID(d_type, direction)
            sc.mids.append(mid)
    return sc
