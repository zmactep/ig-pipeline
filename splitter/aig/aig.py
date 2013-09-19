import copy
import json
from common.enums import SeqType, DomainType


class AIG(object):
    def __init__(self):
        self.name = ""
        self.s_type = SeqType.undefined
        self.d_type = DomainType.undefined
        self.seq = ""
        self.v_gene = (0, 0)
        self.d_gene = (0, 0)
        self.j_gene = (0, 0)
        self.regions = {}
        self.germlines = []
        self.homologs = []

    def dumps(self):
        a = copy.copy(self.__dict__)
        a['s_type'] = a['s_type'].name
        a['d_type'] = a['d_type'].name
        return json.dumps(a)

    def loads(self, s):
        self.__dict__ = json.loads(s)
        self.s_type = getattr(SeqType, self.s_type)
        self.d_type = getattr(SeqType, self.d_type)
