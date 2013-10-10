from collections import Counter
from config import config as cfg
from classifier.dataset import DomainTypeDataset
from classifier.classifier import DomainTypeClassifier
from common.enums import DomainType, Direction
from common.alignment_algo import cut_record_by_dir

mdd_list = [(dt, dr) for dt in DomainType for dr in Direction]


def generate_train(v):
    result = []
    for d_type in v:
        f = []
        r = []
        for rec in v[d_type]:
            if "(RC)" in rec.id:
                if len(r) == cfg.svm_train_size:
                    continue
                else:
                    r.append(str(rec.seq))
            else:
                if len(f) == cfg.svm_train_size:
                    continue
                else:
                    f.append(str(rec.seq))
            if len(f) == len(r) == cfg.svm_train_size:
                break
        result.append((d_type, Direction.forward, f))
        result.append((d_type, Direction.reversed, r))
    return result


def svm_create(v):
    dataset = DomainTypeDataset()
    for d_type, direction, train in generate_train(v):
        dataset.setData(mdd_list.index((d_type, direction)), train)
    cl = DomainTypeClassifier(cfg.svm_train_radius, True)
    print("Training on test dataset")
    cl.train(dataset)
    print("Training end")
    return cl


def svm_split(classifier, config, unsplitted):
    v2 = {}
    res = []
    for rec in unsplitted:
        d_type, direct = mdd_list[classifier.predict(str(rec.seq))]
        res.append((d_type, direct))
        if d_type not in v2:
            v2[d_type] = []
        v2[d_type].append(cut_record_by_dir(rec, direct, config[d_type]))
    return v2
