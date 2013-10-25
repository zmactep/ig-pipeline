import numpy as np
# from sklearn import svm
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from collections import Counter
from classifier.dataset import encode, create_region


class DomainTypeClassifier(object):
    def __init__(self, radius, window_mode=False):
        self.classifier = AdaBoostClassifier(
            DecisionTreeClassifier(max_depth=2),
            n_estimators=20,
            learning_rate=1,
            algorithm="SAMME")
        # svm.SVC(kernel='rbf')
        self.radius = radius
        self.window_mode = window_mode

    def train(self, dataset):
        k = self.radius if not self.window_mode else 2 * self.radius + 1
        rin, rout = dataset.getData(k, self.window_mode)
        print("fitting", len(rin))
        self.classifier.fit(np.asarray(rin, float), np.asarray(rout, float))

    def predict(self, ns):
        k = self.radius if not self.window_mode else 2 * self.radius + 1
        to_predict = []
        for i in range(len(ns)):
            if not self.window_mode:
                to_predict.append(encode(create_region(ns, i, k)))
            else:
                if i > len(ns) - k:
                    break
                to_predict.append(encode(ns[i:i+k]))
        return int(Counter(self.classifier.predict(
            np.asarray(to_predict, float))).most_common(1)[0][0])
