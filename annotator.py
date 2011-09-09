'''Handle the annotation of auxillary informations for equations
'''

class Annotator(object):

    def __init__(self):
        self.annotating = False
        self.states = [[None, None]]

    def annotate(self, eq):
        if self.annotating:
            self.states[-1][1].append(eq)

    def push(self, label=None, annotate=True):
        self.states.append([label, [] if annotate else None])
        self.annotating = annotate
        return label

    def pop(self, label=None):
        assert len(self.states) >= 2
        _label, acc = self.states.pop()
        assert label is None or label == _label
        self.annotating = self.states[-1][1] is not None
        return acc

annotator = Annotator()
