class IRs:
    def __init__(self, sequence, ira, irb):
        self.sequence = sequence  # Bio.SeqRecord
        self.ira = ira            # zci_bio.utils.feature.Feature
        self.irb = irb

    def diffs(self):
        # Returns None if IRs are exactly the same, and something complicated if not.
        pass
