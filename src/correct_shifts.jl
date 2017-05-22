"""rifraf-style fast frameshift correction"""
function correct_shifts(consensus::DNASeq,
                        reference::DNASeq;
                        bandwidth=10)
