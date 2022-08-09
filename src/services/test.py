from ensemblrest import EnsemblRest

def fetch_transcript(transcript_id):
    ens_rest = EnsemblRest()
    transcript = ens_rest.getSequenceById(id=transcript_id,
                                          expand_5prime=5000,
                                          expand_3prime=5000)
    return f"{transcript}"


print(fetch_transcript("ENST00000315713"))