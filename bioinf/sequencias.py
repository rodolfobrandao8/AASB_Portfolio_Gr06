# =========================
# Validação de sequências
# =========================

def validar_dna(seq):
    if not seq:
        return False
    seq = seq.upper()
    return set(seq).issubset({'A', 'C', 'G', 'T'})


def validar_rna(seq):
    if not seq:
        return False
    seq = seq.upper()
    return set(seq).issubset({'A', 'C', 'G', 'U'})


def validar_proteina(seq):
    if not seq:
        return False
    seq = seq.upper()
    aminoacidos = {
        'A','C','D','E','F','G','H','I','K','L',
        'M','N','P','Q','R','S','T','V','W','Y'
    }
    return set(seq).issubset(aminoacidos)


# =========================
# Transcrição
# =========================

def transcricao(dna_seq):
    if not validar_dna(dna_seq):
        return ""
    return dna_seq.upper().replace('T', 'U')


# =========================
# Operações com DNA
# =========================

def complemento(dna_seq):
    if not validar_dna(dna_seq):
        return ""
        
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemento = ""
    for base in dna_seq.upper():
        complemento += mapa[base]
    return complemento


def reverso(seq):
    if not seq:
        return ""
    return seq[::-1]


def complemento_inverso(dna_seq):
    if not validar_dna(dna_seq):
        return ""
        
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemento = ""
    for base in dna_seq.upper():
        complemento += mapa[base]
    return complemento[::-1]


dna = "ATGCCGTA"
print(validar_dna(dna))