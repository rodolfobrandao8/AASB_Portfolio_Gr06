def validar_dna(seq):
    """
    Valida se uma sequência contém apenas nucleótidos de DNA (A, C, G, T).
    Retorna True se válido, False caso contrário.
    """
    if not seq:
        return False
    seq = seq.upper()
    return set(seq).issubset({'A', 'C', 'G', 'T'})


def validar_rna(seq):
    """
    Valida se uma sequência contém apenas nucleótidos de RNA (A, C, G, U).
    Retorna True se válido, False caso contrário.
    """
    if not seq:
        return False
    seq = seq.upper()
    return set(seq).issubset({'A', 'C', 'G', 'U'})


def validar_proteina(seq):
    """
    Valida se uma sequência contém apenas aminoácidos válidos (20 aminoácidos padrão).
    Retorna True se válido, False caso contrário.
    """
    if not seq:
        return False
    seq = seq.upper()
    aminoacidos = {
        'A','C','D','E','F','G','H','I','K','L',
        'M','N','P','Q','R','S','T','V','W','Y'
    }
    return set(seq).issubset(aminoacidos)


def transcricao(dna_seq):
    """
    Transcreve uma sequência de DNA para RNA (T → U).
    Retorna a sequência de RNA correspondente, ou string vazia se não for DNA válido.
    """
    if not validar_dna(dna_seq):
        return ""
    return dna_seq.upper().replace('T', 'U')


def complemento(dna_seq):
    """
    Retorna o complemento de uma sequência de DNA.
    Retorna string vazia se não for DNA válido.
    """
    if not validar_dna(dna_seq):
        return ""
        
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemento = ""
    for base in dna_seq.upper():
        complemento += mapa[base]
    return complemento


def reverso(seq):
    """
    Retorna a sequência invertida.
    Retorna string vazia se a sequência estiver vazia.
    """
    if not seq:
        return ""
    return seq[::-1]


def complemento_inverso(dna_seq):
    """
    Retorna o complemento reverso de uma sequência de DNA (5' → 3').
    Retorna string vazia se não for DNA válido.
    """
    if not validar_dna(dna_seq):
        return ""
        
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemento = ""
    for base in dna_seq.upper():
        complemento += mapa[base]
    return complemento[::-1]
