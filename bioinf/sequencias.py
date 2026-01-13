def validar_dna(seq):
    if not seq:  # Se a string estiver vazia
        return False
        
    seq = seq.upper()
    validos = {'A', 'C', 'G', 'T'}
    return set(seq).issubset(validos)

def transcricao(dna_seq):
    if not dna_seq:
        return ""
    
    # Simplesmente troca T por U e garante maiúsculas
    return dna_seq.upper().replace('T', 'U')

def complemento_inverso(dna_seq):
    if not dna_seq:
        return ""
        
    dna_seq = dna_seq.upper()
    # Dicionário de complementos
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # 1. Calcular o complemento (base a base)
    complemento = ""
    for base in dna_seq:
        # Se a base não estiver no mapa (ex: N), mantém a original ou ignora
        complemento += mapa.get(base, base)
        
    # 2. Inverter a string (slice [::-1] é o truque do Python para inverter)
    return complemento[::-1]