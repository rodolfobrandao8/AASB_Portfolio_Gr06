def blast_simplificado(query, seq_alvo, w=3, match=2, mismatch=-1):
    """
    BLAST simplificado:
    1. Query map (k-mers)
    2. Identificação de hits
    3. Extensão
    4. Melhor alinhamento local (HSP)
    """
    mapa_query = construir_mapa(query, w)
    hits = encontrar_hits(seq_alvo, mapa_query, w)

    if not hits:
        return "", "", 0, -1

    melhor_hsp = (0, 0, 0, 0)  # score, q_start, t_start, length

    for hit in hits:
        hsp = estender_hit(query, seq_alvo, hit, w, match, mismatch)
        if hsp[0] > melhor_hsp[0]:
            melhor_hsp = hsp

    score, q_start, t_start, length = melhor_hsp
    sub_q = query[q_start:q_start + length]
    sub_t = seq_alvo[t_start:t_start + length]

    return sub_q, sub_t, score, t_start


def construir_mapa(query, w):
    """Cria o query map: {k-mer: [posições]}"""
    mapa = {}
    for i in range(len(query) - w + 1):
        palavra = query[i:i+w]
        mapa.setdefault(palavra, []).append(i)
    return mapa


def encontrar_hits(seq, mapa, w):
    """Identifica hits entre query map e sequência alvo"""
    hits = []
    for i in range(len(seq) - w + 1):
        palavra = seq[i:i+w]
        if palavra in mapa:
            for pos_q in mapa[palavra]:
                hits.append((pos_q, i))
    return hits


def estender_hit(query, seq, hit, w, match, mismatch):
    """
    Extensão de um hit para esquerda e direita (HSP)
    com critério de drop-off simples
    """
    q_start, t_start = hit
    q_end = q_start + w
    t_end = t_start + w

    score_atual = w * match
    max_score = score_atual

    best_q_start = q_start
    best_t_start = t_start
    best_len = w

    # Extensão à esquerda
    i, j = q_start - 1, t_start - 1
    temp_score = score_atual

    while i >= 0 and j >= 0:
        temp_score += match if query[i] == seq[j] else mismatch

        if temp_score > max_score:
            max_score = temp_score
            best_q_start = i
            best_t_start = j
