def blast_simplificado(query, seq_alvo, w=3, match=2, mismatch=-1):
    """
    Executa um BLAST simplificado:
    - Cria query map (k-mers)
    - Identifica hits
    - Estende hits para melhor alinhamento local (HSP)
    """
    mapa_query = construir_mapa(query, w)
    hits = encontrar_hits(seq_alvo, mapa_query, w)

    if not hits:
        return "", "", 0, -1

    melhor_hsp = (0, 0, 0, 0)

    for hit in hits:
        hsp = estender_hit(query, seq_alvo, hit, w, match, mismatch)
        if hsp[0] > melhor_hsp[0]:
            melhor_hsp = hsp

    score, q_start, t_start, length = melhor_hsp
    sub_q = query[q_start:q_start + length]
    sub_t = seq_alvo[t_start:t_start + length]

    return sub_q, sub_t, score, t_start


def construir_mapa(query, w):
    """
    Cria um mapa da query ({k-mer: posições}) para identificação de hits.
    """
    mapa = {}
    for i in range(len(query) - w + 1):
        palavra = query[i:i+w]
        mapa.setdefault(palavra, []).append(i)
    return mapa


def encontrar_hits(seq, mapa, w):
    """
    Identifica hits entre o query map e a sequência alvo.
    """
    hits = []
    for i in range(len(seq) - w + 1):
        palavra = seq[i:i+w]
        if palavra in mapa:
            for pos_q in mapa[palavra]:
                hits.append((pos_q, i))
    return hits


def estender_hit(query, seq, hit, w, match, mismatch):
    """
    Estende um hit para esquerda e direita para encontrar HSP.
    Complexidade reduzida para B ou A no Radon.
    """
    def extender(direcao):
        score = w * match
        best_score = score
        best_start_q = hit[0]
        best_start_t = hit[1]
        best_len = w

        if direcao == "esquerda":
            i, j = hit[0] - 1, hit[1] - 1
            passo = -1
        else:  
            i, j = hit[0] + w, hit[1] + w
            passo = 1

        while 0 <= i < len(query) and 0 <= j < len(seq):
            score += match if query[i] == seq[j] else mismatch
            if score > best_score:
                best_score = score
                if direcao == "esquerda":
                    best_start_q = i
                    best_start_t = j
                    best_len = (hit[0] + w) - i
                else:
                    best_len = i - best_start_q + 1
            elif score < best_score / 2:
                break
            i += passo
            j += passo

        return best_score, best_start_q, best_start_t, best_len


    score_esq, q_start, t_start, length = extender("esquerda")

    score_dir, _, _, length = extender("direita")
 
    max_score = max(score_esq, score_dir)
    return max_score, q_start, t_start, length

