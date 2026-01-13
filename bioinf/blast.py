def blast_simplificado(query, seq_alvo, w=3, match=2, mismatch=-1, gap=-1):
    """
    Executa uma pesquisa BLAST simplificada.
    1. Cria mapa de palavras (seeds) da query.
    2. Encontra hits na sequência alvo.
    3. Estende os hits para encontrar o melhor alinhamento local (HSP).
    
    Args:
        query (str): Sequência de procura.
        seq_alvo (str): Sequência onde procurar (Database).
        w (int): Tamanho da palavra (word size). Default 3 (comum para proteínas/DNA curto).
        
    Returns:
        tuple: (Melhor_Subsequencia_Query, Melhor_Subsequencia_Alvo, Score, Indice_Inicio)
    """
    # Passo 1: Mapear a query (k-mers)
    mapa_query = _construir_mapa(query, w)
    
    # Passo 2: Encontrar Hits (Onde as sementes coincidem)
    hits = _encontrar_hits(seq_alvo, mapa_query, w)
    
    if not hits:
        return "", "", 0, -1
        
    # Passo 3: Estender Hits (HSP)
    melhor_hsp = (0, 0, 0, 0) # score, q_start, t_start, length
    
    for h in hits:
        # h é (pos_query, pos_target)
        hsp = _estender_hit(query, seq_alvo, h, w, match, mismatch)
        score_atual = hsp[0]
        
        if score_atual > melhor_hsp[0]:
            melhor_hsp = hsp
            
    # Reconstruir as strings do melhor resultado
    score, q_start, t_start, length = melhor_hsp
    sub_q = query[q_start : q_start + length]
    sub_t = seq_alvo[t_start : t_start + length]
    
    return sub_q, sub_t, score, t_start


def _construir_mapa(query, w):
    """Cria um dicionário {palavra: [lista_posicoes]}."""
    mapa = {}
    for i in range(len(query) - w + 1):
        palavra = query[i : i + w]
        if palavra not in mapa:
            mapa[palavra] = []
        mapa[palavra].append(i)
    return mapa


def _encontrar_hits(seq, mapa, w):
    """Procura as palavras do mapa na sequência alvo."""
    hits = []
    for i in range(len(seq) - w + 1):
        palavra = seq[i : i + w]
        if palavra in mapa:
            # Se encontrou, guarda o par (pos_query, pos_target)
            for pos_q in mapa[palavra]:
                hits.append((pos_q, i))
    return hits


def _estender_hit(query, seq, hit, w, match, mismatch):
    """
    Tenta estender o alinhamento para a esquerda e direita.
    Para quando o score começa a descer (Drop-off simples).
    """
    q_start, t_start = hit
    q_end = q_start + w
    t_end = t_start + w
    
    current_score = w * match # Score inicial da semente (assumindo match perfeito)
    max_score = current_score
    
    # Limites do melhor alinhamento encontrado
    best_q_start = q_start
    best_t_start = t_start
    best_len = w
    
    # 1. Estender para a ESQUERDA
    i = q_start - 1
    j = t_start - 1
    
    temp_score = current_score
    while i >= 0 and j >= 0:
        if query[i] == seq[j]:
            temp_score += match
        else:
            temp_score += mismatch
            
        if temp_score > max_score:
            max_score = temp_score
            best_q_start = i
            best_t_start = j
            best_len = (q_end - i) # Atualiza comprimento total
        elif temp_score < (max_score / 2): # Critério de paragem (drop-off)
            break
            
        i -= 1
        j -= 1
        
    # 2. Estender para a DIREITA
    # (Nota: simplificação - num BLAST real faríamos esquerda e direita separadas)
    # Aqui vamos assumir que a extensão para a esquerda definiu o novo 'start'
    # e vamos tentar crescer para a direita a partir do fim da semente original.
    
    i = q_end
    j = t_end
    # Recuperar o score máximo da esquerda para continuar
    current_run_score = max_score 
    
    while i < len(query) and j < len(seq):
        if query[i] == seq[j]:
            current_run_score += match
        else:
            current_run_score += mismatch
            
        if current_run_score > max_score:
            max_score = current_run_score
            # O comprimento aumenta
            best_len = (i - best_q_start) + 1
        elif current_run_score < (max_score / 2):
            break
            
        i += 1
        j += 1
        
    return max_score, best_q_start, best_t_start, best_len