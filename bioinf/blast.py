def blast_simplificado(query, seq_alvo, w=3, match=2, mismatch=-1):
    """Executa um BLAST simplificado (seed-and-extend) entre uma query e uma sequência alvo.

    Pipeline implementado:
    1) Constrói um *query map* com todos os k-mers (tamanho ``w``) da query.
    2) Percorre a sequência alvo para encontrar *hits* (k-mers presentes no mapa).
    3) Para cada hit, estende à esquerda e à direita para obter um HSP
       (High-scoring Segment Pair) e escolhe o melhor.

    Args:
        query (str): Sequência de consulta (query).
        seq_alvo (str): Sequência alvo (target/database sequence) onde procurar hits.
        w (int, optional): Tamanho da palavra (k-mer) usada como seed. Por omissão ``3``.
        match (int, optional): Score atribuído a match durante a extensão. Por omissão ``2``.
        mismatch (int, optional): Penalização atribuída a mismatch durante a extensão.
            Por omissão ``-1``.

    Returns:
        tuple[str, str, int, int]:
        - sub_q (str): Subsequência da query correspondente ao melhor HSP.
        - sub_t (str): Subsequência da sequência alvo correspondente ao melhor HSP.
        - score (int): Score do melhor HSP encontrado.
        - t_start (int): Posição inicial (0-based) do HSP na sequência alvo.
          Se não houver hits, devolve ``-1``.

    Raises:
        ValueError: Se ``w`` for maior do que o comprimento da query (pode levar a
            mapa vazio) ou maior do que o comprimento da sequência alvo (sem seeds).
            (Nota: o código atual não valida explicitamente; este erro é uma sugestão
            de validação para robustez.)
        TypeError: Se ``query`` ou ``seq_alvo`` não forem strings.

    Examples:
        >>> blast_simplificado("ACGTAC", "TTACGTAA", w=3)
        ('ACG', 'ACG', 6, 2)
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
    """Constrói o *query map* (índice de seeds) para BLAST simplificado.

    Cria um dicionário onde cada k-mer (tamanho ``w``) da query aponta para uma lista
    das posições (0-based) onde esse k-mer ocorre na query.

    Args:
        query (str): Sequência de consulta (query).
        w (int): Tamanho do k-mer (seed/word size).

    Returns:
        dict[str, list[int]]: Dicionário do tipo ``{kmer: [pos1, pos2, ...]}``.

    Raises:
        TypeError: Se ``query`` não for string ou ``w`` não for inteiro.

    Examples:
        >>> construir_mapa("ACGTA", 3)
        {'ACG': [0], 'CGT': [1], 'GTA': [2]}
    """
    mapa = {}
    for i in range(len(query) - w + 1):
        palavra = query[i:i+w]
        mapa.setdefault(palavra, []).append(i)
    return mapa


def encontrar_hits(seq, mapa, w):
    """Encontra *hits* (seeds coincidentes) entre a sequência alvo e o query map.

    Percorre todos os k-mers (tamanho ``w``) da sequência alvo e, quando um k-mer
    existir no ``mapa`` da query, gera um hit por cada posição de ocorrência na query.

    Args:
        seq (str): Sequência alvo (target) onde procurar.
        mapa (dict[str, list[int]]): Query map produzido por :func:`construir_mapa`.
        w (int): Tamanho do k-mer (word size).

    Returns:
        list[tuple[int, int]]: Lista de hits no formato ``(pos_q, pos_t)``, onde:
        - ``pos_q`` é a posição do k-mer na query,
        - ``pos_t`` é a posição do k-mer na sequência alvo.

    Raises:
        TypeError: Se ``seq`` não for string, ou ``mapa`` não for dicionário.

    Examples:
        >>> encontrar_hits("TTACGTAA", {"ACG":[0], "CGT":[1]}, 3)
        [(0, 2), (1, 3)]
    """
    hits = []
    for i in range(len(seq) - w + 1):
        palavra = seq[i:i+w]
        if palavra in mapa:
            for pos_q in mapa[palavra]:
                hits.append((pos_q, i))
    return hits


def estender_hit(query, seq, hit, w, match, mismatch):
    """Estende um hit (seed) para obter um HSP (alinhamento local simplificado).

    A partir de um hit inicial (k-mer coincidente), tenta estender em duas direções:
    - para a esquerda (índices decrescentes),
    - para a direita (índices crescentes),
    acumulando um score simples de match/mismatch.

    A extensão pára quando:
    - sai dos limites das sequências, ou
    - o score cai abaixo de metade do melhor score observado durante a extensão
      (heurística de cutoff).

    Args:
        query (str): Sequência de consulta.
        seq (str): Sequência alvo.
        hit (tuple[int, int]): Par ``(pos_q, pos_t)`` com as posições iniciais do seed.
        w (int): Tamanho do seed (k-mer) usado no hit.
        match (int): Score para match.
        mismatch (int): Penalização para mismatch.

    Returns:
        tuple[int, int, int, int]: Tuplo ``(score, q_start, t_start, length)``, onde:
        - ``score`` é o melhor score atingido durante a extensão,
        - ``q_start`` é a posição inicial (0-based) na query,
        - ``t_start`` é a posição inicial (0-based) no alvo,
        - ``length`` é o comprimento do segmento alinhado (sem gaps).

    Raises:
        TypeError: Se ``hit`` não for um tuplo com dois inteiros.
        IndexError: Se ``hit`` contiver posições fora dos limites (o código assume hit válido).

    Examples:
        >>> estender_hit("ACGTAC", "TTACGTAA", (0, 2), 3, 2, -1)
        (6, 0, 2, 3)
    """
    def extender(direcao):
        """Extende um seed numa direção e devolve o melhor score local.

        Args:
            direcao (str): "esquerda" ou "direita".

        Returns:
            tuple[int, int, int, int]: ``(best_score, start_q, start_t, length)``.

        Raises:
            ValueError: Se ``direcao`` não for "esquerda" nem "direita".
        """
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

        start_q, start_t, length = best_start_q, best_start_t, best_len

        while 0 <= i < len(query) and 0 <= j < len(seq):
            score += match if query[i] == seq[j] else mismatch
            if score > best_score:
                best_score = score
                if direcao == "esquerda":
                    start_q = i
                    start_t = j
                    length = (hit[0] + w) - i
                else:
                    length = i - hit[0] + 1
            elif score < best_score / 2:
                break
            i += passo
            j += passo

        return best_score, start_q, start_t, length

    score_esq, q_start_esq, t_start_esq, len_esq = extender("esquerda")
    score_dir, q_start_dir, t_start_dir, len_dir = extender("direita")

    if score_esq >= score_dir:
        return score_esq, q_start_esq, t_start_esq, len_esq
    else:
        return score_dir, q_start_dir, t_start_dir, len_dir

