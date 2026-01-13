BLOSUM62 = {
    'A': {'A': 4,  'C': 0,  'G': 0,  'T': 0},
    'C': {'A': 0,  'C': 9,  'G': -3, 'T': -1},
    'G': {'A': 0,  'C': -3, 'G': 6,  'T': -2},
    'T': {'A': 0,  'C': -1, 'G': -2, 'T': 5}
}


def score_subst(a, b, matriz):
    """Retorna o score de substituição entre dois caracteres."""
    return matriz[a][b]


def dot_plot(seq1, seq2):
    """Cria matriz de pontos (1 = match, 0 = mismatch)."""
    return [[1 if a == b else 0 for b in seq2] for a in seq1]


def needleman_wunsch(seq1, seq2, matriz_subst=BLOSUM62, gap=-5):
    """Alinhamento global (Needleman-Wunsch)."""
    rows, cols = len(seq1)+1, len(seq2)+1
    m = [[0]*cols for _ in range(rows)]

    for i in range(1, rows):
        m[i][0] = m[i-1][0] + gap
    for j in range(1, cols):
        m[0][j] = m[0][j-1] + gap

    for i in range(1, rows):
        for j in range(1, cols):
            diag = m[i-1][j-1] + score_subst(seq1[i-1], seq2[j-1], matriz_subst)
            cima = m[i-1][j] + gap
            esq  = m[i][j-1] + gap
            m[i][j] = max(diag, cima, esq)

    return _traceback_nw(seq1, seq2, m, matriz_subst, gap)


def _traceback_nw(seq1, seq2, m, matriz_subst, gap):
    """Reconstrói alinhamento global a partir da matriz de scores."""
    i, j = len(seq1), len(seq2)
    a1, a2 = "", ""

    while i > 0 or j > 0:
        if i > 0 and j > 0 and m[i][j] == m[i-1][j-1] + score_subst(seq1[i-1], seq2[j-1], matriz_subst):
            a1 += seq1[i-1]
            a2 += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0 and m[i][j] == m[i-1][j] + gap:
            a1 += seq1[i-1]
            a2 += "-"
            i -= 1
        else:
            a1 += "-"
            a2 += seq2[j-1]
            j -= 1

    return a1[::-1], a2[::-1], m[len(seq1)][len(seq2)]


def smith_waterman(seq1, seq2, matriz_subst=BLOSUM62, gap=-5):
    """Alinhamento local (Smith-Waterman)."""
    rows, cols = len(seq1)+1, len(seq2)+1
    m = [[0]*cols for _ in range(rows)]
    max_score, max_pos = 0, (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            diag = m[i-1][j-1] + score_subst(seq1[i-1], seq2[j-1], matriz_subst)
            cima = m[i-1][j] + gap
            esq  = m[i][j-1] + gap
            m[i][j] = max(0, diag, cima, esq)
            if m[i][j] > max_score:
                max_score = m[i][j]
                max_pos = (i, j)

    return _traceback_sw(seq1, seq2, m, max_pos, matriz_subst, gap)


def _traceback_sw(seq1, seq2, m, start, matriz_subst, gap):
    """Reconstrói alinhamento local a partir da matriz de scores."""
    i, j = start
    a1, a2 = "", ""
    while m[i][j] != 0:
        if m[i][j] == m[i-1][j-1] + score_subst(seq1[i-1], seq2[j-1], matriz_subst):
            a1 += seq1[i-1]
            a2 += seq2[j-1]
            i -= 1
            j -= 1
        elif m[i][j] == m[i-1][j] + gap:
            a1 += seq1[i-1]
            a2 += "-"
            i -= 1
        else:
            a1 += "-"
            a2 += seq2[j-1]
            j -= 1

    return a1[::-1], a2[::-1], m[start[0]][start[1]]


def consenso_multiplas(alinhamento):
    """Gera sequência consenso; desempate pelo primeiro elemento da coluna."""
    consenso = ""
    for col in zip(*alinhamento):
        freq = {c: col.count(c) for c in set(col)}
        max_freq = max(freq.values())
        candidatos = [c for c in col if freq[c] == max_freq]
        consenso += candidatos[0]
    return consenso


def alinhamento_multiplo(seqs, matriz_subst=BLOSUM62, gap=-5):
    """Alinhamento múltiplo progressivo usando Needleman-Wunsch."""
    if len(seqs) == 1:
        return [[seqs[0]]], seqs[0]

    alinhamentos = [[s] for s in seqs]

    while len(alinhamentos) > 1:
        melhor_score = -1e9
        melhor_par = None
        melhor_alin = None

        for i in range(len(alinhamentos)):
            for j in range(i+1, len(alinhamentos)):
                a1, a2, score = needleman_wunsch(alinhamentos[i][0], alinhamentos[j][0], matriz_subst, gap)
                if score > melhor_score:
                    melhor_score = score
                    melhor_par = (i, j)
                    melhor_alin = (a1, a2)

        i, j = melhor_par
        alinhamentos.pop(j)
        alinhamentos.pop(i)
        alinhamentos.append([melhor_alin[0], melhor_alin[1]])

    return alinhamentos[0], consenso_multiplas(alinhamentos[0])
