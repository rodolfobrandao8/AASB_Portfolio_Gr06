BLOSUM62 = {
    'A': {'A': 4,  'C': 0,  'G': 0,  'T': 0},
    'C': {'A': 0,  'C': 9,  'G': -3, 'T': -1},
    'G': {'A': 0,  'C': -3, 'G': 6,  'T': -2},
    'T': {'A': 0,  'C': -1, 'G': -2, 'T': 5}
}
"""
dict[str, dict[str, int]]: Matriz de substituição usada por omissão nas funções de alinhamento.

Nota:
    Apesar do nome ``BLOSUM62``, esta matriz é definida apenas para o alfabeto
    ``A,C,G,T`` (DNA). Se forem usados caracteres fora deste alfabeto, as
    funções que acedem a esta matriz podem levantar ``KeyError``.
"""


def score_subst(a, b, matriz):
    """Obtém o score de substituição entre dois símbolos.

    Esta função faz um acesso direto a ``matriz[a][b]`` e assume que a matriz
    está completa para os símbolos fornecidos.

    Args:
        a (str): Símbolo da sequência 1 (tipicamente um carácter).
        b (str): Símbolo da sequência 2 (tipicamente um carácter).
        matriz (dict): Matriz de substituição no formato
            ``dict[símbolo][símbolo] -> score`` (ex.: PAM/BLOSUM ou matriz custom).

    Returns:
        int | float: Score de substituição definido na matriz.

    Raises:
        KeyError: Se ``a`` ou ``b`` não existir na matriz fornecida.

    Examples:
        >>> m = {'A': {'A': 1, 'C': -1}, 'C': {'A': -1, 'C': 1}}
        >>> score_subst('A', 'C', m)
        -1
    """
    return matriz[a][b]


def dot_plot(seq1, seq2):
    """Cria uma matriz de pontos (dot plot) binária entre duas sequências.

    Produz uma matriz ``len(seq1) x len(seq2)`` em que cada célula vale:
    - ``1`` se os caracteres forem iguais (match),
    - ``0`` caso contrário (mismatch).

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.

    Returns:
        list[list[int]]: Matriz de inteiros (0/1) com dimensão
        ``len(seq1) x len(seq2)``.

    Examples:
        >>> dot_plot("AC", "AGC")
        [[1, 0, 0], [0, 0, 1]]
    """
    return [[1 if a == b else 0 for b in seq2] for a in seq1]


def needleman_wunsch(seq1, seq2, matriz_subst=BLOSUM62, gap=-5):
    """Executa alinhamento global (Needleman–Wunsch) entre duas sequências.

    Implementa a fase de preenchimento da matriz de scores para alinhamento global
    e, no final, reconstrói um alinhamento ótimo via traceback.

    Args:
        seq1 (str): Primeira sequência a alinhar.
        seq2 (str): Segunda sequência a alinhar.
        matriz_subst (dict, optional): Matriz de substituição (ex.: BLOSUM/PAM ou
            matriz custom). Por omissão usa ``BLOSUM62``.
        gap (int, optional): Penalização linear de gap (inserção/deleção).
            Por omissão ``-5``.

    Returns:
        tuple[str, str, int | float]: Triplo ``(seq1_alinhada, seq2_alinhada, score)``,
        onde ``seq1_alinhada`` e ``seq2_alinhada`` incluem gaps ``'-'``, e
        ``score`` é o score ótimo do alinhamento global.

    Raises:
        KeyError: Se a matriz de substituição não contiver algum símbolo presente
            nas sequências.

    Examples:
        >>> a1, a2, sc = needleman_wunsch("AC", "AGC")
        >>> len(a1) == len(a2)
        True
    """
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
    """Reconstrói um alinhamento global ótimo a partir da matriz de scores.

    Esta função percorre a matriz ``m`` do canto inferior direito para o superior
    esquerdo, escolhendo transições compatíveis com o score ótimo (diagonal,
    cima, esquerda), gerando duas strings alinhadas (com gaps).

    Args:
        seq1 (str): Primeira sequência original.
        seq2 (str): Segunda sequência original.
        m (list[list[int | float]]): Matriz de programação dinâmica já preenchida.
        matriz_subst (dict): Matriz de substituição usada no preenchimento.
        gap (int): Penalização linear de gap usada no preenchimento.

    Returns:
        tuple[str, str, int | float]: Triplo ``(seq1_alinhada, seq2_alinhada, score_final)``.

    Raises:
        KeyError: Se a matriz de substituição não contiver símbolos necessários
            para comparar durante o traceback.
        IndexError: Se for chamado com ``m`` incompatível com os comprimentos de
            ``seq1`` e ``seq2``.

    Notes:
        Em empates, a preferência de caminho é:
        1) diagonal (match/mismatch),
        2) cima (gap na seq2),
        3) esquerda (gap na seq1).

    Examples:
        >>> a1, a2, sc = needleman_wunsch("A", "A")
        >>> (a1, a2, sc)
        ('A', 'A', 4)
    """
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
    """Executa alinhamento local (Smith–Waterman) entre duas sequências.

    Implementa a fase de preenchimento da matriz de scores para alinhamento local,
    aplicando a regra de truncar a 0 (não permitir scores negativos), e depois
    reconstrói o melhor alinhamento local via traceback a partir da célula de
    score máximo.

    Args:
        seq1 (str): Primeira sequência a alinhar.
        seq2 (str): Segunda sequência a alinhar.
        matriz_subst (dict, optional): Matriz de substituição (ex.: BLOSUM/PAM ou
            matriz custom). Por omissão usa ``BLOSUM62``.
        gap (int, optional): Penalização linear de gap (inserção/deleção).
            Por omissão ``-5``.

    Returns:
        tuple[str, str, int | float]: Triplo ``(subseq1_alinhada, subseq2_alinhada, score)``,
        onde as subsequências alinhadas correspondem ao melhor alinhamento local.

    Raises:
        KeyError: Se a matriz de substituição não contiver algum símbolo presente
            nas sequências.

    Examples:
        >>> a1, a2, sc = smith_waterman("ACGT", "TACG")
        >>> sc >= 0
        True
    """
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
    """Reconstrói o melhor alinhamento local a partir da matriz de scores.

    Começa na posição ``start`` (tipicamente a célula de score máximo) e segue
    o caminho de traceback até encontrar uma célula com score 0.

    Args:
        seq1 (str): Primeira sequência original.
        seq2 (str): Segunda sequência original.
        m (list[list[int | float]]): Matriz de programação dinâmica já preenchida.
        start (tuple[int, int]): Coordenadas (i, j) do ponto inicial do traceback.
        matriz_subst (dict): Matriz de substituição usada no preenchimento.
        gap (int): Penalização linear de gap usada no preenchimento.

    Returns:
        tuple[str, str, int | float]: Triplo ``(subseq1_alinhada, subseq2_alinhada, score_local)``.

    Raises:
        KeyError: Se a matriz de substituição não contiver símbolos necessários
            para comparar durante o traceback.
        IndexError: Se ``start`` estiver fora dos limites de ``m``.

    Notes:
        Tal como no traceback global, em empates a preferência de caminho é:
        1) diagonal, 2) cima, 3) esquerda.

    Examples:
        >>> a1, a2, sc = smith_waterman("A", "T")
        >>> sc
        0
    """
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
    """Calcula a sequência consenso a partir de um alinhamento múltiplo.

    Para cada coluna do alinhamento, escolhe o carácter mais frequente.
    Em caso de empate, escolhe o primeiro carácter que aparece nessa coluna.

    Args:
        alinhamento (list[str]): Lista de sequências já alinhadas (todas com o
            mesmo comprimento), tipicamente contendo gaps ``'-'``.

    Returns:
        str: Sequência consenso (inclui gaps se forem o símbolo mais frequente).

    Raises:
        TypeError: Se ``alinhamento`` não for iterável.
        ValueError: Se as sequências em ``alinhamento`` não tiverem o mesmo comprimento.

    Examples:
        >>> consenso_multiplas(["A-C", "ACC", "A-C"])
        'A-C'
    """
    consenso = ""
    for col in zip(*alinhamento):
        freq = {c: col.count(c) for c in set(col)}
        max_freq = max(freq.values())
        
        for c in col:
            if freq[c] == max_freq:
                consenso += c
                break
    return consenso


def alinhamento_multiplo(seqs, matriz_subst=BLOSUM62, gap=-5):
    """Realiza um alinhamento múltiplo “progressivo” simplificado.

    Estratégia implementada:
        - Inicia cada sequência como um alinhamento “unitário”.
        - Repetidamente, escolhe o par com melhor score de Needleman–Wunsch
          (aplicado às strings representativas atuais) e junta-o num alinhamento
          de 2 sequências.
        - No final devolve o alinhamento final e a sequência consenso.

    Args:
        seqs (list[str]): Lista de sequências a alinhar.
        matriz_subst (dict, optional): Matriz de substituição usada no NW.
            Por omissão usa ``BLOSUM62``.
        gap (int, optional): Penalização linear de gap. Por omissão ``-5``.

    Returns:
        tuple[list[str], str]: Par ``(alinhamento, consenso)`` onde:
        - ``alinhamento`` é uma lista de strings alinhadas (com gaps).
        - ``consenso`` é a sequência consenso calculada por ``consenso_multiplas``.

    Raises:
        IndexError: Se ``seqs`` estiver vazio (o código assume pelo menos 1 sequência).
        KeyError: Se a matriz de substituição não contiver símbolos usados.

    Examples:
        >>> alin, cons = alinhamento_multiplo(["AC", "AGC"])
        >>> len(alin)
        2
        >>> len(alin[0]) == len(alin[1])
        True
    """
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
    
        combinacao = [melhor_alin[0], melhor_alin[1]]
        alinhamentos.pop(j)
        alinhamentos.pop(i)
        alinhamentos.append(combinacao)

    return alinhamentos[0], consenso_multiplas(alinhamentos[0])
