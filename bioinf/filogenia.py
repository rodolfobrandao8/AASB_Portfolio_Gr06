def distancia_levenshtein(seq1, seq2):
    """Calcula a distância de edição (Levenshtein) entre duas strings.

    A distância de Levenshtein é o número mínimo de operações para transformar
    ``seq1`` em ``seq2``, onde as operações permitidas são:
    - inserção,
    - deleção,
    - substituição.

    Implementação por programação dinâmica com matriz de dimensão
    ``(len(seq1)+1) x (len(seq2)+1)``.

    Args:
        seq1 (str): Primeira string/sequência.
        seq2 (str): Segunda string/sequência.

    Returns:
        int: Distância de Levenshtein entre ``seq1`` e ``seq2``.

    Raises:
        TypeError: Se ``seq1`` ou ``seq2`` não forem strings (ou não suportarem ``len`` e indexação).

    Examples:
        >>> distancia_levenshtein("GATTACA", "GATTTCA")
        1
        >>> distancia_levenshtein("", "ABC")
        3
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matriz = [[0]*cols for _ in range(rows)]

    for i in range(rows):
        matriz[i][0] = i
    for j in range(cols):
        matriz[0][j] = j

    for i in range(1, rows):
        for j in range(1, cols):
            cost = 0 if seq1[i-1] == seq2[j-1] else 1
            matriz[i][j] = min(
                matriz[i-1][j] + 1,
                matriz[i][j-1] + 1,
                matriz[i-1][j-1] + cost
            )

    return matriz[-1][-1]


def matriz_distancias(seqs):
    """Constrói uma matriz de distâncias par-a-par para um conjunto de sequências.

    A distância usada é a distância de Levenshtein calculada por
    :func:`distancia_levenshtein`.

    A matriz é representada como um dicionário com chaves ``(seq_i, seq_j)``.
    O resultado inclui entradas simétricas, isto é, guarda (i,j) e (j,i).

    Args:
        seqs (list[str]): Lista de sequências/strings.

    Returns:
        dict[tuple[str, str], int]: Dicionário com as distâncias par-a-par.

    Raises:
        TypeError: Se ``seqs`` não for iterável.
        ValueError: Se ``seqs`` tiver menos de 2 sequências (a matriz ficará vazia).

    Examples:
        >>> m = matriz_distancias(["AA", "AB", "BB"])
        >>> m[("AA", "AB")]
        1
        >>> m[("AB", "AA")]
        1
    """
    matriz = {}
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            d = distancia_levenshtein(seqs[i], seqs[j])
            matriz[(seqs[i], seqs[j])] = d
            matriz[(seqs[j], seqs[i])] = d
    return matriz


def upgma(seqs):
    """Constrói uma árvore filogenética simplificada usando UPGMA.

    UPGMA (Unweighted Pair Group Method with Arithmetic Mean) é um método de
    clustering hierárquico que, a cada iteração:
    - escolhe o par de clusters com menor distância média,
    - funde-os num novo cluster,
    - repete até existir um único cluster.

    Nesta implementação:
    - a distância base entre sequências é Levenshtein,
    - a distância entre clusters é a média das distâncias entre todas as pares
      de sequências (um de cada cluster),
    - a árvore devolvida é uma estrutura de tuplos aninhados (ex.: ``('A', ('B','C'))``),
      não um objeto com comprimentos de ramos.

    Args:
        seqs (list[str]): Lista de sequências/strings a agrupar.

    Returns:
        object: Raiz da árvore (cluster final), representada como:
        - uma string (se só houver uma sequência), ou
        - um tuplo ``(cluster1, cluster2)`` recursivamente.

    Raises:
        IndexError: Se ``seqs`` estiver vazio (o código assume pelo menos 1 sequência).
        TypeError: Se ``seqs`` não for iterável.

    Examples:
        >>> upgma(["AA", "AB", "BB"])
        (('AA', 'AB'), 'BB')
    """
    clusters = {s: [s] for s in seqs}
    dist = matriz_distancias(seqs)

    while len(clusters) > 1:
        menor_dist = float('inf')
        par = (None, None)

        keys = list(clusters.keys())

        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                c1, c2 = keys[i], keys[j]

                dists = [
                    dist[(s1, s2)]
                    for s1 in clusters[c1]
                    for s2 in clusters[c2]
                ]
                media = sum(dists) / len(dists)

                if media < menor_dist:
                    menor_dist = media
                    par = (c1, c2)

        c1, c2 = par
        novo_cluster = (c1, c2)

        clusters[novo_cluster] = clusters[c1] + clusters[c2]
        del clusters[c1]
        del clusters[c2]

    return list(clusters.keys())[0]

