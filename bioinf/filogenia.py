def distancia_levenshtein(seq1, seq2):
    """
    Calcula a distância de edição (Levenshtein) entre duas strings.
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
    """
    Constrói a matriz de distâncias par-a-par entre sequências.
    """
    matriz = {}
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            d = distancia_levenshtein(seqs[i], seqs[j])
            matriz[(seqs[i], seqs[j])] = d
            matriz[(seqs[j], seqs[i])] = d
    return matriz


def upgma(seqs):
    """
    Constrói uma árvore filogenética simplificada usando UPGMA.
    A distância entre clusters é a média das distâncias entre as sequências.
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
