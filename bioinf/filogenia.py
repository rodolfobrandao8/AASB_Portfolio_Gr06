def distancia_levenshtein(seq1, seq2):
    """Calcula a distância de edição entre duas strings."""
    if seq1 == seq2: return 0
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matriz = [[0 for _ in range(cols)] for _ in range(rows)]
    
    for i in range(1, rows): matriz[i][0] = i
    for j in range(1, cols): matriz[0][j] = j
            
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i-1] == seq2[j-1]:
                cost = 0
            else:
                cost = 1
            matriz[i][j] = min(matriz[i-1][j]+1, matriz[i][j-1]+1, matriz[i-1][j-1]+cost)
    return matriz[rows-1][cols-1]

def upgma_simples(lista_seqs):
    # Copiar a lista para não estragar a original
    clusters = list(lista_seqs)
    
    while len(clusters) > 1:
        melhor_dist = float('inf')
        melhor_par = (-1, -1)
        
        # Comparar todos contra todos para achar o par mais próximo
        # Nota: Isto recalcula distâncias a cada passo (menos eficiente, mas mais fácil de ler/implementar)
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                
                # Se for tuplo, temos de "imaginar" uma string representativa ou usar média
                # Simplificação: Usamos a representação string para calcular distância
                # (Isto não é o UPGMA matemático perfeito das médias, mas serve para gerar a topologia)
                s1 = str(clusters[i]) if not isinstance(clusters[i], str) else clusters[i]
                s2 = str(clusters[j]) if not isinstance(clusters[j], str) else clusters[j]
                
                # Se forem tuplos, extrair apenas uma string "amostra" para comparar
                # (Truque para o código não crashar)
                if isinstance(clusters[i], tuple): s1 = _extrair_seq(clusters[i])
                if isinstance(clusters[j], tuple): s2 = _extrair_seq(clusters[j])
                
                d = distancia_levenshtein(s1, s2)
                
                if d < melhor_dist:
                    melhor_dist = d
                    melhor_par = (i, j)
        
        # Agrupar
        i, j = melhor_par
        novo_cluster = (clusters[i], clusters[j])
        
        # Remover (do maior para o menor índice)
        clusters.pop(max(i, j))
        clusters.pop(min(i, j))
        
        # Adicionar novo
        clusters.append(novo_cluster)
        
    return clusters[0]

def _extrair_seq(cluster):
    """Função auxiliar para tirar uma string de dentro de tuplos aninhados."""
    if isinstance(cluster, str): return cluster
    return _extrair_seq(cluster[0]) # Vai sempre pela esquerda