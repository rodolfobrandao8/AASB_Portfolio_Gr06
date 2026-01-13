def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-1):

    # 1. Calcular a matriz de scores
    matriz = _criar_matriz_pontuacao(seq1, seq2, match, mismatch, gap)
    
    # 2. Obter o score final (canto inferior direito)
    score = matriz[len(seq1)][len(seq2)]
    
    # 3. Reconstruir o alinhamento (Traceback)
    alin1, alin2 = _traceback(seq1, seq2, matriz, match, mismatch, gap)
    
    return alin1, alin2, score


def _criar_matriz_pontuacao(seq1, seq2, match, mismatch, gap):
    """Função auxiliar para preencher a matriz de programação dinâmica."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    
    # Criar matriz preenchida com zeros
    matriz = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Inicializar primeira linha e coluna com gaps acumulados
    for i in range(1, rows):
        matriz[i][0] = matriz[i-1][0] + gap
    for j in range(1, cols):
        matriz[0][j] = matriz[0][j-1] + gap
        
    # Preencher o resto da matriz
    for i in range(1, rows):
        for j in range(1, cols):
            # Calcular pontuação diagonal (Match ou Mismatch?)
            if seq1[i-1] == seq2[j-1]:
                score_diagonal = matriz[i-1][j-1] + match
            else:
                score_diagonal = matriz[i-1][j-1] + mismatch
                
            # Calcular pontuação vinda de cima (Gap na seq2)
            score_cima = matriz[i-1][j] + gap
            
            # Calcular pontuação vinda da esquerda (Gap na seq1)
            score_esquerda = matriz[i][j-1] + gap
            
            # Escolher o máximo
            matriz[i][j] = max(score_diagonal, score_cima, score_esquerda)
            
    return matriz


def _traceback(seq1, seq2, matriz, match, mismatch, gap):
    """Função auxiliar para reconstruir o caminho ótimo (do fim para o início)."""
    alin1 = ""
    alin2 = ""
    
    # Começar no canto inferior direito
    i = len(seq1)
    j = len(seq2)
    
    while i > 0 or j > 0:
        score_atual = matriz[i][j]
        
        # Verificar se veio da Diagonal (Match/Mismatch)
        # Nota: i>0 e j>0 garante que não estamos na borda
        score_diag = matriz[i-1][j-1] if (i > 0 and j > 0) else -99999
        match_score = match if (i > 0 and j > 0 and seq1[i-1] == seq2[j-1]) else mismatch
        
        # Verificar se veio de Cima ou Esquerda
        score_cima = matriz[i-1][j] if i > 0 else -99999
        score_esq = matriz[i][j-1] if j > 0 else -99999
        
        # Lógica de prioridade: Diagonal > Cima > Esquerda (podes ajustar)
        if i > 0 and j > 0 and score_atual == score_diag + match_score:
            alin1 += seq1[i-1]
            alin2 += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0 and score_atual == score_cima + gap:
            alin1 += seq1[i-1]
            alin2 += "-"
            i -= 1
        else: # Veio da esquerda
            alin1 += "-"
            alin2 += seq2[j-1]
            j -= 1
            
    # Inverter as strings porque construímos do fim para o início
    return alin1[::-1], alin2[::-1]

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1):
   
    # 1. Calcular a matriz (com a regra dos zeros)
    matriz, max_i, max_j = _criar_matriz_sw(seq1, seq2, match, mismatch, gap)
    
    # 2. O score final é o valor máximo encontrado na matriz inteira
    score = matriz[max_i][max_j]
    
    # 3. Reconstruir o alinhamento a partir desse máximo
    alin1, alin2 = _traceback_sw(seq1, seq2, matriz, max_i, max_j, match, mismatch, gap)
    
    return alin1, alin2, score


def _criar_matriz_sw(seq1, seq2, match, mismatch, gap):
    """Gera a matriz para Smith-Waterman e guarda onde está o valor máximo."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matriz = [[0 for _ in range(cols)] for _ in range(rows)]
    
    max_score = 0
    max_i, max_j = 0, 0
    
    # Nota: A primeira linha e coluna já são zeros, não precisamos de fazer nada (diferente do NW)
    
    for i in range(1, rows):
        for j in range(1, cols):
            # Calcular scores
            if seq1[i-1] == seq2[j-1]:
                score_diag = matriz[i-1][j-1] + match
            else:
                score_diag = matriz[i-1][j-1] + mismatch
                
            score_cima = matriz[i-1][j] + gap
            score_esq = matriz[i][j-1] + gap
            
            # A regra de ouro: O score nunca pode ser negativo (mínimo é 0)
            matriz[i][j] = max(0, score_diag, score_cima, score_esq)
            
            # Guardar o recorde (onde começa o melhor alinhamento local)
            if matriz[i][j] > max_score:
                max_score = matriz[i][j]
                max_i, max_j = i, j
                
    return matriz, max_i, max_j


def _traceback_sw(seq1, seq2, matriz, start_i, start_j, match, mismatch, gap):
    """Reconstrói o caminho parando ao encontrar zero."""
    alin1 = ""
    alin2 = ""
    i, j = start_i, start_j
    
    # Loop continua enquanto não encontrarmos um zero na matriz
    while matriz[i][j] != 0:
        score_atual = matriz[i][j]
        
        # Recalcular de onde poderia ter vindo para decidir o caminho
        score_diag = matriz[i-1][j-1]
        match_score = match if seq1[i-1] == seq2[j-1] else mismatch
        
        score_cima = matriz[i-1][j]
        score_esq = matriz[i][j-1]
        
        if score_atual == score_diag + match_score:
            alin1 += seq1[i-1]
            alin2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_atual == score_cima + gap:
            alin1 += seq1[i-1]
            alin2 += "-"
            i -= 1
        else: # Esquerda
            alin1 += "-"
            alin2 += seq2[j-1]
            j -= 1
            
    return alin1[::-1], alin2[::-1]

def gerar_consenso(alin1, alin2):
    """
    Gera uma sequência de consenso a partir de duas sequências alinhadas.
    Regra simples: Se iguais mantém, se diferentes escolhe o caráter não-gap ou o primeiro.
    """
    consenso = ""
    for c1, c2 in zip(alin1, alin2):
        if c1 == c2:
            consenso += c1
        elif c1 == '-':
            consenso += c2
        elif c2 == '-':
            consenso += c1
        else:
            # Em caso de conflito (ex: A vs T), escolhemos arbitrariamente o primeiro
            # Numa implementação avançada, usaríamos contagem de frequências
            consenso += c1 
    return consenso

def alinhamento_multiplo(lista_seqs, match=2, mismatch=-1, gap=-1):
    """
    Realiza alinhamento múltiplo progressivo (Heurística Gulosa).
    1. Procura o par de sequências com melhor score.
    2. Alinha-as e gera um consenso.
    3. Substitui o par pelo consenso na lista.
    4. Repete até sobrar apenas uma sequência (o consenso global).
    
    Retorna:
        str: A sequência de consenso final.
    """
    # Copiar lista para não estragar a original
    seqs_atuais = lista_seqs[:]
    
    while len(seqs_atuais) > 1:
        melhor_score = -float('inf')
        melhor_par = (-1, -1)
        melhor_alinhamento = ("", "")
        
        # 1. Comparar todos contra todos para achar o par mais compatível
        for i in range(len(seqs_atuais)):
            for j in range(i + 1, len(seqs_atuais)):
                # Usamos o NW que já criaste para ver o score
                s1, s2, score = needleman_wunsch(seqs_atuais[i], seqs_atuais[j], match, mismatch, gap)
                
                if score > melhor_score:
                    melhor_score = score
                    melhor_par = (i, j)
                    melhor_alinhamento = (s1, s2)
        
        # 2. Gerar consenso do melhor par
        s1_alin, s2_alin = melhor_alinhamento
        novo_consenso = gerar_consenso(s1_alin, s2_alin)
        
        # 3. Remover as sequências originais e adicionar o consenso
        i, j = melhor_par
        # Remover do maior índice para o menor para não afetar a ordem
        seqs_atuais.pop(max(i, j))
        seqs_atuais.pop(min(i, j))
        
        seqs_atuais.append(novo_consenso)
        
    # No final, sobra apenas uma string na lista: o consenso global
    return seqs_atuais[0]