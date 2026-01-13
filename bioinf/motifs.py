import re

def prosite_para_regex(padrao_prosite):
    # 1. Remover o '-' que separa os elementos
    regex = padrao_prosite.replace('-', '')
    
    # 2. 'x' vira '.' (qualquer caracter)
    regex = regex.replace('x', '.')
    
    # 3. '{AB}' vira '[^AB]' (qualquer coisa EXCETO A e B)
    regex = regex.replace('{', '[^').replace('}', ']')
    
    # 4. '(3)' vira '{3}' (repetição)
    regex = regex.replace('(', '{').replace(')', '}')
    
    return regex

def procurar_motifs(sequencia, padrao_prosite):

    regex = prosite_para_regex(padrao_prosite)
    
    # Usamos re.finditer para achar todas as posições
    # (?=...) é um lookahead para permitir overlaps, mas para simplificar
    # vamos usar a busca padrão que é mais comum em exercícios básicos.
    matches = re.finditer(regex, sequencia)
    
    posicoes = []
    for m in matches:
        posicoes.append(m.start())
        
    return posicoes

def criar_pwm(lista_seqs, pseudocount=1):
   
    if not lista_seqs: return []
    
    comp = len(lista_seqs[0])
    num_seqs = len(lista_seqs)
    alfabeto = "ACGT"
    
    pwm = []
    
    for i in range(comp):
        # 1. Contar
        coluna = [s[i] for s in lista_seqs]
        contagens = {base: 0 for base in alfabeto}
        for base in coluna:
            if base in contagens: contagens[base] += 1
            
        # 2. Converter para Probabilidade (Laplace)
        # Prob = (Count + 1) / (Total + 4)
        probs = {}
        denominador = num_seqs + (4 * pseudocount)
        
        for base in alfabeto:
            probs[base] = (contagens[base] + pseudocount) / denominador
            
        pwm.append(probs)
        
    return pwm

def probabilidade_seq_pwm(pwm, sequencia):
    """Calcula a probabilidade de uma sequência gerar a PWM (multiplicação)."""
    prob = 1.0
    for i, base in enumerate(sequencia):
        if i < len(pwm) and base in pwm[i]:
            prob *= pwm[i][base]
        else:
            prob = 0 # Se a base não estiver na PWM ou tamanho errado
    return prob

def subsequencia_mais_provavel(pwm, sequencia_alvo):
    """
    Desliza a PWM pela sequência alvo para achar a posição mais provável.
    Retorna (Posição_Inicial, Subsequência, Probabilidade).
    """
    tamanho_motif = len(pwm)
    melhor_prob = -1
    melhor_pos = -1
    melhor_sub = ""
    
    # Janela deslizante
    for i in range(len(sequencia_alvo) - tamanho_motif + 1):
        subseq = sequencia_alvo[i : i + tamanho_motif]
        prob = probabilidade_seq_pwm(pwm, subseq)
        
        if prob > melhor_prob:
            melhor_prob = prob
            melhor_pos = i
            melhor_sub = subseq
            
    return melhor_pos, melhor_sub, melhor_prob