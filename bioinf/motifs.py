import re
import math


def prosite_para_regex(padrao):
    """
    Converte um padrão PROSITE para expressão regular.
    Suporta x, { }, (n), (n,m), < >
    """
    regex = padrao.replace('-', '')
    regex = regex.replace('x', '.')
    regex = regex.replace('{', '[^').replace('}', ']')

    regex = re.sub(r'\((\d+),(\d+)\)', r'{\1,\2}', regex)
    regex = re.sub(r'\((\d+)\)', r'{\1}', regex)

    if regex.startswith('<'):
        regex = '^' + regex[1:]
    if regex.endswith('>'):
        regex = regex[:-1] + '$'

    return regex


def procurar_motifs(sequencia, padrao_prosite):
    """
    Procura posições de motifs na sequência usando padrão PROSITE.
    """
    regex = prosite_para_regex(padrao_prosite)
    return [m.start() for m in re.finditer(regex, sequencia)]


def enzima_para_regex(sitio):
    """
    Converte sítio de enzima (ex: G^AATTC) em sequência e posição de corte.
    """
    corte = sitio.index('^')
    seq = sitio.replace('^', '')
    return seq, corte

def fragmentar_dna(sequencia, sitio_enzima):
    """
    Fragmenta uma sequência de DNA com base num sítio de enzima de restrição.
    """
    padrao, corte = enzima_para_regex(sitio_enzima)
    posicoes = [m.start() for m in re.finditer(padrao, sequencia)]

    cortes = [p + corte for p in posicoes]
    fragmentos = []
    inicio = 0

    for c in cortes:
        fragmentos.append(sequencia[inicio:c])
        inicio = c

    # Adiciona o fragmento restante, mesmo que vazio
    fragmentos.append(sequencia[inicio:])

    return fragmentos, cortes



def criar_pwm(lista_seqs, pseudocount=1):
    """
    Cria uma PWM a partir de uma lista de sequências.
    """
    if not lista_seqs:
        return []

    alfabeto = "ACGT"
    comprimento = len(lista_seqs[0])
    pwm = []

    for i in range(comprimento):
        contagens = {b: pseudocount for b in alfabeto}
        for seq in lista_seqs:
            contagens[seq[i]] += 1

        total = sum(contagens.values())
        pwm.append({b: contagens[b] / total for b in alfabeto})

    return pwm


def probabilidade_seq_pwm(pwm, sequencia):
    """
    Calcula a probabilidade de uma sequência gerar a PWM.
    """
    prob = 1.0
    for i, base in enumerate(sequencia):
        prob *= pwm[i].get(base, 0)
    return prob


def subsequencia_mais_provavel(pwm, sequencia_alvo):
    """
    Encontra a subsequência mais provável na sequência alvo com base na PWM.
    """
    tamanho = len(pwm)
    melhor_prob = -1
    melhor_pos = -1
    melhor_sub = ""

    for i in range(len(sequencia_alvo) - tamanho + 1):
        sub = sequencia_alvo[i:i+tamanho]
        p = probabilidade_seq_pwm(pwm, sub)
        if p > melhor_prob:
            melhor_prob = p
            melhor_pos = i
            melhor_sub = sub

    return melhor_pos, melhor_sub, melhor_prob


def pwm_para_pssm(pwm, bg=None):
    """
    Converte PWM em PSSM usando log2 odds.
    """
    if bg is None:
        bg = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}

    pssm = []
    for coluna in pwm:
        pssm.append({
            base: math.log2(coluna[base] / bg[base])
            for base in coluna
        })
    return pssm


def score_seq_pssm(pssm, sequencia):
    """
    Calcula o score de uma sequência usando PSSM.
    """
    score = 0
    for i, base in enumerate(sequencia):
        score += pssm[i].get(base, -float('inf'))
    return score
