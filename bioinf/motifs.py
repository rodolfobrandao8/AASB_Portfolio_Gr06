import re
import math


def prosite_para_regex(padrao):
    """Converte um padrão PROSITE para uma expressão regular (regex) equivalente.

    Implementa uma conversão simples para suportar construções comuns de PROSITE:
    - ``-`` (separador) é removido;
    - ``x`` vira ``.`` (qualquer carácter);
    - ``{...}`` vira classe negada ``[^...]``;
    - ``(n)`` e ``(n,m)`` viram quantificadores regex ``{n}`` e ``{n,m}``;
    - ``<`` no início ancora ao início da string (``^``);
    - ``>`` no fim ancora ao fim da string (``$``).

    Args:
        padrao (str): Padrão PROSITE (ex.: ``"C-x(2)-C"``).

    Returns:
        str: Expressão regular compatível com o módulo :mod:`re`.

    Raises:
        TypeError: Se ``padrao`` não for string.

    Examples:
        >>> prosite_para_regex("C-x(2)-C")
        'C.{2}C'
        >>> prosite_para_regex("<A-x(3)-G>")
        '^A.{3}G$'
        >>> prosite_para_regex("C-{GP}-x-C")
        'C[^GP].C'
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
    """Procura ocorrências de um motif PROSITE numa sequência.

    Converte o padrão PROSITE para regex (via :func:`prosite_para_regex`) e usa
    :func:`re.finditer` para obter as posições iniciais (0-based) de cada match.

    Args:
        sequencia (str): Sequência onde procurar (DNA/RNA/proteína, dependendo do padrão).
        padrao_prosite (str): Padrão no formato PROSITE.

    Returns:
        list[int]: Lista de posições iniciais (0-based) onde o motif ocorre.

    Raises:
        re.error: Se o padrão convertido gerar uma regex inválida.
        TypeError: Se os argumentos não forem strings.

    Examples:
        >>> procurar_motifs("ACCCAC", "C-x(2)-C")
        [1]
    """
    regex = prosite_para_regex(padrao_prosite)
    return [m.start() for m in re.finditer(regex, sequencia)]


def enzima_para_regex(sitio):
    """Interpreta um sítio de restrição com marca de corte.

    Recebe uma string com o símbolo ``^`` a indicar a posição de corte
    (ex.: ``"G^AATTC"``) e devolve:
    - a sequência do sítio sem ``^``,
    - a posição de corte relativa ao início do sítio (0-based).

    Nota:
        Apesar do nome, esta função **não** converte códigos de ambiguidade para
        regex; devolve apenas a sequência literal e o índice do corte.

    Args:
        sitio (str): Sítio de enzima com ``^`` (ex.: ``"G^AATTC"``).

    Returns:
        tuple[str, int]: ``(seq, corte)``, onde:
        - ``seq`` é o sítio sem ``^``,
        - ``corte`` é o índice (0-based) em ``seq`` onde ocorre o corte.

    Raises:
        ValueError: Se ``'^'`` não existir em ``sitio`` (por causa de ``index``).
        TypeError: Se ``sitio`` não for string.

    Examples:
        >>> enzima_para_regex("G^AATTC")
        ('GAATTC', 1)
    """
    corte = sitio.index('^')
    seq = sitio.replace('^', '')
    return seq, corte


def fragmentar_dna(sequencia, sitio_enzima):
    """Fragmenta uma sequência de DNA com base num sítio de restrição.

    O processo é:
    1) Extrair o padrão do sítio e a posição de corte (via :func:`enzima_para_regex`).
    2) Encontrar todas as ocorrências do padrão na sequência (como *match literal*).
    3) Gerar as posições de corte (offset + corte) e devolver os fragmentos.

    Args:
        sequencia (str): Sequência de DNA a fragmentar.
        sitio_enzima (str): Sítio com marca de corte ``^`` (ex.: ``"G^AATTC"``).

    Returns:
        tuple[list[str], list[int]]: ``(fragmentos, cortes)``, onde:
        - ``fragmentos`` é a lista de fragmentos resultantes (strings),
        - ``cortes`` é a lista de posições de corte (0-based) na sequência original.

    Raises:
        ValueError: Se ``sitio_enzima`` não contiver ``^``.
        re.error: Se o padrão produzir um erro de regex (pouco provável aqui porque é literal).
        TypeError: Se os argumentos não forem strings.

    Examples:
        >>> fragmentar_dna("TTGAATTCAA", "G^AATTC")
        (['TTG', 'AATTCAA'], [3])
    """
    padrao, corte = enzima_para_regex(sitio_enzima)
    posicoes = [m.start() for m in re.finditer(padrao, sequencia)]
    cortes = [p + corte for p in posicoes]

    fragmentos = []
    inicio = 0

    for c in cortes:
        fragmentos.append(sequencia[inicio:c])
        inicio = c

    fragmentos.append(sequencia[inicio:] if inicio < len(sequencia) else '')

    return fragmentos, cortes


def criar_pwm(lista_seqs, pseudocount=1):
    """Cria uma PWM (Position Weight Matrix) a partir de uma lista de sequências.

    Esta implementação assume:
    - alfabeto fixo ``"ACGT"``;
    - todas as sequências têm o mesmo comprimento (usa o comprimento da primeira);
    - aplica pseudocounts para evitar probabilidades zero.

    Args:
        lista_seqs (list[str]): Lista de sequências (mesmo comprimento).
        pseudocount (int, optional): Valor inicial somado às contagens de cada base
            em cada coluna. Por omissão ``1``.

    Returns:
        list[dict[str, float]]: PWM representada como lista de colunas; cada coluna
        é um dicionário ``{'A':pA, 'C':pC, 'G':pG, 'T':pT}``.
        Se ``lista_seqs`` estiver vazia, devolve ``[]``.

    Raises:
        IndexError: Se existirem sequências com comprimento menor do que a primeira.
        KeyError: Se aparecerem símbolos fora de A/C/G/T (ao indexar contagens).
        TypeError: Se ``lista_seqs`` não for iterável.

    Examples:
        >>> criar_pwm(["AAA", "AAT"], pseudocount=1)[2]["A"] > criar_pwm(["AAA", "AAT"], 1)[2]["T"]
        True
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
    """Calcula a probabilidade de uma sequência segundo uma PWM.

    A probabilidade é o produto das probabilidades por posição:
    ``prod_i PWM[i][base_i]``.
    Se uma base não existir na coluna, considera-se probabilidade 0.

    Args:
        pwm (list[dict[str, float]]): PWM no formato devolvido por :func:`criar_pwm`.
        sequencia (str): Sequência (tipicamente com comprimento igual ao da PWM).

    Returns:
        float: Probabilidade (pode ser 0.0).

    Raises:
        IndexError: Se ``sequencia`` for maior do que o comprimento do PWM.
        TypeError: Se ``pwm`` não for uma lista de dicionários.

    Examples:
        >>> pwm = criar_pwm(["AAAA", "AAAT"], pseudocount=1)
        >>> probabilidade_seq_pwm(pwm, "AAAA") >= probabilidade_seq_pwm(pwm, "TTTT")
        True
    """
    prob = 1.0
    for i, base in enumerate(sequencia):
        prob *= pwm[i].get(base, 0)
    return prob


def subsequencia_mais_provavel(pwm, sequencia_alvo):
    """Encontra a subsequência mais provável numa sequência alvo dada uma PWM.

    Varre todas as janelas contíguas de tamanho ``len(pwm)`` e calcula a
    probabilidade de cada janela com :func:`probabilidade_seq_pwm`, devolvendo
    a melhor.

    Args:
        pwm (list[dict[str, float]]): PWM (lista de colunas).
        sequencia_alvo (str): Sequência onde procurar (target).

    Returns:
        tuple[int, str, float]: ``(pos, subseq, prob)``, onde:
        - ``pos`` é a posição inicial (0-based) da melhor subsequência,
        - ``subseq`` é a subsequência encontrada,
        - ``prob`` é a probabilidade dessa subsequência.

    Raises:
        ValueError: Se ``pwm`` estiver vazio (o código assume tamanho > 0).
        IndexError: Se ``sequencia_alvo`` for menor do que ``len(pwm)``.
        TypeError: Se os argumentos forem de tipos inválidos.

    Examples:
        >>> pwm = criar_pwm(["AAA", "AAT"], pseudocount=1)
        >>> subsequencia_mais_provavel(pwm, "TTTAAATTT")[0]
        3
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
    """Converte uma PWM numa PSSM usando log2-odds.

    Para cada posição e base calcula:
    ``log2(P(base|posição) / P(base|background))``.

    Args:
        pwm (list[dict[str, float]]): PWM (lista de colunas).
        bg (dict[str, float] | None, optional): Distribuição de background.
            Se ``None``, assume uniforme (0.25 para A/C/G/T).

    Returns:
        list[dict[str, float]]: PSSM no mesmo formato (lista de colunas), mas com
        scores log-odds em vez de probabilidades.

    Raises:
        KeyError: Se ``bg`` não tiver uma base presente no PWM.
        ValueError: Se alguma probabilidade for 0 (pode gerar divisão por zero).
            (Nota: com pseudocounts > 0 isto é menos provável.)
        TypeError: Se ``pwm``/``bg`` forem de tipos inválidos.

    Examples:
        >>> pwm = criar_pwm(["AAAA", "AAAT"], pseudocount=1)
        >>> pssm = pwm_para_pssm(pwm)
        >>> isinstance(pssm[0]["A"], float)
        True
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
    """Calcula o score total de uma sequência segundo uma PSSM.

    Soma os scores por posição. Se uma base não existir numa coluna da PSSM,
    atribui ``-inf`` a essa posição (via ``-float('inf')``), tornando o score final
    muito baixo.

    Args:
        pssm (list[dict[str, float]]): PSSM no formato devolvido por :func:`pwm_para_pssm`.
        sequencia (str): Sequência a pontuar (tipicamente com comprimento igual ao da PSSM).

    Returns:
        float: Score total (pode ser ``-inf`` se houver bases desconhecidas).

    Raises:
        IndexError: Se ``sequencia`` for maior do que o comprimento da PSSM.
        TypeError: Se ``pssm`` não for uma lista de dicionários.

    Examples:
        >>> pwm = criar_pwm(["AAAA", "AAAT"], pseudocount=1)
        >>> pssm = pwm_para_pssm(pwm)
        >>> score_seq_pssm(pssm, "AAAA") >= score_seq_pssm(pssm, "TTTT")
        True
    """
    score = 0
    for i, base in enumerate(sequencia):
        score += pssm[i].get(base, -float('inf'))
    return score
